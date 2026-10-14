"""
Test suite for the Universal LST pipeline.

These tests use synthetic rasters generated in-memory so they don't depend
on external data. They verify every bug fix from the code review.

Run:  pytest tests/ -v
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest
import rasterio
from rasterio.transform import from_origin

# Make the source module importable
REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

import lst_universal_research_safe as lst  # noqa: E402


# ------------------------------------------------------------------
# Fixtures
# ------------------------------------------------------------------
def _write_geotiff(path, arr, tags=None, nodata=None, dtype=None):
    """Write a 2D array as a single-band GeoTIFF."""
    if dtype is None:
        dtype = arr.dtype
    profile = {
        "driver": "GTiff",
        "height": arr.shape[0],
        "width": arr.shape[1],
        "count": 1,
        "dtype": dtype,
        "crs": "EPSG:32645",
        "transform": from_origin(600000, 2370000, 30, 30),
    }
    if nodata is not None:
        profile["nodata"] = nodata
    with rasterio.open(path, "w", **profile) as dst:
        dst.write(arr.astype(dtype), 1)
        if tags:
            dst.update_tags(**tags)


@pytest.fixture
def landsat8_thermal(tmp_path):
    """A synthetic Landsat 8 B10 GeoTIFF with realistic DNs."""
    H, W = 64, 64
    y, x = np.mgrid[0:H, 0:W]
    dn = 30 + 6 * np.sin(y / 20) * np.cos(x / 25)
    dn = np.clip(dn, 15, 60).astype(np.uint16)
    # Embed fill pixels (DN=0) AFTER clipping so they survive
    dn[10:15, 10:15] = 0
    path = tmp_path / "LE08_B10.tif"
    _write_geotiff(path, dn, tags={"SATELLITE": "LANDSAT_8", "PLATFORM": "LC08"},
                   nodata=0)
    return path, dn


@pytest.fixture
def red_nir_bands(tmp_path):
    """Synthetic red and NIR bands for the NDVI emissivity test."""
    H, W = 64, 64
    red = np.full((H, W), 0.05, dtype=np.float32)
    nir = np.full((H, W), 0.45, dtype=np.float32)
    # Make one half bare soil (low NDVI), one half vegetation (high NDVI)
    red[:, :W // 2] = 0.30
    nir[:, :W // 2] = 0.35
    red_p = tmp_path / "red.tif"
    nir_p = tmp_path / "nir.tif"
    _write_geotiff(red_p, red, dtype="float32")
    _write_geotiff(nir_p, nir, dtype="float32")
    return red_p, nir_p


@pytest.fixture
def modis_lst(tmp_path):
    """Synthetic MOD11A1 GeoTIFF with embedded scale_factor tag."""
    H, W = 64, 64
    target_k = 305 + np.random.normal(0, 1, (H, W))
    dn = (target_k / 0.02).astype(np.uint16)
    dn = np.clip(dn, 7500, 65535)
    path = tmp_path / "MOD11A1.tif"
    _write_geotiff(path, dn,
                   tags={"SATELLITE": "MODIS", "PRODUCT": "MOD11A1",
                         "scale_factor": "0.02", "add_offset": "0.0",
                         "UNITS": "Kelvin"},
                   nodata=0)
    return path, dn


# ------------------------------------------------------------------
# Bug-fix tests
# ------------------------------------------------------------------
class TestBugFixes:
    """Each test verifies one specific bug fix from the code review."""

    def test_bug_1_fill_pixels_masked(self, landsat8_thermal):
        """DN=0 fill pixels must NOT produce spurious LST around -125 °C."""
        path, dn = landsat8_thermal
        with rasterio.open(path) as src:
            arr = src.read(1)
        out = lst.calculate_raw_lst(
            arr.astype(np.float32),
            emissivity=0.985,
            sensor_name="LANDSAT_8",
            sensor_cfg=lst.SENSOR_REGISTRY["LANDSAT_8"],
            nodata=0,
        )
        # The fill region must be NaN, not -125 °C
        fill_region = out[10:15, 10:15]
        assert np.all(np.isnan(fill_region)), \
            f"Fill pixels leaked through: min={np.nanmin(fill_region)}"
        # And no pixel should be below -100 °C (sanity check on the rest)
        valid = out[~np.isnan(out)]
        assert np.nanmin(valid) > -50, \
            f"Implausibly low LST detected: {np.nanmin(valid)}"

    def test_bug_4_invalid_emissivity_rejected(self, landsat8_thermal):
        """ε <= 0 or ε > 1 must raise, not silently produce NaN/inf."""
        path, dn = landsat8_thermal
        with pytest.raises(ValueError):
            lst.brightness_temperature_to_lst(
                np.array([[300.0]]), emissivity=0.0, wavelength_um=10.9
            )
        with pytest.raises(ValueError):
            lst.brightness_temperature_to_lst(
                np.array([[300.0]]), emissivity=1.5, wavelength_um=10.9
            )

    def test_bug_5_sanity_filter(self):
        """Output outside [lo, hi] °C must be NaN."""
        arr = np.array([[-100.0, 50.0, 200.0]])
        filtered = lst.apply_sanity_filter(arr, lo=-50, hi=90)
        assert np.isnan(filtered[0, 0])
        assert filtered[0, 1] == 50.0
        assert np.isnan(filtered[0, 2])

    def test_bug_6_method_label_honest(self, landsat8_thermal, tmp_path):
        """METHOD tag must say 'Emissivity-Corrected BT', not 'SCA (RAW)'."""
        path, _ = landsat8_thermal
        args = lst.build_arg_parser().parse_args([
            "--input", str(path),
            "--sensor", "LANDSAT_8",
            "--emissivity", "0.985",
            "--non-interactive",
            "--output-dir", str(tmp_path),
        ])
        lst.setup_logging(quiet=True)
        stats = lst.process_file(path, args)
        assert stats is not None
        # Find the output file
        out_files = list(tmp_path.glob("LST_LANDSAT_8_*.tif"))
        assert len(out_files) == 1
        with rasterio.open(out_files[0]) as src:
            assert src.tags().get("METHOD") == "Emissivity-Corrected BT"
            assert src.tags().get("LST_UNITS") == "Celsius"

    def test_bug_8_landsat_wavelength(self):
        """Default wavelength for Landsat 8 must be 10.9, not 10.8."""
        assert lst.SENSOR_REGISTRY["LANDSAT_8"]["default_wavelength"] == 10.9
        assert lst.SENSOR_REGISTRY["LANDSAT_9"]["default_wavelength"] == 10.9

    def test_bug_9_sensor_detection_no_master_clash(self, tmp_path):
        """'MASTER' must not match 'ASTER'."""
        # Write a dummy file with MASTER tag
        path = tmp_path / "MASTER_scene.tif"
        _write_geotiff(path, np.zeros((8, 8), dtype=np.uint16),
                       tags={"SATELLITE": "MASTER"})
        with rasterio.open(path) as src:
            sensor = lst.detect_sensor(src, path)
        # Should not return ASTER for MASTER
        assert sensor != "ASTER", "ASTER matched MASTER — word-boundary bug"

    def test_bug_9_filename_detection(self, tmp_path):
        """Landsat Collection 2 filename LC08_... must be detected."""
        path = tmp_path / "LC08_L1TP_137045_20240115_B10.tif"
        _write_geotiff(path, np.zeros((8, 8), dtype=np.uint16),
                       tags={})  # no sensor tags
        with rasterio.open(path) as src:
            sensor = lst.detect_sensor(src, path)
        assert sensor == "LANDSAT_8"

    def test_bug_10_scale_from_rasterio_field(self, tmp_path):
        """Scale must be read from src.scales[0] before tags."""
        H, W = 8, 8
        arr = (300.0 / 0.02 * np.ones((H, W))).astype(np.uint16)
        arr = np.clip(arr, 7500, 65535)
        path = tmp_path / "MOD02.tif"
        profile = {
            "driver": "GTiff", "height": H, "width": W, "count": 1,
            "dtype": "uint16", "crs": "EPSG:32645",
            "transform": from_origin(600000, 2370000, 30, 30),
            "nodata": 0,
        }
        with rasterio.open(path, "w", **profile) as dst:
            dst.write(arr, 1)
            dst.scales = (0.02,)    # proper API
            dst.offsets = (0.0,)
            dst.update_tags(SATELLITE="MODIS", PRODUCT="MOD11A1",
                            UNITS="Kelvin")
        with rasterio.open(path) as src:
            assert src.scales[0] == 0.02, "Test setup: scale not persisted"
            result = lst.calculate_product_lst(
                src, "MODIS", lst.SENSOR_REGISTRY["MODIS"]
            )
        # 300 K -> 26.85 C
        assert abs(np.nanmean(result) - 26.85) < 0.5

    def test_bug_12_unique_output_filename(self, landsat8_thermal, tmp_path):
        """Two runs must not overwrite each other."""
        path, _ = landsat8_thermal
        args = lst.build_arg_parser().parse_args([
            "--input", str(path),
            "--sensor", "LANDSAT_8",
            "--emissivity", "0.985",
            "--non-interactive",
            "--output-dir", str(tmp_path),
        ])
        lst.setup_logging(quiet=True)
        lst.process_file(path, args)
        lst.process_file(path, args)
        out_files = list(tmp_path.glob("LST_LANDSAT_8_*.tif"))
        assert len(out_files) == 2, \
            f"Expected 2 outputs, got {len(out_files)} (silent overwrite)"

    def test_bug_2_ecostress_units_check(self, tmp_path):
        """ECOSTRESS file already in Celsius must not have 273.15 subtracted."""
        H, W = 8, 8
        # Pretend this ECOSTRESS file already stores LST in Celsius, DN=30°C
        arr = np.full((H, W), 30, dtype=np.uint16)
        path = tmp_path / "ECO_celsius.tif"
        _write_geotiff(path, arr,
                       tags={"SATELLITE": "ECOSTRESS",
                             "scale_factor": "1.0",
                             "UNITS": "Celsius"})
        with rasterio.open(path) as src:
            result = lst.calculate_product_lst(
                src, "ECOSTRESS", lst.SENSOR_REGISTRY["ECOSTRESS"]
            )
        assert abs(np.nanmean(result) - 30.0) < 0.1


# ------------------------------------------------------------------
# Functionality tests
# ------------------------------------------------------------------
class TestNewFunctionality:

    def test_ndvi_emissivity(self, red_nir_bands):
        """NDVI-based emissivity must produce sensible per-pixel values."""
        red_p, nir_p = red_nir_bands
        with rasterio.open(red_p) as rds:
            red = rds.read(1).astype(np.float32)
        with rasterio.open(nir_p) as nds:
            nir = nds.read(1).astype(np.float32)
        eps, ndvi = lst.compute_ndvi_emissivity(red, nir)
        assert eps.shape == red.shape
        # All values in [0.9, 1.0]
        valid = eps[~np.isnan(eps)]
        assert np.all(valid >= 0.9) and np.all(valid <= 1.0)
        # Vegetation side has higher ε than bare-soil side
        veg_eps = np.nanmean(eps[:, eps.shape[1] // 2:])
        soil_eps = np.nanmean(eps[:, :eps.shape[1] // 2])
        assert veg_eps > soil_eps

    def test_mtl_parser(self, tmp_path):
        """MTL parser must extract RADIANCE_MULT, K1, K2 etc."""
        mtl_text = """
        GROUP = L1_METADATA_FILE
        RADIANCE_MULT_BAND_10 = 3.3420E-01
        RADIANCE_ADD_BAND_10 = 0.10000
        K1_CONSTANT_BAND_10 = 774.8853
        K2_CONSTANT_BAND_10 = 1321.0789
        END_GROUP = L1_METADATA_FILE
        """
        mtl_path = tmp_path / "LE08_MTL.txt"
        mtl_path.write_text(mtl_text)
        mtl = lst.parse_mtl(mtl_path)
        consts = lst.landsat_constants_from_mtl(mtl, band=10)
        assert abs(consts["rad_mult"] - 0.33420) < 1e-6
        assert abs(consts["k1"] - 774.8853) < 1e-4
        assert abs(consts["k2"] - 1321.0789) < 1e-4

    def test_cli_args_parsed(self):
        """argparse must accept all expected flags."""
        args = lst.build_arg_parser().parse_args([
            "--input", "x.tif", "--sensor", "LANDSAT_8",
            "--emissivity", "0.985", "--mtl", "x_MTL.txt",
            "--red", "b4.tif", "--nir", "b5.tif",
            "--plot", "--csv", "stats.csv", "--cog",
        ])
        assert args.sensor == "LANDSAT_8"
        assert args.emissivity == 0.985
        assert args.plot is True
        assert args.cog is True

    def test_batch_mode(self, tmp_path):
        """Batch mode must process all .tif files in a folder."""
        # Create 3 dummy Landsat scenes
        for i in range(3):
            H, W = 16, 16
            dn = np.full((H, W), 30 + i, dtype=np.uint16)
            path = tmp_path / f"LE08_B10_{i}.tif"
            _write_geotiff(path, dn, tags={"SATELLITE": "LANDSAT_8",
                                            "PLATFORM": "LC08"})
        out_dir = tmp_path / "out"
        args = lst.build_arg_parser().parse_args([
            "--batch", str(tmp_path),
            "--sensor", "LANDSAT_8",
            "--emissivity", "0.985",
            "--non-interactive",
            "--output-dir", str(out_dir),
        ])
        lst.setup_logging(quiet=True)
        rc = lst.main([
            "--batch", str(tmp_path),
            "--sensor", "LANDSAT_8",
            "--emissivity", "0.985",
            "--non-interactive",
            "--output-dir", str(out_dir),
        ])
        assert rc == 0
        outputs = list(out_dir.glob("LST_LANDSAT_8_*.tif"))
        assert len(outputs) == 3

    def test_end_to_end_landsat_with_ndvi(self, landsat8_thermal,
                                          red_nir_bands, tmp_path):
        """Full Landsat 8 pipeline with NDVI emissivity must produce LST
        in the plausible 0-50 °C range."""
        thermal_path, _ = landsat8_thermal
        red_p, nir_p = red_nir_bands
        out_dir = tmp_path / "out"
        rc = lst.main([
            "--input", str(thermal_path),
            "--sensor", "LANDSAT_8",
            "--red", str(red_p),
            "--nir", str(nir_p),
            "--non-interactive",
            "--output-dir", str(out_dir),
            "--quiet",
        ])
        assert rc == 0
        out_files = list(out_dir.glob("LST_LANDSAT_8_*.tif"))
        assert len(out_files) == 1
        with rasterio.open(out_files[0]) as src:
            arr = src.read(1).astype(np.float32)
            arr[arr == -9999] = np.nan
        valid = arr[~np.isnan(arr)]
        assert 0 < np.nanmean(valid) < 50, \
            f"LST mean {np.nanmean(valid)} outside plausible range"

    def test_end_to_end_modis(self, modis_lst, tmp_path):
        """Full MODIS pipeline must produce LST around 30 °C."""
        path, _ = modis_lst
        out_dir = tmp_path / "out"
        rc = lst.main([
            "--input", str(path),
            "--sensor", "MODIS",
            "--non-interactive",
            "--output-dir", str(out_dir),
            "--quiet",
        ])
        assert rc == 0
        out_files = list(out_dir.glob("LST_MODIS_*.tif"))
        assert len(out_files) == 1
        with rasterio.open(out_files[0]) as src:
            arr = src.read(1).astype(np.float32)
            arr[arr == -9999] = np.nan
        valid = arr[~np.isnan(arr)]
        # 305 K -> 31.85 C
        assert 25 < np.nanmean(valid) < 40, \
            f"MODIS LST mean {np.nanmean(valid)} outside plausible range"
