"""Tests for lst_universal_research_safe.py. Run:  pytest -q"""
import csv
from pathlib import Path

import numpy as np
import pytest
import rasterio
from rasterio.transform import from_origin

import lst_universal_research_safe as lst

# Real Landsat 8 band-10 constants (as printed in Collection 2 MTL files)
K1, K2, MULT, ADD = 774.8853, 1321.0789, 3.3420e-4, 0.1
SUN_ELEV = 60.0


def dn_for_bt(bt_k: float) -> int:
    """DN a real Landsat 8 band-10 pixel would have for a given brightness temperature."""
    radiance = K1 / (np.exp(K2 / bt_k) - 1.0)
    return int(round((radiance - ADD) / MULT))


def expected_lst_c(bt_k: float, eps: float, wl: float = 10.904) -> float:
    return bt_k / (1 + (wl * bt_k / 14380.0) * np.log(eps)) - 273.15


def write_tif(path: Path, arr, nodata=None):
    arr = np.asarray(arr)
    with rasterio.open(path, "w", driver="GTiff", height=arr.shape[0], width=arr.shape[1],
                       count=1, dtype=arr.dtype, crs="EPSG:32646",
                       transform=from_origin(600000, 2800000, 30, 30), nodata=nodata) as d:
        d.write(arr, 1)


def read_tif(path: Path):
    with rasterio.open(path) as s:
        return s.read(1), s.nodata


MTL_TEXT = f"""GROUP = L1_METADATA_FILE
  SPACECRAFT_ID = "LANDSAT_8"
  SUN_ELEVATION = {SUN_ELEV}
  RADIANCE_MULT_BAND_10 = 3.3420E-04
  RADIANCE_ADD_BAND_10 = 0.10000
  K1_CONSTANT_BAND_10 = 774.8853
  K2_CONSTANT_BAND_10 = 1321.0789
  REFLECTANCE_MULT_BAND_4 = 2.0000E-05
  REFLECTANCE_ADD_BAND_4 = -0.100000
  REFLECTANCE_MULT_BAND_5 = 2.0000E-05
  REFLECTANCE_ADD_BAND_5 = -0.100000
END_GROUP = L1_METADATA_FILE
"""


def refl_to_dn(rho: float) -> int:
    return int(round((rho * np.sin(np.radians(SUN_ELEV)) + 0.1) / 2e-5))


def scene_dn(bt_k=305.0, size=20):
    dn = np.full((size, size), dn_for_bt(bt_k), dtype=np.uint16)
    dn[:3, :] = 0                      # fill border like a real scene
    return dn


# ---------------------------------------------------------------- Landsat RAW
def test_landsat8_realistic_dn_with_mtl(tmp_path):
    write_tif(tmp_path / "LC08_X_B10.TIF", scene_dn(), nodata=0)
    (tmp_path / "LC08_X_MTL.txt").write_text(MTL_TEXT)
    rc = lst.main(["-i", str(tmp_path / "LC08_X_B10.TIF"), "--mtl", str(tmp_path / "LC08_X_MTL.txt"),
                   "--emissivity", "0.98", "-o", str(tmp_path / "out"), "-q"])
    assert rc == 0
    out, nodata = read_tif(next((tmp_path / "out").glob("LST_LANDSAT_8_*.tif")))
    valid = out[out != nodata]
    assert valid.size == 17 * 20                      # fill border removed
    assert valid.mean() == pytest.approx(expected_lst_c(305.0, 0.98), abs=0.1)
    assert 30 < valid.mean() < 36


def test_landsat8_builtin_defaults_work_on_realistic_dn(tmp_path):
    """Regression: the old RADIANCE_MULT default was 1000x too big and dropped every pixel."""
    write_tif(tmp_path / "LC08_X_B10.TIF", scene_dn(), nodata=0)
    rc = lst.main(["-i", str(tmp_path / "LC08_X_B10.TIF"), "--emissivity", "0.98",
                   "--non-interactive", "-o", str(tmp_path / "out"), "-q"])
    assert rc == 0


def test_landsat9_without_mtl_fails_instead_of_guessing(tmp_path):
    write_tif(tmp_path / "LC09_X_B10.TIF", scene_dn(), nodata=0)
    rc = lst.main(["-i", str(tmp_path / "LC09_X_B10.TIF"), "--emissivity", "0.98",
                   "--non-interactive", "-o", str(tmp_path / "out"), "-q"])
    assert rc == 1
    assert not (tmp_path / "out").exists() or not list((tmp_path / "out").glob("*.tif"))


def test_ndvi_emissivity_from_dn_uses_reflectance(tmp_path):
    n = 20
    write_tif(tmp_path / "LC08_X_B10.TIF", scene_dn(size=n), nodata=0)
    write_tif(tmp_path / "LC08_X_B4.TIF", np.full((n, n), refl_to_dn(0.05), np.uint16))
    write_tif(tmp_path / "LC08_X_B5.TIF", np.full((n, n), refl_to_dn(0.40), np.uint16))
    (tmp_path / "LC08_X_MTL.txt").write_text(MTL_TEXT)
    rc = lst.main(["-i", str(tmp_path / "LC08_X_B10.TIF"), "--mtl", str(tmp_path / "LC08_X_MTL.txt"),
                   "--red", str(tmp_path / "LC08_X_B4.TIF"), "--nir", str(tmp_path / "LC08_X_B5.TIF"),
                   "-o", str(tmp_path / "out"), "-q"])
    assert rc == 0
    out, nodata = read_tif(next((tmp_path / "out").glob("*.tif")))
    # NDVI = 0.78 > 0.5 -> eps = 0.987
    assert out[out != nodata].mean() == pytest.approx(expected_lst_c(305.0, 0.987), abs=0.1)


def test_raw_red_nir_without_mtl_is_rejected(tmp_path):
    write_tif(tmp_path / "LC08_X_B10.TIF", scene_dn(), nodata=0)
    write_tif(tmp_path / "r.tif", np.full((20, 20), 9000, np.uint16))
    write_tif(tmp_path / "n.tif", np.full((20, 20), 9000, np.uint16))
    rc = lst.main(["-i", str(tmp_path / "LC08_X_B10.TIF"), "--red", str(tmp_path / "r.tif"),
                   "--nir", str(tmp_path / "n.tif"), "--non-interactive", "-q"])
    assert rc == 1


def test_emissivity_function_is_nearly_continuous_and_in_range():
    red = np.full(200, 0.1, np.float32)
    ndvi = np.linspace(-0.2, 0.9, 200, dtype=np.float32)
    nir = red * (1 + ndvi) / (1 - ndvi)          # invert NDVI for a fixed red
    eps, nd = lst.ndvi_emissivity(red, nir)
    assert np.isfinite(eps).all()
    assert eps.min() > 0.96 and eps.max() <= 0.987 + 1e-6
    assert np.abs(np.diff(eps)).max() < 0.005    # no 0.05-sized jumps (old code had ~0.05)


def test_invalid_emissivity_rejected(tmp_path):
    write_tif(tmp_path / "LC08_X_B10.TIF", scene_dn(), nodata=0)
    assert lst.main(["-i", str(tmp_path / "LC08_X_B10.TIF"), "--emissivity", "1.5",
                     "--non-interactive", "-q"]) == 1


# ------------------------------------------------------------------ PRODUCTS
def test_modis_uint16_scale_and_valid_range(tmp_path):
    dn = np.full((10, 10), 15000, np.uint16)       # 300 K
    dn[0, :] = 0                                   # fill
    dn[1, :] = 5000                                # below valid range 7500
    write_tif(tmp_path / "MOD11A1.A2020001.h28v07_LST_Day_1km.tif", dn)
    rc = lst.main(["-i", str(tmp_path / "MOD11A1.A2020001.h28v07_LST_Day_1km.tif"),
                   "-o", str(tmp_path / "out"), "-q"])
    assert rc == 0
    out, nodata = read_tif(next((tmp_path / "out").glob("LST_MODIS_*.tif")))
    valid = out[out != nodata]
    assert valid.size == 80
    assert valid.mean() == pytest.approx(300 - 273.15, abs=1e-3)


def test_ecostress_float_kelvin_is_not_rescaled(tmp_path):
    """Regression: a float Kelvin file must not be multiplied by a default scale."""
    arr = np.full((10, 10), 300.0, np.float32)
    arr[0, :] = np.nan
    write_tif(tmp_path / "ECOv002_L2T_LST_01234_001_X.tif", arr)
    rc = lst.main(["-i", str(tmp_path / "ECOv002_L2T_LST_01234_001_X.tif"),
                   "--sensor", "ECOSTRESS", "-o", str(tmp_path / "out"), "-q"])
    assert rc == 0
    out, nodata = read_tif(next((tmp_path / "out").glob("*.tif")))
    assert out[out != nodata].mean() == pytest.approx(26.85, abs=1e-3)


def test_landsat_l2_st_scale_offset_and_detection(tmp_path):
    dn = int(round((300.0 - 149.0) / 0.00341802))
    f = tmp_path / "LC08_L2SP_137043_20200101_20200101_02_T1_ST_B10.TIF"
    write_tif(f, np.full((10, 10), dn, np.uint16), nodata=0)
    with rasterio.open(f) as s:
        assert lst.detect_sensor(s, f) == "LANDSAT_L2ST"   # not mistaken for raw LANDSAT_8
    assert lst.main(["-i", str(f), "-o", str(tmp_path / "out"), "-q"]) == 0
    out, nodata = read_tif(next((tmp_path / "out").glob("*.tif")))
    assert out[out != nodata].mean() == pytest.approx(26.85, abs=0.01)


# --------------------------------------------------------------------- BATCH
def make_scene(folder: Path, scene: str, bt: float):
    folder.mkdir(parents=True)
    write_tif(folder / f"{scene}_B10.TIF", scene_dn(bt), nodata=0)
    write_tif(folder / f"{scene}_B4.TIF", np.full((20, 20), refl_to_dn(0.05), np.uint16))
    write_tif(folder / f"{scene}_B5.TIF", np.full((20, 20), refl_to_dn(0.40), np.uint16))
    (folder / f"{scene}_MTL.txt").write_text(MTL_TEXT)


def test_batch_uses_only_thermal_band_writes_elsewhere_and_labels_csv(tmp_path):
    make_scene(tmp_path / "s1", "LC08_L1TP_A", 300.0)
    make_scene(tmp_path / "s2", "LC08_L1TP_B", 310.0)
    assert lst.main(["--batch", str(tmp_path), "-q"]) == 0
    out_dir = tmp_path / "lst_output"
    assert len(list(out_dir.glob("*.tif"))) == 2          # not 6: B4/B5 are not thermal
    rows = list(csv.DictReader(open(out_dir / "lst_stats.csv")))
    assert {r["scene"] for r in rows} == {"LC08_L1TP_A_B10", "LC08_L1TP_B_B10"}
    assert all(r["sensor"] == "LANDSAT_8" for r in rows)
    means = {r["scene"]: float(r["mean"]) for r in rows}
    assert means["LC08_L1TP_A_B10"] == pytest.approx(expected_lst_c(300.0, 0.987), abs=0.1)
    assert means["LC08_L1TP_B_B10"] == pytest.approx(expected_lst_c(310.0, 0.987), abs=0.1)
    # second run must not feed the outputs back in
    assert lst.main(["--batch", str(tmp_path), "-q"]) == 0
    names = [p.name for p in out_dir.glob("*.tif")]
    assert len(names) == 4 and not any("LST_LANDSAT_8_LST_" in n for n in names)


# ------------------------------------------------------------- ROBUSTNESS
def test_closed_stdin_exits_quickly(monkeypatch):
    def eof(_):
        raise EOFError
    monkeypatch.setattr("builtins.input", eof)
    with pytest.raises(SystemExit):
        lst.ask("number: ", float)


def test_bad_input_gives_up_after_few_tries(monkeypatch):
    monkeypatch.setattr("builtins.input", lambda _: "abc")
    with pytest.raises(SystemExit):
        lst.ask("number: ", float, tries=3)


def test_cog_output_opens_and_has_tags(tmp_path):
    write_tif(tmp_path / "LC08_X_B10.TIF", scene_dn(), nodata=0)
    assert lst.main(["-i", str(tmp_path / "LC08_X_B10.TIF"), "--emissivity", "0.98", "--cog",
                     "-o", str(tmp_path / "out"), "-q"]) == 0
    out = next((tmp_path / "out").glob("*.tif"))
    with rasterio.open(out) as s:
        assert s.tags()["LST_UNITS"] == "Celsius"
        assert "no atmospheric correction" in s.tags()["METHOD"]


def test_toa_reflectance_applies_sun_elevation(tmp_path):
    (tmp_path / "m_MTL.txt").write_text(MTL_TEXT)
    mtl = lst.parse_mtl(tmp_path / "m_MTL.txt")
    dn = np.array([[refl_to_dn(0.05), refl_to_dn(0.40), 0]], dtype=np.uint16)
    refl = lst.toa_reflectance(dn, mtl, 4)
    assert refl[0, 0] == pytest.approx(0.05, abs=2e-4)
    assert refl[0, 1] == pytest.approx(0.40, abs=2e-4)
    assert np.isnan(refl[0, 2])                       # fill stays masked
