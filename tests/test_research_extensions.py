"""Regression tests for safer retrieval and optional research workflows."""
from argparse import Namespace
from pathlib import Path
import numpy as np
import pytest
import rasterio
from rasterio.transform import from_origin
import lst_universal_research_safe as lst
import lst_research_extensions as research


def raster(path, x, *, scale=30.0):
    with rasterio.open(path, "w", driver="GTiff", width=x.shape[1], height=x.shape[0],
                       count=1, dtype="float32", crs="EPSG:32646",
                       transform=from_origin(400000, 2700000, scale, scale),
                       nodata=-9999.0) as dst:
        dst.write(np.asarray(x, dtype="float32"), 1)


def test_water_vapor_extrapolation_requires_opt_in(tmp_path):
    args = Namespace(water_vapor=5.0, water_vapor_raster=None, non_interactive=True,
                     allow_wv_extrapolation=False)
    with pytest.raises(ValueError, match="exceeds"):
        lst.get_water_vapor(None, args)
    args.allow_wv_extrapolation = True
    assert lst.get_water_vapor(None, args) == 5.0


def test_water_vapor_negative_and_nonfinite():
    for value in (-1.0, float("nan")):
        args = Namespace(water_vapor=value, water_vapor_raster=None,
                         non_interactive=True, allow_wv_extrapolation=False)
        with pytest.raises(ValueError):
            lst.get_water_vapor(None, args)


def test_odd_station_window_required(tmp_path):
    path = tmp_path / "surface.tif"
    raster(path, np.ones((10, 10)) * 27)
    obs = tmp_path / "obs.csv"
    obs.write_text("lon,lat,temp_c\n91,24,27\n")
    with rasterio.open(path) as src:
        with pytest.raises(ValueError, match="odd"):
            lst.validate_insitu(np.ones((10, 10)), src, obs, 4,
                               tmp_path / "points.csv")


def test_trend_constant_slope(tmp_path):
    rasters = []
    for year, value in [(2020, 25), (2021, 27), (2022, 29)]:
        p = tmp_path / f"scene_{year}.tif"
        raster(p, np.full((4, 4), value, dtype=float))
        rasters.append(p)
    args = Namespace(rasters=rasters, dates=["2020-01-01", "2021-01-01", "2022-01-01"],
                     min_observations=3, output_dir=tmp_path / "out")
    assert research.trend(args) == 0
    with rasterio.open(args.output_dir / "lst_theil_sen_c_per_year.tif") as src:
        values = src.read(1)
        assert np.allclose(values, 2.0, atol=0.01)


def test_downscale_spatial_holdout(tmp_path):
    sklearn = pytest.importorskip("sklearn")
    yy, xx = np.indices((30, 30))
    ndvi = (xx + yy) / 100
    ndbi = xx / 100
    temp = 25 + 8 * ndvi - 2 * ndbi
    coarse_files = []
    fine_files = []
    for name, array in (("ndvi", ndvi), ("ndbi", ndbi)):
        coarse = tmp_path / f"coarse_{name}.tif"
        fine = tmp_path / f"fine_{name}.tif"
        raster(coarse, array)
        raster(fine, array)
        coarse_files.append(coarse)
        fine_files.append(fine)
    target = tmp_path / "coarse_lst.tif"
    raster(target, temp)
    out = tmp_path / "prediction.tif"
    args = Namespace(lst=target, coarse=coarse_files, fine=fine_files,
                     output=out, model="ridge", block_size=5, test_size=0.25,
                     seed=42, residual_correction=False, model_output=None)
    assert research.downscale(args) == 0
    with rasterio.open(out) as src:
        prediction = src.read(1)
        assert np.max(np.abs(prediction - temp)) < 1.0
        assert src.tags()["PRODUCT"] == "MODELLED_DOWNSCALED_LST"
