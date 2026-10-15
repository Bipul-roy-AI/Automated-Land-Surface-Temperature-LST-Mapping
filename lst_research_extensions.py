"""Optional scientific LST research extensions.

Inspired by the workflows of pansism/downscale-satelliteLST,
leonsnill/lst_landsat and wqlin08/LST-Downscaling-Project.
No third-party source code is copied.

All rasters must be co-registered before use. Downscaled maps are predictions,
not directly observed fine-resolution thermal measurements.
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.warp import reproject


def _raster(path):
    with rasterio.open(path) as src:
        x = src.read(1, masked=True).astype(np.float64).filled(np.nan)
        profile = src.profile.copy()
        grid = (src.crs, src.transform, src.shape)
    return x, profile, grid


def _write(path, data, profile, tags):
    dst_profile = profile.copy()
    for key in ("blockxsize", "blockysize", "tiled", "photometric", "predictor"):
        dst_profile.pop(key, None)
    dst_profile.update(driver="GTiff", count=1, dtype="float32",
                       nodata=-9999.0, compress="deflate")
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(path, "w", **dst_profile) as dst:
        dst.write(np.where(np.isfinite(data), data, -9999.0).astype("float32"), 1)
        dst.update_tags(**tags)


def _align(path, grid, resampling=Resampling.average):
    crs, transform, shape = grid
    with rasterio.open(path) as src:
        arr = src.read(1, masked=True).astype(np.float64).filled(np.nan)
        if src.crs is None or crs is None:
            raise ValueError("Input raster has no coordinate reference system")
        if src.crs == crs and src.shape == shape and src.transform.almost_equals(transform):
            return arr
        dst = np.full(shape, np.nan, dtype=np.float64)
        reproject(arr, dst, src_transform=src.transform, src_crs=src.crs,
                  dst_transform=transform, dst_crs=crs,
                  src_nodata=np.nan, dst_nodata=np.nan,
                  resampling=resampling)
    return dst


def _metrics(y, prediction):
    from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
    return dict(n=int(y.size), rmse=float(np.sqrt(mean_squared_error(y, prediction))),
                mae=float(mean_absolute_error(y, prediction)),
                bias=float(np.mean(prediction - y)),
                r2=float(r2_score(y, prediction)))


def downscale(args):
    """Train on coarse predictors, evaluate on held-out spatial blocks, predict at fine scale."""
    from sklearn.ensemble import RandomForestRegressor, StackingRegressor, HistGradientBoostingRegressor
    from sklearn.linear_model import Ridge
    from sklearn.model_selection import GroupShuffleSplit
    import joblib

    coarse, _, coarse_grid = _raster(args.lst)
    fine_first, fine_profile, fine_grid = _raster(args.fine[0])
    if not args.coarse or len(args.coarse) != len(args.fine):
        raise ValueError("Provide matching --coarse and --fine predictor lists (same order).")
    if any(Path(p) == Path(args.lst) for p in args.coarse):
        raise ValueError("Target LST must not be used as a predictor.")
    Xcoarse = np.stack([_align(p, coarse_grid) for p in args.coarse], axis=-1)
    Xfine = np.stack([_align(p, fine_grid, Resampling.bilinear) for p in args.fine], axis=-1)
    valid = np.isfinite(coarse) & np.all(np.isfinite(Xcoarse), axis=-1)
    rr, cc = np.indices(coarse.shape)
    group = (rr // args.block_size) * (1 + coarse.shape[1] // args.block_size) + cc // args.block_size
    X = Xcoarse[valid]
    y = coarse[valid]
    groups = group[valid]
    if X.shape[0] < 30 or len(np.unique(groups)) < 3:
        raise ValueError("Insufficient valid samples or spatial blocks for evaluation.")
    tr, te = next(GroupShuffleSplit(n_splits=1, test_size=args.test_size,
                                   random_state=args.seed).split(X, y, groups))
    models = {"ridge": Ridge(alpha=1.0),
              "rf": RandomForestRegressor(n_estimators=200, min_samples_leaf=3,
                                          n_jobs=-1, random_state=args.seed),
              "boost": HistGradientBoostingRegressor(random_state=args.seed)}
    if args.model == "stack":
        model = StackingRegressor(
            estimators=[("ridge", Ridge(alpha=1.0)),
                        ("rf", RandomForestRegressor(n_estimators=100,
                                                     min_samples_leaf=3,
                                                     random_state=args.seed, n_jobs=-1))],
            final_estimator=Ridge(alpha=1.0), cv=5)
    else:
        model = models[args.model]
    model.fit(X[tr], y[tr])
    report = _metrics(y[te], model.predict(X[te]))
    report.update(model=args.model, random_seed=args.seed,
                  train_pixels=int(len(tr)), test_pixels=int(len(te)),
                  block_size_pixels=args.block_size,
                  note="Spatial-block holdout at the input LST scale; not independent fine-resolution validation.")
    model.fit(X, y)
    fine_valid = np.all(np.isfinite(Xfine), axis=-1)
    fine_pred = np.full(fine_grid[2], np.nan, dtype=np.float64)
    fine_pred[fine_valid] = model.predict(Xfine[fine_valid])
    if args.residual_correction:
        # Add coarse-scale model errors, interpolated back to fine grid.
        coarse_pred = np.full(coarse.shape, np.nan, dtype=np.float64)
        coarse_pred[valid] = model.predict(Xcoarse[valid])
        residual = coarse - coarse_pred
        src_crs, src_transform, _ = coarse_grid
        fine_crs, fine_transform, _ = fine_grid
        correction = np.full(fine_grid[2], np.nan, dtype=np.float64)
        reproject(residual, correction, src_transform=src_transform, src_crs=src_crs,
                  dst_transform=fine_transform, dst_crs=fine_crs,
                  src_nodata=np.nan, dst_nodata=np.nan, resampling=Resampling.bilinear)
        use = fine_valid & np.isfinite(correction)
        fine_pred[use] += correction[use]
        report["residual_correction"] = True
    else:
        report["residual_correction"] = False
    _write(args.output, fine_pred, fine_profile,
           {"PRODUCT": "MODELLED_DOWNSCALED_LST", "MODEL": args.model,
            "REFERENCE_LST": Path(args.lst).name,
            "LIMITATION": "Predicted fine-scale LST; not directly measured"})
    with open(str(args.output) + ".json", "w") as f:
        json.dump(report, f, indent=2)
    if args.model_output:
        joblib.dump(model, args.model_output)
    print(json.dumps(report, indent=2))
    return 0


def trend(args):
    """Theil-Sen LST slope and temporal summaries per pixel, using dates as elapsed years."""
    from scipy.stats import theilslopes
    from datetime import date
    if len(args.rasters) != len(args.dates) or len(args.rasters) < 3:
        raise ValueError("Provide at least three matching --rasters and --dates.")
    dates = [date.fromisoformat(d) for d in args.dates]
    t = np.array([(d - dates[0]).days / 365.2425 for d in dates], dtype=float)
    if len(set(dates)) != len(dates):
        raise ValueError("Dates must be unique.")
    first, profile, grid = _raster(args.rasters[0])
    layers = np.stack([first] + [_align(p, grid, Resampling.nearest)
                                 for p in args.rasters[1:]])
    n = np.sum(np.isfinite(layers), axis=0)
    slope = np.full(first.shape, np.nan, dtype=np.float64)
    for r, c in zip(*np.where(n >= args.min_observations)):
        good = np.isfinite(layers[:, r, c])
        slope[r, c] = theilslopes(layers[good, r, c], t[good])[0]
    out = Path(args.output_dir)
    _write(out / "lst_theil_sen_c_per_year.tif", slope, profile,
           {"METHOD": "Theil-Sen", "UNITS": "degC/year",
            "CAUTION": "Acquisition season and time-of-day differences affect trend interpretation"})
    with np.errstate(all="ignore"):
        mean = np.nanmean(layers, axis=0)
    mean[n == 0] = np.nan
    _write(out / "lst_temporal_mean_c.tif", mean, profile,
           {"METHOD": "Temporal mean", "UNITS": "degC"})
    _write(out / "lst_observation_count.tif", n.astype(float), profile,
           {"METHOD": "Valid observation count", "UNITS": "count"})
    print(f"Saved temporal layers to {out}")
    return 0


def main(argv=None):
    p = argparse.ArgumentParser(description="Optional LST research workflows")
    sub = p.add_subparsers(dest="command", required=True)
    d = sub.add_parser("downscale", help="Spatially evaluated coarse-to-fine ML LST prediction")
    d.add_argument("--lst", type=Path, required=True, help="Coarse-resolution reference LST")
    d.add_argument("--coarse", nargs="+", type=Path, required=True, help="Predictors on coarse grid")
    d.add_argument("--fine", nargs="+", type=Path, required=True, help="Matching fine-scale predictors")
    d.add_argument("--output", type=Path, required=True)
    d.add_argument("--model", choices=["ridge", "rf", "boost", "stack"], default="rf")
    d.add_argument("--block-size", type=int, default=20)
    d.add_argument("--test-size", type=float, default=0.25)
    d.add_argument("--seed", type=int, default=42)
    d.add_argument("--residual-correction", action="store_true")
    d.add_argument("--model-output", type=Path)
    d.set_defaults(func=downscale)
    t = sub.add_parser("trend", help="Temporal statistics and Theil-Sen slope")
    t.add_argument("--rasters", nargs="+", type=Path, required=True)
    t.add_argument("--dates", nargs="+", required=True, help="ISO dates, one per raster")
    t.add_argument("--min-observations", type=int, default=3)
    t.add_argument("--output-dir", type=Path, required=True)
    t.set_defaults(func=trend)
    args = p.parse_args(argv)
    if getattr(args, "block_size", 1) < 1:
        p.error("--block-size must be positive")
    if getattr(args, "test_size", 0.25) <= 0 or getattr(args, "test_size", 0.25) >= 1:
        p.error("--test-size must be strictly between 0 and 1")
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
