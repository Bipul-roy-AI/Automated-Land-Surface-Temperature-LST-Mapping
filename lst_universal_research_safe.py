"""
Universal Satellite LST Calculator (v3.0)
=========================================

Land Surface Temperature (LST, in degrees Celsius) from Landsat 8/9, ASTER,
MODIS and ECOSTRESS GeoTIFFs.

Two processing paths
  RAW      (Landsat 8/9 Level-1 band 10, ASTER L1T):
           DN -> radiance -> brightness temperature (Planck inversion)
           -> emissivity correction.  NO atmospheric correction.
  PRODUCT  (Landsat Collection 2 Level-2 ST_B10, MODIS MOD11/MYD11, ECOSTRESS):
           DN * scale + offset -> Kelvin -> Celsius, with valid-range masking.

Constants policy: only values verified against primary documentation are built
in.  Landsat 9 and ASTER have NO built-in radiometric constants -- pass --mtl
(Landsat) or the explicit values (ASTER).  Scene-specific MTL values always
beat registry defaults.

Author: Bipul Roy (original).  License: MIT
"""

from __future__ import annotations

import argparse
import csv
import logging
import re
import sys
from pathlib import Path
from typing import Optional

import numpy as np
import rasterio

logger = logging.getLogger("lst")

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
RHO_UM_K = 14380.0          # hc/k_B in micrometre-kelvin (second radiation constant)
LST_MIN_C, LST_MAX_C = -50.0, 90.0   # sanity window for the output (deg C)
NODATA_OUT = -9999.0
KELVIN = 273.15

# Registry. "keys" are matched (upper-case, alphanumeric word boundaries)
# against file name and GeoTIFF tags. Order matters: more specific first.
SENSORS: dict[str, dict] = {
    # Landsat 8/9 Collection 2 Level-2 Surface Temperature (scale/offset from
    # the USGS Level-2 Science Product Guide: 0.00341802 * DN + 149.0, Kelvin).
    "LANDSAT_L2ST": {
        "type": "PRODUCT", "keys": ["ST_B10"],
        "scale": 0.00341802, "offset": 149.0, "dn_range": (293, 61440),
        "file_filter": "ST_B10",
    },
    # Landsat 8 band 10: constants as printed in real Collection 2 MTL files.
    "LANDSAT_8": {
        "type": "RAW", "keys": ["LANDSAT_8", "LANDSAT08", "LC08", "LO08", "LT08"],
        "wavelength": 10.904,
        "k1": 774.8853, "k2": 1321.0789, "rad_mult": 3.3420e-4, "rad_add": 0.1,
    },
    # Landsat 9: constants differ between scenes/processing dates -> use --mtl.
    "LANDSAT_9": {
        "type": "RAW", "keys": ["LANDSAT_9", "LANDSAT09", "LC09", "LO09", "LT09"],
        "wavelength": 10.9,
    },
    # ASTER: no defaults (band-specific); user supplies ucc, k1, k2, wavelength.
    "ASTER": {"type": "RAW", "keys": ["AST_L1T", "AST_L1A", "ASTER"]},
    # MODIS LST: uint16, scale 0.02 K, valid range 7500-65535, fill 0.
    "MODIS": {
        "type": "PRODUCT",
        "keys": ["MOD11A1", "MYD11A1", "MOD11A2", "MYD11A2", "MOD11B1", "MYD11B1",
                 "MOD11B2", "MYD11B2", "MOD11C1", "MYD11C1", "MOD11C2", "MYD11C2",
                 "MOD11_L2", "MYD11_L2"],
        "scale": 0.02, "offset": 0.0, "dn_range": (7500, 65535),
    },
    # ECOSTRESS LST: uint16 (scale 0.02) in some products, float32 Kelvin in
    # others (e.g. L2T GeoTIFF). Float input is treated as already physical.
    "ECOSTRESS": {
        "type": "PRODUCT",
        "keys": ["ECOSTRESS", "ECO_L2T_LSTE", "ECO_L2G_LSTE", "ECO_L2_LSTE", "L2T_LST"],
        "scale": 0.02, "offset": 0.0, "dn_range": (7500, 65535),
    },
}


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------
def setup_logging(verbose: bool = False, quiet: bool = False) -> None:
    level = logging.WARNING if quiet else logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format="%(asctime)s [%(levelname)s] %(message)s",
                        datefmt="%H:%M:%S")


def ask(prompt: str, cast=str, tries: int = 5):
    """Prompt the user; give up cleanly on closed stdin or repeated bad input."""
    for _ in range(tries):
        try:
            return cast(input(prompt).strip())
        except EOFError:
            raise SystemExit("No input available (stdin is closed). Pass the values as options.")
        except ValueError:
            print("Invalid input, please try again.")
    raise SystemExit("Too many invalid inputs.")


def detect_sensor(src, path: Path) -> Optional[str]:
    """Detect sensor from file name, dataset tags and band-1 tags."""
    parts = [path.name.upper()]
    for tags in (src.tags(), src.tags(1)):
        parts.append(" ".join(str(v) for v in tags.values()).upper())
    text = " | ".join(parts)
    for name, cfg in SENSORS.items():
        for key in cfg["keys"]:
            if re.search(rf"(?<![A-Z0-9]){re.escape(key)}(?![A-Z0-9])", text):
                return name
    return None


# ---------------------------------------------------------------------------
# Landsat MTL
# ---------------------------------------------------------------------------
def parse_mtl(path: Path) -> dict:
    if not path.exists():
        raise FileNotFoundError(f"MTL file not found: {path}")
    params: dict[str, str] = {}
    for line in path.read_text(errors="replace").splitlines():
        m = re.match(r"^\s*([A-Z_0-9]+)\s*=\s*(.+?)\s*$", line)
        if m:
            params[m.group(1)] = m.group(2).strip('"')
    return params


def landsat_thermal_constants(mtl: dict, band: int = 10) -> dict:
    out = {}
    for label, key in (("rad_mult", f"RADIANCE_MULT_BAND_{band}"),
                       ("rad_add", f"RADIANCE_ADD_BAND_{band}"),
                       ("k1", f"K1_CONSTANT_BAND_{band}"),
                       ("k2", f"K2_CONSTANT_BAND_{band}")):
        if key not in mtl:
            raise KeyError(f"{key} missing in MTL file")
        out[label] = float(mtl[key])
    return out


def toa_reflectance(dn: np.ndarray, mtl: dict, band: int) -> np.ndarray:
    """Landsat DN -> TOA reflectance corrected for sun elevation (USGS formula)."""
    mult = float(mtl[f"REFLECTANCE_MULT_BAND_{band}"])
    add = float(mtl[f"REFLECTANCE_ADD_BAND_{band}"])
    sun = np.radians(float(mtl["SUN_ELEVATION"]))
    dn = dn.astype(np.float32)
    refl = (mult * dn + add) / np.sin(sun)
    refl[dn == 0] = np.nan                     # Landsat fill value
    return refl.astype(np.float32)


# ---------------------------------------------------------------------------
# Emissivity
# ---------------------------------------------------------------------------
def ndvi_emissivity(red: np.ndarray, nir: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    NDVI-threshold emissivity for Landsat 8/9 band 10 (reflectance inputs, 0-1).
    Coefficients: Skokovic et al. (2014) as tabulated in Remote Sensing 12(2):294,
    built on the Sobrino et al. (2004) NDVI-threshold method (NDVImin=0.2, NDVImax=0.5):

      NDVI < 0.2            eps = 0.979 - 0.046 * red
      0.2 <= NDVI <= 0.5    eps = 0.987*Pv + 0.971*(1-Pv),   Pv = ((NDVI-0.2)/0.3)^2
      NDVI > 0.5            eps = 0.987

    The cavity/roughness term (d_eps) is NOT included. The method is slightly
    discontinuous at NDVI = 0.2 (about 0.003 in eps, roughly 0.2 K).
    Returns (emissivity, ndvi). Pixels with non-finite input stay NaN.
    """
    red = red.astype(np.float32)
    nir = nir.astype(np.float32)
    with np.errstate(divide="ignore", invalid="ignore"):
        ndvi = (nir - red) / (nir + red)
    pv = np.clip(((ndvi - 0.2) / 0.3) ** 2, 0.0, 1.0)
    eps = np.full(ndvi.shape, np.nan, dtype=np.float32)
    soil = ndvi < 0.2
    eps[soil] = 0.979 - 0.046 * np.clip(red[soil], 0.0, 1.0)
    mixed = (ndvi >= 0.2) & (ndvi <= 0.5)
    eps[mixed] = 0.987 * pv[mixed] + 0.971 * (1.0 - pv[mixed])
    eps[ndvi > 0.5] = 0.987
    return eps, ndvi


def get_emissivity(shape, args, mtl: Optional[dict], red: Optional[Path], nir: Optional[Path]):
    """Return a scalar or per-pixel emissivity array."""
    if args.emissivity_raster:
        with rasterio.open(args.emissivity_raster) as er:
            eps = er.read(1).astype(np.float32)
        if eps.shape != shape:
            raise ValueError(f"Emissivity raster shape {eps.shape} != thermal shape {shape}; resample first.")
        return eps
    if red and nir:
        with rasterio.open(red) as r, rasterio.open(nir) as n:
            red_dn, nir_dn = r.read(1), n.read(1)
        if red_dn.shape != shape or nir_dn.shape != shape:
            raise ValueError("Red/NIR shape differs from thermal band shape; resample first.")
        if args.refl_input:
            red_r, nir_r = red_dn.astype(np.float32), nir_dn.astype(np.float32)
        else:
            if mtl is None:
                raise ValueError("--red/--nir are raw DN: pass --mtl to convert them to reflectance, "
                                 "or --refl-input if they already hold reflectance (0-1).")
            red_r = toa_reflectance(red_dn, mtl, 4)
            nir_r = toa_reflectance(nir_dn, mtl, 5)
        eps, _ = ndvi_emissivity(red_r, nir_r)
        logger.info("NDVI emissivity: mean %.4f (min %.4f, max %.4f)",
                    np.nanmean(eps), np.nanmin(eps), np.nanmax(eps))
        return eps
    if args.emissivity is not None:
        eps = args.emissivity
    elif args.non_interactive:
        raise ValueError("No emissivity source: pass --emissivity, --emissivity-raster, or --red + --nir.")
    else:
        eps = ask("Surface emissivity (0-1): ", float)
    if not (0 < eps <= 1):
        raise ValueError(f"Invalid emissivity {eps}; must be in (0, 1].")
    return eps


# ---------------------------------------------------------------------------
# Physics
# ---------------------------------------------------------------------------
def bt_to_lst_c(bt_k: np.ndarray, eps, wavelength_um: float) -> np.ndarray:
    """Emissivity correction of brightness temperature (Artis & Carnahan, 1982).

    LST = BT / (1 + (lambda * BT / rho) * ln(eps)),   rho = 14380 um K.
    No atmospheric correction (transmittance assumed 1, no path radiance).
    """
    if isinstance(eps, np.ndarray):
        eps = np.where((eps <= 0) | (eps > 1), np.nan, eps)
    with np.errstate(divide="ignore", invalid="ignore"):
        lst_k = bt_k / (1.0 + (wavelength_um * bt_k / RHO_UM_K) * np.log(eps))
    return lst_k - KELVIN


def resolve_constants(sensor: str, cfg: dict, args, mtl: Optional[dict]) -> dict:
    """Priority: command-line override > MTL file > registry default."""
    keys = ["rad_mult", "rad_add", "k1", "k2", "wavelength"] if sensor != "ASTER" \
        else ["ucc", "k1", "k2", "wavelength"]
    c = {k: cfg.get(k) for k in keys}
    used_defaults = mtl is None and sensor == "LANDSAT_8"
    if mtl is not None and sensor != "ASTER":
        c.update(landsat_thermal_constants(mtl, args.band))
    for k in keys:
        v = getattr(args, k, None)
        if v is not None:
            c[k] = v
            used_defaults = False
    missing = [k for k in keys if c.get(k) is None]
    if missing:
        if args.non_interactive:
            hint = "pass --mtl" if sensor.startswith("LANDSAT") else "pass them as options"
            raise ValueError(f"{sensor}: missing {', '.join(missing)} (no built-in default). {hint}.")
        for k in missing:
            c[k] = ask(f"{sensor}: enter {k}: ", float)
    if used_defaults:
        logger.warning("Using built-in Landsat 8 defaults; pass --mtl for scene-specific constants.")
    return c


def raw_lst(dn: np.ndarray, eps, sensor: str, c: dict, nodata=None) -> np.ndarray:
    dn = dn.astype(np.float32)
    dn[dn == 0] = np.nan                      # fill value
    if nodata is not None and not np.isnan(nodata):
        dn[dn == nodata] = np.nan
    rad = (dn - 1.0) * c["ucc"] if sensor == "ASTER" else dn * c["rad_mult"] + c["rad_add"]
    with np.errstate(invalid="ignore", divide="ignore"):
        rad[rad <= 0] = np.nan
        bt = c["k2"] / np.log(c["k1"] / rad + 1.0)
    return bt_to_lst_c(bt, eps, c["wavelength"])


def product_lst(src, cfg: dict, scale=None, offset=None) -> np.ndarray:
    """DN * scale + offset -> Kelvin -> Celsius, with fill and valid-range masking."""
    arr = src.read(1)
    is_float = np.issubdtype(arr.dtype, np.floating)
    dn = arr.astype(np.float32)
    tags = src.tags()
    if src.nodata is not None and not np.isnan(src.nodata):
        dn[dn == src.nodata] = np.nan
    dn[dn == 0] = np.nan                      # 0 is the fill value of all these products
    if not is_float and cfg.get("dn_range"):
        lo, hi = cfg["dn_range"]
        dn[(dn < lo) | (dn > hi)] = np.nan    # DN-domain valid range (integers only)

    if scale is None:
        if "scale_factor" in tags:
            scale = float(tags["scale_factor"])
        elif src.scales and src.scales[0] not in (None, 1.0):
            scale = float(src.scales[0])
        elif is_float:
            scale = 1.0                        # float rasters hold physical values
        else:
            scale = cfg.get("scale")
        if scale is None:
            raise ValueError("Cannot determine scale factor; pass --scale.")
    if offset is None:
        if "add_offset" in tags:
            offset = float(tags["add_offset"])
        elif src.offsets and src.offsets[0] not in (None, 0.0):
            offset = float(src.offsets[0])
        elif is_float:
            offset = 0.0
        else:
            offset = cfg.get("offset", 0.0)
    logger.info("Product scaling: scale=%s offset=%s (%s input)", scale, offset,
                "float" if is_float else "integer")
    val = dn * scale + offset
    units = str(tags.get("UNITS", tags.get("units", "Kelvin"))).upper()
    return val if units.startswith("C") else val - KELVIN


# ---------------------------------------------------------------------------
# Stats and output
# ---------------------------------------------------------------------------
def sanity_filter(lst: np.ndarray, lo: float, hi: float) -> np.ndarray:
    out = lst.copy()
    with np.errstate(invalid="ignore"):
        bad = (out < lo) | (out > hi)
    if bad.any():
        logger.warning("Sanity filter removed %d pixels outside [%.0f, %.0f] C", int(bad.sum()), lo, hi)
    out[bad] = np.nan
    return out


def compute_stats(lst: np.ndarray) -> dict:
    v = lst[np.isfinite(lst)]
    if v.size == 0:
        return {"n_valid": 0}
    return {"n_valid": int(v.size), "mean": float(v.mean()), "std": float(v.std()),
            "min": float(v.min()), "max": float(v.max()),
            "p05": float(np.percentile(v, 5)), "p50": float(np.percentile(v, 50)),
            "p95": float(np.percentile(v, 95))}


def unique_path(out_dir: Path, sensor: str, scene: str, suffix: str = ".tif") -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    p = out_dir / f"LST_{sensor}_{scene}{suffix}"
    i = 1
    while p.exists():
        p = out_dir / f"LST_{sensor}_{scene}_{i}{suffix}"
        i += 1
    return p


def write_geotiff(lst: np.ndarray, profile: dict, path: Path, tags: dict, cog: bool) -> None:
    prof = dict(profile)
    for k in ("blockxsize", "blockysize", "tiled", "compress", "predictor",
              "photometric", "interleave"):
        prof.pop(k, None)
    prof.update(driver="GTiff", dtype="float32", count=1, nodata=NODATA_OUT, compress="deflate")
    arr = np.where(np.isfinite(lst), lst, NODATA_OUT).astype("float32")
    tmp = path.with_name(path.stem + "_tmp.tif") if cog else path
    with rasterio.open(tmp, "w", **prof) as dst:
        dst.write(arr, 1)
        dst.update_tags(**tags)
    if cog:
        try:
            from rasterio.shutil import copy as rio_copy
            rio_copy(tmp, path, driver="COG", compress="deflate")
            tmp.unlink()
        except Exception as e:  # GDAL < 3.1 has no COG driver
            logger.warning("COG driver unavailable (%s); wrote a plain GeoTIFF.", e)
            tmp.replace(path)
    logger.info("Wrote %s", path)


def write_quicklook(lst: np.ndarray, path: Path, title: str, stats: dict) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        logger.warning("matplotlib not installed; skipping quick-look PNG.")
        return
    fig, (a0, a1) = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
    im = a0.imshow(lst, cmap="inferno", vmin=stats["p05"], vmax=stats["p95"])
    a0.set_title(title)
    a0.axis("off")
    fig.colorbar(im, ax=a0, label="LST (°C)", shrink=0.85)
    v = lst[np.isfinite(lst)]
    a1.hist(v, bins=60, color="#d9622b")
    a1.axvline(stats["mean"], ls="--", color="k", label=f"mean {stats['mean']:.1f} °C")
    a1.set_xlabel("LST (°C)")
    a1.set_ylabel("pixels")
    a1.legend()
    fig.savefig(path, dpi=130)
    plt.close(fig)


CSV_FIELDS = ["scene", "sensor", "input_file", "n_valid", "mean", "std", "min", "max",
              "p05", "p50", "p95"]


def append_csv(row: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    new = not path.exists() or path.stat().st_size == 0
    with open(path, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=CSV_FIELDS)
        if new:
            w.writeheader()
        w.writerow({k: (f"{v:.4f}" if isinstance(v, float) else v) for k, v in row.items()})


# ---------------------------------------------------------------------------
# One scene
# ---------------------------------------------------------------------------
def process(thermal: Path, args, *, mtl_path: Optional[Path] = None,
            red: Optional[Path] = None, nir: Optional[Path] = None,
            sensor: Optional[str] = None, out_dir: Optional[Path] = None,
            csv_path: Optional[Path] = None) -> Optional[dict]:
    """Process one thermal GeoTIFF. Returns the stats dict, or None on failure."""
    logger.info("Processing %s", thermal.name)
    try:
        with rasterio.open(thermal) as src:
            sensor = sensor or args.sensor or detect_sensor(src, thermal)
            if sensor is None:
                if args.non_interactive:
                    raise ValueError("Cannot detect sensor; pass --sensor.")
                names = list(SENSORS)
                for i, n in enumerate(names, 1):
                    print(f"  {i}. {n}")
                sensor = names[ask("Sensor number: ", int) - 1]
            cfg = SENSORS[sensor]
            logger.info("Sensor: %s", sensor)
            mtl = parse_mtl(mtl_path) if mtl_path else None

            if cfg["type"] == "RAW":
                eps = get_emissivity(src.shape, args, mtl, red, nir)
                consts = resolve_constants(sensor, cfg, args, mtl)
                lst = raw_lst(src.read(1), eps, sensor, consts, src.nodata)
                method = "Emissivity-corrected brightness temperature (no atmospheric correction)"
            else:
                lst = product_lst(src, cfg, args.scale, args.offset)
                method = "Product scaling"

            lst = sanity_filter(lst, args.lst_min, args.lst_max)
            stats = compute_stats(lst)
            if stats["n_valid"] == 0:
                raise ValueError("No valid pixels left after masking. Check inputs and constants.")

            scene = thermal.stem
            path = unique_path(out_dir or thermal.parent, sensor, scene)
            write_geotiff(lst, src.profile, path, {
                "SENSOR": sensor, "LST_UNITS": "Celsius", "METHOD": method,
                "SOURCE_FILE": thermal.name, "GENERATED_BY": "lst_universal_research_safe v3.0",
            }, args.cog)
            if args.plot:
                write_quicklook(lst, path.with_suffix(".png"), f"{sensor} LST", stats)
            row = {"scene": scene, "sensor": sensor, "input_file": thermal.name, **stats}
            if csv_path:
                append_csv(row, csv_path)
            print(f"\n{scene}  [{sensor}]  valid={stats['n_valid']:,}  "
                  f"mean={stats['mean']:.2f} C  std={stats['std']:.2f}  "
                  f"min={stats['min']:.2f}  max={stats['max']:.2f}")
            return row
    except Exception as e:
        logger.error("Failed on %s: %s", thermal.name, e)
        logger.debug("Traceback", exc_info=True)
        return None


# ---------------------------------------------------------------------------
# Batch mode
# ---------------------------------------------------------------------------
def landsat_jobs(folder: Path, args) -> list[dict]:
    """One job per Landsat scene, found through its *_MTL.txt (searched recursively)."""
    jobs = []
    mtls = {p for p in folder.rglob("*") if p.is_file() and p.name.upper().endswith("_MTL.TXT")}
    for mtl_path in sorted(mtls):
        prefix = mtl_path.name[:-len("_MTL.txt")]
        files = {p.name.lower(): p for p in mtl_path.parent.iterdir() if p.is_file()}
        craft = parse_mtl(mtl_path).get("SPACECRAFT_ID", "").upper()
        if craft not in ("LANDSAT_8", "LANDSAT_9"):
            logger.warning("Skipping %s: SPACECRAFT_ID '%s' is not Landsat 8/9", mtl_path.name, craft)
            continue
        if args.sensor and args.sensor != craft:
            continue
        thermal = files.get(f"{prefix}_b{args.band}.tif".lower())
        if thermal is None:
            logger.warning("Skipping %s: no band %d GeoTIFF next to it", mtl_path.name, args.band)
            continue
        jobs.append({"thermal": thermal, "mtl": mtl_path, "sensor": craft,
                     "red": files.get(f"{prefix}_b4.tif".lower()),
                     "nir": files.get(f"{prefix}_b5.tif".lower())})
    return jobs


def product_jobs(folder: Path, args) -> list[dict]:
    cfg = SENSORS[args.sensor]
    flt = cfg.get("file_filter")
    files = sorted({p for p in folder.iterdir()
                    if p.is_file() and p.suffix.lower() == ".tif"
                    and not p.name.upper().startswith("LST_")
                    and (flt is None or flt in p.name.upper())})
    return [{"thermal": p, "mtl": None, "sensor": args.sensor, "red": None, "nir": None}
            for p in files]


def run_batch(args) -> int:
    folder = args.batch
    if not folder.is_dir():
        logger.error("Batch folder not found: %s", folder)
        return 2
    if args.sensor == "ASTER":
        logger.error("Batch mode does not support ASTER (needs per-band constants). Use --input.")
        return 2
    if args.sensor in (None, "LANDSAT_8", "LANDSAT_9"):
        jobs = landsat_jobs(folder, args)
        if not jobs:
            logger.error("No Landsat scenes (*_MTL.txt with band %d) found. For MODIS, ECOSTRESS "
                         "or Landsat L2 ST folders pass --sensor.", args.band)
            return 2
    else:
        jobs = product_jobs(folder, args)
        if not jobs:
            logger.error("No matching .tif files in %s", folder)
            return 2
    out_dir = args.output_dir or folder / "lst_output"
    csv_path = args.csv or out_dir / "lst_stats.csv"
    ok = 0
    for j in jobs:
        r = process(j["thermal"], args, mtl_path=j["mtl"], red=j["red"], nir=j["nir"],
                    sensor=j["sensor"], out_dir=out_dir, csv_path=csv_path)
        ok += r is not None
    logger.info("Batch finished: %d/%d scenes succeeded. Output: %s", ok, len(jobs), out_dir)
    return 0 if ok else 1


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="lst_universal_research_safe.py",
        description="Land Surface Temperature (deg C) from Landsat 8/9, ASTER, MODIS, ECOSTRESS.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
Examples:
  # Landsat 8/9 Level-1: constants from MTL, NDVI emissivity from bands 4 and 5
  python lst_universal_research_safe.py -i LC08_..._B10.TIF --mtl LC08_..._MTL.txt \\
      --red LC08_..._B4.TIF --nir LC08_..._B5.TIF --plot

  # Same, with a fixed emissivity
  python lst_universal_research_safe.py -i LC08_..._B10.TIF --mtl LC08_..._MTL.txt --emissivity 0.98

  # Landsat Collection 2 Level-2 surface temperature
  python lst_universal_research_safe.py -i LC08_..._ST_B10.TIF

  # MODIS or ECOSTRESS
  python lst_universal_research_safe.py -i MOD11A1_LST_Day_1km.tif --sensor MODIS

  # Batch: a folder of Landsat scene folders (needs *_MTL.txt); results in <folder>/lst_output
  python lst_universal_research_safe.py --batch ./scenes --non-interactive
""")
    p.add_argument("--input", "-i", type=Path, help="Thermal GeoTIFF")
    p.add_argument("--batch", "-b", type=Path, help="Folder of scenes (see examples)")
    p.add_argument("--sensor", choices=list(SENSORS), help="Force sensor (skip auto-detect)")
    p.add_argument("--mtl", type=Path, help="Landsat *_MTL.txt (single-file mode)")
    p.add_argument("--band", type=int, default=10, help="Landsat thermal band (default 10)")
    g = p.add_argument_group("emissivity (RAW sensors only)")
    g.add_argument("--emissivity", type=float, help="Scalar emissivity, 0 < e <= 1")
    g.add_argument("--emissivity-raster", type=Path, help="Per-pixel emissivity GeoTIFF (same grid)")
    g.add_argument("--red", type=Path, help="Red band GeoTIFF (Landsat band 4) for NDVI emissivity")
    g.add_argument("--nir", type=Path, help="NIR band GeoTIFF (Landsat band 5) for NDVI emissivity")
    g.add_argument("--refl-input", action="store_true",
                   help="--red/--nir already contain reflectance (0-1), not raw DN")
    c = p.add_argument_group("constant overrides (beat MTL and registry)")
    c.add_argument("--rad-mult", type=float)
    c.add_argument("--rad-add", type=float)
    c.add_argument("--k1", type=float)
    c.add_argument("--k2", type=float)
    c.add_argument("--wavelength", type=float, help="Effective wavelength (um)")
    c.add_argument("--ucc", type=float, help="ASTER unit conversion coefficient")
    c.add_argument("--scale", type=float, help="Product scale factor override")
    c.add_argument("--offset", type=float, help="Product offset override")
    o = p.add_argument_group("output")
    o.add_argument("--output-dir", "-o", type=Path, help="Default: next to the input "
                   "(batch: <folder>/lst_output)")
    o.add_argument("--plot", action="store_true", help="Also write a quick-look PNG")
    o.add_argument("--csv", type=Path, help="Append per-scene statistics to this CSV")
    o.add_argument("--cog", action="store_true", help="Write a Cloud-Optimized GeoTIFF")
    o.add_argument("--lst-min", type=float, default=LST_MIN_C)
    o.add_argument("--lst-max", type=float, default=LST_MAX_C)
    p.add_argument("--non-interactive", action="store_true", help="Never prompt; fail if input is missing")
    p.add_argument("--verbose", "-v", action="store_true")
    p.add_argument("--quiet", "-q", action="store_true")
    return p


def main(argv: Optional[list[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    setup_logging(args.verbose, args.quiet)

    if args.batch is not None:
        args.non_interactive = True          # never prompt in batch mode
        return run_batch(args)

    if args.input is None:
        if args.non_interactive:
            parser.error("--input or --batch is required with --non-interactive")
        print("Interactive mode (use --help for all options).")
        args.input = Path(ask("Thermal GeoTIFF path: ").strip("\"'"))
        mtl = ask("Landsat MTL file (blank to skip): ").strip("\"'")
        args.mtl = Path(mtl) if mtl else None
        if args.mtl:
            red = ask("Red band file for NDVI emissivity (blank to skip): ").strip("\"'")
            nir = ask("NIR band file (blank to skip): ").strip("\"'") if red else ""
            args.red, args.nir = (Path(red), Path(nir)) if red and nir else (None, None)

    if not args.input.exists():
        logger.error("Input file not found: %s", args.input)
        return 2
    row = process(args.input, args, mtl_path=args.mtl, red=args.red, nir=args.nir,
                  out_dir=args.output_dir, csv_path=args.csv)
    return 0 if row else 1


if __name__ == "__main__":
    sys.exit(main())
