"""
Universal Satellite LST Calculator (v3.1)
=========================================

Land Surface Temperature (LST, in degrees Celsius) from Landsat 8/9, ASTER,
MODIS and ECOSTRESS GeoTIFFs.

Processing paths
  RAW      (Landsat 8/9 Level-1 band 10, ASTER L1T)
           DN -> radiance -> brightness temperature, then one of
             --atmo none  emissivity correction only (default; no atmosphere)
             --atmo sc    single-channel algorithm, Jimenez-Munoz et al. (2014)
             --atmo sw    split-window algorithm (bands 10+11), same paper
           sc and sw need total column water vapour (g/cm2) from YOU
           (--water-vapor or --water-vapor-raster); this tool does not fetch it.
  PRODUCT  (Landsat Collection 2 Level-2 ST_B10, MODIS MOD11/MYD11, ECOSTRESS)
           DN * scale + offset -> Kelvin -> Celsius, with valid-range masking.

Extras
  --qa            Landsat QA_PIXEL cloud / shadow / cirrus masking
  --compare-st    compare against the USGS Level-2 ST_B10 product
  --insitu-csv    compare against your own station measurements

Constants policy: only values verified against primary documentation are built
in.  Landsat 9 and ASTER have NO built-in radiometric constants -- pass --mtl
(Landsat) or explicit values (ASTER).  MTL values always beat registry defaults.

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
from rasterio.enums import Resampling
from rasterio.warp import reproject
from rasterio.warp import transform as warp_transform

logger = logging.getLogger("lst")

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
RHO_UM_K = 14380.0           # hc/k_B (um K), used in the emissivity correction
C1 = 1.19104e8               # Planck c1, W um^4 m^-2 sr^-1
C2_UM_K = 1.43877e4          # Planck c2, um K  (K2 = c2 / wavelength)
LST_MIN_C, LST_MAX_C = -50.0, 90.0
NODATA_OUT = -9999.0
KELVIN = 273.15
WV_TABULATED_MAX = 4.5       # g/cm2: upper end of the SC atmospheric-function table

# Landsat Collection 2 QA_PIXEL bits (USGS): 0 fill, 1 dilated cloud, 2 cirrus,
# 3 cloud, 4 cloud shadow, 5 snow, 6 clear, 7 water
QA_CLOUD_BITS = (0, 1, 2, 3, 4)
QA_SNOW_BIT, QA_WATER_BIT = 5, 7

# NDVI-threshold emissivity coefficients, Skokovic et al. (2014) as tabulated in
# Remote Sensing 12(2):294 (Landsat 8 TIRS band 10 and 11).
EMIS_COEF = {
    10: {"a": 0.979, "b": 0.046, "ev": 0.987, "es": 0.971},
    11: {"a": 0.982, "b": 0.027, "ev": 0.989, "es": 0.977},
}

# Split-window coefficients, Jimenez-Munoz et al. (2014), Landsat 8 bands 10/11
SW_COEF = (-0.268, 1.378, 0.183, 54.30, -2.238, -129.20, 16.40)

SENSORS: dict[str, dict] = {
    # Landsat 8/9 Collection 2 Level-2 Surface Temperature (USGS Level-2 guide:
    # 0.00341802 * DN + 149.0 Kelvin, valid DN 293-61440).
    "LANDSAT_L2ST": {
        "type": "PRODUCT", "keys": ["ST_B10"],
        "scale": 0.00341802, "offset": 149.0, "dn_range": (293, 61440),
        "file_filter": "ST_B10",
    },
    # Landsat 8 band 10: constants as printed in real Collection 2 MTL files.
    # The effective wavelength is derived from K2 (lambda = c2 / K2).
    "LANDSAT_8": {
        "type": "RAW", "keys": ["LANDSAT_8", "LANDSAT08", "LC08", "LO08", "LT08"],
        "k1": 774.8853, "k2": 1321.0789, "rad_mult": 3.3420e-4, "rad_add": 0.1,
    },
    # Landsat 9: constants differ between scenes / processing dates -> use --mtl.
    "LANDSAT_9": {
        "type": "RAW", "keys": ["LANDSAT_9", "LANDSAT09", "LC09", "LO09", "LT09"],
    },
    "ASTER": {"type": "RAW", "keys": ["AST_L1T", "AST_L1A", "ASTER"]},
    # MODIS LST: uint16, scale 0.02 K, valid range 7500-65535, fill 0.
    "MODIS": {
        "type": "PRODUCT",
        "keys": ["MOD11A1", "MYD11A1", "MOD11A2", "MYD11A2", "MOD11B1", "MYD11B1",
                 "MOD11B2", "MYD11B2", "MOD11C1", "MYD11C1", "MOD11C2", "MYD11C2",
                 "MOD11_L2", "MYD11_L2"],
        "scale": 0.02, "offset": 0.0, "dn_range": (7500, 65535),
    },
    # ECOSTRESS: uint16 (scale 0.02) in some products, float32 Kelvin in others
    # (e.g. L2T GeoTIFF). Float input is treated as already physical.
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
# Grid alignment
# ---------------------------------------------------------------------------
def same_grid(s, ref) -> bool:
    return s.crs == ref.crs and s.shape == ref.shape and s.transform.almost_equals(ref.transform)


def to_grid_array(arr: np.ndarray, s, ref, resampling) -> np.ndarray:
    """Resample an array that belongs to dataset `s` onto the grid of `ref`."""
    if same_grid(s, ref):
        return arr
    if s.crs is None or ref.crs is None:
        raise ValueError("Cannot align rasters without a CRS on both.")
    dst = np.full(ref.shape, np.nan, dtype=np.float32)
    reproject(source=arr.astype(np.float32), destination=dst,
              src_transform=s.transform, src_crs=s.crs,
              dst_transform=ref.transform, dst_crs=ref.crs,
              resampling=resampling, src_nodata=np.nan, dst_nodata=np.nan)
    return dst


def read_on_grid(path: Path, ref, resampling=Resampling.bilinear,
                 zero_is_fill: bool = False) -> np.ndarray:
    """Read band 1 of `path` on the grid of `ref` (resampling only if grids differ)."""
    with rasterio.open(path) as s:
        arr = s.read(1).astype(np.float32)
        if s.nodata is not None and not np.isnan(s.nodata):
            arr[arr == s.nodata] = np.nan
        if same_grid(s, ref):
            return arr
        logger.info("Resampling %s onto the thermal grid (%s)", Path(path).name, resampling.name)
        if zero_is_fill:
            arr[arr == 0] = np.nan             # keep fill out of the interpolation
        return to_grid_array(arr, s, ref, resampling)


# ---------------------------------------------------------------------------
# Landsat MTL and QA
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


def qa_bit(qa: np.ndarray, bit: int) -> np.ndarray:
    """Boolean array: is `bit` set in QA_PIXEL? Missing values count as fill."""
    q = np.nan_to_num(qa, nan=1.0).astype(np.uint32)
    return ((q >> bit) & 1).astype(bool)


def qa_mask(qa: np.ndarray, mask_snow: bool = False, mask_water: bool = False) -> np.ndarray:
    bits = list(QA_CLOUD_BITS) + ([QA_SNOW_BIT] if mask_snow else []) \
        + ([QA_WATER_BIT] if mask_water else [])
    m = np.zeros(qa.shape, dtype=bool)
    for b in bits:
        m |= qa_bit(qa, b)
    return m


# ---------------------------------------------------------------------------
# Emissivity
# ---------------------------------------------------------------------------
def ndvi_emissivity(red: np.ndarray, nir: np.ndarray, band: int = 10):
    """
    NDVI-threshold emissivity from TOA reflectance (0-1), Landsat 8/9 band 10 or 11.
    Method: Sobrino et al. (2004) with Skokovic et al. (2014) coefficients (a, b, ev, es):

      NDVI < 0.2            eps = a - b * red
      0.2 <= NDVI <= 0.5    eps = ev*Pv + es*(1-Pv),   Pv = ((NDVI-0.2)/0.3)^2
      NDVI > 0.5            eps = ev

    No roughness (cavity) term. Slightly discontinuous at NDVI = 0.2 (~0.003).
    Returns (emissivity, ndvi); pixels with non-finite input stay NaN.
    """
    k = EMIS_COEF[band]
    red = red.astype(np.float32)
    nir = nir.astype(np.float32)
    with np.errstate(divide="ignore", invalid="ignore"):
        ndvi = (nir - red) / (nir + red)
    pv = np.clip(((ndvi - 0.2) / 0.3) ** 2, 0.0, 1.0)
    eps = np.full(ndvi.shape, np.nan, dtype=np.float32)
    soil = ndvi < 0.2
    eps[soil] = k["a"] - k["b"] * np.clip(red[soil], 0.0, 1.0)
    mixed = (ndvi >= 0.2) & (ndvi <= 0.5)
    eps[mixed] = k["ev"] * pv[mixed] + k["es"] * (1.0 - pv[mixed])
    eps[ndvi > 0.5] = k["ev"]
    return eps, ndvi


def get_emissivity(src, args, mtl: Optional[dict], red, nir, water: Optional[np.ndarray]):
    """Return (eps_band10, eps_band11); each is a scalar or an array on the thermal grid."""
    if args.emissivity_raster:
        eps10 = read_on_grid(args.emissivity_raster, src, Resampling.bilinear)
        eps11 = eps10
    elif red and nir:
        r = read_on_grid(red, src, Resampling.bilinear, zero_is_fill=not args.refl_input)
        n = read_on_grid(nir, src, Resampling.bilinear, zero_is_fill=not args.refl_input)
        if not args.refl_input:
            if mtl is None:
                raise ValueError("--red/--nir are raw DN: pass --mtl to convert them to reflectance, "
                                 "or --refl-input if they already hold reflectance (0-1).")
            r, n = toa_reflectance(r, mtl, 4), toa_reflectance(n, mtl, 5)
        eps10, _ = ndvi_emissivity(r, n, 10)
        eps11, _ = ndvi_emissivity(r, n, 11)
        logger.info("NDVI emissivity (band 10): mean %.4f, min %.4f, max %.4f",
                    np.nanmean(eps10), np.nanmin(eps10), np.nanmax(eps10))
    else:
        if args.emissivity is not None:
            eps = args.emissivity
        elif args.non_interactive:
            raise ValueError("No emissivity source: pass --emissivity, --emissivity-raster, or --red + --nir.")
        else:
            eps = ask("Surface emissivity (0-1): ", float)
        if not (0 < eps <= 1):
            raise ValueError(f"Invalid emissivity {eps}; must be in (0, 1].")
        eps10 = eps11 = eps
    for name, e in (("band 10", eps10), ("band 11", eps11)):
        v = np.asarray(e, dtype=np.float64)
        if np.any(np.isfinite(v) & ((v <= 0) | (v > 1))):
            raise ValueError(f"Invalid emissivity values for {name}; require 0 < e <= 1")
    if water is not None and args.water_emissivity is not None:
        we = args.water_emissivity
        if not (0 < we <= 1):
            raise ValueError(f"Invalid --water-emissivity {we}.")
        eps10 = np.where(water, we, np.broadcast_to(eps10, water.shape)).astype(np.float32)
        eps11 = np.where(water, we, np.broadcast_to(eps11, water.shape)).astype(np.float32)
    return eps10, eps11


# ---------------------------------------------------------------------------
# Physics
# ---------------------------------------------------------------------------
def dn_to_rad_bt(dn: np.ndarray, sensor: str, c: dict, nodata=None):
    """DN -> (radiance, brightness temperature in K). Fill pixels become NaN."""
    dn = dn.astype(np.float32)
    dn[dn == 0] = np.nan
    if nodata is not None and not np.isnan(nodata):
        dn[dn == nodata] = np.nan
    rad = (dn - 1.0) * c["ucc"] if sensor == "ASTER" else dn * c["rad_mult"] + c["rad_add"]
    with np.errstate(invalid="ignore", divide="ignore"):
        rad[rad <= 0] = np.nan
        bt = c["k2"] / np.log(c["k1"] / rad + 1.0)
    return rad, bt


def bt_to_lst_c(bt_k, eps, wavelength_um: float):
    """Emissivity correction of brightness temperature (Artis & Carnahan, 1982):
    LST = BT / (1 + (lambda * BT / rho) * ln(eps)).  No atmospheric correction."""
    if isinstance(eps, np.ndarray):
        eps = np.where((eps <= 0) | (eps > 1), np.nan, eps)
    with np.errstate(divide="ignore", invalid="ignore"):
        lst_k = bt_k / (1.0 + (wavelength_um * bt_k / RHO_UM_K) * np.log(eps))
    return lst_k - KELVIN


def single_channel_lst_c(rad, bt_k, eps, w, wavelength_um: float):
    """Single-channel algorithm, Jimenez-Munoz et al. (2014), Landsat 8 band 10.

    LST = gamma * [ (psi1 * L + psi2) / eps + psi3 ] + delta
    gamma = c1 * lam * T^2 / ( c2 * L * (lam^5 * L + c1) ),   delta = T - gamma * L
    psi1..3 are quadratic in column water vapour w (g/cm2).  L = at-sensor radiance
    (W m-2 sr-1 um-1), T = at-sensor brightness temperature (K).
    """
    L = np.asarray(rad, dtype=np.float64)
    T = np.asarray(bt_k, dtype=np.float64)
    w = np.asarray(w, dtype=np.float64)
    lam = wavelength_um
    with np.errstate(divide="ignore", invalid="ignore"):
        gamma = C1 * lam * T ** 2 / (C2_UM_K * L * (lam ** 5 * L + C1))
        delta = T - gamma * L
        psi1 = 0.04019 * w ** 2 + 0.02916 * w + 1.01523
        psi2 = -0.38333 * w ** 2 - 1.50294 * w + 0.20324
        psi3 = 0.00918 * w ** 2 + 1.36072 * w - 0.27514
        lst_k = gamma * ((psi1 * L + psi2) / eps + psi3) + delta
    return (lst_k - KELVIN).astype(np.float32)


def split_window_lst_c(bt10, bt11, eps10, eps11, w):
    """Split-window algorithm, Jimenez-Munoz et al. (2014), Landsat 8 bands 10 + 11.

    LST = T10 + c0 + c1 (T10-T11) + c2 (T10-T11)^2
          + (c3 + c4 w)(1 - eps_mean) + (c5 + c6 w) d_eps
    eps_mean = (eps10 + eps11)/2,  d_eps = eps10 - eps11.  T in K, w in g/cm2.
    """
    c0, c1, c2, c3, c4, c5, c6 = SW_COEF
    T10 = np.asarray(bt10, dtype=np.float64)
    T11 = np.asarray(bt11, dtype=np.float64)
    d = T10 - T11
    em = (np.asarray(eps10, dtype=np.float64) + np.asarray(eps11, dtype=np.float64)) / 2.0
    de = np.asarray(eps10, dtype=np.float64) - np.asarray(eps11, dtype=np.float64)
    w = np.asarray(w, dtype=np.float64)
    lst_k = T10 + c0 + c1 * d + c2 * d ** 2 + (c3 + c4 * w) * (1 - em) + (c5 + c6 * w) * de
    return (lst_k - KELVIN).astype(np.float32)


def get_water_vapor(src, args):
    """Column water vapour in g/cm2: scalar or array on the thermal grid."""
    if args.water_vapor_raster:
        w = read_on_grid(args.water_vapor_raster, src, Resampling.bilinear)
    elif args.water_vapor is not None:
        w = args.water_vapor
    elif args.non_interactive:
        raise ValueError("--atmo sc/sw needs water vapour: pass --water-vapor (g/cm2) or "
                         "--water-vapor-raster. (ERA5 total column water vapour is in kg/m2; "
                         "divide by 10 to get g/cm2.)")
    else:
        w = ask("Column water vapour in g/cm2 (ERA5 kg/m2 / 10): ", float)
    wmax = float(np.nanmax(w))
    if not np.all(np.isfinite(w)):
        raise ValueError("Water vapour contains missing or nonfinite values.")
    if float(np.nanmin(w)) < 0:
        raise ValueError("Water vapour must be >= 0 g/cm2.")
    if wmax > WV_TABULATED_MAX:
        if not args.allow_wv_extrapolation:
            raise ValueError(f"Water vapour {wmax:.2f} g/cm2 exceeds the documented SC range "
                             f"({WV_TABULATED_MAX:.1f}); use --allow-wv-extrapolation "
                             "only with independent validation.")
        logger.warning("EXTRAPOLATING atmospheric functions beyond %.1f g/cm2", WV_TABULATED_MAX)
    return w


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
    if c.get("wavelength") is None and c.get("k2") and sensor.startswith("LANDSAT"):
        c["wavelength"] = C2_UM_K / c["k2"]          # effective wavelength from K2
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
        dn[(dn < lo) | (dn > hi)] = np.nan

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
# Stats, validation, output
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


def _errors(pred: np.ndarray, ref: np.ndarray) -> dict:
    d = pred - ref
    out = {"n": int(d.size), "bias": float(d.mean()), "rmse": float(np.sqrt((d ** 2).mean())),
           "mae": float(np.abs(d).mean()), "r": float("nan")}
    if d.size >= 3 and pred.std() > 0 and ref.std() > 0:
        out["r"] = float(np.corrcoef(pred, ref)[0, 1])
    return out


def compare_with_st(lst: np.ndarray, ref_src, st_path: Path) -> dict:
    """Compare against the USGS Level-2 ST_B10 product (a model product, not ground truth)."""
    with rasterio.open(st_path) as s:
        st = to_grid_array(product_lst(s, SENSORS["LANDSAT_L2ST"]), s, ref_src, Resampling.bilinear)
    both = np.isfinite(lst) & np.isfinite(st)
    if both.sum() < 2:
        logger.warning("Fewer than 2 co-valid pixels against %s", st_path.name)
        return {"st_n": int(both.sum())}
    e = _errors(lst[both], st[both])
    logger.info("vs L2 ST_B10: n=%d  bias=%.2f K  RMSE=%.2f K  r=%.3f", e["n"], e["bias"], e["rmse"], e["r"])
    return {"st_n": e["n"], "st_bias": e["bias"], "st_rmse": e["rmse"], "st_r": e["r"]}


def validate_insitu(lst: np.ndarray, ref_src, csv_path: Path, window: int, out_csv: Path) -> dict:
    """Compare with station data. CSV columns: lon, lat, temp_c (WGS84, degrees C)."""
    if ref_src.crs is None:
        raise ValueError("Raster has no CRS; cannot locate in-situ points.")
    with open(csv_path, newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise ValueError(f"{csv_path.name} has no rows.")
    cols = {k.strip().lower(): k for k in rows[0]}
    for need in ("lon", "lat", "temp_c"):
        if need not in cols:
            raise ValueError(f"{csv_path.name} needs columns lon, lat, temp_c (found {list(rows[0])}).")
    lons = [float(r[cols["lon"]]) for r in rows]
    lats = [float(r[cols["lat"]]) for r in rows]
    obs = [float(r[cols["temp_c"]]) for r in rows]
    xs, ys = warp_transform("EPSG:4326", ref_src.crs, lons, lats)
    if window < 1 or window % 2 != 1:
        raise ValueError("--insitu-window must be a positive odd integer.")
    half = window // 2
    h, w = lst.shape
    kept = []
    for lon, lat, o, x, y in zip(lons, lats, obs, xs, ys):
        r, c = ref_src.index(x, y)
        if not (0 <= r < h and 0 <= c < w):
            continue
        patch = lst[max(r - half, 0):r + half + 1, max(c - half, 0):c + half + 1]
        v = patch[np.isfinite(patch)]
        if v.size:
            kept.append((lon, lat, o, float(v.mean()), int(v.size)))
    logger.info("In-situ: %d of %d points inside the raster with valid pixels", len(kept), len(rows))
    if len(kept) < 1:
        return {"insitu_n": 0}
    with open(out_csv, "w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(["lon", "lat", "obs_c", "sat_c", "n_pixels"])
        wr.writerows([(a, b, o, round(s, 3), n) for a, b, o, s, n in kept])
    e = _errors(np.array([k[3] for k in kept]), np.array([k[2] for k in kept]))
    logger.info("vs in-situ: n=%d  bias=%.2f K  RMSE=%.2f K  MAE=%.2f K", e["n"], e["bias"], e["rmse"], e["mae"])
    return {"insitu_n": e["n"], "insitu_bias": e["bias"], "insitu_rmse": e["rmse"], "insitu_mae": e["mae"]}


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


CSV_FIELDS = ["scene", "sensor", "input_file", "atmo", "pct_masked", "n_valid", "mean", "std",
              "min", "max", "p05", "p50", "p95", "st_n", "st_bias", "st_rmse", "st_r",
              "insitu_n", "insitu_bias", "insitu_rmse", "insitu_mae"]


def append_csv(row: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and path.stat().st_size > 0:
        with open(path, newline="") as f:
            header = next(csv.reader(f), [])
        if header != CSV_FIELDS:
            alt = path.with_name(path.stem + "_v31" + path.suffix)
            logger.warning("%s has an older column layout; writing to %s instead.", path.name, alt.name)
            path = alt
    new = not path.exists() or path.stat().st_size == 0

    def fmt(v):
        if isinstance(v, float):
            return "" if np.isnan(v) else f"{v:.4f}"
        return v

    with open(path, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=CSV_FIELDS, restval="")
        if new:
            w.writeheader()
        w.writerow({k: fmt(v) for k, v in row.items()})


# ---------------------------------------------------------------------------
# One scene
# ---------------------------------------------------------------------------
METHOD_LABEL = {
    "none": "Emissivity-corrected brightness temperature (no atmospheric correction)",
    "sc": "Single-channel algorithm, Jimenez-Munoz et al. 2014 (water vapour user-supplied)",
    "sw": "Split-window algorithm, Jimenez-Munoz et al. 2014 (water vapour user-supplied)",
}


def process(thermal: Path, args, *, mtl_path: Optional[Path] = None,
            red: Optional[Path] = None, nir: Optional[Path] = None,
            b11: Optional[Path] = None, qa: Optional[Path] = None,
            sensor: Optional[str] = None, out_dir: Optional[Path] = None,
            csv_path: Optional[Path] = None) -> Optional[dict]:
    """Process one thermal GeoTIFF. Returns the stats row, or None on failure."""
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
            qa_arr = read_on_grid(qa, src, Resampling.nearest) if qa else None
            tags = {}
            atmo = "none"
            w_used = None

            if cfg["type"] == "RAW":
                wanted = args.water_vapor is not None or args.water_vapor_raster is not None
                atmo = args.atmo or ("sc" if wanted else "none")
                if atmo != "none" and sensor not in ("LANDSAT_8", "LANDSAT_9"):
                    raise ValueError("--atmo sc/sw is implemented for Landsat 8/9 band 10 only.")
                if atmo != "none" and args.band != 10:
                    raise ValueError("--atmo sc/sw works with band 10 only.")
                water = qa_bit(qa_arr, QA_WATER_BIT) if qa_arr is not None else None
                eps10, eps11 = get_emissivity(src, args, mtl, red, nir, water)
                c10 = resolve_constants(sensor, cfg, args, mtl)
                rad10, bt10 = dn_to_rad_bt(src.read(1), sensor, c10, src.nodata)
                if atmo == "none":
                    lst = bt_to_lst_c(bt10, eps10, c10["wavelength"])
                else:
                    w_used = get_water_vapor(src, args)
                    if atmo == "sc":
                        lst = single_channel_lst_c(rad10, bt10, eps10, w_used, c10["wavelength"])
                    else:
                        if mtl is None or b11 is None:
                            raise ValueError("--atmo sw needs --mtl and --b11 (band 11 GeoTIFF).")
                        c11 = landsat_thermal_constants(mtl, 11)
                        dn11 = read_on_grid(b11, src, Resampling.bilinear, zero_is_fill=True)
                        _, bt11 = dn_to_rad_bt(dn11, sensor, c11)
                        lst = split_window_lst_c(bt10, bt11, eps10, eps11, w_used)
                method = METHOD_LABEL[atmo]
            else:
                lst = product_lst(src, cfg, args.scale, args.offset)
                method = "Product scaling"

            pct_masked = float("nan")
            if qa_arr is not None:
                m = qa_mask(qa_arr, args.mask_snow, args.mask_water)
                lst = np.where(m, np.nan, lst).astype(np.float32)
                pct_masked = 100.0 * float(m.mean())
                logger.info("QA_PIXEL mask removed %.1f%% of pixels", pct_masked)

            lst = sanity_filter(lst, args.lst_min, args.lst_max)
            stats = compute_stats(lst)
            if stats["n_valid"] == 0:
                raise ValueError("No valid pixels left after masking. Check inputs and constants.")

            path = unique_path(out_dir or thermal.parent, sensor, thermal.stem)
            val: dict = {}
            if args.compare_st:
                val.update(compare_with_st(lst, src, args.compare_st))
            if args.insitu_csv:
                val.update(validate_insitu(lst, src, args.insitu_csv, args.insitu_window,
                                           path.with_name(path.stem + "_insitu.csv")))
            tags.update({"SENSOR": sensor, "LST_UNITS": "Celsius", "METHOD": method,
                         "ATMOSPHERIC_CORRECTION": atmo, "SOURCE_FILE": thermal.name,
                         "GENERATED_BY": "lst_universal_research_safe v3.1"})
            if w_used is not None:
                tags["WATER_VAPOR_G_CM2"] = f"{float(np.nanmean(w_used)):.3f}"
            write_geotiff(lst, src.profile, path, tags, args.cog)
            if args.plot:
                write_quicklook(lst, path.with_suffix(".png"), f"{sensor} LST", stats)
            row = {"scene": thermal.stem, "sensor": sensor, "input_file": thermal.name,
                   "atmo": atmo, "pct_masked": pct_masked, **stats, **val}
            if csv_path:
                append_csv(row, csv_path)
            print(f"\n{thermal.stem}  [{sensor}, atmo={atmo}]  valid={stats['n_valid']:,}  "
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
                     "nir": files.get(f"{prefix}_b5.tif".lower()),
                     "b11": files.get(f"{prefix}_b11.tif".lower()),
                     "qa": None if args.no_qa_mask else files.get(f"{prefix}_qa_pixel.tif".lower())})
    return jobs


def product_jobs(folder: Path, args) -> list[dict]:
    cfg = SENSORS[args.sensor]
    flt = cfg.get("file_filter")
    files = {p.name.lower(): p for p in folder.iterdir() if p.is_file()}
    picked = sorted({p for p in files.values()
                     if p.suffix.lower() == ".tif" and not p.name.upper().startswith("LST_")
                     and (flt is None or flt in p.name.upper())})
    jobs = []
    for p in picked:
        qa = None
        if args.sensor == "LANDSAT_L2ST" and not args.no_qa_mask:
            qa = files.get(p.name.upper().replace("_ST_B10", "_QA_PIXEL").lower())
        jobs.append({"thermal": p, "mtl": None, "sensor": args.sensor,
                     "red": None, "nir": None, "b11": None, "qa": qa})
    return jobs


def run_batch(args) -> int:
    folder = args.batch
    if not folder.is_dir():
        logger.error("Batch folder not found: %s", folder)
        return 2
    if args.sensor == "ASTER":
        logger.error("Batch mode does not support ASTER (needs per-band constants). Use --input.")
        return 2
    if args.compare_st:
        logger.warning("--compare-st is for single scenes; ignored in batch mode.")
        args.compare_st = None
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
                    b11=j["b11"], qa=j["qa"], sensor=j["sensor"], out_dir=out_dir,
                    csv_path=csv_path)
        ok += r is not None
    logger.info("Batch finished: %d/%d scenes succeeded. Output: %s", ok, len(jobs), out_dir)
    return 0 if ok == len(jobs) else 1


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
  # Landsat 8/9 Level-1: MTL constants, NDVI emissivity, cloud mask
  python lst_universal_research_safe.py -i LC08_..._B10.TIF --mtl LC08_..._MTL.txt \\
      --red LC08_..._B4.TIF --nir LC08_..._B5.TIF --qa LC08_..._QA_PIXEL.TIF --plot

  # Same with single-channel atmospheric correction (water vapour in g/cm2)
  python lst_universal_research_safe.py -i LC08_..._B10.TIF --mtl LC08_..._MTL.txt \\
      --red ... --nir ... --qa ... --atmo sc --water-vapor 3.2

  # Split-window (needs band 11)
  python lst_universal_research_safe.py -i ..._B10.TIF --b11 ..._B11.TIF --mtl ... \\
      --red ... --nir ... --atmo sw --water-vapor 3.2

  # Check against USGS Level-2 ST_B10, or against your own stations (lon,lat,temp_c)
  python lst_universal_research_safe.py ... --compare-st LC08_..._ST_B10.TIF
  python lst_universal_research_safe.py ... --insitu-csv stations.csv

  # Landsat Level-2 surface temperature / MODIS / ECOSTRESS
  python lst_universal_research_safe.py -i LC08_..._ST_B10.TIF --qa LC08_..._QA_PIXEL.TIF
  python lst_universal_research_safe.py -i MOD11A1_LST_Day_1km.tif --sensor MODIS

  # Batch: folder of Landsat scenes (needs *_MTL.txt); results in <folder>/lst_output
  python lst_universal_research_safe.py --batch ./scenes --atmo sc --water-vapor 3.2
""")
    p.add_argument("--input", "-i", type=Path, help="Thermal GeoTIFF")
    p.add_argument("--insitu-variable", choices=["surface", "air"], default="surface",
                   help="Station temperature type; air observations are rejected for direct validation")
    p.add_argument("--batch", "-b", type=Path, help="Folder of scenes (see examples)")
    p.add_argument("--sensor", choices=list(SENSORS), help="Force sensor (skip auto-detect)")
    p.add_argument("--mtl", type=Path, help="Landsat *_MTL.txt (single-file mode)")
    p.add_argument("--band", type=int, default=10, help="Landsat thermal band (default 10)")
    g = p.add_argument_group("emissivity (RAW sensors only)")
    g.add_argument("--emissivity", type=float, help="Scalar emissivity, 0 < e <= 1")
    g.add_argument("--emissivity-raster", type=Path, help="Per-pixel emissivity GeoTIFF")
    g.add_argument("--red", type=Path, help="Red band GeoTIFF (Landsat band 4) for NDVI emissivity")
    g.add_argument("--nir", type=Path, help="NIR band GeoTIFF (Landsat band 5) for NDVI emissivity")
    g.add_argument("--refl-input", action="store_true",
                   help="--red/--nir already contain reflectance (0-1), not raw DN")
    g.add_argument("--water-emissivity", type=float,
                   help="Emissivity to use on QA water pixels (needs --qa; no built-in default)")
    a = p.add_argument_group("atmospheric correction (Landsat 8/9 band 10)")
    a.add_argument("--atmo", choices=["none", "sc", "sw"],
                   help="none | sc single-channel | sw split-window (default: sc if water "
                        "vapour is given, else none)")
    a.add_argument("--water-vapor", type=float, help="Column water vapour, g/cm2")
    a.add_argument("--water-vapor-raster", type=Path, help="Water vapour raster, g/cm2")
    a.add_argument("--allow-wv-extrapolation", action="store_true",
                   help="Explicitly allow extrapolation beyond atmospheric-function calibration range")
    a.add_argument("--b11", type=Path, help="Band 11 GeoTIFF (for --atmo sw)")
    q = p.add_argument_group("masking")
    q.add_argument("--qa", type=Path, help="Landsat *_QA_PIXEL.TIF: mask cloud, shadow, cirrus, fill")
    q.add_argument("--mask-snow", action="store_true", help="Also mask snow (QA bit 5)")
    q.add_argument("--mask-water", action="store_true", help="Also mask water (QA bit 7)")
    q.add_argument("--no-qa-mask", action="store_true", help="Batch: ignore QA_PIXEL files")
    v = p.add_argument_group("validation")
    v.add_argument("--compare-st", type=Path, help="USGS Level-2 ST_B10 GeoTIFF to compare against")
    v.add_argument("--insitu-csv", type=Path, help="CSV with lon, lat, temp_c columns")
    v.add_argument("--insitu-window", type=int, default=3,
                   help="Pixel window (odd) averaged around each station (default 3)")
    c = p.add_argument_group("constant overrides (beat MTL and registry)")
    c.add_argument("--rad-mult", type=float)
    c.add_argument("--rad-add", type=float)
    c.add_argument("--k1", type=float)
    c.add_argument("--k2", type=float)
    c.add_argument("--wavelength", type=float, help="Effective wavelength (um); default c2/K2")
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
    actual_argv = list(sys.argv[1:] if argv is None else argv)
    if actual_argv and actual_argv[0] == "research":
        from lst_research_extensions import main as research_main
        return research_main(actual_argv[1:])
    parser = build_parser()
    args = parser.parse_args(argv)
    setup_logging(args.verbose, args.quiet)
    if args.insitu_csv and args.insitu_variable != "surface":
        parser.error("Air temperature is not a direct LST validation reference; use radiometric surface temperature.")

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
                  b11=args.b11, qa=args.qa, out_dir=args.output_dir, csv_path=args.csv)
    return 0 if row else 1


if __name__ == "__main__":
    sys.exit(main())
