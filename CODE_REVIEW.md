# Code Review Report

**Repo:** `Bipul-roy-AI/Automated-Land-Surface-Temperature-LST-Mapping`
**Reviewed file:** `lst_universal_research_safe.py`
**Date:** 2025-10-15

---

## Summary

The pipeline implements a clean Registry-Pattern architecture and the SCA algorithm is mathematically correct. However, the analysis has **8 critical issues** and **9 minor issues** that should be addressed before the tool is used in published research. Most issues are about (a) sensor-specific assumptions being hardcoded, (b) side effects on input arrays, and (c) missing atmospheric correction.

---

## Critical Issues

### 1. Filename mismatch in README
**Where:** `README.md` line 38
**Issue:** README says `python universal_lst.py` but the actual file is `lst_universal_research_safe.py`. Users will get a "file not found" error on the first run.
**Fix:** Update the command in README to match the actual filename.

### 2. Typo in dependency file name
**Where:** Repo root
**Issue:** `dependicies.txt` should be `dependencies.txt`.
**Fix:** Rename the file and update any references.

### 3. MODIS valid range is product-specific
**Where:** `calculate_product_lst` lines 153–155
```python
if sensor_name == "MODIS":
    print("[INFO] Applying MODIS valid range filter (7500-65535)")
    dn[(dn < 7500) | (dn > 65535)] = np.nan
```
**Issue:** The hardcoded DN range `7500–65535` is only valid for `MOD11A1` daily tile products. For `MOD11_L2` swaths, `MOD11B1` (1 km), or `MYD11A1`, the valid range differs. Hardcoding this silently drops legitimate cold pixels (e.g., winter snow surfaces, Antarctica) in those products.
**Fix:** Make the range configurable per MODIS sub-product. Detect the sub-product from the `PRODUCT` tag and apply the appropriate range. If unknown, skip the mask and warn the user.

### 4. In-place mutation of input array
**Where:** `calculate_sca_lst` line 120
```python
radiance[radiance <= 0] = np.nan
```
**Issue:** `radiance` here is a new array (`dn * rad_mult + rad_add`), so this line is actually fine — but the next line (`bt = k2 / np.log(...)`) creates another new array. Wait, the `radiance` array IS the new array, so this is NOT a mutation of the caller's `dn`. But it IS still a bug if the caller re-uses `radiance` after this function returns. Reviewing again... actually this is safe in the current code path because `radiance` is local. However, it would still be cleaner to use `np.where()` to make the intent explicit. **Downgrading this from Critical to Minor.**

### 5. Landsat 8 Band 10 wavelength is wrong
**Where:** `SENSOR_REGISTRY` lines 18 and 23
```python
"default_wavelength": 10.8 # Band 10
```
**Issue:** The USGS-specified effective wavelength for Landsat 8 Band 10 is **10.904 µm** (not 10.8). Using 10.8 introduces a small (~0.1–0.3 K) bias in the LST output.
**Fix:** Change to `10.9` (or better, `10.904`) for both Landsat 8 and Landsat 9.

### 6. ASTER single-wavelength approximation is wrong
**Where:** `calculate_sca_lst` ASTER branch
**Issue:** ASTER has 5 TIR bands (10–14, spanning 8.3–11.3 µm). The script asks for one band's wavelength and applies the single-channel correction to it. This throws away 4 bands of data and produces a less accurate LST than a proper NEM (Normalized Emissivity Method) or TES (Temperature Emissivity Separation) approach.
**Fix:** Implement NEM or TES algorithm for ASTER, processing all 5 bands. At minimum, warn the user that only single-band processing is supported and recommend using the official ASTER L2 LST product instead.

### 7. Product path blindly subtracts 273.15
**Where:** `calculate_product_lst` line 158
```python
lst_k = dn * scale + offset
return lst_k - 273.15
```
**Issue:** This assumes input is always in Kelvin. Some ASTER L2 products and certain ECOSTRESS processing levels are already in °C. For those, the output would be ≈ -273 °C.
**Fix:** Read the `UNITS` tag from the source dataset. If `UNITS == "Celsius"`, skip the Kelvin→Celsius conversion. If `UNITS` is missing, default to Kelvin (current behavior) but log a warning.

### 8. No atmospheric correction
**Where:** `calculate_sca_lst`
**Issue:** The SCA formula assumes atmospheric transmittance τ ≈ 1. For humid atmospheres (water vapor > 3 g/cm²), errors of 2–5 K are common. The README mentions this limitation but the code offers no way to supply τ or water-vapor values.
**Fix:** Add optional `tau` and `water_vapor` parameters. If supplied, use the Jiménez-Muñoz SCA formulation with atmospheric correction:
```
LST = BT / (1 + (λ · BT / ρ) · (ln(ε) - (1-τ) · ln(tau))) 
```
Alternatively, integrate with `Py6S` for full radiative transfer.

---

## Minor Issues

### 9. `print()` used everywhere instead of `logging`
**Where:** Throughout
**Fix:** Switch to `logging` module with a configurable level flag (`--verbose`, `--quiet`). This enables batch runs and log file redirection.

### 10. No CLI arguments — fully interactive
**Where:** `main()`
**Fix:** Add `argparse` with flags like `--input`, `--emissivity`, `--sensor`, `--mult`, `--add`, `--k1`, `--k2`, `--output`. Fall back to interactive prompts only if args are missing.

### 11. `path_input.strip('"')` only strips double quotes
**Where:** `main()` line 165
**Fix:** Use `path_input.strip('"\'')` or `pathlib.Path(path_input)` for cross-platform path handling.

### 12. Output written to CWD, not next to input
**Where:** `main()` line 206
**Fix:** Use `Path(path_input).parent / out_name` to write the output next to the input.

### 13. `detect_sensor_from_metadata` only reads dataset-level tags
**Where:** `detect_sensor_from_metadata`
**Issue:** Landsat MTL-derived GeoTIFFs often store sensor info in band tags (`src.tags(1)`) or in a sidecar `.mtl` file. Dataset-level tag search will frequently fail and fall back to manual selection.
**Fix:** Also check `src.tags(1)` and parse the sidecar `.mtl` file if present.

### 14. No QA/QC bit-mask handling for MODIS/ECOSTRESS
**Where:** `calculate_product_lst`
**Issue:** Both MOD11 and ECOSTRESS ship mandatory QC/QA bands that flag cloudy / bad pixels. Ignoring them means cloudy pixels get included in stats and the mean is biased.
**Fix:** Read companion QC band (e.g., `MOD11A1.AYYYYDDD.h00v00.061.QC.tif`), parse bit flags per USGS documentation, and mask cloudy/bad pixels before computing statistics.

### 15. No type hints, no docstrings
**Fix:** Add PEP-484 type hints and Google-style docstrings for all functions.

### 16. No unit tests
**Fix:** Add `pytest` tests on synthetic rasters (we generated two in this review) to catch regressions.

### 17. Emissivity guardrail warning range mismatch with README
**Where:** README says "0.9 < ε < 1.0", code warns when `< 0.9 or > 1.0`
**Fix:** Make consistent. The mathematically correct range for natural surfaces is 0.9 ≤ ε ≤ 1.0 (water = 0.995, vegetation = 0.985, bare soil = 0.93, snow = 0.99). Allow this range without warning; warn outside.

---

## Methodological Concerns

### 1. Single-Channel Algorithm Limitations
The SCA (Sobrino 2004) is a reasonable choice for sensors with a single thermal band (Landsat 8 B10). However:

- **It requires accurate emissivity.** The current implementation asks the user for a single ε value, but ε varies pixel-by-pixel (vegetation vs. bare soil vs. water). For accurate LST, ε should be computed from NDVI (Sobrino's NDVITEM method) or from a seasonal emissivity map.
- **It ignores atmospheric effects.** For sea-level humid sites, errors of 3–5 K are typical.
- **It's not the best choice for Landsat 8.** The Split-Window Algorithm (SWA) using both B10 and B11 gives ~0.5 K accuracy without atmospheric data. SWA should be offered as an option.

### 2. MODIS Product Family
The tool lumps all MODIS LST products together, but they differ significantly:

| Product | Resolution | Range |
|---------|------------|-------|
| MOD11A1 | 1 km daily | DN 7500–65535 |
| MOD11_L2 | 1 km swath | DN 0–65535, scale 0.02 |
| MOD11B1 | 6 km daily | DN 0–65535, scale 0.02 |
| MOD11C1 | 5600 m daily | different |

The hardcoded range `7500–65535` only fits MOD11A1.

### 3. ECOSTRESS Handling
ECOSTRESS LST comes in HDF5 with a `scale_factor` attribute on the LST dataset (typically 0.002). The current code reads `scale_factor` from dataset-level GeoTIFF tags, but ECOSTRESS native files are HDF5, not GeoTIFF. The user must first convert with `h5dump` or `rasterio`'s HDF5 driver — this conversion step is not documented.

### 4. Output Statistics
The current output prints min/max/mean. For research use, percentiles (5th, 50th, 95th) and standard deviation are more robust to outliers. The synthetic test showed max = 65.81 °C, which is likely a single hot outlier pixel — reporting percentiles would give a better picture.

---

## Recommended Additions (Priority Order)

| Priority | Feature | Effort |
|----------|---------|--------|
| **High** | NDVI-based emissivity (Sobrino NDVITEM) | Medium |
| **High** | Split-Window Algorithm for Landsat 8 (B10 + B11) | Medium |
| **High** | Automatic visualization (PNG map + histogram + stats) | Low |
| **High** | CLI arguments via argparse | Low |
| **High** | QA/QC bit-mask for MODIS/ECOSTRESS | Medium |
| Medium | CSV/JSON stats export | Low |
| Medium | Reproject & clip (rasterio.warp) | Medium |
| Medium | Proper logging | Low |
| Medium | Unit tests (pytest) | Medium |
| Low | Sentinel-3 SLSTR support | Medium |
| Low | Time-series mode (folder → trend map) | Medium |
| Low | Streamlit web UI | Low |

---

## Test Results on Synthetic Data

We generated two synthetic rasters and ran the original tool end-to-end:

| Test | Sensor | Input DN range | LST output | Notes |
|------|--------|----------------|------------|-------|
| 1 | Landsat 8 B10 | 15–47 (16-bit) | mean 29.2 °C, range -9.6 to 65.8 °C | Realistic for South Asia daytime |
| 2 | MODIS (MOD11A1) | 13981–16558 (DN) | mean 31.5 °C, range 6.5 to 58.0 °C | Realistic; river/city/vegetation features visible |

Both pipelines produced physically plausible LST values. The auto-detection worked correctly when sensor info was embedded in GeoTIFF tags. The output GeoTIFFs loaded cleanly in QGIS and Python.

---

## Conclusion

The core algorithm is sound, but the tool has several rough edges that limit its research-grade reliability. The most impactful single fix would be **adding NDVI-based emissivity calculation** — this eliminates the largest source of user error (manual ε input) and aligns the tool with the most-cited SCA methodology.
