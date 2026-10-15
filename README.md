# Universal Satellite Land Surface Temperature (LST) Calculator

Turns a thermal GeoTIFF into a Land Surface Temperature map in °C. Supports Landsat 8/9, ASTER, MODIS and ECOSTRESS, with optional atmospheric correction, cloud masking, and validation against the USGS Level-2 product or your own weather stations.

```
pip install -r requirements.txt
python lst_universal_research_safe.py -i LC08_..._B10.TIF --mtl LC08_..._MTL.txt \
    --red LC08_..._B4.TIF --nir LC08_..._B5.TIF --qa LC08_..._QA_PIXEL.TIF \
    --atmo sc --water-vapor 3.0 --plot
```

![LST quick-look](screenshots/screenshot_lst_sc.png)

---

## What this tool does

| Input | What happens |
|---|---|
| **Landsat 8/9 Level-1**, band 10 (raw DN) | DN → radiance → brightness temperature, then one of:<br>• `--atmo none` — emissivity correction only (default)<br>• `--atmo sc` — single-channel algorithm (Jiménez-Muñoz et al. 2014), needs water vapour<br>• `--atmo sw` — split-window algorithm (same paper), needs band 11 + water vapour |
| **Landsat 8/9 Collection 2 Level-2** `ST_B10` | scale + offset → °C (USGS surface-temperature product, already atmospherically corrected) |
| **MODIS** MOD11/MYD11 and **ECOSTRESS** LST | scale + offset, fill and valid-range masking → °C |
| **ASTER** L1T thermal bands | same as Landsat Level-1, but you must supply the constants |

**Emissivity** (Level-1 inputs only): one fixed value (`--emissivity 0.98`), your own raster (`--emissivity-raster`), or computed from NDVI using Landsat red and NIR bands (`--red`, `--nir`, with `--mtl`).

**Cloud masking**: pass `--qa path/to/QA_PIXEL.TIF` and the tool masks fill, dilated cloud, cirrus, cloud, and cloud shadow (USGS Collection 2 bit layout). Add `--mask-snow` or `--mask-water` to also mask those classes.

**Validation**: `--compare-st ST_B10.tif` compares your result against the USGS Level-2 surface-temperature product (bias, RMSE, r). `--insitu-csv stations.csv` compares against ground station data (CSV columns `lon, lat, temp_c`).

**Batch**: `--batch <folder>` processes every Landsat scene it finds (via `*_MTL.txt`), or every MODIS / ECOSTRESS / Level-2 ST file with `--sensor`. Results go to `<folder>/lst_output`, with one CSV row per scene.

---

## Files you need

### Landsat 8 / 9 Level-1 (most common case)

Download from [USGS EarthExplorer](https://earthexplorer.usgs.gov/) — select *Landsat Collection 2 Level-1* product. One scene unzips to ~12 files; you need these five:

| File | Purpose | Required? |
|------|---------|-----------|
| `LC08_..._B10.TIF` | Thermal band 10 (the input the tool computes LST from) | **Yes** |
| `LC08_..._MTL.txt` | Metadata file with radiometric constants (RADIANCE_MULT, K1, K2, etc.) | **Yes** for Landsat 9; strongly recommended for Landsat 8 |
| `LC08_..._B4.TIF` | Red band, for NDVI-based emissivity | For NDVI ε |
| `LC08_..._B5.TIF` | NIR band, for NDVI-based emissivity | For NDVI ε |
| `LC08_..._QA_PIXEL.TIF` | Quality band, for cloud / shadow / fill masking | For cloud mask |
| `LC08_..._B11.TIF` | Thermal band 11, only for `--atmo sw` (split-window) | For SW only |

If you skip `--red`/`--nir`, pass `--emissivity 0.98` (or your best estimate). If you skip `--qa`, fill pixels (DN=0) are still masked but clouds are not.

### Landsat 8/9 Collection 2 Level-2 ST_B10 (already atmospherically corrected)

| File | Purpose |
|------|---------|
| `LC08_..._ST_B10.TIF` | USGS surface-temperature product (scale 0.00341802, offset 149.0 K) |
| `LC08_..._QA_PIXEL.TIF` | Optional, for cloud masking |

Run: `python lst_universal_research_safe.py -i LC08_..._ST_B10.TIF` — the sensor is auto-detected from the filename.

### MODIS (MOD11A1 / MYD11A1 / MOD11_L2 etc.)

| File | Purpose |
|------|---------|
| `MOD11A1.AYYYYDDD.hVVvVV.061.LST_Day_1km.tif` | LST day or night band (uint16, scale 0.02 K) |

Convert the HDF to GeoTIFF first with `gdal_translate -sds MOD11A1.hdf LST_Day_1km.tif`. The tool auto-detects the sensor from the filename and reads `scale_factor` from the GeoTIFF tags.

### ECOSTRESS

| File | Purpose |
|------|---------|
| `ECOv002_L2T_LST_*.tif` | LST GeoTIFF (float32 Kelvin, or uint16 with scale 0.02) |

The tool detects the input dtype: float rasters are treated as already-physical (Kelvin), integer rasters are scaled.

### ASTER L1T

| File | Purpose |
|------|---------|
| `AST_L1T_..._Band_10.tif` (or 11–14) | Thermal band GeoTIFF |

ASTER has no built-in constants — you must pass `--ucc`, `--k1`, `--k2`, `--wavelength` for the specific band (these are in the ASTER HDF header).

---

## How to use it

### Quick start — Landsat 8 with everything

```bash
python lst_universal_research_safe.py \
    -i LC08_L1TP_137045_20240115_B10.TIF \
    --mtl LC08_L1TP_137045_20240115_MTL.txt \
    --red LC08_L1TP_137045_20240115_B4.TIF \
    --nir LC08_L1TP_137045_20240115_B5.TIF \
    --qa  LC08_L1TP_137045_20240115_QA_PIXEL.TIF \
    --atmo sc --water-vapor 3.0 \
    --plot --csv stats.csv
```

This runs the single-channel atmospheric correction with 3.0 g/cm² column water vapour, computes per-pixel emissivity from NDVI, masks clouds via QA_PIXEL, writes a quick-look PNG and appends a stats row to `stats.csv`.

### Terminal output

![Run 1 — Basic](screenshots/screenshot_run1_basic.png)

![Run 2 — Single-channel atmospheric correction](screenshots/screenshot_run2_sc.png)

![Run 3 — MODIS product](screenshots/screenshot_run3_modis.png)

### Output LST quick-look PNG (`--plot`)

The tool writes a two-panel PNG next to each output GeoTIFF: an LST map (inferno colourmap, 5th–95th percentile stretch) and a histogram with the mean marked.

![LST map — single-channel run](screenshots/screenshot_lst_sc.png)

![LST map — basic run](screenshots/screenshot_lst_basic.png)

![LST map — MODIS](screenshots/screenshot_lst_modis.png)

### All command-line options

```
python lst_universal_research_safe.py --help
```

| Flag | Purpose |
|------|---------|
| `-i, --input` | Thermal GeoTIFF |
| `-b, --batch` | Folder of scenes (Landsat: needs `*_MTL.txt`) |
| `--sensor` | Force sensor (`LANDSAT_8`, `LANDSAT_9`, `LANDSAT_L2ST`, `ASTER`, `MODIS`, `ECOSTRESS`) |
| `--mtl` | Landsat `*_MTL.txt` (reads RADIANCE_MULT, K1, K2, REFLECTANCE_* automatically) |
| `--band` | Landsat thermal band number (default 10) |
| **Emissivity** | |
| `--emissivity` | Scalar emissivity (0 < ε ≤ 1) |
| `--emissivity-raster` | Per-pixel emissivity GeoTIFF |
| `--red` / `--nir` | Red (B4) and NIR (B5) GeoTIFFs for NDVI emissivity |
| `--refl-input` | `--red`/`--nir` already hold reflectance (0–1), not raw DN |
| `--water-emissivity` | Emissivity for water pixels (needs `--qa`) |
| **Atmospheric correction** (Landsat 8/9 band 10) | |
| `--atmo` | `none` \| `sc` (single-channel) \| `sw` (split-window). Default: `sc` if water vapour given, else `none` |
| `--water-vapor` | Column water vapour, g/cm² (ERA5 kg/m² ÷ 10) |
| `--water-vapor-raster` | Water vapour raster, g/cm² |
| `--b11` | Band 11 GeoTIFF (for `--atmo sw`) |
| **Masking** | |
| `--qa` | Landsat `*_QA_PIXEL.TIF`: mask cloud, shadow, cirrus, fill |
| `--mask-snow` | Also mask snow (QA bit 5) |
| `--mask-water` | Also mask water (QA bit 7) |
| `--no-qa-mask` | Batch mode: ignore QA_PIXEL files |
| **Validation** | |
| `--compare-st` | USGS Level-2 ST_B10 GeoTIFF to compare against |
| `--insitu-csv` | CSV with `lon, lat, temp_c` columns |
| `--insitu-window` | Pixel window (odd) averaged around each station (default 3) |
| **Constant overrides** (beat MTL and registry) | |
| `--rad-mult` / `--rad-add` | RADIANCE_MULT / ADD override |
| `--k1` / `--k2` | Thermal constants override |
| `--wavelength` | Effective wavelength (µm); default c2/K2 |
| `--ucc` | ASTER unit conversion coefficient |
| `--scale` / `--offset` | Product scale / offset override |
| **Output** | |
| `-o, --output-dir` | Default: next to input (batch: `<folder>/lst_output`) |
| `--plot` | Also write a quick-look PNG |
| `--csv` | Append per-scene statistics to this CSV |
| `--cog` | Write a Cloud-Optimized GeoTIFF |
| `--lst-min` / `--lst-max` | Sanity filter range (default -50 / 90 °C) |
| `--non-interactive` | Never prompt; fail if input is missing |
| `-v, --verbose` / `-q, --quiet` | Logging level |

---

## Where to get water vapour for `--atmo sc` / `--atmo sw`

The single-channel and split-window algorithms need total column water vapour in g/cm². Free sources:

- **ERA5** (Copernicus, free with registration): download `Total column water vapour` (variable `tcwv`) in kg/m² — **divide by 10** to get g/cm².
- **MODIS MOD05_L2** (water vapour product): scale 0.001 g/cm².
- **NCEP/NCAR reanalysis**: less precise but global daily coverage.
- **AERONET** stations: sun-photometer measurements, the gold standard but sparse.

Example: ERA5 says 30 kg/m² → `--water-vapor 3.0`.

---

## Output

The tool writes one GeoTIFF per scene: `LST_<sensor>_<scene>_<n>.tif`, with:

- **dtype:** `float32`, **nodata:** `-9999`
- **units:** degrees Celsius
- **CRS & transform:** inherited from the input
- **tags:** `SENSOR`, `LST_UNITS`, `METHOD`, `ATMOSPHERIC_CORRECTION`, `SOURCE_FILE`, `GENERATED_BY`, and `WATER_VAPOR_G_CM2` (when applicable)

The METHOD tag reflects what was actually run:

| `--atmo` | METHOD tag |
|----------|------------|
| `none` | `Emissivity-corrected brightness temperature (no atmospheric correction)` |
| `sc` | `Single-channel algorithm, Jimenez-Muñoz et al. 2014 (water vapour user-supplied)` |
| `sw` | `Split-window algorithm, Jimenez-Muñoz et al. 2014 (water vapour user-supplied)` |
| (products) | `Product scaling` |

### Reading the output in Python

```python
import rasterio, numpy as np
with rasterio.open("LST_LANDSAT_8_LC08_..._B10_1.tif") as src:
    lst = src.read(1).astype(np.float32)
    lst[lst == -9999] = np.nan
    print(f"Sensor: {src.tags()['SENSOR']}")
    print(f"Method: {src.tags()['METHOD']}")
    print(f"Mean LST: {np.nanmean(lst):.2f} °C")
```

### CSV columns (`--csv`)

`scene, sensor, input_file, atmo, pct_masked, n_valid, mean, std, min, max, p05, p50, p95, st_n, st_bias, st_rmse, st_r, insitu_n, insitu_bias, insitu_rmse, insitu_mae`

---

## Limitations

- **Atmospheric correction is optional and approximate.** With `--atmo none` (the default), water vapour and path radiance are ignored, so results can be several kelvin off in humid air. With `--atmo sc` / `--atmo sw`, the Jiménez-Muñoz et al. (2014) atmospheric functions are used — these are tabulated only up to 4.5 g/cm² water vapour, so the tool extrapolates beyond that and warns you. For best accuracy, prefer the USGS Level-2 `ST_B10` product (`LANDSAT_L2ST`), which includes a full atmospheric correction.
- **Cloud masking requires the QA_PIXEL band.** Without `--qa`, only fill pixels (DN = 0) are removed — clouds and shadows stay in the data and bias the mean. Pass `--qa LC08_..._QA_PIXEL.TIF` to mask them.
- **NDVI emissivity is approximate.** It is an NDVI-threshold method (Sobrino et al. 2004) with Landsat 8 band-10 coefficients from Skokovic et al. (2014) (also used for Landsat 9). It has no roughness (cavity) term, no separate treatment of water or built-up areas, and it is slightly discontinuous at NDVI = 0.2 (about 0.2 K). For water bodies, use `--qa` + `--water-emissivity 0.995`.
- **Landsat thermal pixels are 100 m** natively (resampled to 30 m in the product). Red, NIR, QA and band 11 are auto-resampled onto the thermal grid if their grids differ.
- **Built-in constants are verified for Landsat 8 band 10 only.** The effective wavelength is derived from K2 via `λ = c2 / K2` (≈10.9 µm). Landsat 9 and ASTER have no built-in radiometric constants — always pass `--mtl` (Landsat 9) or the explicit values (ASTER).
- **ASTER** needs `--ucc`, `--k1`, `--k2` and `--wavelength` for the band you use, and is not supported in batch mode.
- **Split-window uses Landsat 8 band 11**, which USGS warns has known stray-light errors. Prefer `--atmo sc` unless you have a specific reason to use SW.
- **Output values** outside -50 to 90 °C are removed (change with `--lst-min` / `--lst-max`). No-data in the output is -9999.
- **Validation against in-situ data** (`--insitu-csv`) compares a satellite pixel window to a point station measurement — this conflates spatial representativeness errors with retrieval errors. Use bias/RMSE values as indicative, not definitive.

---

## Testing

```bash
pip install -r requirements.txt
pytest tests/ -v
```

The 15 tests use realistic Landsat 8 DN values (~30,000, as in real Collection 2 scenes) and real MTL constants. They verify: fill-pixel masking, NDVI emissivity from reflectance, the single-channel arithmetic, MODIS / ECOSTRESS / Level-2 ST scaling, batch mode, closed-stdin handling, and COG output.

---

## References

- **Sobrino, J. A., Jiménez-Muñoz, J. C., & Paolini, L. (2004).** *Land surface temperature retrieval from LANDSAT TM 5.* Remote Sensing of Environment, 90(4), 434–440.
- **Jiménez-Muñoz, J. C., Sobrino, J. A., Skoković, D., Mattar, C., & Cristóbal, J. (2014).** *Land surface temperature retrieval methods from Landsat-8 thermal infrared sensor data.* IEEE Geoscience and Remote Sensing Letters, 11(10), 1840–1843.
- **Skoković, D., Sobrino, J. A., Jiménez-Muñoz, J. C., Sòria, G., Julien, Y., Mattar, C., & Cristóbal, J. (2014).** *Calibration and validation of land surface temperature for Landsat-8 TIRS.* Remote Sensing, 12(2):294. (Note: the correct citation is *Remote Sensing of Environment*, 2014 — see the code for the coefficient values.)
- **USGS (2023).** *Landsat 8-9 Collection 2 Level 1 / Level 2 Science Product Guide.*
- **MODIS Land Surface Temperature Team.** *MOD11A1 MODIS/Terra LST/Emissivity Daily L3 Global 1 km.* https://lpdaac.usgs.gov/products/mod11a1v061/

---

## License

MIT — see [LICENSE](LICENSE).


---

## Optional research workflows (development branch)

The default single-scene physical retrieval is unchanged. New research features are
explicitly opt-in, because statistical downscaling is not directly measured LST.

**Spatially evaluated LST downscaling** (inspired by the workflows in
pansism/downscale-satelliteLST and wqlin08/LST-Downscaling-Project):

```bash
pip install scikit-learn scipy joblib
python lst_universal_research_safe.py research downscale \
  --lst landsat_lst_c.tif \
  --coarse coarse_ndvi.tif coarse_ndbi.tif coarse_mndwi.tif \
  --fine sentinel_ndvi.tif sentinel_ndbi.tif sentinel_mndwi.tif \
  --model rf --block-size 20 --output modelled_fine_lst_c.tif
```

Supported models: `ridge`, `rf`, `boost` (histogram gradient boosting), and
`stack` (Ridge and Random Forest ensemble). `--residual-correction` is
available but off by default. Predictors must describe compatible surface variables,
use matching band order, units and acquisition conditions; inputs should be
preprocessed, spatially co-registered and quality masked. Coarse and fine predictor
rasters must use comparable physical definitions. A JSON report records spatial-block
holdout metrics, but **does not establish fine-resolution accuracy**. Use independent
fine-resolution thermal measurements for that claim. No cloud filling is inferred to
be valid and no resolution gain is claimed as a thermal sensor measurement.

**Multi-date temporal analysis** (inspired by leonsnill/lst_landsat):

```bash
python lst_universal_research_safe.py research trend \
  --rasters lst_2021.tif lst_2022.tif lst_2023.tif \
  --dates 2021-04-01 2022-04-01 2023-04-01 \
  --output-dir trends
```

Outputs: mean LST (Celsius), valid observation count, Theil-Sen slope
(Celsius/year). All inputs must have compatible units, QA screening,
acquisition season and overpass time; results are descriptive, not
a significance test. Trend slopes calculated from seasonal mixtures may
be misleading. The code uses the first input grid and aligns other scenes
using nearest-neighbour resampling.

**Stricter safeguards:** positive odd in-situ windows, rejection of
explicit air-temperature reference data for direct LST validation,
emissivity range checks, and nonzero return code when a batch only partially
succeeds. Water-vapour values above the supported single-channel
range now require explicit `--allow-wv-extrapolation`.

**Scientific scope:** satellite LST retrieval, statistical downscaling, and
temporal trend estimation are separate products. The new modules are
research prototypes pending independent quantitative benchmarking.

**Acknowledgements:** Workflow concepts were informed by
[pansism/downscale-satelliteLST](https://github.com/pansism/downscale-satelliteLST),
[leonsnill/lst_landsat](https://github.com/leonsnill/lst_landsat),
and [wqlin08/LST-Downscaling-Project](https://github.com/wqlin08/LST-Downscaling-Project).
The optional module is a new implementation rather than copied upstream code.
