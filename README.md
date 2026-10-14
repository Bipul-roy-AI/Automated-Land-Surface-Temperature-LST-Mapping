# Universal Satellite Land Surface Temperature (LST) Calculator

## What this tool does

Turns a thermal GeoTIFF into a Land Surface Temperature map in °C (GeoTIFF, with optional quick-look PNG, statistics CSV and Cloud-Optimized GeoTIFF).

| Input | What happens |
|---|---|
| **Landsat 8/9 Level-1**, band 10 (raw DN) | DN → radiance → brightness temperature → emissivity correction |
| **Landsat 8/9 Collection 2 Level-2** `ST_B10` | scale + offset → °C (this is the USGS surface-temperature product) |
| **MODIS** MOD11/MYD11 and **ECOSTRESS** LST | scale + offset, fill and valid-range masking → °C |
| **ASTER** L1T thermal bands | same as Landsat Level-1, but you must supply the constants |

**Emissivity** (Level-1 inputs only): one fixed value (`--emissivity 0.98`), your own raster (`--emissivity-raster`), or computed from NDVI using Landsat red and NIR bands (`--red`, `--nir`, with `--mtl`).

**Constants:** Landsat 8 has built-in constants. Landsat 9 does not; always pass the scene's `*_MTL.txt` with `--mtl`. Constants from the MTL file are used in preference to built-in defaults.

**Batch:** `--batch <folder>` processes every Landsat scene it finds (via `*_MTL.txt`), or every MODIS / ECOSTRESS / Level-2 ST file with `--sensor`. Results go to `<folder>/lst_output`, with one CSV row per scene.

```
pip install -r requirements.txt
python lst_universal_research_safe.py -i LC08_..._B10.TIF --mtl LC08_..._MTL.txt \
    --red LC08_..._B4.TIF --nir LC08_..._B5.TIF --plot
python lst_universal_research_safe.py --help
pytest -q
```

## Limitations

- **No atmospheric correction for Level-1 inputs.** Water vapour and path radiance are ignored, so results can be several kelvin off, most of all in humid air. For research use, prefer the Level-2 `ST_B10` product (`LANDSAT_L2ST`), which includes an atmospheric correction.
- **No cloud, shadow or water masking.** Only fill pixels (DN = 0) are removed. Mask clouds yourself, for example with the Landsat `QA_PIXEL` band.
- **NDVI emissivity is approximate.** It is an NDVI-threshold method with Landsat 8 band-10 coefficients (also used for Landsat 9). It has no roughness term, no separate treatment of water or built-up areas, and it jumps slightly at NDVI = 0.2 (about 0.2 K).
- **Landsat thermal pixels are 100 m** (resampled to 30 m). Red and NIR must be on the same grid as the thermal band.
- **Built-in constants are verified for Landsat 8 band 10 only.** The thermal effective wavelength is 10.904 µm; other sources quote 10.8–10.9 µm, which changes LST by less than 0.05 K. USGS has recalibrated Landsat 9 thermal data, so always use the MTL file that came with your scene.
- **ASTER** needs `--ucc`, `--k1`, `--k2` and `--wavelength` for the band you use, and is not supported in batch mode.
- **Output values** outside -50 to 90 °C are removed (change with `--lst-min` / `--lst-max`). No-data in the output is -9999.
- **Not validated against in-situ measurements.** Tests check the arithmetic, not the physical accuracy of the method.
