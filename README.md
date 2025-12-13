# Automated-Land-Surface-Temperature-LST-Mapping
Universal Satellite LST Calculator: A unified Python pipeline for retrieving Land Surface Temperature from Landsat 8/9, ASTER, MODIS, and ECOSTRESS. Implements physics-based SCA retrieval for raw sensors and automated scaling for L2 products
# Universal Satellite Land Surface Temperature (LST) Calculator

## Overview
This repository hosts a **sensor-agnostic Python pipeline** designed to standardize the retrieval of Land Surface Temperature (LST) across multi-resolution satellite missions.

Unlike single-sensor scripts, this tool employs a **Registry Pattern** architecture to handle diverse data maturity levels:
1.  **Physics-Based Retrieval (Level 1):** Implements the **Simplified Single-Channel Algorithm (SCA)** for raw sensors (Landsat 8/9, ASTER), performing dynamic radiometric calibration and Planck’s Law inversion.
2.  **Product Scaling (Level 2):** Automates metadata parsing, unit conversion, and valid-range masking for operational products (MODIS, ECOSTRESS).

## Key Features
* **Multi-Sensor Support:** Single interface for processing **Landsat 8/9 (TIRS)**, **ASTER (TIR)**, **MODIS**, and **ECOSTRESS**.
* **Scientific Rigor:**
    * Uses sensor-specific center wavelengths ($\lambda$) and calibration constants ($K_1, K_2$).
    * Enforces physical guardrails (e.g., Emissivity sanity checks: $0.9 < \epsilon < 1.0$).
    * Standardized physical constants ($\rho = 14380 \mu m \cdot K$) based on Sobrino et al.
* **Automated Metadata Detection:** Auto-detects sensor type from GeoTIFF tags to select the appropriate processing algorithm.
* **Defensible Output:** Tags output GeoTIFFs with processing metadata (`METHOD`, `SENSOR`, `UNITS`), ensuring reproducibility in research workflows.

## Methodology

### 1. Physics-Based Retrieval (Landsat/ASTER)
For raw Digital Numbers (DN), the pipeline reverses the sensor transfer function:
1.  **Radiance ($L_\lambda$):** $DN \times Gain + Offset$ (or Unit Conversion Coeffs for ASTER).
2.  **Brightness Temperature ($BT$):** Inversion of Planck's Law.
3.  **Emissivity Correction:**
    $$LST = \frac{BT}{1 + (\lambda \times BT / \rho) \times \ln(\epsilon)}$$
    *(Note: This implementation assumes atmospheric transmittance $\tau \approx 1$. For hyperspectral accuracy, users should integrate MODTRAN parameters.)*

### 2. Product-Based Retrieval (MODIS/ECOSTRESS)
Parses internal metadata to apply scale factors and offsets automatically:
$$LST_{Celsius} = (DN \times Scale + Offset) - 273.15$$
*Includes automated masking for valid data ranges (e.g., MODIS 7500-65535).*

## Usage
Run the main driver script and follow the interactive prompts:
```bash
python universal_lst.py
