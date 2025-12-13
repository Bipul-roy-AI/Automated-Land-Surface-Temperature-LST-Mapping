import numpy as np
import rasterio
import sys
import warnings

# ==============================
# SCIENTIFIC CONSTANTS
# ==============================
RHO_MICRON_K = 14380.0  # µm·K (Sobrino et al.)

# ==============================
# SENSOR REGISTRY
# ==============================
SENSOR_REGISTRY = {
    "LANDSAT_8": {
        "type": "RAW",
        "platform_keys": ["LANDSAT_8", "LC08"],
        "default_wavelength": 10.8 # Band 10
    },
    "LANDSAT_9": {
        "type": "RAW",
        "platform_keys": ["LANDSAT_9", "LC09"],
        "default_wavelength": 10.8 # Band 10
    },
    "ASTER": {
        "type": "RAW",
        "platform_keys": ["ASTER"],
        "default_wavelength": None # Force user input (Multiband)
    },
    "MODIS": {
        "type": "PRODUCT",
        "platform_keys": ["MODIS", "MOD11", "MYD11"]
    },
    "ECOSTRESS": {
        "type": "PRODUCT",
        "platform_keys": ["ECOSTRESS"]
    },
    # Removed Sentinel-3 to ensure GeoTIFF compatibility only
}

# ==============================
# HELPER FUNCTIONS
# ==============================
def safe_float(prompt):
    while True:
        try:
            return float(input(prompt))
        except ValueError:
            print("Invalid numeric input. Please enter a number (e.g., 0.98).")

def detect_sensor_from_metadata(src):
    try:
        tags = src.tags()
        # Combine all metadata values into one search string
        text = " ".join([str(v) for v in tags.values()]).upper()
        
        for sensor, cfg in SENSOR_REGISTRY.items():
            for key in cfg["platform_keys"]:
                if key in text:
                    return sensor
    except Exception:
        pass # Fail silently and fall back to manual
    return None

def ask_user_sensor():
    print("\n[INPUT REQUIRED] Could not auto-detect sensor from metadata.")
    print("Select Satellite Sensor:")
    sensors = list(SENSOR_REGISTRY.keys())
    for i, s in enumerate(sensors):
        print(f"  {i+1}. {s}")
    
    while True:
        try:
            idx = int(input("Enter number: ")) - 1
            if 0 <= idx < len(sensors):
                return sensors[idx]
        except ValueError:
            pass
        print("Invalid selection.")

# ==============================
# PHYSICS-BASED SCA (Refined)
# ==============================
def calculate_sca_lst(dn, emissivity, sensor_name, default_wavelength=None):
    print(f"\n--- Calibrating {sensor_name} ---")

    # 1. Radiance Conversion
    if sensor_name == "ASTER":
        print("ASTER requires Unit Conversion Coefficient (UCC).")
        ucc = safe_float("Enter UCC (W/m²·sr·µm) from header: ")
        # ASTER formula: (DN - 1) * UCC
        radiance = (dn - 1) * ucc
    else:
        # Landsat formula: DN * Mult + Add
        print("Enter Radiometric Rescaling Constants (from MTL file):")
        rad_mult = safe_float("  RADIANCE_MULT_BAND_x: ")
        rad_add  = safe_float("  RADIANCE_ADD_BAND_x: ")
        radiance = dn * rad_mult + rad_add

    # 2. Brightness Temperature
    print("Enter Thermal Constants (from metadata):")
    k1 = safe_float("  K1 constant: ")
    k2 = safe_float("  K2 constant: ")
    
    # IMPROVEMENT: Wavelength handling
    if default_wavelength:
        use_def = input(f"Use default center wavelength {default_wavelength} µm? (y/n): ")
        if use_def.lower() == 'y':
            wavelength = default_wavelength
        else:
            wavelength = safe_float("Enter central wavelength (µm): ")
    else:
        # Essential for ASTER (Band 10-14 vary from 8.3 to 11.3 µm)
        wavelength = safe_float("Enter central wavelength (µm) for this band: ")

    # Calculate BT
    # Suppress divide by zero for background pixels
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        radiance[radiance <= 0] = np.nan
        bt = k2 / np.log((k1 / radiance) + 1)

    # 3. LST Calculation
    # Formula: BT / (1 + (lambda * BT / rho) * ln(epsilon))
    lst_k = bt / (1 + (wavelength * bt / RHO_MICRON_K) * np.log(emissivity))
    return lst_k - 273.15

# ==============================
# PRODUCT-BASED LST
# ==============================
def calculate_product_lst(src, sensor_name):
    dn = src.read(1).astype(np.float32)
    tags = src.tags()
    
    # IMPROVEMENT: Try-catch for scale parsing
    scale = tags.get("scale_factor")
    offset = tags.get("add_offset", 0.0)

    if scale is None:
        print(f"\n[WARNING] 'scale_factor' missing in {sensor_name} metadata.")
        scale = safe_float("Enter scale factor manually (e.g., 0.02): ")
    else:
        scale = float(scale)
        print(f"[INFO] Auto-detected scale factor: {scale}")

    offset = float(offset)

    # Handle NoData
    if src.nodata is not None:
        dn[dn == src.nodata] = np.nan

    # Sensor-specific masking
    if sensor_name == "MODIS":
        print("[INFO] Applying MODIS valid range filter (7500-65535)")
        dn[(dn < 7500) | (dn > 65535)] = np.nan

    lst_k = dn * scale + offset
    return lst_k - 273.15

# ==============================
# MAIN DRIVER
# ==============================
def main():
    print("Universal Satellite LST Calculator (Registry-Based)")
    path_input = input("Enter thermal image path: ").strip('"')

    try:
        with rasterio.open(path_input) as src:
            dn = src.read(1).astype(np.float32)

            # 1. Detection Logic
            detected = detect_sensor_from_metadata(src)
            if detected:
                print(f"\n[SUCCESS] Auto-detected sensor: {detected}")
                sensor = detected
            else:
                sensor = ask_user_sensor()

            sensor_cfg = SENSOR_REGISTRY[sensor]

            # 2. Calculation Dispatcher
            if sensor_cfg["type"] == "RAW":
                # IMPROVEMENT: Emissivity Safety Check
                emissivity = safe_float("Enter surface emissivity (0.0 - 1.0): ")
                if emissivity < 0.9 or emissivity > 1.0:
                    print(f"[WARNING] Emissivity {emissivity} is unusual. Verify your input.")
                
                lst = calculate_sca_lst(dn, emissivity, sensor, sensor_cfg.get("default_wavelength"))

            else:
                lst = calculate_product_lst(src, sensor)

            # 3. Stats & Output
            valid = lst[~np.isnan(lst)]
            if valid.size == 0:
                print("[ERROR] Output is empty. Check your inputs.")
                return

            print("\n--- Results ---")
            print(f"Valid pixels: {len(valid)}")
            print(f"Mean: {np.mean(valid):.2f} °C")
            print(f"Min : {np.min(valid):.2f} °C")
            print(f"Max : {np.max(valid):.2f} °C")

            # Save Output
            out_name = f"LST_{sensor}_Calculated.tif"
            profile = src.profile.copy()
            profile.update(dtype=rasterio.float32, count=1, nodata=-9999)

            with rasterio.open(out_name, "w", **profile) as dst:
                dst.write(np.nan_to_num(lst, nan=-9999).astype(np.float32), 1)
                dst.update_tags(
                    SENSOR=sensor,
                    LST_UNITS="Celsius",
                    METHOD="SCA (RAW)" if sensor_cfg["type"] == "RAW" else "Product Scaling",
                    GENERATED_BY="Universal LST Pipeline"
                )

            print(f"\n[SUCCESS] Output saved: {out_name}")

    except Exception as e:
        print(f"\n[CRITICAL ERROR] {e}")

if __name__ == "__main__":
    main()
