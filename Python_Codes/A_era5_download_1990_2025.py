import cdsapi
import os

# CDS API 
client = cdsapi.Client()

dataset = "reanalysis-era5-single-levels"

# Output folder
output_folder = "era5_reanalysis"
os.makedirs(output_folder, exist_ok=True)

# Output file
target = os.path.join(output_folder, "era5_mumbai_1940_2024.grib")

# Years and months
years = [str(y) for y in range(1940, 2026)]  # 1940-2025
months = [f"{m:02d}" for m in range(1, 13)]

# Variables
variables = [
    "10m_u_component_of_wind",
    "10m_v_component_of_wind",
    "2m_dewpoint_temperature",
    "2m_temperature",
    "mean_sea_level_pressure",
    "sea_surface_temperature"
]

# Days and hours
days = [f"{d:02d}" for d in range(1, 32)]
hours = [f"{h:02d}:00" for h in range(24)]

# Request
request = {
    "product_type": ["reanalysis"],
    "variable": variables,
    "year": years,
    "month": months,
    "day": days,
    "time": hours,
    "data_format": "grib",
    "download_format": "unarchived",
    "area": [20, 70, 17, 74],   # North, West, South, East
    "grid": [1.0, 1.0]
}

# Download
client.retrieve(dataset, request).download(target)

print(f"\nDownload complete.\nSaved to: {target}")