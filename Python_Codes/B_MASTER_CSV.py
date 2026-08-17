import os
import re
import numpy as np
import pandas as pd
import xarray as xr

folder = "era5_reanalysis"      # Folder containing GRIB files
records = []

for file in sorted(os.listdir(folder)):

    if not file.endswith(".grib"):
        continue

    path = os.path.join(folder, file)

    # Extract year and month  
    m = re.search(r"(\d{4}).*?(\d{2})", file)
    year = int(m.group(1)) if m else None
    month = int(m.group(2)) if m else None

    ds = xr.open_dataset(path, engine="cfgrib")

    sst = ds["sst"].mean(dim=["latitude","longitude"]).values
    t2m = ds["t2m"].mean(dim=["latitude","longitude"]).values
    d2m = ds["d2m"].mean(dim=["latitude","longitude"]).values
    u10 = ds["u10"].mean(dim=["latitude","longitude"]).values
    v10 = ds["v10"].mean(dim=["latitude","longitude"]).values

    wind_speed = np.sqrt(u10**2 + v10**2)
    deltaT = t2m - sst

    records.append({
    "year": year,
    "month": month,
    "format": "grib",
    "sst_mean": float(np.nanmean(sst)),
    "t2m_mean": float(np.nanmean(t2m)),
    "d2m_mean": float(np.nanmean(d2m)),
    "wind_speed_mean": float(np.nanmean(wind_speed)),
    "deltaT_mean": float(np.nanmean(deltaT)),
    "n_hours": len(sst)
})

df = pd.DataFrame(records)
df = df.sort_values(["year", "month"])
df.to_csv("master_1980_2025.csv", index=False)

print(df.head())
print("\nSaved as master_1980_2025.csv")

