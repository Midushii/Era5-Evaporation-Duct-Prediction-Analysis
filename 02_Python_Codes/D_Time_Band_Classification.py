import os
import pandas as pd
import xarray as xr

folder = "era5_reanalysis"
records = []


def assign_band(hour):
    if 0 <= hour < 6:
        return "00-06"
    elif 6 <= hour < 12:
        return "06-12"
    elif 12 <= hour < 18:
        return "12-18"
    else:
        return "18-24"


for file in sorted(os.listdir(folder)):

    if not file.endswith(".grib"):
        continue

    path = os.path.join(folder, file)

    ds = xr.open_dataset(
        path,
        engine="cfgrib"
    )

    # ERA5 time(UTC) converting to IST
    times = (
        pd.to_datetime(ds["time"].values)
        + pd.Timedelta(hours=5, minutes=30)
    )

    temp = pd.DataFrame({
        "valid_time_ist": times
    })

    temp["hour_ist"] = temp["valid_time_ist"].dt.hour
    temp["time_band"] = temp["hour_ist"].apply(assign_band)

    records.append(temp)

    ds.close()


master = pd.concat(records, ignore_index=True)

master.to_csv("time_bands.csv", index=False)

print(master.head())
print("\nTime-band counts:")
print(master["time_band"].value_counts().sort_index())