import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import re

LAT, LON = 19.0, 71.0
MONTHS = [f"{m:02d}" for m in range(1, 13)]
TIME_BANDS = ["00-06", "06-12", "12-18", "18-24"]

month_names = {
    "01": "January", "02": "February", "03": "March", "04": "April",
    "05": "May", "06": "June", "07": "July", "08": "August",
    "09": "September", "10": "October", "11": "November", "12": "December"
}

def saturation_vapor_pressure(T_K):
    T_C = T_K - 273.15
    return 6.1094 * np.exp(17.625 * T_C / (T_C + 243.04))

def assign_band(h):
    if 0 <= h < 6: return "00-06"
    elif 6 <= h < 12: return "06-12"
    elif 12 <= h < 18: return "12-18"
    else: return "18-24"

def plot_histogram(values, var_label, unit, bins, color, band, month_name, n_years, filename):
    values = np.array(values)
    values = values[~np.isnan(values)]
    if len(values) == 0:
        return
    counts, edges = np.histogram(values, bins=bins)
    pct = 100 * counts / counts.sum()
    bin_labels = [f"{int(edges[i])}-{int(edges[i+1])}" for i in range(len(edges)-1)]
    mean_val = np.mean(values)

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.barh(bin_labels, pct, color=color, edgecolor="black", height=0.6)
    ax.set_xlabel("Percent occurrence (%)")
    ax.set_ylabel(f"{var_label} ({unit})")
    ax.set_title(f"{month_name} Climatology (1980-2025)  Latitude: {LAT}N  Longitude: {LON}E\nTime: {band} IST (pooled across {n_years} years)", fontsize=11)
    fig.text(0.15, 0.02, f"■ {band} IST average: {mean_val:.1f} {unit}   n={len(values)}", color=color, fontsize=10)

    plt.tight_layout(rect=[0, 0.05, 1, 1])
    plt.savefig(filename, dpi=110)
    plt.close(fig)


pattern = re.compile(r"era5_(\d{4})_(\d{2})\.(nc|grib)$")
all_files = os.listdir("era5_raw")

base_out_dir = "Histograms/duct63_style/year_climatology"

for MONTH in MONTHS:
    month_name = month_names[MONTH]
    print(f"\n{'='*60}\nMONTH: {month_name}\n{'='*60}")

    year_file_map = {}
    for f in all_files:
        m = pattern.match(f)
        if not m:
            continue
        year, month, ext = m.groups()
        if month != MONTH or int(year) < 1980:
            continue
        if year not in year_file_map or (year_file_map[year][1] == "grib" and ext == "nc"):
            year_file_map[year] = (f, ext)

    print(f"Found {len(year_file_map)} years for {month_name}")

    pooled = {band: {"wind_speed": [], "abs_humidity": [], "deltaT": []} for band in TIME_BANDS}

    for year, (fname, ext) in sorted(year_file_map.items()):
        filepath = os.path.join("era5_raw", fname)
        try:
            ds = xr.open_dataset(filepath) if ext == "nc" else xr.open_dataset(filepath, engine="cfgrib")
            point = ds.sel(latitude=LAT, longitude=LON)

            times_utc = pd.to_datetime(point.valid_time.values)
            t2m = point.t2m.values
            d2m = point.d2m.values
            sst = point.sst.values
            u10 = point.u10.values
            v10 = point.v10.values
            wind_speed = np.sqrt(u10**2 + v10**2)
            e_hPa = saturation_vapor_pressure(d2m)
            abs_humidity = 216.7 * e_hPa / t2m
            deltaT = t2m - sst

            time_ist = times_utc + pd.Timedelta(hours=5, minutes=30)
            hour_ist = time_ist.hour
            bands = np.array([assign_band(h) for h in hour_ist])
            ds.close()

            for band in TIME_BANDS:
                mask = bands == band
                pooled[band]["wind_speed"].extend(wind_speed[mask])
                pooled[band]["abs_humidity"].extend(abs_humidity[mask])
                pooled[band]["deltaT"].extend(deltaT[mask])

            print(f"  {year}: processed")
        except Exception as e:
            print(f"  {year}: FAILED - {e}")

    duct_pooled = {band: [] for band in TIME_BANDS}
    for year in year_file_map.keys():
        duct_path = f"Tables/duct_climatology/duct_{year}_{MONTH}.csv"
        if not os.path.exists(duct_path):
            continue
        duct_df = pd.read_csv(duct_path, parse_dates=["valid_time"])
        duct_df = duct_df[(duct_df["latitude"] == LAT) & (duct_df["longitude"] == LON)]
        duct_df["time_ist"] = duct_df["valid_time"] + pd.Timedelta(hours=5, minutes=30)
        duct_df["hour_ist"] = duct_df["time_ist"].dt.hour
        duct_df["time_band"] = duct_df["hour_ist"].apply(assign_band)
        for band in TIME_BANDS:
            vals = duct_df[duct_df["time_band"] == band]["duct_height_m"].dropna().values
            duct_pooled[band].extend(vals)

    month_out_dir = os.path.join(base_out_dir, MONTH)
    os.makedirs(month_out_dir, exist_ok=True)

    n_years = len(year_file_map)
    for band in TIME_BANDS:
        plot_histogram(duct_pooled[band], "Evaporation Duct Height", "m",
                       np.arange(0, 82, 4), "tab:red", band, month_name, n_years,
                       f"{month_out_dir}/duct_height_{band}.png")
        plot_histogram(pooled[band]["wind_speed"], "Surface Wind Speed", "m/s",
                       np.arange(0, 22, 2), "tab:blue", band, month_name, n_years,
                       f"{month_out_dir}/wind_speed_{band}.png")
        plot_histogram(pooled[band]["abs_humidity"], "Absolute Humidity", "g/m3",
                       np.arange(0, 32, 2), "tab:green", band, month_name, n_years,
                       f"{month_out_dir}/abs_humidity_{band}.png")
        plot_histogram(pooled[band]["deltaT"], "Air/Sea Temp Difference", "K",
                       np.arange(-6, 6, 1), "tab:orange", band, month_name, n_years,
                       f"{month_out_dir}/deltaT_{band}.png")

    print(f"  16 climatology plots saved to {month_out_dir}/")

print("\nDONE.")