import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import re

LAT, LON = 19.0, 71.0
TIME_BANDS = ["00-06", "06-12", "12-18", "18-24"]

month_names = {
    "01": "January", "02": "February", "03": "March", "04": "April",
    "05": "May", "06": "June", "07": "July", "08": "August",
    "09": "September", "10": "October", "11": "November", "12": "December"
}
month_order_top_to_bottom = list(reversed(list(month_names.items())))  

def saturation_vapor_pressure(T_K):
    T_C = T_K - 273.15
    return 6.1094 * np.exp(17.625 * T_C / (T_C + 243.04))

def assign_band(h):
    if 0 <= h < 6: return "00-06"
    elif 6 <= h < 12: return "06-12"
    elif 12 <= h < 18: return "12-18"
    else: return "18-24"

def extract_year_month_data(year, month, lat, lon):
    """Returns dict: {band: {'wind_speed':arr, 'abs_humidity':arr, 'deltaT':arr}}, and duct_by_band"""
    candidates = [f"era5_raw/era5_{year}_{month}.nc", f"era5_raw/era5_{year}_{month}.grib"]
    filepath = next((f for f in candidates if os.path.exists(f)), None)
    if filepath is None:
        return None, None

    ext = "nc" if filepath.endswith(".nc") else "grib"
    ds = xr.open_dataset(filepath) if ext == "nc" else xr.open_dataset(filepath, engine="cfgrib")
    point = ds.sel(latitude=lat, longitude=lon)

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

    met_by_band = {}
    for band in TIME_BANDS:
        mask = bands == band
        met_by_band[band] = {
            "wind_speed": wind_speed[mask],
            "abs_humidity": abs_humidity[mask],
            "deltaT": deltaT[mask]
        }

    duct_by_band = {band: np.array([]) for band in TIME_BANDS}
    duct_path = f"Tables/duct_climatology/duct_{year}_{month}.csv"
    if os.path.exists(duct_path):
        duct_df = pd.read_csv(duct_path, parse_dates=["valid_time"])
        duct_df = duct_df[(duct_df["latitude"] == lat) & (duct_df["longitude"] == lon)]
        duct_df["time_ist"] = duct_df["valid_time"] + pd.Timedelta(hours=5, minutes=30)
        duct_df["hour_ist"] = duct_df["time_ist"].dt.hour
        duct_df["time_band"] = duct_df["hour_ist"].apply(assign_band)
        for band in TIME_BANDS:
            vals = duct_df[duct_df["time_band"] == band]["duct_height_m"].dropna().values
            duct_by_band[band] = vals

    return met_by_band, duct_by_band


def plot_yearly_summary(month_stats, var_label, unit, color, band, year, filename):
    """month_stats: list of (month_name, mean_val, n) in top-to-bottom display order"""
    months_present = [m for m in month_stats if m[2] > 0]
    if not months_present:
        return

    labels = [m[0] for m in months_present]
    means = [m[1] for m in months_present]
    ns = [m[2] for m in months_present]

    total_n = sum(ns)
    lo, hi = min(means), max(means)
    pad = max((hi - lo) * 0.15, 0.1)
    xlim_lo, xlim_hi = lo - pad, hi + pad

    fig, ax = plt.subplots(figsize=(9, 7))
    bars = ax.barh(labels, means, color=color, edgecolor="black", height=0.6, left=xlim_lo)
    

    for bar, n in zip(bars, ns):
        ax.text(bar.get_x() + bar.get_width()*0.02, bar.get_y() + bar.get_height()/2,
                f"# Obs: {n}", va="center", ha="left", fontsize=8, fontweight="bold", color="black")

    ax.set_xlim(xlim_lo, xlim_hi)
    ax.set_xlabel(f"Average {var_label} ({unit})")
    ax.set_ylabel("Months")
    ax.set_title(f"{year}   Total Obs: {total_n}\nTime: {band} IST", fontsize=12)

    plt.tight_layout()
    plt.savefig(filename, dpi=110)
    plt.close(fig)


var_configs = [
    ("duct_height", "Evaporation Duct Height", "m", "tab:red"),
    ("wind_speed", "Surface Wind Speed", "m/s", "tab:blue"),
    ("abs_humidity", "Absolute Humidity", "g/m3", "tab:green"),
    ("deltaT", "Air/Sea Temp Difference", "K", "tab:orange"),
]

years = [str(y) for y in range(1980, 2026)]
base_out_dir = "Cumulative"

for year in years:
    print(f"Processing {year}...")
    year_out_dir = os.path.join(base_out_dir, year)
    os.makedirs(year_out_dir, exist_ok=True)

    
    all_month_data = {} 
    for month in month_names.keys():
        met, duct = extract_year_month_data(year, month, LAT, LON)
        if met is not None:
            all_month_data[month] = (met, duct)

    for band in TIME_BANDS:
        for var_key, var_label, unit, color in var_configs:
            month_stats = []
            for month, month_name in month_order_top_to_bottom:
                if month not in all_month_data:
                    continue
                met, duct = all_month_data[month]
                if var_key == "duct_height":
                    vals = duct[band]
                else:
                    vals = met[band][var_key]
                vals = vals[~np.isnan(vals)] if len(vals) > 0 else vals
                if len(vals) > 0:
                    month_stats.append((month_name, np.mean(vals), len(vals)))
                else:
                    month_stats.append((month_name, 0, 0))

            plot_yearly_summary(
                month_stats, var_label, unit, color, band, year,
                f"{year_out_dir}/{var_key}_{band}_yearly_summary.png"
            )

    print(f"  Saved 16 yearly-summary plots for {year}")

print("\nDONE.")