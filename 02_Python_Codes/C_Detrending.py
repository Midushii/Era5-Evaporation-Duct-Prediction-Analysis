import pandas as pd
import numpy as np
import pyhomogeneity as hg

df = pd.read_csv("master_1980_2025.csv")
df = df.groupby("year")[["sst_mean","t2m_mean","d2m_mean","deltaT_mean","wind_speed_mean"]].mean().reset_index()

years = df["year"].values

for var in ["sst_mean", "t2m_mean", "d2m_mean", "deltaT_mean", "wind_speed_mean"]:

    series = df[var].values
    coeffs = np.polyfit(years, series, 1)
    trend_line = np.polyval(coeffs, years)
    detrended = series - trend_line

    r = hg.pettitt_test(detrended)

    print(f"\n{var}")
    print("Change Point:", years[r.cp] if r.cp is not None else "None")
    print(f"p-value: {r.p:.4f}")
    print("Break Persists:", r.p < 0.05)


    