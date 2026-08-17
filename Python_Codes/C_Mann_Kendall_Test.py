import pandas as pd
import pymannkendall as mk

df = pd.read_csv("master_1980_2025.csv")
df = df.groupby("year")[["sst_mean","t2m_mean","d2m_mean","deltaT_mean","wind_speed_mean"]].mean().reset_index()

for var in ["sst_mean", "t2m_mean", "d2m_mean", "deltaT_mean", "wind_speed_mean"]:
    r = mk.original_test(df[var])

    print(f"\n{var}")
    print("Trend:", r.trend)
    print(f"p-value: {r.p:.4f}")
    print(f"Slope: {r.slope:.5f} per year")
    print("Significant:", r.p < 0.05)



    