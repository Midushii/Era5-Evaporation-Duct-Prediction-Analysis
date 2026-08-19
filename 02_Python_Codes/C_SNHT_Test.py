import pandas as pd
import pyhomogeneity as hg

df = pd.read_csv("master_1980_2025.csv")
df = df.groupby("year")[["sst_mean","t2m_mean","d2m_mean","deltaT_mean","wind_speed_mean"]].mean().reset_index()

for var in ["sst_mean", "t2m_mean", "d2m_mean", "deltaT_mean", "wind_speed_mean"]:
    r = hg.snht_test(df[var])

    print(f"\n{var}")
    print("Change Point:", df.loc[r.cp, "year"] if r.cp is not None else "None")
    print("p-value:", r.p)
    print("Significant:", r.p < 0.05)

    