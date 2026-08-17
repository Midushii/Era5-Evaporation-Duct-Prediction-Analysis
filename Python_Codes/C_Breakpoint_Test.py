import pandas as pd
import ruptures as rpt

df = pd.read_csv("master_1980_2025.csv")
df = df.groupby("year")[["sst_mean","t2m_mean","d2m_mean","deltaT_mean","wind_speed_mean"]].mean().reset_index()

for var in ["sst_mean","t2m_mean","d2m_mean","deltaT_mean","wind_speed_mean"]:
    algo = rpt.Pelt(model="rbf").fit(df[var].values)
    bp = algo.predict(pen=3)
    years = [df.loc[i, "year"] for i in bp[:-1]]
    print(f"\n{var}")
    print("Breakpoints:", years if years else "None")


    