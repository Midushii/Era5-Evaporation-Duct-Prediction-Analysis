import pandas as pd
import numpy as np
import pyhomogeneity as hg
import pymannkendall as mk
from scipy.stats import norm

df = pd.read_csv("master_1980_2025.csv")
df = df.groupby("year")[["sst_mean","t2m_mean","d2m_mean","deltaT_mean","wind_speed_mean"]].mean().reset_index()

def von_neumann_ratio_test(series):
    n = len(series)
    diffs = np.diff(series)
    N = np.sum(diffs**2) / np.sum((series - np.mean(series))**2)
    varN = 4*(n-2)/((n+1)*(n-1))
    z = (N-2)/np.sqrt(varN)
    p = 2*(1-norm.cdf(abs(z)))
    return N, p

for var in ["sst_mean","t2m_mean","d2m_mean","deltaT_mean","wind_speed_mean"]:

    years = df["year"].values
    series = df[var].values

    pettitt = hg.pettitt_test(series)
    snht = hg.snht_test(series)
    buishand = hg.buishand_range_test(series)
    vn, vn_p = von_neumann_ratio_test(series)
    mk_result = mk.original_test(series)

    coeffs = np.polyfit(years, series, 1)
    detrended = series - np.polyval(coeffs, years)
    det = hg.pettitt_test(detrended)

    print(f"\n{var}")
    print("Pettitt:", pettitt.p)
    print("SNHT:", snht.p)
    print("Buishand:", buishand.p)
    print("Von Neumann:", vn_p)
    print("Mann-Kendall:", mk_result.p)
    print("Detrended Pettitt:", det.p)

    