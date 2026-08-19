import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import re

LAT, LON = 19.0, 71.0
k = 0.4
g = 9.81
z_air = 2.0
z_wind = 10.0
HEIGHT_MAX = 80
N_LEVELS = 400
heights_global = np.linspace(0.05, HEIGHT_MAX, N_LEVELS)
BOUNDARY_ZONE_FRACTION = 0.10
boundary_zone_start_idx = int(N_LEVELS * (1 - BOUNDARY_ZONE_FRACTION))

def saturation_vapor_pressure(T_K):
    T_C = T_K - 273.15
    return 6.1094 * np.exp(17.625 * T_C / (T_C + 243.04))

def specific_humidity(e_hPa, P_hPa):
    return 0.622 * e_hPa / (P_hPa - 0.378 * e_hPa)

def psi_m(zeta):
    zeta = np.atleast_1d(zeta).astype(float)
    x = np.full_like(zeta, np.nan)
    unstable_mask = zeta < 0
    x[unstable_mask] = (1 - 16*zeta[unstable_mask])**0.25
    unstable = np.zeros_like(zeta)
    unstable[unstable_mask] = (
        2*np.log((1+x[unstable_mask])/2)
        + np.log((1+x[unstable_mask]**2)/2)
        - 2*np.arctan(x[unstable_mask]) + np.pi/2
    )
    stable = -5*zeta
    result = np.where(zeta < 0, unstable, stable)
    return result if result.size > 1 else result.item()

def psi_h(zeta):
    zeta = np.atleast_1d(zeta).astype(float)
    x = np.full_like(zeta, np.nan)
    unstable_mask = zeta < 0
    x[unstable_mask] = (1 - 16*zeta[unstable_mask])**0.25
    unstable = np.zeros_like(zeta)
    unstable[unstable_mask] = 2*np.log((1+x[unstable_mask]**2)/2)
    stable = -5*zeta
    result = np.where(zeta < 0, unstable, stable)
    return result if result.size > 1 else result.item()

def refractivity(T_K, q_kgkg, P_hPa):
    e_hPa = q_kgkg * P_hPa / (0.622 + 0.378*q_kgkg)
    N = 77.6*(P_hPa/T_K) + 3.73e5*(e_hPa/T_K**2)
    return N

def get_M_profile(T2m_K, Td2m_K, SST_K, MSLP_hPa, U10):
    try:
        e_air = saturation_vapor_pressure(Td2m_K)
        q_air = specific_humidity(e_air, MSLP_hPa)
        e_sea = saturation_vapor_pressure(SST_K)
        q_sea = 0.98 * specific_humidity(e_sea, MSLP_hPa)

        Ta_virtual = T2m_K * (1 + 0.61*q_air)
        delta_T = T2m_K - SST_K
        delta_q = q_air - q_sea

        z0 = 0.0002
        L = 1e6

        if U10 < 0.5:
            U10 = 0.5

        for _ in range(20):
            u_star = k * U10 / (
                np.log(z_wind/z0) - psi_m(z_wind/L)
            )

            T_star = k * delta_T / (
                np.log(z_air/z0) - psi_h(z_air/L)
            )

            q_star = k * delta_q / (
                np.log(z_air/z0) - psi_h(z_air/L)
            )

            if u_star <= 0 or not np.isfinite(u_star):
                return None

            L_new = (
                u_star**2 * Ta_virtual
            ) / (
                k * g * (T_star + 0.61*Ta_virtual*q_star)
            )

            z0_new = (
                0.011 * u_star**2 / g
                + 0.11 * 1.5e-5 / u_star
            )

            if not np.isfinite(L_new) or not np.isfinite(z0_new):
                return None

            if abs(L_new - L) < 0.01 and abs(z0_new - z0) < 1e-6:
                L, z0 = L_new, z0_new
                break

            L, z0 = L_new, z0_new

        T_profile = SST_K + (T_star/k) * (
            np.log(heights_global/z0)
            - psi_h(heights_global/L)
        )

        q_profile = q_sea + (q_star/k) * (
            np.log(heights_global/z0)
            - psi_h(heights_global/L)
        )

        P_profile = MSLP_hPa * np.exp(
            -heights_global * g / (287.05 * T2m_K)
        )

        N_profile = refractivity(
            T_profile,
            q_profile,
            P_profile
        )

        return N_profile + 0.157 * heights_global

    except Exception:
        return None

def plot_duct_profile(mean_M, duct_height, title, filename, n_hours):
    fig, ax = plt.subplots(figsize=(7, 8))

    ax.plot(
        mean_M,
        heights_global,
        color="black",
        linewidth=1.8
    )

    if duct_height is not None:
        ax.axhline(
            duct_height,
            color="red",
            linestyle="--",
            linewidth=1.2,
            label=f"Duct Height = {duct_height:.2f} m"
        )

    ax.set_xlabel("Modified Refractivity, M (M-units)")
    ax.set_ylabel("Height (m)")
    ax.set_title(
        f"{title}\nn = {n_hours} hours"
    )

    ax.set_ylim(0, HEIGHT_MAX)
    ax.grid(alpha=0.3)

    if duct_height is not None:
        ax.legend()

    plt.tight_layout()
    plt.savefig(filename, dpi=110)
    plt.close(fig)

month_names = {
    "01": "January",
    "02": "February",
    "03": "March",
    "04": "April",
    "05": "May",
    "06": "June",
    "07": "July",
    "08": "August",
    "09": "September",
    "10": "October",
    "11": "November",
    "12": "December"
}

pattern = re.compile(r"era5_(\d{4})_(\d{2})\.(nc|grib)$")
all_files = os.listdir("era5_raw")

base_out_dir = "Histograms/duct63_style/profiles_climatology"
os.makedirs(base_out_dir, exist_ok=True)

for MONTH, month_name in month_names.items():

    out_path = os.path.join(
        base_out_dir,
        f"profile_climatology_{MONTH}.png"
    )

    if os.path.exists(out_path):
        print(f"{month_name}: already exists, skipping")
        continue

    print(f"Processing {month_name} (climatology, all years)...")

    year_file_map = {}

    for f in all_files:
        m = pattern.match(f)

        if not m:
            continue

        year, month, ext = m.groups()

        if month != MONTH or int(year) < 1980:
            continue

        if (
            year not in year_file_map
            or (
                year_file_map[year][1] == "grib"
                and ext == "nc"
            )
        ):
            year_file_map[year] = (f, ext)

    all_profiles = []

    for year, (fname, ext) in sorted(year_file_map.items()):

        filepath = os.path.join("era5_raw", fname)

        try:
            ds = (
                xr.open_dataset(filepath)
                if ext == "nc"
                else xr.open_dataset(filepath, engine="cfgrib")
            )

            point = ds.sel(
                latitude=LAT,
                longitude=LON
            )

            T2m_arr = point.t2m.values
            Td2m_arr = point.d2m.values
            SST_arr = point.sst.values
            MSLP_arr = point.msl.values / 100.0
            u10_arr = point.u10.values
            v10_arr = point.v10.values

            for i in range(len(T2m_arr)):

                U10 = np.sqrt(
                    u10_arr[i]**2 +
                    v10_arr[i]**2
                )

                prof = get_M_profile(
                    T2m_arr[i],
                    Td2m_arr[i],
                    SST_arr[i],
                    MSLP_arr[i],
                    U10
                )

                if prof is not None:
                    all_profiles.append(prof)

            ds.close()

        except Exception as e:
            print(f"    {year}: ERROR - {e}")

    if len(all_profiles) == 0:
        print(f"  {month_name}: no valid profiles, skipping")
        continue

    mean_M = np.mean(
        all_profiles,
        axis=0
    )

    min_idx = np.argmin(mean_M)

    if min_idx < boundary_zone_start_idx:
        duct_height = heights_global[min_idx]
    else:
        duct_height = None

    plot_duct_profile(
        mean_M,
        duct_height,
        f"{month_name} Climatology (1980-2025)\nMean Evaporation Duct Profile ({LAT}N, {LON}E)",
        out_path,
        len(all_profiles)
    )

    print(
        f"  {month_name}: done, "
        f"n={len(all_profiles)} hours "
        f"({len(year_file_map)} years), "
        f"duct={duct_height}"
    )

print(
    "\nDONE.",
    base_out_dir
)