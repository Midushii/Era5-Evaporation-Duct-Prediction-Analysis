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
boundary_zone_start_idx = int(
    N_LEVELS * (1 - BOUNDARY_ZONE_FRACTION)
)


def saturation_vapor_pressure(T_K):
    T_C = T_K - 273.15
    return 6.1094 * np.exp(
        17.625 * T_C / (T_C + 243.04)
    )


def specific_humidity(e_hPa, P_hPa):
    return 0.622 * e_hPa / (
        P_hPa - 0.378 * e_hPa
    )


def psi_m(zeta):
    zeta = np.atleast_1d(zeta).astype(float)

    x = np.full_like(zeta, np.nan)
    unstable_mask = zeta < 0

    x[unstable_mask] = (
        1 - 16 * zeta[unstable_mask]
    ) ** 0.25

    unstable = np.zeros_like(zeta)

    unstable[unstable_mask] = (
        2 * np.log(
            (1 + x[unstable_mask]) / 2
        )
        + np.log(
            (1 + x[unstable_mask] ** 2) / 2
        )
        - 2 * np.arctan(
            x[unstable_mask]
        )
        + np.pi / 2
    )

    stable = -5 * zeta

    result = np.where(
        zeta < 0,
        unstable,
        stable
    )

    return result if result.size > 1 else result.item()


def psi_h(zeta):
    zeta = np.atleast_1d(zeta).astype(float)

    x = np.full_like(zeta, np.nan)
    unstable_mask = zeta < 0

    x[unstable_mask] = (
        1 - 16 * zeta[unstable_mask]
    ) ** 0.25

    unstable = np.zeros_like(zeta)

    unstable[unstable_mask] = (
        2 * np.log(
            (1 + x[unstable_mask] ** 2) / 2
        )
    )

    stable = -5 * zeta

    result = np.where(
        zeta < 0,
        unstable,
        stable
    )

    return result if result.size > 1 else result.item()


def refractivity(T_K, q_kgkg, P_hPa):
    e_hPa = (
        q_kgkg * P_hPa
        / (0.622 + 0.378 * q_kgkg)
    )

    N = (
        77.6 * (P_hPa / T_K)
        + 3.73e5 * (e_hPa / T_K**2)
    )

    return N


def get_M_profile(
    T2m_K,
    Td2m_K,
    SST_K,
    MSLP_hPa,
    U10
):
    try:
        e_air = saturation_vapor_pressure(Td2m_K)
        q_air = specific_humidity(
            e_air,
            MSLP_hPa
        )

        e_sea = saturation_vapor_pressure(SST_K)
        q_sea = (
            0.98
            * specific_humidity(
                e_sea,
                MSLP_hPa
            )
        )

        Ta_virtual = (
            T2m_K
            * (1 + 0.61 * q_air)
        )

        delta_T = T2m_K - SST_K
        delta_q = q_air - q_sea

        z0 = 0.0002
        L = 1e6

        if U10 < 0.5:
            U10 = 0.5

        for _ in range(20):

            u_star = (
                k * U10
                /
                (
                    np.log(z_wind / z0)
                    - psi_m(z_wind / L)
                )
            )

            T_star = (
                k * delta_T
                /
                (
                    np.log(z_air / z0)
                    - psi_h(z_air / L)
                )
            )

            q_star = (
                k * delta_q
                /
                (
                    np.log(z_air / z0)
                    - psi_h(z_air / L)
                )
            )

            if (
                u_star <= 0
                or not np.isfinite(u_star)
            ):
                return None

            L_new = (
                u_star**2
                * Ta_virtual
                /
                (
                    k * g
                    * (
                        T_star
                        + 0.61
                        * Ta_virtual
                        * q_star
                    )
                )
            )

            z0_new = (
                0.011 * u_star**2 / g
                + 0.11 * 1.5e-5 / u_star
            )

            if (
                not np.isfinite(L_new)
                or not np.isfinite(z0_new)
            ):
                return None

            if (
                abs(L_new - L) < 0.01
                and
                abs(z0_new - z0) < 1e-6
            ):
                L = L_new
                z0 = z0_new
                break

            L = L_new
            z0 = z0_new

        T_profile = (
            SST_K
            + (T_star / k)
            * (
                np.log(
                    heights_global / z0
                )
                - psi_h(
                    heights_global / L
                )
            )
        )

        q_profile = (
            q_sea
            + (q_star / k)
            * (
                np.log(
                    heights_global / z0
                )
                - psi_h(
                    heights_global / L
                )
            )
        )

        P_profile = (
            MSLP_hPa
            * np.exp(
                -heights_global * g
                /
                (287.05 * T2m_K)
            )
        )

        N_profile = refractivity(
            T_profile,
            q_profile,
            P_profile
        )

        M_profile = (
            N_profile
            + 0.157 * heights_global
        )

        return M_profile

    except Exception:
        return None


def plot_duct_profile(
    mean_M,
    duct_height,
    year,
    n_hours,
    filename
):
    fig, ax = plt.subplots(
        figsize=(8, 7)
    )

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
            linewidth=1.2
        )

        ax.text(
            mean_M.min(),
            duct_height + 1,
            f"Duct height = {duct_height:.2f} m",
            fontsize=10
        )

    ax.set_xlabel(
        "Modified Refractivity, M (M-units)"
    )

    ax.set_ylabel(
        "Height (m)"
    )

    ax.set_ylim(
        0,
        HEIGHT_MAX
    )

    ax.set_title(
        f"{year} Mean Evaporation Duct Profile\n"
        f"All Months | {LAT}°N, {LON}°E | "
        f"n = {n_hours} hours"
    )

    ax.grid(alpha=0.3)

    plt.tight_layout()

    plt.savefig(
        filename,
        dpi=150,
        bbox_inches="tight"
    )

    plt.close(fig)


pattern = re.compile(
    r"era5_(\d{4})_(\d{2})\.(nc|grib)$"
)

all_files = os.listdir("era5_raw")
year_month_map = {}

for f in all_files:

    m = pattern.match(f)

    if not m:
        continue

    year, month, ext = m.groups()

    if int(year) < 1980:
        continue

    if year not in year_month_map:
        year_month_map[year] = {}

    if (
        month not in year_month_map[year]
        or
        (
            year_month_map[year][month][0] == "grib"
            and ext == "nc"
        )
    ):
        year_month_map[year][month] = (
            ext,
            f
        )


years = sorted(
    year_month_map.keys()
)

print(
    f"Years found: {len(years)}"
)

base_out_dir = (
    "Histograms/duct63_style/"
    "profiles_per_year"
)

os.makedirs(
    base_out_dir,
    exist_ok=True
)


for yidx, year in enumerate(
    years,
    1
):

    out_path = os.path.join(
        base_out_dir,
        f"profile_{year}.png"
    )

    if os.path.exists(out_path):

        print(
            f"[{yidx}/{len(years)}] "
            f"{year}: already exists, skipping"
        )

        continue

    all_profiles = []
    months_in_year = year_month_map[year]

    for month, (ext, fname) in sorted(
        months_in_year.items()
    ):

        filepath = os.path.join(
            "era5_raw",
            fname
        )

        try:

            if ext == "nc":
                ds = xr.open_dataset(
                    filepath
                )
            else:
                ds = xr.open_dataset(
                    filepath,
                    engine="cfgrib"
                )

            point = ds.sel(
                latitude=LAT,
                longitude=LON
            )

            T2m_arr = point.t2m.values
            Td2m_arr = point.d2m.values
            SST_arr = point.sst.values

            MSLP_arr = (
                point.msl.values / 100.0
            )

            u10_arr = point.u10.values
            v10_arr = point.v10.values

            for i in range(
                len(T2m_arr)
            ):

                U10 = np.sqrt(
                    u10_arr[i]**2
                    + v10_arr[i]**2
                )

                prof = get_M_profile(
                    T2m_arr[i],
                    Td2m_arr[i],
                    SST_arr[i],
                    MSLP_arr[i],
                    U10
                )

                if prof is not None:
                    all_profiles.append(
                        prof
                    )

            ds.close()

        except Exception as e:

            print(
                f"    {year}-{month}: "
                f"ERROR - {e}"
            )

    if len(all_profiles) == 0:

        print(
            f"[{yidx}/{len(years)}] "
            f"{year}: no valid profiles"
        )

        continue

    mean_M = np.mean(
        all_profiles,
        axis=0
    )

    min_idx = np.argmin(
        mean_M
    )

    if (
        min_idx < boundary_zone_start_idx
        and min_idx != 0
    ):
        duct_height = (
            heights_global[min_idx]
        )
    else:
        duct_height = None

    plot_duct_profile(
        mean_M,
        duct_height,
        year,
        len(all_profiles),
        out_path
    )

    print(
        f"[{yidx}/{len(years)}] "
        f"{year}: done | "
        f"n={len(all_profiles)} hours | "
        f"duct={duct_height}"
    )


print("\nDONE.")
print(
    "Yearly duct profile plots saved to:",
    base_out_dir
)