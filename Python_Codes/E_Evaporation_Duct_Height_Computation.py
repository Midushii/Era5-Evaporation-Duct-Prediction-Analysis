import numpy as np
import xarray as xr
k = 0.4          # von Karman constant
g = 9.81         # gravity (m/s^2)
z_air = 2.0      # height of temp/humidity measurement (m)
z_wind = 10.0    # height of wind measurement (m)
LSM_THRESHOLD = 0.05
HEIGHT_MAX = 80
N_LEVELS = 400
heights = np.linspace(0.05, HEIGHT_MAX, N_LEVELS)
BOUNDARY_ZONE_FRACTION = 0.10
boundary_zone_start_idx = int(N_LEVELS * (1 - BOUNDARY_ZONE_FRACTION))
# select valid ocean grid points
static = xr.open_dataset("era5_static.nc")
lsm = static.lsm.isel(valid_time=0)

ocean_points = []
for lat_val in lsm.latitude.values:
    for lon_val in lsm.longitude.values:
        lsm_val = lsm.sel(latitude=lat_val, longitude=lon_val).values.item()
        if lsm_val < LSM_THRESHOLD:
            ocean_points.append((lat_val, lon_val))
# 13 ocean grid points 
def saturation_vapor_pressure(T_K):
    T_C = T_K - 273.15
    return 6.1094 * np.exp(17.625 * T_C / (T_C + 243.04))

def specific_humidity(e_hPa, P_hPa):
    return 0.622 * e_hPa / (P_hPa - 0.378 * e_hPa)

def psi_m(zeta):  # momentum stability correction
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

def psi_h(zeta):  # heat/humidity stability correction
    zeta = np.atleast_1d(zeta).astype(float)
    x = np.full_like(zeta, np.nan)
    unstable_mask = zeta < 0
    x[unstable_mask] = (1 - 16*zeta[unstable_mask])**0.25
    unstable = np.zeros_like(zeta)
    unstable[unstable_mask] = 2*np.log((1+x[unstable_mask]**2)/2)
    stable = -5*zeta
    result = np.where(zeta < 0, unstable, stable)
    return result if result.size > 1 else result.item()

# ITU-R P.453 REFRACTIVITY 
def refractivity(T_K, q_kgkg, P_hPa):
    e_hPa = q_kgkg * P_hPa / (0.622 + 0.378*q_kgkg)
    N = 77.6*(P_hPa/T_K) + 3.73e5*(e_hPa/T_K**2)
    return N
# Paulus-Jeske profile reconstruction 
def compute_duct_height(T2m_K, Td2m_K, SST_K, MSLP_hPa, U10):
    try:
        
        e_air = saturation_vapor_pressure(Td2m_K)
        q_air = specific_humidity(e_air, MSLP_hPa)
        e_sea = saturation_vapor_pressure(SST_K)
        q_sea = 0.98 * specific_humidity(e_sea, MSLP_hPa)  # 0.98 :- salinity correction

        Ta_virtual = T2m_K * (1 + 0.61*q_air)
        delta_T = T2m_K - SST_K
        delta_q = q_air - q_sea

        z0 = 0.0002
        L = 1e6
        if U10 < 0.5:
            U10 = 0.5

        
        for _ in range(20):
            u_star = k * U10 / (np.log(z_wind/z0) - psi_m(z_wind/L))
            T_star = k * delta_T / (np.log(z_air/z0) - psi_h(z_air/L))
            q_star = k * delta_q / (np.log(z_air/z0) - psi_h(z_air/L))

            if u_star <= 0 or not np.isfinite(u_star):
                return np.nan, "invalid_ustar"
            L_new = (u_star**2 * Ta_virtual) / (k * g * (T_star + 0.61*Ta_virtual*q_star))
            z0_new = 0.011 * u_star**2 / g + 0.11 * 1.5e-5 / u_star  # Charnock relation
            if not np.isfinite(L_new) or not np.isfinite(z0_new):
                return np.nan, "invalid_L_or_z0"
            if abs(L_new - L) < 0.01 and abs(z0_new - z0) < 1e-6:
                L, z0 = L_new, z0_new
                break
            L, z0 = L_new, z0_new
        #Vertical profile from 0.05 m to 80 m
        T_profile = SST_K + (T_star/k) * (np.log(heights/z0) - psi_h(heights/L))
        q_profile = q_sea + (q_star/k) * (np.log(heights/z0) - psi_h(heights/L))
        P_profile = MSLP_hPa * np.exp(-heights * g / (287.05 * T2m_K))       
        N_profile = refractivity(T_profile, q_profile, P_profile)
        M_profile = N_profile + 0.157 * heights   # modified refractivity       
        min_idx = np.argmin(M_profile)
        if min_idx >= boundary_zone_start_idx:  # last 10% of search range
            return np.nan, "boundary_hit"
        if min_idx == 0:
            return np.nan, "surface_hit"
        return heights[min_idx], None
    except Exception:
        return np.nan, "exception"


    