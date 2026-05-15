# prosail_run.py
import numpy as np
from pypro4sail import pypro4sail, four_sail, prospect

def run_prosail_single(LAI, skyl, solar_zenith, solar_azimuth, view_zenith=10, view_azimuth=180):
    params = {
        "N": 1.5,
        "chloro": 40,
        "caroten": 8,
        "brown": 0.0,
        "EWT": 0.01,
        "LMA": 0.009,
        "Ant": 1.0,
        "hot_spot": 0.01,
        "solar_zenith": solar_zenith,
        "solar_azimuth": solar_azimuth,
        "view_zenith": view_zenith,
        "view_azimuth": view_azimuth,
        "LIDF": (-0.35, -0.15),
        "skyl": skyl,
        "soilType": pypro4sail.DEFAULT_SOIL,
        "LAI": LAI
    }

    soil_path = four_sail.SOIL_LIBRARY / params["soilType"]
    rsoil = np.genfromtxt(soil_path)[:, 1]

    wl, rho_leaf, tau_leaf = prospect.prospectd(
        params["N"], params["chloro"], params["caroten"],
        params["brown"], params["EWT"], params["LMA"], params["Ant"]
    )
    alpha_leaf = 1 - rho_leaf - tau_leaf

    lidf = four_sail.calc_lidf_verhoef(*params["LIDF"])
    psi = abs(params["solar_azimuth"] - params["view_azimuth"])

    sail_out = four_sail.foursail(
        LAI, params["hot_spot"], lidf,
        solar_zenith, view_zenith, psi,
        rho_leaf, tau_leaf, rsoil
    )

    rdot = sail_out[14]
    rsot = sail_out[17]
    tau_canopy = sail_out[15] * skyl + rsot * (1 - skyl)
    rho_canopy = rsot * (1 - skyl) + rdot * skyl
    alpha_canopy = 1 - rho_canopy - tau_canopy

    return {
        "Wavelength_nm": wl,
        "Leaf_Reflectance": rho_leaf,
        "Leaf_Transmittance": tau_leaf,
        "Leaf_Absorptance": alpha_leaf,
        "Canopy_Reflectance": rho_canopy,
        "Canopy_Transmittance": tau_canopy,
        "Canopy_Absorptance": alpha_canopy
    }
