"""
I2EM Double-Bounce Water-Wall Scattering Model for Urban Flood Detection.

Implements the dihedral (double-bounce) scattering mechanism where radar
energy reflects off a rough water surface and then a building wall (or
vice versa). This is the dominant scattering mechanism for detecting
urban flooding in SAR imagery (e.g., ICEYE X-band).

Physics:
    Path 1: Radar → water surface (forward scatter) → wall → back to radar
    Path 2: Radar → wall → water surface (forward scatter) → back to radar
    Total σ⁰_db = coherent sum of both paths (they are reciprocal)

The water surface component uses the I2EM bistatic model at the specular
geometry. The wall component uses Fresnel reflection coefficients for the
wall material. Wall orientation relative to the radar look direction
modulates the effective cross-section.

Reference:
    - Ulaby & Long (2014), Microwave Radar and Radiometric Remote Sensing
    - Hong & Wdowinski (2014), Double-bounce component in cross-pol SAR
    - Ferro et al. (2011), Flood mapping with SAR: dihedral geometry

Units: GHz, meters, degrees throughout (consistent with pyi2em).
"""

import numpy as np

try:
    import pyi2em
except ImportError:
    pyi2em = None


# Common wall materials: complex permittivity at X-band (~9.65 GHz)
WALL_MATERIALS = {
    "concrete": complex(6.0, 0.3),      # dry concrete
    "wet_concrete": complex(9.0, 1.5),   # wet concrete
    "brick": complex(4.5, 0.2),          # dry brick
    "metal": complex(1e6, 1e7),          # metal (effectively PEC)
    "wood": complex(2.0, 0.4),           # dry wood
    "glass": complex(6.5, 0.1),          # glass
}

# Default water roughness parameters for calm/light-wind conditions
WATER_ROUGHNESS = {
    "calm": {"rms_height_m": 0.001, "corr_length_m": 0.05},
    "light_wind": {"rms_height_m": 0.003, "corr_length_m": 0.08},
    "moderate_wind": {"rms_height_m": 0.008, "corr_length_m": 0.12},
    "rough": {"rms_height_m": 0.015, "corr_length_m": 0.15},
}


def fresnel_reflection(theta_i_rad, er):
    """
    Fresnel reflection coefficients for a planar dielectric boundary.

    Parameters
    ----------
    theta_i_rad : float or ndarray
        Incidence angle in radians.
    er : complex
        Complex relative permittivity of the reflecting medium.

    Returns
    -------
    R_h, R_v : complex
        Horizontal and vertical Fresnel reflection coefficients.
    """
    cos_i = np.cos(theta_i_rad)
    sin_i = np.sin(theta_i_rad)
    sqrt_term = np.sqrt(er - sin_i**2)

    R_h = (cos_i - sqrt_term) / (cos_i + sqrt_term)
    R_v = (er * cos_i - sqrt_term) / (er * cos_i + sqrt_term)
    return R_h, R_v


def water_permittivity(freq_ghz, temperature_c=20.0, salinity_ppt=0.0):
    """
    Debye model for water permittivity.

    Simple single-relaxation Debye model for fresh/brackish water.

    Parameters
    ----------
    freq_ghz : float
        Frequency in GHz.
    temperature_c : float
        Water temperature in Celsius.
    salinity_ppt : float
        Salinity in parts per thousand (0 for fresh water).

    Returns
    -------
    er : complex
        Complex relative permittivity of water.
    """
    # Static permittivity (Stogryn, 1971 approx)
    eps_s = 87.74 - 0.4008 * temperature_c + 9.398e-4 * temperature_c**2
    eps_inf = 4.9  # high-frequency limit
    # Relaxation frequency (GHz)
    tau_ns = 1.0 / (2 * np.pi * (17.0 - 0.07 * temperature_c))  # ns
    f_rel_ghz = 1.0 / (2 * np.pi * tau_ns)

    # Ionic conductivity from salinity (S/m)
    sigma = 0.18 * salinity_ppt  # rough linear approximation

    eps_real = eps_inf + (eps_s - eps_inf) / (1 + (freq_ghz / f_rel_ghz) ** 2)
    eps_imag = ((eps_s - eps_inf) * (freq_ghz / f_rel_ghz)
                / (1 + (freq_ghz / f_rel_ghz) ** 2)
                + sigma / (2 * np.pi * freq_ghz * 1e9 * 8.854e-12))

    return complex(eps_real, eps_imag)


def _orientation_factor(wall_azimuth_deg, look_azimuth_deg):
    """
    Compute the orientation modulation factor for double-bounce.

    The double-bounce return is maximized when the wall normal is parallel
    to the radar look direction (wall face perpendicular to look). It falls
    off as cos²(delta_azimuth) because both the incident and reflected
    projections onto the wall are reduced.

    Parameters
    ----------
    wall_azimuth_deg : float or ndarray
        Wall normal azimuth in degrees (direction the wall faces).
    look_azimuth_deg : float or ndarray
        Radar look direction azimuth in degrees.

    Returns
    -------
    factor : float or ndarray
        Orientation modulation factor in [0, 1].
        1.0 when wall faces the radar, 0.0 when wall is parallel to look.
    """
    delta = np.deg2rad(wall_azimuth_deg - look_azimuth_deg)
    # cos^2 because both legs of the double-bounce are projected
    factor = np.cos(delta) ** 2
    # Clip: wall facing away from radar has no double-bounce
    return np.clip(factor, 0.0, 1.0)


def sigma0_double_bounce(
    freq_ghz=9.65,
    theta_deg=30.0,
    wall_material="concrete",
    wall_er=None,
    water_er=None,
    water_temperature_c=20.0,
    water_salinity_ppt=0.0,
    rms_height_m=0.003,
    corr_length_m=0.08,
    water_condition=None,
    correl="exponential",
    wall_azimuth_deg=0.0,
    look_azimuth_deg=0.0,
    return_db=True,
    include_components=False,
):
    """
    Double-bounce backscatter from water-wall dihedral.

    Computes σ⁰ for the double-bounce (dihedral) scattering path:
    radar → rough water → building wall → radar.

    Parameters
    ----------
    freq_ghz : float
        Radar frequency in GHz. Default 9.65 (ICEYE X-band).
    theta_deg : float or array_like
        Incidence angle in degrees.
    wall_material : str
        Wall material name (key in WALL_MATERIALS). Ignored if wall_er given.
    wall_er : complex, optional
        Complex permittivity of wall. Overrides wall_material.
    water_er : complex, optional
        Complex permittivity of water. If None, computed from Debye model.
    water_temperature_c : float
        Water temperature for Debye model.
    water_salinity_ppt : float
        Salinity for Debye model.
    rms_height_m : float
        Water surface RMS roughness height in meters.
    corr_length_m : float
        Water surface correlation length in meters.
    water_condition : str, optional
        Key in WATER_ROUGHNESS to set rms_height_m and corr_length_m.
    correl : str
        Correlation function ('exponential', 'gaussian', etc.).
    wall_azimuth_deg : float or array_like
        Wall normal azimuth in degrees.
    look_azimuth_deg : float or array_like
        Radar look direction azimuth in degrees.
    return_db : bool
        If True, return σ⁰ in dB. Otherwise linear.
    include_components : bool
        If True, also return water and wall reflection components.

    Returns
    -------
    dict
        Keys: 'hh', 'vv' (and optionally 'water_hh', 'water_vv',
        'wall_Rh', 'wall_Rv', 'orientation_factor').
        Values in dB (default) or linear.
    """
    if pyi2em is None:
        raise ImportError("pyi2em C extension required for I2EM water surface model")

    theta_deg = np.atleast_1d(np.asarray(theta_deg, dtype=float))

    # Water permittivity
    if water_er is None:
        water_er = water_permittivity(freq_ghz, water_temperature_c, water_salinity_ppt)

    # Wall permittivity
    if wall_er is None:
        if wall_material not in WALL_MATERIALS:
            raise ValueError(
                f"Unknown wall_material '{wall_material}'. "
                f"Options: {list(WALL_MATERIALS.keys())}"
            )
        wall_er = WALL_MATERIALS[wall_material]

    # Water roughness from condition preset
    if water_condition is not None:
        if water_condition not in WATER_ROUGHNESS:
            raise ValueError(
                f"Unknown water_condition '{water_condition}'. "
                f"Options: {list(WATER_ROUGHNESS.keys())}"
            )
        rms_height_m = WATER_ROUGHNESS[water_condition]["rms_height_m"]
        corr_length_m = WATER_ROUGHNESS[water_condition]["corr_length_m"]

    # Orientation modulation
    orient_factor = _orientation_factor(wall_azimuth_deg, look_azimuth_deg)

    # Compute per-angle
    hh_out = np.zeros_like(theta_deg)
    vv_out = np.zeros_like(theta_deg)
    water_hh_arr = np.zeros_like(theta_deg)
    water_vv_arr = np.zeros_like(theta_deg)

    for i, th in enumerate(theta_deg):
        if th < 0.1 or th > 89.9:
            hh_out[i] = -999.0 if return_db else 0.0
            vv_out[i] = -999.0 if return_db else 0.0
            continue

        # --- Water surface: I2EM bistatic at specular geometry ---
        # Double-bounce geometry: the water surface scatters from
        # incidence angle θ to the complementary angle (90° - θ) for
        # the wall reflection. For a vertical wall, the specular path
        # is: incident at θ → forward scatter at θ (specular off water)
        # → wall reflects back at θ.
        #
        # The key bistatic angle is ths = thi (specular forward scatter
        # off water towards the wall), with phs = 0° (forward direction).
        thi = th
        ths = th  # specular forward scatter
        phs = 0.0  # forward scatter azimuth

        water_sig = pyi2em.sigma0_bistatic(
            freq_ghz=freq_ghz,
            rms_height_m=rms_height_m,
            corr_length_m=corr_length_m,
            thi_deg=thi,
            ths_deg=ths,
            phs_deg=phs,
            er_complex=water_er,
            correl=correl,
            xcoeff=1.0,
            return_db=False,  # need linear for multiplication
        )

        water_hh_lin = water_sig["hh"]
        water_vv_lin = water_sig["vv"]

        # --- Wall reflection: Fresnel at complementary angle ---
        # After scattering off water at angle θ, the signal hits the
        # vertical wall at angle (90° - θ) from the wall normal.
        wall_theta_rad = np.deg2rad(90.0 - th)
        R_h, R_v = fresnel_reflection(wall_theta_rad, wall_er)

        wall_Rh2 = np.abs(R_h) ** 2  # power reflection coefficient
        wall_Rv2 = np.abs(R_v) ** 2

        # --- Double-bounce σ⁰ ---
        # σ⁰_db = σ⁰_water(bistatic) × |R_wall|² × orientation_factor
        # Factor of 2 accounts for both paths (water→wall and wall→water)
        # which are coherent and reciprocal.
        sigma_hh_lin = 2.0 * water_hh_lin * wall_Rh2 * orient_factor
        sigma_vv_lin = 2.0 * water_vv_lin * wall_Rv2 * orient_factor

        if return_db:
            hh_out[i] = 10.0 * np.log10(max(sigma_hh_lin, 1e-30))
            vv_out[i] = 10.0 * np.log10(max(sigma_vv_lin, 1e-30))
        else:
            hh_out[i] = sigma_hh_lin
            vv_out[i] = sigma_vv_lin

        water_hh_arr[i] = water_hh_lin
        water_vv_arr[i] = water_vv_lin

    # Squeeze scalar input back to scalar output
    if hh_out.size == 1:
        hh_out = float(hh_out[0])
        vv_out = float(vv_out[0])

    result = {"hh": hh_out, "vv": vv_out}

    if include_components:
        result["water_hh_linear"] = water_hh_arr if water_hh_arr.size > 1 else float(water_hh_arr[0])
        result["water_vv_linear"] = water_vv_arr if water_vv_arr.size > 1 else float(water_vv_arr[0])
        result["wall_er"] = wall_er
        result["water_er"] = water_er
        result["orientation_factor"] = orient_factor

    return result


def sigma0_double_bounce_orientation_curve(
    freq_ghz=9.65,
    theta_deg=30.0,
    wall_material="concrete",
    wall_er=None,
    water_er=None,
    rms_height_m=0.003,
    corr_length_m=0.08,
    water_condition=None,
    correl="exponential",
    look_azimuth_deg=0.0,
    wall_azimuths_deg=None,
    return_db=True,
):
    """
    Compute double-bounce σ⁰ as a function of wall orientation.

    Sweeps wall_azimuth_deg from 0° to 360° (or user-specified values)
    and returns the orientation-dependent backscatter curve.

    Parameters
    ----------
    wall_azimuths_deg : array_like, optional
        Wall normal azimuths to evaluate. Default: 0° to 360° in 1° steps.
    (other parameters as in sigma0_double_bounce)

    Returns
    -------
    dict
        'wall_azimuths_deg': ndarray, 'hh': ndarray, 'vv': ndarray
    """
    if wall_azimuths_deg is None:
        wall_azimuths_deg = np.arange(0, 361, 1.0)
    else:
        wall_azimuths_deg = np.atleast_1d(np.asarray(wall_azimuths_deg, dtype=float))

    hh_curve = np.zeros(len(wall_azimuths_deg))
    vv_curve = np.zeros(len(wall_azimuths_deg))

    for i, waz in enumerate(wall_azimuths_deg):
        result = sigma0_double_bounce(
            freq_ghz=freq_ghz,
            theta_deg=theta_deg,
            wall_material=wall_material,
            wall_er=wall_er,
            water_er=water_er,
            rms_height_m=rms_height_m,
            corr_length_m=corr_length_m,
            water_condition=water_condition,
            correl=correl,
            wall_azimuth_deg=waz,
            look_azimuth_deg=look_azimuth_deg,
            return_db=return_db,
        )
        hh_curve[i] = result["hh"]
        vv_curve[i] = result["vv"]

    return {
        "wall_azimuths_deg": wall_azimuths_deg,
        "hh": hh_curve,
        "vv": vv_curve,
    }


def sigma0_double_bounce_incidence_curve(
    freq_ghz=9.65,
    theta_min_deg=15.0,
    theta_max_deg=60.0,
    theta_step_deg=1.0,
    wall_material="concrete",
    wall_er=None,
    water_er=None,
    rms_height_m=0.003,
    corr_length_m=0.08,
    water_condition=None,
    correl="exponential",
    wall_azimuth_deg=0.0,
    look_azimuth_deg=0.0,
    return_db=True,
):
    """
    Compute double-bounce σ⁰ as a function of incidence angle.

    Parameters
    ----------
    theta_min_deg, theta_max_deg, theta_step_deg : float
        Incidence angle sweep range.
    (other parameters as in sigma0_double_bounce)

    Returns
    -------
    dict
        'theta_deg': ndarray, 'hh': ndarray, 'vv': ndarray
    """
    thetas = np.arange(theta_min_deg, theta_max_deg + theta_step_deg / 2, theta_step_deg)

    result = sigma0_double_bounce(
        freq_ghz=freq_ghz,
        theta_deg=thetas,
        wall_material=wall_material,
        wall_er=wall_er,
        water_er=water_er,
        rms_height_m=rms_height_m,
        corr_length_m=corr_length_m,
        water_condition=water_condition,
        correl=correl,
        wall_azimuth_deg=wall_azimuth_deg,
        look_azimuth_deg=look_azimuth_deg,
        return_db=return_db,
    )

    return {
        "theta_deg": thetas,
        "hh": result["hh"],
        "vv": result["vv"],
    }
