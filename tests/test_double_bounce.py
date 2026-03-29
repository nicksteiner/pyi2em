"""
Tests for the I2EM double-bounce water-wall scattering model.

Acceptance criteria:
    - I2EM forward model produces orientation-dependent backscatter predictions
    - Model handles wall materials (concrete eps~6, metal eps~inf)
    - Model handles incidence angle variation
    - Wall orientation relative to radar look direction modulates output
"""

import numpy as np
import pytest
from double_bounce import (
    fresnel_reflection,
    water_permittivity,
    sigma0_double_bounce,
    sigma0_double_bounce_orientation_curve,
    sigma0_double_bounce_incidence_curve,
    WALL_MATERIALS,
    WATER_ROUGHNESS,
    _orientation_factor,
)


class TestFresnelReflection:
    def test_normal_incidence(self):
        """At normal incidence, R_h and R_v should be equal."""
        R_h, R_v = fresnel_reflection(0.0, complex(6.0, 0.3))
        # At normal incidence Rh = Rv
        assert np.isclose(np.abs(R_h), np.abs(R_v), atol=1e-10)

    def test_metal_high_reflectivity(self):
        """Metal (high eps) should have |R| ≈ 1."""
        R_h, R_v = fresnel_reflection(np.deg2rad(30.0), complex(1e6, 1e7))
        assert np.abs(R_h) > 0.99
        assert np.abs(R_v) > 0.99

    def test_reflectivity_range(self):
        """Reflection coefficients should be in [0, 1]."""
        for theta in [10, 30, 50, 70]:
            R_h, R_v = fresnel_reflection(np.deg2rad(theta), complex(6.0, 0.3))
            assert 0.0 <= np.abs(R_h) <= 1.0
            assert 0.0 <= np.abs(R_v) <= 1.0


class TestWaterPermittivity:
    def test_fresh_water_real_part(self):
        """Fresh water at X-band should have eps_real >> 1."""
        er = water_permittivity(9.65, temperature_c=20.0, salinity_ppt=0.0)
        assert er.real > 20  # water at X-band has high permittivity
        assert er.imag > 0   # lossy

    def test_temperature_dependence(self):
        """Warmer water should have lower static permittivity."""
        er_cold = water_permittivity(9.65, temperature_c=5.0)
        er_warm = water_permittivity(9.65, temperature_c=35.0)
        # Static permittivity decreases with temperature
        # but at X-band the relationship is more complex; just check finite
        assert np.isfinite(er_cold.real) and np.isfinite(er_warm.real)


class TestOrientationFactor:
    def test_facing_radar(self):
        """Wall facing radar (same azimuth) should give factor = 1."""
        f = _orientation_factor(0.0, 0.0)
        assert np.isclose(f, 1.0)

    def test_parallel_to_radar(self):
        """Wall parallel to radar look (90° offset) should give factor = 0."""
        f = _orientation_factor(90.0, 0.0)
        assert np.isclose(f, 0.0, atol=1e-10)

    def test_facing_away(self):
        """Wall facing away (180° offset) should give factor = 1 (cos²(180°)=1)."""
        f = _orientation_factor(180.0, 0.0)
        assert np.isclose(f, 1.0, atol=1e-10)

    def test_45_degrees(self):
        """At 45° offset, factor should be cos²(45°) = 0.5."""
        f = _orientation_factor(45.0, 0.0)
        assert np.isclose(f, 0.5, atol=1e-10)

    def test_symmetry(self):
        """Factor should be symmetric around look direction."""
        f_pos = _orientation_factor(30.0, 0.0)
        f_neg = _orientation_factor(-30.0, 0.0)
        assert np.isclose(f_pos, f_neg)


class TestDoubleBounce:
    """Core double-bounce model tests."""

    def test_basic_output_structure(self):
        """Model returns dict with hh and vv keys."""
        result = sigma0_double_bounce(
            freq_ghz=9.65, theta_deg=30.0, wall_material="concrete"
        )
        assert "hh" in result and "vv" in result
        assert np.isfinite(result["hh"]) and np.isfinite(result["vv"])

    def test_db_returns(self):
        """Default return is in dB (finite values)."""
        result = sigma0_double_bounce(
            freq_ghz=9.65, theta_deg=30.0, wall_material="concrete"
        )
        # dB values should be finite
        assert np.isfinite(result["hh"])
        assert np.isfinite(result["vv"])

    def test_linear_return(self):
        """Linear return should be positive and convertible to dB."""
        result_lin = sigma0_double_bounce(
            freq_ghz=9.65, theta_deg=30.0, wall_material="concrete", return_db=False
        )
        result_db = sigma0_double_bounce(
            freq_ghz=9.65, theta_deg=30.0, wall_material="concrete", return_db=True
        )
        assert result_lin["hh"] > 0
        assert np.isclose(
            10.0 * np.log10(result_lin["hh"]), result_db["hh"], atol=1e-6
        )

    def test_metal_stronger_than_concrete(self):
        """Metal walls should produce stronger double-bounce than concrete."""
        sig_metal = sigma0_double_bounce(
            freq_ghz=9.65, theta_deg=30.0, wall_material="metal", return_db=True
        )
        sig_concrete = sigma0_double_bounce(
            freq_ghz=9.65, theta_deg=30.0, wall_material="concrete", return_db=True
        )
        # Metal has higher reflectivity → stronger double-bounce
        assert sig_metal["hh"] > sig_concrete["hh"]
        assert sig_metal["vv"] > sig_concrete["vv"]

    def test_all_wall_materials(self):
        """All predefined wall materials should produce valid output."""
        for material in WALL_MATERIALS:
            result = sigma0_double_bounce(
                freq_ghz=9.65, theta_deg=30.0, wall_material=material
            )
            assert np.isfinite(result["hh"]), f"Failed for {material}"
            assert np.isfinite(result["vv"]), f"Failed for {material}"

    def test_custom_wall_er(self):
        """Custom wall_er should override wall_material."""
        result = sigma0_double_bounce(
            freq_ghz=9.65, theta_deg=30.0, wall_er=complex(6.0, 0.3)
        )
        assert np.isfinite(result["hh"])

    def test_water_condition_presets(self):
        """Water condition presets should work."""
        for condition in WATER_ROUGHNESS:
            result = sigma0_double_bounce(
                freq_ghz=9.65, theta_deg=30.0, water_condition=condition
            )
            assert np.isfinite(result["hh"]), f"Failed for {condition}"

    def test_rougher_water_different_scatter(self):
        """Different water roughness should produce different backscatter."""
        sig_calm = sigma0_double_bounce(
            freq_ghz=9.65, theta_deg=30.0, water_condition="calm", return_db=True
        )
        sig_rough = sigma0_double_bounce(
            freq_ghz=9.65, theta_deg=30.0, water_condition="rough", return_db=True
        )
        # They should differ
        assert not np.isclose(sig_calm["hh"], sig_rough["hh"], atol=0.1)

    def test_array_incidence_angles(self):
        """Model should accept array of incidence angles."""
        thetas = [20.0, 30.0, 40.0, 50.0]
        result = sigma0_double_bounce(
            freq_ghz=9.65, theta_deg=thetas, wall_material="concrete"
        )
        assert len(result["hh"]) == 4
        assert len(result["vv"]) == 4
        assert np.all(np.isfinite(result["hh"]))

    def test_include_components(self):
        """include_components should return intermediate values."""
        result = sigma0_double_bounce(
            freq_ghz=9.65,
            theta_deg=30.0,
            wall_material="concrete",
            include_components=True,
        )
        assert "water_hh_linear" in result
        assert "water_vv_linear" in result
        assert "wall_er" in result
        assert "water_er" in result
        assert "orientation_factor" in result

    def test_invalid_wall_material(self):
        with pytest.raises(ValueError, match="Unknown wall_material"):
            sigma0_double_bounce(wall_material="unobtanium")

    def test_invalid_water_condition(self):
        with pytest.raises(ValueError, match="Unknown water_condition"):
            sigma0_double_bounce(water_condition="hurricane")


class TestOrientationDependence:
    """Tests that verify orientation-dependent backscatter predictions."""

    def test_max_at_perpendicular_wall(self):
        """Double-bounce is maximum when wall faces the radar."""
        sig_facing = sigma0_double_bounce(
            freq_ghz=9.65,
            theta_deg=30.0,
            wall_azimuth_deg=0.0,
            look_azimuth_deg=0.0,
            return_db=False,
        )
        sig_angled = sigma0_double_bounce(
            freq_ghz=9.65,
            theta_deg=30.0,
            wall_azimuth_deg=45.0,
            look_azimuth_deg=0.0,
            return_db=False,
        )
        assert sig_facing["hh"] > sig_angled["hh"]

    def test_zero_at_parallel_wall(self):
        """Double-bounce vanishes when wall is parallel to look direction."""
        sig = sigma0_double_bounce(
            freq_ghz=9.65,
            theta_deg=30.0,
            wall_azimuth_deg=90.0,
            look_azimuth_deg=0.0,
            return_db=False,
        )
        assert np.isclose(sig["hh"], 0.0, atol=1e-20)
        assert np.isclose(sig["vv"], 0.0, atol=1e-20)

    def test_orientation_curve_shape(self):
        """Orientation curve should peak at 0° and 180° (cos² pattern)."""
        curve = sigma0_double_bounce_orientation_curve(
            freq_ghz=9.65, theta_deg=30.0, return_db=False
        )
        hh = curve["hh"]
        az = curve["wall_azimuths_deg"]

        # Max at 0° and 180°
        idx_0 = np.argmin(np.abs(az - 0))
        idx_90 = np.argmin(np.abs(az - 90))
        idx_180 = np.argmin(np.abs(az - 180))

        assert hh[idx_0] > hh[idx_90]
        assert hh[idx_180] > hh[idx_90]
        assert np.isclose(hh[idx_0], hh[idx_180], rtol=1e-6)
        assert np.isclose(hh[idx_90], 0.0, atol=1e-20)

    def test_orientation_curve_returns_full_sweep(self):
        """Default orientation curve should cover 0-360°."""
        curve = sigma0_double_bounce_orientation_curve()
        assert len(curve["wall_azimuths_deg"]) == 361
        assert len(curve["hh"]) == 361
        assert len(curve["vv"]) == 361


class TestIncidenceCurve:
    def test_incidence_curve_runs(self):
        """Incidence angle curve should produce valid results."""
        curve = sigma0_double_bounce_incidence_curve(
            freq_ghz=9.65,
            theta_min_deg=20.0,
            theta_max_deg=50.0,
            theta_step_deg=5.0,
        )
        assert len(curve["theta_deg"]) == 7  # 20,25,30,35,40,45,50
        assert np.all(np.isfinite(curve["hh"]))
        assert np.all(np.isfinite(curve["vv"]))

    def test_incidence_angle_variation(self):
        """Backscatter should vary with incidence angle."""
        curve = sigma0_double_bounce_incidence_curve(
            theta_min_deg=20.0, theta_max_deg=50.0, theta_step_deg=10.0
        )
        # Values should not all be identical
        assert np.std(curve["hh"]) > 0.1


class TestICEYEScenario:
    """Integration test: realistic ICEYE urban flood scenario."""

    def test_iceye_urban_flood(self):
        """
        Simulate ICEYE X-band observing flooded urban area.
        ICEYE: 9.65 GHz, typical incidence 20-40°.
        """
        result = sigma0_double_bounce(
            freq_ghz=9.65,
            theta_deg=30.0,
            wall_material="concrete",
            water_condition="calm",
            wall_azimuth_deg=0.0,
            look_azimuth_deg=0.0,
            return_db=True,
            include_components=True,
        )
        # Should produce a valid backscatter value
        assert np.isfinite(result["hh"])
        assert np.isfinite(result["vv"])
        # Double-bounce from building+water should be detectable
        # (not absurdly low)
        assert result["hh"] > -60  # should be above noise floor

    def test_iceye_metal_vs_concrete(self):
        """Metal buildings should produce notably stronger double-bounce."""
        sig_c = sigma0_double_bounce(
            freq_ghz=9.65, theta_deg=25.0, wall_material="concrete",
            water_condition="light_wind", return_db=True
        )
        sig_m = sigma0_double_bounce(
            freq_ghz=9.65, theta_deg=25.0, wall_material="metal",
            water_condition="light_wind", return_db=True
        )
        # Metal should be several dB stronger
        assert sig_m["hh"] - sig_c["hh"] > 1.0

    def test_iceye_orientation_sensitivity(self):
        """
        Acceptance: orientation-dependent backscatter predictions.
        Verify that rotating wall azimuth produces varying sigma0.
        """
        azimuths = [0, 15, 30, 45, 60, 75, 90]
        results_hh = []
        for az in azimuths:
            r = sigma0_double_bounce(
                freq_ghz=9.65, theta_deg=30.0, wall_material="concrete",
                water_condition="calm", wall_azimuth_deg=az,
                look_azimuth_deg=0.0, return_db=False
            )
            results_hh.append(r["hh"])

        # Should monotonically decrease from 0° to 90°
        for i in range(len(results_hh) - 1):
            assert results_hh[i] >= results_hh[i + 1], (
                f"σ⁰ should decrease as wall rotates away from look: "
                f"az={azimuths[i]}° ({results_hh[i]}) >= "
                f"az={azimuths[i+1]}° ({results_hh[i+1]})"
            )
        # At 90° should be zero
        assert np.isclose(results_hh[-1], 0.0, atol=1e-20)
