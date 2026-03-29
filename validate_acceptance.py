#!/usr/bin/env python3
"""
Acceptance validation for P3-I2EM:
  "I2EM forward model produces orientation-dependent backscatter predictions"

Checks:
  1. I2EM forward model runs (water surface via C extension + wall Fresnel)
  2. Wall materials: concrete (eps~6) and metal (eps~inf) produce valid output
  3. Orientation-dependent: wall azimuth modulates backscatter as cos^2

Exit 0 on success, 1 on failure.
"""

import sys
import numpy as np
from double_bounce import (
    sigma0_double_bounce,
    sigma0_double_bounce_orientation_curve,
    WALL_MATERIALS,
)


def check_forward_model():
    """Check 1: I2EM forward model produces backscatter predictions."""
    r = sigma0_double_bounce(
        freq_ghz=9.65, theta_deg=30.0, wall_material="concrete",
        return_db=True, include_components=True,
    )
    assert np.isfinite(r["hh"]) and np.isfinite(r["vv"]), "Non-finite output"
    assert r["water_hh_linear"] > 0, "Water I2EM component must be positive"
    assert r["water_vv_linear"] > 0, "Water I2EM component must be positive"
    print(f"  Forward model: HH={r['hh']:.2f} dB, VV={r['vv']:.2f} dB")
    return True


def check_wall_materials():
    """Check 2: Concrete (eps~6) and metal (eps~inf) produce valid, distinct output."""
    sig_c = sigma0_double_bounce(freq_ghz=9.65, theta_deg=30.0, wall_material="concrete", return_db=True)
    sig_m = sigma0_double_bounce(freq_ghz=9.65, theta_deg=30.0, wall_material="metal", return_db=True)

    assert np.isfinite(sig_c["hh"]) and np.isfinite(sig_m["hh"])
    # Metal (|R|~1) must be stronger than concrete
    assert sig_m["hh"] > sig_c["hh"], f"Metal ({sig_m['hh']:.1f}) should exceed concrete ({sig_c['hh']:.1f})"
    # Verify eps values
    assert np.isclose(WALL_MATERIALS["concrete"].real, 6.0, atol=0.5), "Concrete eps~6"
    assert WALL_MATERIALS["metal"].real > 1e4, "Metal eps~inf"
    print(f"  Concrete: HH={sig_c['hh']:.2f} dB | Metal: HH={sig_m['hh']:.2f} dB (delta={sig_m['hh']-sig_c['hh']:.1f} dB)")
    return True


def check_orientation_dependence():
    """Check 3: Wall orientation modulates backscatter (cos^2 pattern)."""
    curve = sigma0_double_bounce_orientation_curve(
        freq_ghz=9.65, theta_deg=30.0, wall_material="concrete",
        look_azimuth_deg=0.0, return_db=False,
    )
    hh = curve["hh"]
    az = curve["wall_azimuths_deg"]

    # Max at 0 and 180, zero at 90 and 270
    idx_0 = np.argmin(np.abs(az - 0))
    idx_90 = np.argmin(np.abs(az - 90))
    idx_180 = np.argmin(np.abs(az - 180))
    idx_270 = np.argmin(np.abs(az - 270))

    assert hh[idx_0] > 0, "Must have positive backscatter at 0 deg"
    assert np.isclose(hh[idx_90], 0.0, atol=1e-15), "Must be zero at 90 deg"
    assert np.isclose(hh[idx_0], hh[idx_180], rtol=1e-6), "Symmetric at 0/180"
    assert np.isclose(hh[idx_90], hh[idx_270], atol=1e-15), "Symmetric at 90/270"

    # Monotonic decrease 0->90
    for i in range(idx_0, idx_90):
        assert hh[i] >= hh[i + 1], f"Not monotonically decreasing at az={az[i]}"

    peak_db = 10 * np.log10(max(hh[idx_0], 1e-30))
    print(f"  Orientation: peak={peak_db:.2f} dB at 0/180 deg, null at 90/270 deg")
    print(f"  Dynamic range: {peak_db - 10*np.log10(max(hh[idx_0+45], 1e-30)):.1f} dB over 45 deg rotation")
    return True


if __name__ == "__main__":
    checks = [
        ("I2EM forward model", check_forward_model),
        ("Wall materials (concrete/metal)", check_wall_materials),
        ("Orientation-dependent predictions", check_orientation_dependence),
    ]
    passed = 0
    for name, fn in checks:
        try:
            fn()
            print(f"  PASS: {name}")
            passed += 1
        except Exception as e:
            print(f"  FAIL: {name} — {e}")

    print(f"\n{passed}/{len(checks)} acceptance checks passed")
    sys.exit(0 if passed == len(checks) else 1)
