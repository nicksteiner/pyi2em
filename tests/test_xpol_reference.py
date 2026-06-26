"""Cross-pol (VH/HV) validation against Fung's published I2EM reference curves.

Reference: Fung, "Microwave Scattering and Emission Models for Users"
(Artech House, 2009), Fig. 3.9 ("epsilon variations with exponential
correlation for cross-pol"). The companion Mathematica notebook
(fung-037: "Fig3.9 cross pol dielectric variation.nb") defines
IEMX[fr, sig, corr, theta, er, sp, xx] -- the same multiple-scattering
cross-pol surface model implemented here as IEMX_model() -- and caches its
evaluated output in a Print cell.

Fig. 3.9 parameters, transcribed and verified against the notebook source:
    f   = 5          (GHz)   -- "Step 2" parameter cell
    sig = 0.3        (cm)    -- "
    L   = 4.0        (cm)    -- "
    sp  = 1                  -- Switch[sp, 1 -> sig/L, ...] == exponential
    xx  = 1.5               (unused for sp=1)
    er  = 3.5, 7.5, 15  (real, no loss)  -- Step 3 Table calls (vh1/vh2/vh3)
    theta = 0.001 .. 70.001 deg, 10 deg steps (astart/aend/astep)
Unit convention confirmed from the notebook: k = 6.283*fr/30 with
"fr in GHz and length in cm". pyi2em takes SI, so sig/L are passed as
0.003 m and 0.04 m. The sp->correlation map (1=exponential, 2=gaussian,
3=x-power) matches the C++ CORREL_* codes exactly.

The REFERENCE_DB values below are the full-precision numbers cached in the
notebook's Print cell (the evaluated {vh1, vh2, vh3}).

OBSERVED DISCREPANCY (facts only; cause not yet diagnosed):
pyi2em's hv does NOT reproduce the reference. At the verified Fig 3.9
parameters the model sits below the reference by a roughly constant offset
that depends on dielectric and flips sign at grazing angles:
  * er=3.5:  ~ -3.9 dB at theta<=50, crossing to about +4.5 dB at theta=70
  * er=7.5:  ~ -3.4 dB across theta<=60, narrowing to ~ -0.8 dB at theta=70
  * er=15:   ~ -2.2 dB at low theta, widening to ~ -3 dB toward grazing
The offset shrinks as dielectric rises, and only the lowest-dielectric /
highest-angle corner overshoots. Curve shapes track the reference. Cause
not yet diagnosed. This test asserts tight agreement so it FAILS now and
the printed per-point diff table is the artifact documenting the gap;
tighten/keep TOL_DB once the model is reconciled with the reference.
"""

import numpy as np
import pytest
import pyi2em as i2em

# Fig 3.9 fixed parameters (verified against notebook)
FREQ_GHZ = 5.0
RMS_M = 0.003   # sig = 0.3 cm
CORR_M = 0.04   # L = 4.0 cm
CORREL = "exponential"  # sp = 1

THETAS_DEG = [0.001, 10.001, 20.001, 30.001, 40.001,
              50.001, 60.001, 70.001]

# {er_real: [sigma0_vh dB at each theta]} -- full-precision notebook output
# (the cached Print of {vh1, vh2, vh3}).
REFERENCE_DB = {
    3.5: [-42.904266772196785, -43.15528512011496, -43.842952557732644,
          -44.92112906686049, -46.4643103768544, -48.70396720281782,
          -52.12150253340404, -57.72356488895055],
    7.5: [-35.998638199333286, -36.239235622464314, -36.9040603792551,
          -37.95927982058538, -39.48976679092703, -41.7370753350591,
          -45.18858812624008, -50.845825963194414],
    15.0: [-32.17767645993628, -32.4108048879079, -33.05824432310883,
           -34.091066012108485, -35.593083741932745, -37.797522523554704,
           -41.16986252327103, -46.66578594795672],
}

# Tight tolerance: we expect the model to reproduce the published reference.
# Set to a few tenths of a dB so the current ~2-4 dB gap is a hard failure.
TOL_DB = 0.3


def _hv_curve(er_real):
    out = i2em.sigma0_backscatter(
        FREQ_GHZ, RMS_M, CORR_M, THETAS_DEG, complex(er_real, 0.0),
        correl=CORREL, include_hv=True, return_db=True,
    )
    return np.asarray(out["hv"], dtype=float)


def test_xpol_matches_fung_fig39():
    """hv must reproduce Fung Fig. 3.9 to within TOL_DB at every point.

    On failure, prints the full side-by-side reference/pyi2em/diff table so
    the discrepancy is recorded explicitly.
    """
    rows = []
    worst = []  # (abs_diff, message) for points exceeding TOL_DB
    header = f"{'er':>5} {'theta':>6} {'ref_dB':>10} {'pyi2em_dB':>10} {'diff_dB':>9}"
    rows.append(header)
    rows.append("-" * len(header))
    for er_real in sorted(REFERENCE_DB):
        hv = _hv_curve(er_real)
        for theta, ref, got in zip(THETAS_DEG, REFERENCE_DB[er_real], hv):
            diff = got - ref
            flag = "  <-- exceeds tol" if abs(diff) > TOL_DB else ""
            rows.append(
                f"{er_real:5.1f} {theta:6.1f} {ref:10.3f} {got:10.3f} {diff:+9.3f}{flag}"
            )
            if abs(diff) > TOL_DB:
                worst.append(abs(diff))

    table = "\n".join(rows)
    if worst:
        pytest.fail(
            f"pyi2em hv disagrees with Fung Fig. 3.9 (tol={TOL_DB} dB); "
            f"max |diff| = {max(worst):.3f} dB over {len(worst)} point(s).\n"
            f"{table}\n"
            "See module docstring: discrepancy is real, cause not yet diagnosed."
        )


# --- Sanity properties that hold legitimately (these pass) -------------------

@pytest.mark.parametrize("er_real", sorted(REFERENCE_DB))
def test_xpol_finite(er_real):
    """hv is finite across the full incidence-angle sweep."""
    hv = _hv_curve(er_real)
    assert np.all(np.isfinite(hv)), f"non-finite hv for er={er_real}"


@pytest.mark.parametrize("er_real", sorted(REFERENCE_DB))
def test_xpol_monotonic_decrease(er_real):
    """hv decreases monotonically with incidence angle, as in the reference."""
    hv = _hv_curve(er_real)
    diffs = np.diff(hv)
    assert np.all(diffs <= 1e-6), (
        f"er={er_real}: hv not monotonically decreasing: {hv}"
    )


def test_xpol_dielectric_ordering():
    """Higher dielectric -> stronger cross-pol return (3.5 < 7.5 < 15), as in Fig 3.9."""
    h35 = _hv_curve(3.5)
    h75 = _hv_curve(7.5)
    h15 = _hv_curve(15.0)
    for i, theta in enumerate(THETAS_DEG):
        assert h35[i] < h75[i] < h15[i], (
            f"dielectric ordering broken at theta={theta:.0f}: "
            f"{h35[i]:.2f} < {h75[i]:.2f} < {h15[i]:.2f}"
        )
