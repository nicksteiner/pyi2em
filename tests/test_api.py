import numpy as np
import pyi2em as i2em


def test_emissivity_runs():
    eh, ev = i2em.emissivity(1.26, 0.02, 0.1, 40.0, 5+1j, "exponential")
    assert 0 <= eh <= 1 and 0 <= ev <= 1


def test_sigma0_backscatter_runs():
    out = i2em.sigma0_backscatter(1.26, 0.02, 0.1, [20, 30, 40], 5+1j)
    assert set(out.keys()) == {"hh", "vv", "hv"}
    assert len(out["vv"]) == 3


def test_sigma0_bistatic_runs():
    out = i2em.sigma0_bistatic(1.26, 0.02, 0.1, 30, 30, 180, 5+1j)
    assert "vv" in out and "hh" in out


def test_backscatter_vs_bistatic_equivalence():
    # Backscatter geometry (ths=thi, phs=180) should match sigma0_backscatter for HH/VV
    freq = 5.0
    rms, L = 0.01, 0.05
    theta = 40.0
    er = 15 + 3j

    back = i2em.sigma0_backscatter(freq, rms, L, theta, er, correl="gaussian", include_hv=True, return_db=True)
    bi = i2em.sigma0_bistatic(freq, rms, L, theta, theta, 180.0, er, correl="gaussian", xcoeff=1.0, return_db=True)

    assert np.isfinite(back["hh"]).all() and np.isfinite(back["vv"]).all()
    assert np.isfinite(bi["hh"]) and np.isfinite(bi["vv"]) 
    # back arrays may be 0-d arrays; convert to float
    hh_back = float(np.asarray(back["hh"]))
    vv_back = float(np.asarray(back["vv"]))
    assert np.isclose(hh_back, bi["hh"], atol=1e-6)
    assert np.isclose(vv_back, bi["vv"], atol=1e-6)


def test_theta_array_and_hv_toggle():
    freq = 5.0
    rms, L = 0.01, 0.05
    thetas = np.array([20.0, 30.0, 40.0])
    er = 11.3 + 1.5j

    out1 = i2em.sigma0_backscatter(freq, rms, L, thetas, er, correl="exponential", include_hv=True, return_db=True)
    out2 = i2em.sigma0_backscatter(freq, rms, L, thetas, er, correl="exponential", include_hv=False, return_db=True)

    # Shapes
    assert out1["hh"].shape == (3,) and out1["vv"].shape == (3,) and out1["hv"].shape == (3,)
    assert out2["hh"].shape == (3,) and out2["vv"].shape == (3,) and out2["hv"].shape == (3,)
    # HV present when requested, zero when disabled
    assert np.all(np.isfinite(out1["hv"]))
    assert np.all(out2["hv"] == 0.0)


def test_return_db_linear_conversion():
    freq = 5.0
    rms, L = 0.01, 0.05
    thetas = [25.0, 35.0]
    er = 12 + 2j

    out_db = i2em.sigma0_backscatter(freq, rms, L, thetas, er, correl="gaussian", include_hv=True, return_db=True)
    out_lin = i2em.sigma0_backscatter(freq, rms, L, thetas, er, correl="gaussian", include_hv=True, return_db=False)

    def to_db(x):
        return 10.0 * np.log10(np.maximum(np.asarray(x), 1e-300))

    for pol in ("hh", "vv", "hv"):
        assert np.allclose(to_db(out_lin[pol]), np.asarray(out_db[pol]), atol=1e-6)


def test_input_validation_errors():
    # Invalid frequency
    with np.testing.assert_raises(Exception):
        i2em.sigma0_backscatter(0.0, 0.01, 0.05, 30.0, 5+1j)
    # Invalid corr length
    with np.testing.assert_raises(Exception):
        i2em.emissivity(5.0, 0.01, 0.0, 30.0, 5+1j)
    # Invalid theta range
    with np.testing.assert_raises(Exception):
        i2em.sigma0_bistatic(5.0, 0.01, 0.05, 95.0, 40.0, 180.0, 5+1j)
    with np.testing.assert_raises(Exception):
        i2em.sigma0_bistatic(5.0, 0.01, 0.05, 40.0, 95.0, 180.0, 5+1j)


def test_correlation_options_and_unknown():
    freq = 3.0
    rms, L = 0.005, 0.05
    theta = 35.0
    er = 12 + 2.5j
    for name in ("exponential", "gaussian", "x_power", "x_exponential"):
        eh, ev = i2em.emissivity(freq, rms, L, theta, er, correl=name)
        assert 0 <= eh <= 1 and 0 <= ev <= 1
    # Unknown correlation must error
    with np.testing.assert_raises(Exception):
        i2em.emissivity(freq, rms, L, theta, er, correl="unknown")
