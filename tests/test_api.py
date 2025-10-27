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
