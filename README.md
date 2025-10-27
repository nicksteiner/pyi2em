![Tests](https://github.com/nicksteiner/pyi2em/actions/workflows/tests.yml/badge.svg?branch=master)

# pyi2em

Python bindings for the I2EM rough-surface scattering and emissivity model.

Reference: Ulaby, F.T. and Long, D.G. (2014). Microwave Radar and Radiometric Remote Sensing. University of Michigan Press.

## Install


```bash
pip install pyi2em
```

Notes on wheels vs source builds

- PyPI wheels bundle required native libraries, including GSL and embedded cubature routines. End users do not need system packages.
- Source builds require GSL available via pkg-config. The cubature integration code is compiled into the extension automatically.

From source

- Ubuntu: `sudo apt-get install -y cmake build-essential libgsl-dev`
- macOS: `brew install cmake gsl`

Then

```bash
pip install .
```

Note: The cubature library (hcubature/pcubature) is fetched automatically during
build; you do not need to install a system package for it. GSL remains a system
dependency for source builds.

## Test
```bash
pytest -q
```

## Usage

The high-level Python API exposes three functions. Units are GHz, meters, and degrees throughout.

1) Emissivity

```python
from pyi2em import emissivity

eh, ev = emissivity(
    freq_ghz=5.3,
    rms_height_m=0.0025,
    corr_length_m=0.10,
    theta_deg=30.0,
    er_complex=complex(11.3, 1.5),
    correl="gaussian",         # 'exponential' | 'gaussian' | 'x_power' | 'x_exponential'
)
print(eh, ev)  # linear emissivities
```

2) Monostatic backscatter (scalar or array incidence)

```python
import numpy as np
from pyi2em import sigma0_backscatter

thetas = np.array([20.0, 30.0, 40.0])
sig = sigma0_backscatter(
    freq_ghz=5.0,
    rms_height_m=0.01,
    corr_length_m=0.05,
    theta_deg=thetas,                 # scalar or array
    er_complex=complex(15.0, 3.0),
    correl="gaussian",
    include_hv=True,
    return_db=True,                   # return dB (default)
)
print(sig["hh"], sig["vv"], sig["hv"])  # arrays in dB
```

3) Bistatic backscatter

```python
from pyi2em import sigma0_bistatic

sig = sigma0_bistatic(
    freq_ghz=5.0,
    rms_height_m=0.01,
    corr_length_m=0.05,
    thi_deg=40.0,
    ths_deg=40.0,
    phs_deg=180.0,
    er_complex=complex(15.0, 3.0),
    correl="gaussian",
    xcoeff=1.0,
    return_db=True,
)
print(sig["hh"], sig["vv"])  # dB
```

Correlation options: 'exponential' | 'gaussian' | 'x_power' | 'x_exponential'.

Notes
- sigma0_* return dB by default; set return_db=False for linear values.
- sigma0_backscatter(theta_deg=...) accepts a scalar or a NumPy array.
- emissivity returns linear values in [0, 1].

## License and Notices

- Wrapper: MIT (see `LICENSE`)
- GSL: GPL-licensed dependency; the `NOTICE` and `COPYING.GPL` files are included.
