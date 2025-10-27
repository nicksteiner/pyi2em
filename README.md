![Tests](https://github.com/nicksteiner/pyi2em/actions/workflows/tests.yml/badge.svg?branch=master)

![Tests](https://github.com/nicksteiner/pyi2em/actions/workflows/tests.yml/badge.svg?branch=master)

# I2EM Python Library

This code calculates the emission from rough surfaces using the I2EM model          

Ulaby, F.T. and Long D.G.(2014), Microwave Radar and Radiometric Remote Sensing, The University of Michigan Press  

## Install

Wheels (recommended)

```bash
pip install pyi2em
```

Notes on wheels vs source builds

- PyPI wheels bundle required native libraries, including GSL and the embedded cubature routines, so end users do not need to install system packages.
- If you build from source, you must have GSL installed on your system and discoverable via pkg-config. The cubature integration code is compiled into the extension automatically.

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

### Emissivity Calculation

```python
    
import pyi2em

#set frequency in ghz
fr = 3.0

# set complex dielectric of the soil
el = 11.3
ei = 1.5

# set correlation length [m]
l = 0.10  # 10 cm

# set standard deviation of the surface height variation (rms) [m]
sig = 0.0025  # .25 cm

# set incidence angle [deg]
theta_d = 30.0

# set type of surface correlation function (1) exponential (2) Gaussian
sp = 2

e_ = pyi2em.I2EM(fr, sig, l, theta_d, el, ei, sp)

print("FREQ: 3.0 [GHZ], CDC: 11.3 + i1.5, CL: 10 [cm], RMS: .25 [cm], INC: 30 [deg], CORRF: Gaussian")

print("EMISSIVITY: {:g} [V], {:g} [H]".format(*e_))
```
```
FREQ: 3.0 [GHZ], CDC: 11.3 + i1.5, CL: 10 [cm], RMS: .25 [cm], INC: 30 [deg], CORRF: Gaussian
Emissivity: 0.646743 [V], 0.75031 [H]
```

### Bistatic Scattering 

```python
import pyi2em

# Set frequency in GHz
freq_ghz = 5.0

# Set complex dielectric constant
el = 15.0  # real part
ei = 3.0   # imaginary part

# Set correlation length [m]
correl_length = 0.05  # 5 cm

# Set standard deviation of surface height variation (rms) [m]
rmsheight = 0.01  # 1 cm

# Set incident and scatter angles [degrees]
thi = 40.0    # incident angle
ths = 40.0    # scatter angle (same as incident for backscatter)
phs = 180.0   # scatter azimuth angle (180° for backscatter)

# Set correlation function type
correl_func = pyi2em.CORREL_GAUSSIAN

# Get bistatic scattering coefficients
sigma0 = pyi2em.I2EM_Bistat(freq_ghz, rmsheight, correl_length,
                            thi, ths, phs, el, ei, correl_func)

print(f"σ⁰_HH = {sigma0[0]:.3f} dB")
print(f"σ⁰_VV = {sigma0[1]:.3f} dB")
```

#### Available Correlation Functions

- `pyi2em.CORREL_EXPONENTIAL` (1) - Exponential correlation function
- `pyi2em.CORREL_GAUSSIAN` (2) - Gaussian correlation function  
- `pyi2em.CORREL_X_POWER` (3) - Power law correlation function
- `pyi2em.CORREL_X_EXPONENTIAL` (4) - Exponential power correlation function

For more examples, see `test_bistat.py`

### Backscatter, Single-Scale Random Surface

```python
import pyi2em

freq_ghz = 5.0
rmsheight = 0.01      # m
correl_length = 0.05  # m
theta = 40.0          # deg
el, ei = 15.0, 3.0

# [HH_dB, VV_dB, HV_dB or NaN]
sigma0 = pyi2em.I2EM_Backscatter_Single(freq_ghz, rmsheight, correl_length,
                                        theta, el, ei,
                                        pyi2em.CORREL_GAUSSIAN, 1.0, False)
```

### Backscatter, Periodic Sinusoidal Surface

```python
import pyi2em

theta_0, phi_0 = 30.0, 0.0  # deg
el, ei = 12.0, 2.0
freq_ghz = 5.3
rmsheight = 0.008    # m
correl_length = 0.05 # m
Gmm = 0.20           # period [m]
A = 0.01             # amplitude [m]

# [VV_dB, HH_dB, VH_dB]
sigma0 = pyi2em.I2EM_Backscatter_Periodic(theta_0, phi_0, el, ei,
                                          freq_ghz, rmsheight, correl_length,
                                          pyi2em.CORREL_GAUSSIAN,
                                          Gmm, A)
```
