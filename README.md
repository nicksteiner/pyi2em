# I2EM Library built using PyBind11

This code calculates the emission from rough surfaces using the I2EM model          

Ulaby, F.T. and Long D.G.(2014), Microwave Radar and Radiometric Remote Sensing, The University of Michigan Press  

## Install  
```git clone https://github.com/nicksteiner/pyi2em```  
```pip install ./pyi2em```

## Test
```bash
cd pyi2em
python3 test_i2em.py
```

## Usage

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

    e_ = pyi2em.tI2EM(fr, sig, l, theta_d, el, ei, sp)

    print("FREQ: 3.0 [GHZ], CDC: 11.3 + i1.5, CL: 10 [cm], RMS: .25 [cm], INC: 30 [deg], CORRF: Gaussian")
    print("Emissivity: {} [V], {} [H]".format(*e_))

```

