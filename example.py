#!/usr/bin/env python3
"""
Example demonstrating both emissivity and bistatic scattering models in pyi2em
"""

import pyi2em

def main():
    print("=" * 70)
    print("pyi2em Example: Emissivity and Bistatic Scattering")
    print("=" * 70)
    print()
    
    # Common parameters
    freq_ghz = 5.3
    el = 15.0  # real part of dielectric constant
    ei = 3.0   # imaginary part
    rmsheight = 0.01  # 1 cm
    correl_length = 0.08  # 8 cm
    theta = 40.0  # degrees
    
    print("Surface Parameters:")
    print(f"  Frequency: {freq_ghz} GHz")
    print(f"  Dielectric constant: {el} + i{ei}")
    print(f"  RMS height: {rmsheight*100} cm")
    print(f"  Correlation length: {correl_length*100} cm")
    print(f"  Angle: {theta}°")
    print()
    
    # 1. Emissivity calculation
    print("-" * 70)
    print("1. EMISSIVITY MODEL")
    print("-" * 70)
    emissivity = pyi2em.I2EM(freq_ghz, rmsheight, correl_length, 
                             theta, el, ei, pyi2em.CORREL_GAUSSIAN)
    print(f"  e_V = {emissivity[0]:.4f}")
    print(f"  e_H = {emissivity[1]:.4f}")
    print()
    
    # 2. Backscatter calculation
    print("-" * 70)
    print("2. BISTATIC MODEL - Backscatter Configuration")
    print("-" * 70)
    # For backscatter: scatter angle = incident angle, azimuth = 180°
    sigma0_back = pyi2em.I2EM_Bistat(freq_ghz, rmsheight, correl_length,
                                     theta, theta, 180.0,
                                     el, ei, pyi2em.CORREL_GAUSSIAN)
    print(f"  σ⁰_HH = {sigma0_back[0]:.3f} dB")
    print(f"  σ⁰_VV = {sigma0_back[1]:.3f} dB")
    print()
    
    # 3. Bistatic scattering
    print("-" * 70)
    print("3. BISTATIC MODEL - Forward Scatter Configuration")
    print("-" * 70)
    ths = 50.0  # scatter angle
    phs = 45.0  # scatter azimuth
    sigma0_bistat = pyi2em.I2EM_Bistat(freq_ghz, rmsheight, correl_length,
                                       theta, ths, phs,
                                       el, ei, pyi2em.CORREL_GAUSSIAN)
    print(f"  Incident angle: {theta}°")
    print(f"  Scatter angle: {ths}°")
    print(f"  Scatter azimuth: {phs}°")
    print(f"  σ⁰_HH = {sigma0_bistat[0]:.3f} dB")
    print(f"  σ⁰_VV = {sigma0_bistat[1]:.3f} dB")
    print()
    
    # 4. Compare different correlation functions
    print("-" * 70)
    print("4. EFFECT OF CORRELATION FUNCTION (Backscatter)")
    print("-" * 70)
    
    corr_funcs = [
        ("Exponential", pyi2em.CORREL_EXPONENTIAL),
        ("Gaussian", pyi2em.CORREL_GAUSSIAN),
    ]
    
    print(f"{'Function':<15} {'σ⁰_HH (dB)':<12} {'σ⁰_VV (dB)':<12}")
    print("-" * 40)
    
    for name, func in corr_funcs:
        sigma0 = pyi2em.I2EM_Bistat(freq_ghz, rmsheight, correl_length,
                                    theta, theta, 180.0,
                                    el, ei, func)
        print(f"{name:<15} {sigma0[0]:<12.3f} {sigma0[1]:<12.3f}")
    
    print()
    print("=" * 70)
    print("Example completed!")
    print("=" * 70)


if __name__ == '__main__':
    main()
