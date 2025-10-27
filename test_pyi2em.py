#!/usr/bin/env python3
"""
Test script for I2EM bistatic scattering model

This demonstrates the use of the bistatic scattering model in pyi2em
"""

import pyi2em
import numpy as np

def test_backscatter():
    """Test backscatter configuration (monostatic)"""
    print("=" * 70)
    print("TEST 1: Backscatter Configuration (Monostatic)")
    print("=" * 70)
    
    # Set frequency in GHz
    freq_ghz = 5.0
    
    # Set complex dielectric of the soil
    el = 15.0  # real part
    ei = 3.0   # imaginary part
    
    # Set correlation length [m]
    correl_length = 0.05  # 5 cm
    
    # Set standard deviation of the surface height variation (rms) [m]
    rmsheight = 0.01  # 1 cm
    
    # Set incidence angle [degrees]
    thi = 40.0
    
    # For backscatter: scattered angle = incident angle, azimuth = 180 degrees
    ths = thi
    phs = 180.0
    
    # Set type of surface correlation function
    correl_func = pyi2em.CORREL_GAUSSIAN
    
    # Correlation coefficient (used for exponential power spectrum)
    xcoeff = 1.0
    
    # Using explicit bistatic call in backscatter geometry
    sigma0 = pyi2em.I2EM_Bistat(freq_ghz, rmsheight, correl_length, 
                                thi, ths, phs, el, ei, 
                                correl_func, xcoeff)
    
    print(f"Frequency: {freq_ghz} GHz")
    print(f"Complex dielectric constant: {el} + i{ei}")
    print(f"RMS height: {rmsheight*100} cm")
    print(f"Correlation length: {correl_length*100} cm")
    print(f"Incidence angle: {thi}°")
    print(f"Correlation function: GAUSSIAN")
    print(f"\nBackscatter coefficients (dB):")
    print(f"  σ⁰_HH = {sigma0[0]:.3f} dB")
    print(f"  σ⁰_VV = {sigma0[1]:.3f} dB")
    print()

    # Using new convenience API for Code 10.1 (single-scale backscatter)
    sigma0_single = pyi2em.I2EM_Backscatter_Single(freq_ghz, rmsheight, correl_length,
                                                   thi, el, ei, correl_func, xcoeff, False)
    print("Backscatter (Code 10.1) via convenience API:")
    print(f"  σ⁰_HH = {sigma0_single[0]:.3f} dB, σ⁰_VV = {sigma0_single[1]:.3f} dB, σ⁰_HV = {sigma0_single[2]}")
    print()


def test_bistatic():
    """Test bistatic configuration"""
    print("=" * 70)
    print("TEST 2: Bistatic Configuration")
    print("=" * 70)
    
    # Set frequency in GHz
    freq_ghz = 10.0
    
    # Set complex dielectric of the soil
    el = 20.0  # real part
    ei = 5.0   # imaginary part
    
    # Set correlation length [m]
    correl_length = 0.08  # 8 cm
    
    # Set standard deviation of the surface height variation (rms) [m]
    rmsheight = 0.015  # 1.5 cm
    
    # Set incident angle [degrees]
    thi = 30.0
    
    # Set scattered angle and azimuth [degrees]
    ths = 45.0
    phs = 90.0
    
    # Set type of surface correlation function
    correl_func = pyi2em.CORREL_EXPONENTIAL
    
    # Correlation coefficient
    xcoeff = 1.0
    
    sigma0 = pyi2em.I2EM_Bistat(freq_ghz, rmsheight, correl_length, 
                                thi, ths, phs, el, ei, 
                                correl_func, xcoeff)
    
    print(f"Frequency: {freq_ghz} GHz")
    print(f"Complex dielectric constant: {el} + i{ei}")
    print(f"RMS height: {rmsheight*100} cm")
    print(f"Correlation length: {correl_length*100} cm")
    print(f"Incidence angle: {thi}°")
    print(f"Scatter angle: {ths}°")
    print(f"Scatter azimuth: {phs}°")
    print(f"Correlation function: EXPONENTIAL")
    print(f"\nBistatic scattering coefficients (dB):")
    print(f"  σ⁰_HH = {sigma0[0]:.3f} dB")
    print(f"  σ⁰_VV = {sigma0[1]:.3f} dB")
    print()


def test_angular_response():
    """Test angular response for backscatter"""
    print("=" * 70)
    print("TEST 3: Angular Response (Backscatter)")
    print("=" * 70)
    
    # Set frequency in GHz
    freq_ghz = 5.3
    
    # Set complex dielectric of the soil
    el = 10.0  # real part
    ei = 2.0   # imaginary part
    
    # Set correlation length [m]
    correl_length = 0.10  # 10 cm
    
    # Set standard deviation of the surface height variation (rms) [m]
    rmsheight = 0.008  # 0.8 cm
    
    # Set type of surface correlation function
    correl_func = pyi2em.CORREL_GAUSSIAN
    xcoeff = 1.0
    
    print(f"Frequency: {freq_ghz} GHz")
    print(f"Complex dielectric constant: {el} + i{ei}")
    print(f"RMS height: {rmsheight*100} cm")
    print(f"Correlation length: {correl_length*100} cm")
    print(f"\nAngular response:")
    print(f"{'Angle':<8} {'σ⁰_HH (dB)':<12} {'σ⁰_VV (dB)':<12}")
    print("-" * 35)
    
    angles = [20, 30, 40, 50, 60]
    
    for thi in angles:
        # Backscatter configuration
        ths = thi
        phs = 180.0

        sigma0 = pyi2em.I2EM_Bistat(freq_ghz, rmsheight, correl_length,
                                    thi, ths, phs, el, ei,
                                    correl_func, xcoeff)

        print(f"{thi:<8}° {sigma0[0]:<12.3f} {sigma0[1]:<12.3f}")
    print()


def test_with_default_params():
    """Test with default parameters"""
    print("=" * 70)
    print("TEST 4: Using Default Parameters")
    print("=" * 70)
    
    # Minimal parameters (correl_func and xcoeff use defaults)
    freq_ghz = 3.0
    rmsheight = 0.005  # 0.5 cm
    correl_length = 0.05  # 5 cm
    thi = 35.0
    ths = 35.0
    phs = 180.0
    el = 12.0
    ei = 2.5
    
    # Call with default correlation function (GAUSSIAN) and xcoeff (1.0)
    sigma0 = pyi2em.I2EM_Bistat(freq_ghz, rmsheight, correl_length, 
                                thi, ths, phs, el, ei)
    
    print(f"Using default correlation function (GAUSSIAN) and xcoeff (1.0)")
    print(f"\nBackscatter coefficients (dB):")
    print(f"  σ⁰_HH = {sigma0[0]:.3f} dB")
    print(f"  σ⁰_VV = {sigma0[1]:.3f} dB")
    print()


def test_periodic_surface():
    """Test Code 10.4: Backscatter from periodic sinusoidal surface"""
    print("=" * 70)
    print("TEST 5: Periodic Sinusoidal Surface (Code 10.4)")
    print("=" * 70)

    # Surface and observation parameters
    theta_0 = 30.0   # deg
    phi_0 = 0.0      # deg
    el, ei = 12.0, 2.0
    freq_ghz = 5.3
    rmsheight = 0.008   # 0.8 cm
    correl_length = 0.05 # 5 cm
    correl_func = pyi2em.CORREL_GAUSSIAN
    Gmm = 0.20  # period [m]
    A = 0.01    # amplitude [m]

    s_vv, s_hh, s_vh = pyi2em.I2EM_Backscatter_Periodic(theta_0, phi_0, el, ei,
                                                         freq_ghz, rmsheight, correl_length,
                                                         correl_func, Gmm, A)
    print(f"σ⁰_VV = {s_vv:.3f} dB, σ⁰_HH = {s_hh:.3f} dB, σ⁰_VH = {s_vh:.3f} dB")
    print()


def main():
    """Run all tests"""
    print("\n")
    print("*" * 70)
    print("I2EM Bistatic Scattering Model Test Suite")
    print("*" * 70)
    print()
    
    # Display available correlation functions
    print("Available correlation functions:")
    print(f"  CORREL_EXPONENTIAL = {pyi2em.CORREL_EXPONENTIAL}")
    print(f"  CORREL_GAUSSIAN = {pyi2em.CORREL_GAUSSIAN}")
    print(f"  CORREL_X_POWER = {pyi2em.CORREL_X_POWER}")
    print(f"  CORREL_X_EXPONENTIAL = {pyi2em.CORREL_X_EXPONENTIAL}")
    print()
    
    test_backscatter()
    test_bistatic()
    test_angular_response()
    test_with_default_params()
    test_periodic_surface()
    
    print("*" * 70)
    print("All tests completed successfully!")
    print("*" * 70)
    print()


if __name__ == '__main__':
    main()
