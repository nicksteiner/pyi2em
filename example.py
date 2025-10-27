#!/usr/bin/env python3
"""
Example demonstrating emissivity and scattering APIs in pyi2em
"""

from pyi2em import emissivity, sigma0_backscatter, sigma0_bistatic
import numpy as np

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
    eh, ev = emissivity(
        freq_ghz=freq_ghz,
        rms_height_m=rmsheight,
        corr_length_m=correl_length,
        theta_deg=theta,
        er_complex=complex(el, ei),
        correl="gaussian",
    )
    print(f"  e_H = {eh:.4f}")
    print(f"  e_V = {ev:.4f}")
    print()
    
    # 2. Backscatter calculation
    print("-" * 70)
    print("2. BISTATIC MODEL - Backscatter Configuration")
    print("-" * 70)
    # For backscatter: scatter angle = incident angle, azimuth = 180°
    sig_back = sigma0_backscatter(
        freq_ghz=freq_ghz,
        rms_height_m=rmsheight,
        corr_length_m=correl_length,
        theta_deg=theta,
        er_complex=complex(el, ei),
        correl="gaussian",
        include_hv=True,
        return_db=True,
    )
    print(f"  σ⁰_HH = {sig_back['hh'].item():.3f} dB")
    print(f"  σ⁰_VV = {sig_back['vv'].item():.3f} dB")
    print(f"  σ⁰_HV = {sig_back['hv'].item():.3f} dB")
    print()
    
    # 3. Bistatic scattering
    print("-" * 70)
    print("3. BISTATIC MODEL - Forward Scatter Configuration")
    print("-" * 70)
    ths = 50.0  # scatter angle
    phs = 45.0  # scatter azimuth
    sig_bi = sigma0_bistatic(
        freq_ghz=freq_ghz,
        rms_height_m=rmsheight,
        corr_length_m=correl_length,
        thi_deg=theta,
        ths_deg=ths,
        phs_deg=phs,
        er_complex=complex(el, ei),
        correl="gaussian",
        xcoeff=1.0,
        return_db=True,
    )
    print(f"  Incident angle: {theta}°")
    print(f"  Scatter angle: {ths}°")
    print(f"  Scatter azimuth: {phs}°")
    print(f"  σ⁰_HH = {sig_bi['hh']:.3f} dB")
    print(f"  σ⁰_VV = {sig_bi['vv']:.3f} dB")
    print()
    
    # 4. Compare different correlation functions
    print("-" * 70)
    print("4. EFFECT OF CORRELATION FUNCTION (Backscatter)")
    print("-" * 70)
    
    corr_funcs = [
        ("Exponential", "exponential"),
        ("Gaussian", "gaussian"),
    ]
    
    print(f"{'Function':<15} {'σ⁰_HH (dB)':<12} {'σ⁰_VV (dB)':<12}")
    print("-" * 40)
    
    for name, func in corr_funcs:
        sig = sigma0_bistatic(
            freq_ghz=freq_ghz,
            rms_height_m=rmsheight,
            corr_length_m=correl_length,
            thi_deg=theta,
            ths_deg=theta,
            phs_deg=180.0,
            er_complex=complex(el, ei),
            correl=func,
        )
        print(f"{name:<15} {sig['hh']:<12.3f} {sig['vv']:<12.3f}")
    
    print()
    print("=" * 70)
    print("Example completed!")
    print("=" * 70)


if __name__ == '__main__':
    main()
