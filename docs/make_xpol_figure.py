#!/usr/bin/env python3
"""Generate the cross-pol validation figures for the README.

Compares pyi2em's hv (I2EM multiple-scattering cross-pol) against the
published reference from Fung, "Microwave Scattering and Emission Models
for Users" (Artech House, 2010), Fig. 3.9 -- dielectric variation, exponential
correlation. Reference values are the full-precision numbers cached in the
companion Mathematica notebook.

Parameters (Fig 3.9): f = 5 GHz, sig = 0.3 cm, L = 4 cm, exponential, eps real.

Outputs:
    docs/figures/xpol_fig39_validation.png   (curves + residual panel)

Run:  PYTHONPATH=<build dir> python3 docs/make_xpol_figure.py
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pyi2em

FREQ_GHZ, RMS_M, CORR_M = 5.0, 0.003, 0.04  # 0.3 cm, 4 cm
CORREL = "exponential"

THETAS_DEG = np.array([0.001, 10.001, 20.001, 30.001, 40.001,
                       50.001, 60.001, 70.001])

# Full-precision reference cached in the Fung Fig 3.9 notebook.
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

COLORS = {3.5: "#1f77b4", 7.5: "#2ca02c", 15.0: "#d62728"}


def pyi2em_hv(er_real, thetas):
    out = pyi2em.sigma0_backscatter(
        FREQ_GHZ, RMS_M, CORR_M, list(thetas), complex(er_real, 0.0),
        correl=CORREL, include_hv=True, return_db=True,
    )
    return np.asarray(out["hv"], float)


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    outdir = os.path.join(here, "figures")
    os.makedirs(outdir, exist_ok=True)

    # smooth model curve + reference markers
    theta_fine = np.linspace(0.001, 70.001, 141)

    fig, (ax, axr) = plt.subplots(
        2, 1, figsize=(7.2, 6.6), height_ratios=[3, 1.1], sharex=True,
        gridspec_kw={"hspace": 0.08},
    )

    max_abs = 0.0
    for er in (3.5, 7.5, 15.0):
        c = COLORS[er]
        model = pyi2em_hv(er, theta_fine)
        ax.plot(theta_fine, model, "-", color=c, lw=1.8,
                label=fr"pyi2em  $\varepsilon_r$={er:g}")
        ax.plot(THETAS_DEG, REFERENCE_DB[er], "o", color=c, ms=6,
                mfc="white", mec=c, mew=1.5,
                label=f"Fung Fig 3.9  ($\\varepsilon_r$={er:g})")

        # residual at reference angles
        model_at_ref = pyi2em_hv(er, THETAS_DEG)
        resid = model_at_ref - np.asarray(REFERENCE_DB[er])
        max_abs = max(max_abs, np.abs(resid).max())
        axr.plot(THETAS_DEG, resid, "o-", color=c, ms=4, lw=1.0)

    ax.set_ylabel(r"$\sigma^0_{vh}$  (dB)")
    ax.set_title("I2EM cross-pol (VH) vs Fung Fig 3.9\n"
                 "f = 5 GHz, $\\sigma$ = 0.3 cm, L = 4 cm, exponential",
                 fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, ncol=3, loc="upper right", framealpha=0.9)

    axr.axhline(0.0, color="0.5", lw=0.8)
    axr.set_ylabel("resid.\n(dB)", fontsize=9)
    axr.set_xlabel(r"incidence angle  $\theta$  (deg)")
    axr.grid(True, alpha=0.3)
    axr.set_ylim(-0.6, 0.6)
    axr.text(0.01, 0.04,
             f"max |pyi2em - Fung| = {max_abs:.2f} dB",
             transform=axr.transAxes, fontsize=8, va="bottom")

    fig.savefig(os.path.join(outdir, "xpol_fig39_validation.png"),
                dpi=140, bbox_inches="tight")
    print(f"wrote {os.path.join(outdir, 'xpol_fig39_validation.png')}")
    print(f"max |residual| = {max_abs:.3f} dB")


if __name__ == "__main__":
    main()
