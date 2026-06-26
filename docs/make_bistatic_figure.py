#!/usr/bin/env python3
"""Bistatic validation figure: pyi2em I2EM vs Fung IEM-M (cross-model).

Compares pyi2em's sigma0_bistatic against the cached bistatic reference from
Fung, "Microwave Scattering and Emission Models for Users" (Artech House, 2010),
Fig. 6.1 (VV) and Fig. 6.2 (HH) -- azimuthal bistatic scattering, exponential
correlation. The reference is computed by Fung's IEM-M single-surface bistatic
program (notebook "Fig 6.1, 6.2.nb"); pyi2em implements the I2EM bistatic model
of Ulaby & Long (2014). These are RELATED BUT DIFFERENT formulations, so this is
a labeled cross-model comparison, not a same-model validation.

Geometry (verified from the notebook's parameter cell):
    f = 10 GHz, theta_i = 20.03 deg, theta_s = 50 deg, eps_r = 3 (real),
    L = 3 cm, exponential correlation (sp = 1), phi_s swept 0.03..180.03 deg.
    rms-height family: sigma = 0.2 / 0.4 / 0.7 cm (Fung curves vv1/2/3, hh1/2/3).

Reference values are extracted verbatim from the notebook's cached vv*/hh* cells.

Outputs: docs/figures/bistatic_fig61_62_validation.png

Run: PYTHONPATH=<build dir> python3 docs/make_bistatic_figure.py
"""
import os
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pyi2em

NB = ("/mnt/data5/fung-037/Ch-6 Figs/chap6 - 1 to 11/Fig 6.1, 6.2.nb")

# Verified geometry / surface parameters
FREQ_GHZ = 10.0
THI_DEG = 20.03
THS_DEG = 50.0
ER = complex(3.0, 0.0)
CORR_M = 0.03            # L = 3 cm
CORREL = "exponential"   # sp = 1
SIGMAS_CM = {"1": 0.2, "2": 0.4, "3": 0.7}   # vv1/hh1 .. vv3/hh3
COLORS = {"1": "#1f77b4", "2": "#2ca02c", "3": "#d62728"}


def extract_curve(text, name):
    """Return [(phi_s_deg, sigma0_dB), ...] for a cached vv*/hh* assignment."""
    m = re.search(r"\b" + name + r"\s*=\s*\{", text)
    j = text.index("{", m.start())
    depth = 0
    for k in range(j, len(text)):
        if text[k] == "{":
            depth += 1
        elif text[k] == "}":
            depth -= 1
            if depth == 0:
                block = text[j:k + 1]
                break
    s = block.replace("`", "").replace("\\(", "").replace("\\)", "")
    s = re.sub(r"\s+", "", s)
    return [(float(a), float(b))
            for a, b in re.findall(r"\{(-?\d+\.?\d*),(-?\d+\.?\d*)\}", s)]


def pyi2em_bistatic(sig_cm, phs_deg):
    vv, hh = [], []
    for phs in phs_deg:
        out = pyi2em.sigma0_bistatic(
            FREQ_GHZ, sig_cm / 100.0, CORR_M, THI_DEG, THS_DEG, float(phs),
            ER, correl=CORREL, xcoeff=1.0, return_db=True,
        )
        vv.append(out["vv"])
        hh.append(out["hh"])
    return np.array(vv), np.array(hh)


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    outdir = os.path.join(here, "figures")
    os.makedirs(outdir, exist_ok=True)
    text = open(NB, "r", errors="replace").read()

    ref = {}
    for i in ("1", "2", "3"):
        ref[("vv", i)] = extract_curve(text, f"vv{i}")
        ref[("hh", i)] = extract_curve(text, f"hh{i}")

    # Away-from-null agreement metric: the specular null is a near-vertical
    # notch where a small null-position shift between the two models produces a
    # huge dB difference -- not representative of model agreement. Report the
    # median |diff| over points that are within 12 dB of each curve's own peak
    # (i.e. excluding the deep-null floor).
    fig, (axv, axh) = plt.subplots(1, 2, figsize=(12, 5.2), sharey=True)
    worst = 0.0
    away_diffs = []
    for ax, pol in ((axv, "vv"), (axh, "hh")):
        for i in ("1", "2", "3"):
            c = COLORS[i]
            sig = SIGMAS_CM[i]
            curve = ref[(pol, i)]
            phs = np.array([p for p, _ in curve])
            refdb = np.array([v for _, v in curve])
            model_vv, model_hh = pyi2em_bistatic(sig, phs)
            model = model_vv if pol == "vv" else model_hh
            diff = np.abs(model - refdb)
            worst = max(worst, diff.max())
            away = refdb > (refdb.max() - 12.0)   # exclude deep-null floor
            away_diffs.extend(diff[away].tolist())
            ax.plot(phs, model, "-", color=c, lw=2.0,
                    label=fr"pyi2em  $\sigma$={sig} cm")
            ax.plot(phs, refdb, "o", color=c, ms=5, mfc="white", mec=c,
                    mew=1.4, label=f"Fung IEM-M ($\\sigma$={sig} cm)")
        ax.set_xlabel(r"bistatic azimuth  $\phi_s$  (deg)")
        ax.set_title(f"{pol.upper()}", fontsize=11)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=7, ncol=2, loc="lower left")
    axv.set_ylabel(r"$\sigma^0$  (dB)")
    away_diffs = np.array(away_diffs)
    med = float(np.median(away_diffs))
    p90 = float(np.percentile(away_diffs, 90))
    fig.suptitle(
        "Bistatic scattering: pyi2em I2EM vs Fung IEM-M (cross-model)\n"
        f"f=10 GHz, $\\theta_i$=20deg, $\\theta_s$=50deg, "
        f"$\\varepsilon_r$=3, L=3 cm, exponential   "
        f"(away from specular null: median |diff|={med:.1f} dB, "
        f"90th pct={p90:.1f} dB)",
        fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    p = os.path.join(outdir, "bistatic_fig61_62_validation.png")
    fig.savefig(p, dpi=140, bbox_inches="tight")
    print(f"wrote {p}")
    print(f"away-from-null: median |diff| = {med:.2f} dB, "
          f"90th pct = {p90:.2f} dB")
    print(f"max |diff| (incl. specular null) = {worst:.2f} dB "
          f"(null-position shift artifact, not model disagreement)")


if __name__ == "__main__":
    main()
