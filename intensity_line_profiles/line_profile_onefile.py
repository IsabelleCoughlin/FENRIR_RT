#!/usr/bin/env python3
"""Make four FENRIR line profiles from matched disk-image/emissivity files.

Usage
-----
python line_profile_onefile.py OUT.png SPIN R_OUT \
    DISK_MDOT03.npy DISK_MDOT02.npy DISK_MDOT01.npy DISK_MDOT00.npy \
    EMISS_MDOT03.npy EMISS_MDOT02.npy EMISS_MDOT01.npy EMISS_MDOT00.npy

The input order is thick to thin: mdot = 0.3, 0.2, 0.1, 0.0.
"""

import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


N_ENERGY_BINS = 150
ENERGY_EDGES = np.linspace(0.0, 2.0, N_ENERGY_BINS + 1)
#ENERGY_EDGES = np.linspace(
#    0.002,
#    2.002,
#    N_ENERGY_BINS + 1,
#)
ENERGY_MIDS = 0.5 * (ENERGY_EDGES[1:] + ENERGY_EDGES[:-1])


def r_isco(a):
    """Kerr ISCO radius in units of r_g for -1 <= a <= 1."""
    if not -1.0 <= a <= 1.0:
        raise ValueError("spin must satisfy -1 <= a <= 1")

    z1 = 1.0 + (1.0 - a * a) ** (1.0 / 3.0) * (
        (1.0 + a) ** (1.0 / 3.0) + (1.0 - a) ** (1.0 / 3.0)
    )
    z2 = np.sqrt(3.0 * a * a + z1 * z1)
    root = np.sqrt((3.0 - z1) * (3.0 + z1 + 2.0 * z2))
    return 3.0 + z2 + root if a < 0.0 else 3.0 + z2 - root


def load_emissivity(path):
    """Return sorted, finite radial emissivity samples from a FENRIR NPY."""
    data = np.load(path, allow_pickle=True)
    rho = np.asarray(data[0][0], dtype=float).ravel()
    flux = np.asarray(data[0][1], dtype=float).ravel()

    good = np.isfinite(rho) & np.isfinite(flux) & (rho > 0.0) & (flux >= 0.0)
    rho = rho[good]
    flux = flux[good]
    if rho.size < 2:
        raise ValueError(f"{path}: fewer than two usable emissivity bins")

    order = np.argsort(rho)
    rho = rho[order]
    flux = flux[order]

    # np.interp requires increasing sample positions.  Duplicate radii are
    # not expected, but retaining only the first occurrence makes this safe.
    rho, unique_index = np.unique(rho, return_index=True)
    return rho, flux[unique_index]


def make_line_profile(disk_path, emissivity_path, spin, r_out):
    """Histogram the observed shifts, weighting every image pixel by g^3 eps."""
    disk = np.load(disk_path, allow_pickle=True)
    if len(disk) < 9:
        raise ValueError(f"{disk_path}: expected the nine disk-image arrays")

    g = np.asarray(disk[2], dtype=float).ravel()
    rho = np.asarray(disk[8], dtype=float).ravel()
    if g.shape != rho.shape:
        raise ValueError(f"{disk_path}: g and cylindrical-radius shapes differ")

    em_rho, em_flux = load_emissivity(emissivity_path)
    inner = max(r_isco(spin), em_rho[0])
    outer = min(r_out, em_rho[-1])
    if outer <= inner:
        raise ValueError(
            f"{emissivity_path}: emissivity range does not overlap "
            f"[{r_isco(spin):g}, {r_out:g})"
        )

    good = (
        np.isfinite(rho)
        & np.isfinite(g)
        & (rho >= inner)
        & (rho < outer)
        & (g > ENERGY_EDGES[0])
        & (g < ENERGY_EDGES[-1])
    )
    g = g[good]
    rho = rho[good]
    if g.size == 0:
        raise ValueError(f"{disk_path}: no usable disk rays after the cuts")

    #local_emissivity = np.interp(rho, em_rho, em_flux)
    indices = np.searchsorted(em_rho, rho, side="left")
    indices = np.clip(indices, 0, len(em_flux) - 1)
    local_emissivity = em_flux[indices]
    # This is the convention used by FENRIR's create_lineprof_lp_auto.py:
    # each equal-area image-plane pixel contributes eps(rho) * g^3.
    # Use g^4 instead only if Phi is explicitly defined as an energy-flux
    # profile of a delta-function line rather than the historical convention.
    weights = local_emissivity * g**3
    profile = np.histogram(g, bins=ENERGY_EDGES, weights=weights)[0]
    return profile, g.size


def main():
    if len(sys.argv) != 12:
        raise SystemExit(
            "Usage: python line_profile_onefile.py OUT.png SPIN R_OUT "
            "DISK03 DISK02 DISK01 DISK00 EMISS03 EMISS02 EMISS01 EMISS00"
        )

    outfile = sys.argv[1]
    spin = float(sys.argv[2])
    r_out = float(sys.argv[3])
    disk_files = sys.argv[4:8]
    emissivity_files = sys.argv[8:12]

    mdots = [0.3, 0.2, 0.1, 0.0]
    colors = ["#0055FF", "tab:red", "0.5", "k"]
    profiles = []

    for mdot, disk_path, emissivity_path in zip(
        mdots, disk_files, emissivity_files
    ):
        profile, ray_count = make_line_profile(
            disk_path, emissivity_path, spin, r_out
        )
        profiles.append(profile)
        print(
            f"mdot={mdot:.1f}: rays={ray_count}, "
            f"sum={profile.sum():.8g}, max={profile.max():.8g}"
        )

    # Match the comparison convention: normalize every curve by the maximum
    # of the thin-disk (mdot=0) profile, preserving relative amplitudes.
    normalization = profiles[-1].max()
    if not np.isfinite(normalization) or normalization <= 0.0:
        raise ValueError("the thin-disk profile has no positive normalization")

    print("Energy bins:", N_ENERGY_BINS)
    print("Bin width:", ENERGY_EDGES[1] - ENERGY_EDGES[0])

    for mdot, profile in zip(mdots, profiles):
        index = np.argmax(profile)

        print(
            f"mdot={mdot:.1f}: "
            f"raw maximum={profile[index]:.8g}, "
            f"peak midpoint={ENERGY_MIDS[index]:.8f}, "
            f"total={np.sum(profile):.8g}"
        )

    print("Thin normalization:", profiles[-1].max())

    '''
    fig, ax = plt.subplots(figsize=(7, 5))
    for mdot, color, profile in reversed(list(zip(mdots, colors, profiles))):
        ax.plot(
            ENERGY_MIDS,
            profile / normalization,
            color=color,
            linewidth=2.0,
            alpha=1.0,
            label=rf"$\dot{{m}}={mdot:.1f}$",
        )

    ax.set_xlabel(r"$g=E_{\rm obs}/E_{\rm em}$")
    ax.set_ylabel(r"$\Phi$ [Arbitrary Units]")
    ax.set_xlim(0.0, 1.5)
    ax.set_ylim(bottom=0.0)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(outfile, dpi=250)
    print(f"Created {outfile}")
    '''
        # Configure global styling for thick borders and inward ticks
    plt.rcParams["axes.linewidth"] = 2.5  # Thicker frame spines

    fig, ax = plt.subplots(figsize=(6, 6))  # Square aspect ratio to match reference

    # Plot the lines
    for mdot, color, profile in reversed(list(zip(mdots, colors, profiles))):
        ax.plot(
            ENERGY_MIDS,
            profile / normalization,
            color=color,
            linewidth=2.5,  # Slightly thicker lines to match the reference style
            alpha=1.0,
            label=rf"$\dot{{m}}={mdot:.1f}$",
        )

    # Replicate the specific axis and tick styles
    ax.set_xlim(0.0, 1.5)
    ax.set_ylim(bottom=0.0)

    # Show ticks on all 4 sides pointing inward
    ax.tick_params(
        axis="both",
        direction="in",
        top=True,
        right=True,
        width=2.5,
        length=8,
        labelsize=28,  # Large clean labels like the reference
        pad=10,
    )

    # Explicitly set the x-ticks to match the reference labels exactly
    ax.set_xticks([0.5, 1.0, 1.5])

    # Remove the y-axis numerical labels completely
    ax.set_yticklabels([])

    # The reference image has no labels or legends visible on the canvas
    # Uncomment these if you need them later, but keep them commented to replicate the look:
    # ax.set_xlabel(r"$g=E_{\rm obs}/E_{\rm em}$", fontsize=24)
    # ax.legend(frameon=False)

    # Adjust bounding box to remove whitespace padding
    fig.savefig(outfile, dpi=250, bbox_inches="tight")
    print(f"Created {outfile}")


if __name__ == "__main__":
    main()
