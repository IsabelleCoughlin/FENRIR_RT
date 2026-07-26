# Create intensity plots like in Figure 5

import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

outfile = sys.argv[1]
disk_files = sys.argv[3:7]
emiss_files = sys.argv[7:11]

a = float(sys.argv[2])

def rIsco(a):
	"""Calculates the inner-most stable circular orbit (or rIsco, in gravitational
	radii Rg) for a given dimensionless spin value. |a| <= 1

	Example: rIsco(0) => 6.0 """
	z1 = 1 + (((1-(a*a))**(1./3.))*(((1+a)**(1./3.))+((1-a)**(1./3.))))
	z2 = ((3*(a*a))+(z1*z1))**0.5
	if (a < 0) and (a >= -1):
		rOut = 3+z2+(((3-z1)*(3+z1+(2.*z2)))**0.5)
	elif (a >= 0) and (a <= 1):
		rOut = 3+z2-(((3-z1)*(3+z1+(2.*z2)))**0.5)
	else:
		print("ERROR: Spin parameter must be within range of -1 to 1, inclusive.")
		sys.exit()
	return rOut



titles = [
    r"$\dot{M}/\dot{M}_{Edd}=0.0$",
    r"$\dot{M}/\dot{M}_{Edd}=0.1$",
    r"$\dot{M}/\dot{M}_{Edd}=0.2$",
    r"$\dot{M}/\dot{M}_{Edd}=0.3$",
]

levels = [
    np.arange(4, 73, 4),
    np.arange(4, 73, 4),
    np.arange(5, 81, 5),
    np.arange(5, 86, 5),
]

fig, axes = plt.subplots(2, 2, figsize=(9, 8), constrained_layout=True)

all_intensities = []
rIsco_val = rIsco(a)

payloads = []
for disk_file, emiss_file in zip(disk_files, emiss_files):
    
    # the disk image tracing photons from disk to observer
    disk = np.load(disk_file)
    # emmisivity profiles from corona to disk
    em = np.load(emiss_file, allow_pickle=True)

    # values from disk redshift and position
    x = np.array(disk[0], dtype=float)
    y = np.array(disk[1], dtype=float)
    g = np.array(disk[2], dtype=float)
    rho = np.array(disk[8], dtype=float)

    em_rho = np.array(em[0][0], dtype=float)
    em_flux = np.array(em[0][1], dtype=float)

    valid_em = np.isfinite(em_rho) & np.isfinite(em_flux) & (em_rho > 0) & (em_flux >= 0)
    em_rho = em_rho[valid_em]
    em_flux = em_flux[valid_em]

    order = np.argsort(em_rho)
    em_rho = em_rho[order]
    em_flux = em_flux[order]
    #emiss = np.interp(rho, em_rho, em_flux, left=0.0, right=0.0)
    emiss = np.interp(
        rho,
        em_rho,
        em_flux,
        left=em_flux[0],
        right=0.0,
    )
    
    intensity = emiss * g**4
    # set anything at or below zero to stay at zero
    intensity = np.where(np.isfinite(intensity) & (intensity > 0), intensity, 0.0)

    good = (
        np.isfinite(x)
        & np.isfinite(y)
        & np.isfinite(g)
        & np.isfinite(rho)
        & np.isfinite(intensity)
        & (rho >= rIsco_val)
        & (rho <= 30.0)
        & (g > 0)
    )

    contour_good = (
        np.isfinite(x)
        & np.isfinite(y)
        & np.isfinite(rho)
        & (rho > 0)
    )
        
    payloads.append((x, y, rho, intensity, good, contour_good))
    all_intensities.append(intensity[good])

order = [3, 2, 1, 0]
payloads = [payloads[i] for i in order]

all_i = np.concatenate(all_intensities)

# scale to match plot
target_max = 10**3.6 - 1.0
EMISS_SCALE = target_max / np.nanmax(all_i)


last_sc = None

for ax, title, payload, level in zip(
    axes.flat, titles, payloads, levels
):
    x, y, rho, intensity, intensity_good, contour_good = payload

    intensity_scaled = EMISS_SCALE * intensity

    logI = np.full_like(intensity_scaled, np.nan, dtype=float)
    logI[intensity_good] = np.log10(
        intensity_scaled[intensity_good] + 1.0
    )

    sc = ax.scatter(
        x[intensity_good],
        y[intensity_good],
        c=logI[intensity_good],
        s=1,
        cmap="gnuplot2",
        vmin=0.0,
        vmax=3.6,
        rasterized=True,
    )

    # Use all valid disk-radius data for contours, including rho > 30.
    cs = ax.tricontour(
        x[contour_good],
        y[contour_good],
        rho[contour_good],
        levels=level,
        colors="lawngreen",
        linewidths=0.6,
        alpha=0.75,
    )

    ax.clabel(
        cs,
        inline=True,
        fontsize=7,
        fmt=lambda value: rf"{value:.0f} $r_g$",
    )

    ax.set_title(title)
    ax.set_aspect("equal")
    ax.set_xlim(-35, 35)
    ax.set_ylim(-35, 35)
    ax.set_xlabel(r"$x\ [r_g]$")
    ax.set_ylabel(r"$y\ [r_g]$")

    # The paper has one colorbar per panel.
    colorbar = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
    colorbar.set_label(
        r"$\log_{10}[I_{\rm obs}(x,y)+1]$"
    )
#fig.colorbar(last_sc, ax=axes.ravel().tolist(), label=r"$\log_{10}(I_{\rm obs}(x,y)+1)$")
fig.savefig(outfile, dpi=250)