import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

outfile = sys.argv[1]
a = float(sys.argv[2])
files = sys.argv[3:]

def rIsco(a):
	z1 = 1 + (((1-(a*a))**(1./3.))*(((1+a)**(1./3.))+((1-a)**(1./3.))))
	z2 = ((3*(a*a))+(z1*z1))**0.5
	rOut = 3+z2-(((3-z1)*(3+z1+(2.*z2)))**0.5)
	return rOut


# plotting extras
titles = [
    r"$\dot{M}/\dot{M}_{Edd}=0.0$",
    r"$\dot{M}/\dot{M}_{Edd}=0.1$",
    r"$\dot{M}/\dot{M}_{Edd}=0.2$",
    r"$\dot{M}/\dot{M}_{Edd}=0.3$",
]


levels = [
    [8, 12, 16, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 72, 74],
    [8, 12, 16, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 72, 74],
    [10, 15, 20, 25, 30, 35],
    [10, 15, 20, 25, 30, 35],
]

# compute risco for limiting inner value of rho
rIsco_val = rIsco(a)

plt.rcParams.update({
    "font.family": "serif",
    "mathtext.fontset": "cm",
    "font.size": 11,
})

fig, axes = plt.subplots(
    2,
    2,
    figsize=(11, 10),
    constrained_layout=True,
)

payloads = []

for filename in files:
    data = np.load(filename)

    x = np.asarray(data[0], dtype=float)
    y = np.asarray(data[1], dtype=float)
    g = np.asarray(data[2], dtype=float)
    rho = np.asarray(data[8], dtype=float)

    # Only calculate log10 for positive g.
    logg = np.full_like(g, np.nan, dtype=float)
    positive_g = np.isfinite(g) & (g > 0)
    logg[positive_g] = np.log10(g[positive_g])

    # Colored emitting disk: rISCO <= rho <= 30 rg.
    color_good = (
        np.isfinite(x)
        & np.isfinite(y)
        & np.isfinite(rho)
        & np.isfinite(logg)
        & (rho >= rIsco_val)
        & (rho <= 30.0)
    )

    # Contours use all valid cylindrical radii, including rho > 30.
    contour_good = (
        np.isfinite(x)
        & np.isfinite(y)
        & np.isfinite(rho)
        & (rho > 0)
        & (rho <= 90.0)
    )

    payloads.append(
        (x, y, logg, rho, color_good, contour_good)
    )

# Input order is [0.3, 0.2, 0.1, 0.0].
payloads = [payloads[i] for i in [3, 2, 1, 0]]

# Contour spacing used to resemble the published figure.
levels = [
    np.arange(4.0, 73.0, 4.0),
    np.arange(4.0, 73.0, 4.0),
    np.arange(5.0, 81.0, 5.0),
    np.arange(5.0, 86.0, 5.0),
]

colorbar_ticks = np.arange(-0.30, 0.101, 0.05)

for ax, title, payload, contour_levels in zip(
    axes.flat,
    titles,
    payloads,
    levels,
):
    x, y, logg, rho, color_good, contour_good = payload

    scatter = ax.scatter(
        x[color_good],
        y[color_good],
        c=logg[color_good],
        s=1,
        cmap="gist_rainbow",
        vmin=-0.30,
        vmax=0.10,
        edgecolors="none",
        rasterized=True,
    )

    contours = ax.tricontour(
        x[contour_good],
        y[contour_good],
        rho[contour_good],
        levels=contour_levels,
        colors="0.2",
        linewidths=0.7,
    )

    ax.clabel(
        contours,
        inline=True,
        fontsize=7,
        fmt="%.3f",
    )

    ax.set_title(title)
    ax.set_aspect("equal", adjustable="box")

    ax.set_xlim(-35, 35)
    ax.set_ylim(-35, 35)

    ax.set_xticks(np.arange(-30, 31, 10))
    ax.set_yticks(np.arange(-30, 31, 10))

    ax.set_xlabel(r"$x\ (r_g)$")
    ax.set_ylabel(r"$y\ (r_g)$")

    ax.tick_params(
        direction="in",
        top=True,
        right=True,
        length=5,
    )

    colorbar = fig.colorbar(
        scatter,
        ax=ax,
        ticks=colorbar_ticks,
        fraction=0.046,
        pad=0.04,
    )

    colorbar.set_label(r"$\log_{10}(g)$")
    colorbar.ax.tick_params(direction="in")

fig.savefig(
    outfile,
    dpi=300,
    bbox_inches="tight",
)