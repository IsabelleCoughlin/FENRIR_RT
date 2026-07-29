import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import matplotlib.ticker as ticker

plt.rcParams.update({
    "font.family": "serif",
    "mathtext.fontset": "cm",
    "axes.linewidth": 2.5,

    "xtick.direction": "in",
    "ytick.direction": "in",

    "xtick.major.width": 2.5,
    "ytick.major.width": 2.5,
    "xtick.minor.width": 2.0,
    "ytick.minor.width": 2.0,

    "xtick.major.size": 9,
    "ytick.major.size": 9,
    "xtick.minor.size": 5,
    "ytick.minor.size": 5,

    "xtick.labelsize": 15,
    "ytick.labelsize": 15,
    "axes.labelsize": 10,
})

# python plot_emissivity_multi.py output_plot.png hist0.npy hist1.npy hist2.npy ...


def rIsco(a):
    z1 = 1.0 + (1.0 - a*a)**(1.0/3.0) * (
        (1.0 + a)**(1.0/3.0) + (1.0 - a)**(1.0/3.0)
    )
    z2 = (3.0*a*a + z1*z1)**0.5
    if a < 0:
        return 3.0 + z2 + ((3.0 - z1) * (3.0 + z1 + 2.0*z2))**0.5
    return 3.0 + z2 - ((3.0 - z1) * (3.0 + z1 + 2.0*z2))**0.5




outfile = sys.argv[1]
a = float(sys.argv[2])
#a = 0.99
infiles = sys.argv[3:]

# hist0 -> Mdot/Medd = 0.3
# hist1 -> Mdot/Medd = 0.2
# hist2 -> Mdot/Medd = 0.1
# hist3 -> Mdot/Medd = 0.0
labels = [
    r"$\dot{M}/\dot{M}_{Edd}=0.3$",
    r"$\dot{M}/\dot{M}_{Edd}=0.2$",
    r"$\dot{M}/\dot{M}_{Edd}=0.1$",
    r"$\dot{M}/\dot{M}_{Edd}=0.0$",
]


colors = ["blue", "red", "gray", "black"]

fig, ax = plt.subplots(figsize=(5, 5))
risco = rIsco(a)

for i, infile in enumerate(infiles):
    data = np.load(infile, allow_pickle=True)
    rho = np.array(data[0][0], dtype=float)
    flux = np.array(data[0][1], dtype=float)


    label = labels[i]
    color = colors[i]
    good = (rho >= risco)
    plt.loglog(rho[good], flux[good], lw=2, color=color, label=label)
    
# vertical line
ax.axvline(risco, color="green", lw=3.0)

ax.set_xlim(1, 50)


ax.set_ylim(1e-3, 1e6)

# Show labels only at 10^-3, 10^0, 10^3, 10^6
ax.yaxis.set_major_locator(
    ticker.FixedLocator([1e-3, 1e0, 1e3, 1e6])
)

ax.yaxis.set_major_formatter(
    ticker.LogFormatterMathtext(base=10)
)



ax.set_xlabel(r"$\rho\ [r_g]$", fontsize=10)
ax.set_ylabel("Irradiation Flux [Arbitrary Units]", fontsize=10)

ax.xaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=10))


ax.xaxis.set_minor_locator(
    ticker.LogLocator(base=10.0, subs=np.arange(2, 10), numticks=100)
)
ax.yaxis.set_minor_locator(
    ticker.LogLocator(base=10.0, subs=np.arange(2, 10))
) 

ax.xaxis.set_minor_formatter(ticker.NullFormatter())
ax.yaxis.set_minor_formatter(ticker.NullFormatter())

ax.tick_params(
    axis="both", which="major",
    direction="in",
    length=10,
    width=3,
    labelsize=10,
    top=False,
    right=False,
)

ax.tick_params(
    axis="both", which="minor",
    direction="in",
    length=6,
    width=2,
    top=False,
    right=False,
)

for spine in ax.spines.values():
    spine.set_linewidth(3.5)

fig.subplots_adjust(left=0.25, bottom=0.18, right=0.97, top=0.98)
fig.savefig(outfile, dpi=200)