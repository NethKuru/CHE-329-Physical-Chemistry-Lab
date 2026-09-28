## Makes one IR spectrum graph for every CSV in NK_AJ_AW ##
import os
import matplotlib.pyplot as plt
import pandas as pd
from scipy.signal import find_peaks

def load_IR_csv(name):
    d = pd.read_csv(f"NK_AJ_AW/{name}.csv", skiprows=1)
    d.columns = ["wn", "T"]
    return d

def fit_labels(ax, labels, pad_pts=4):
    # Lower the bottom of the y axis until every peak label fits inside the plot.
    # Labels are a fixed size in points, so measure how far each one hangs below
    # its dip, then solve for the y-limit that leaves that much room.
    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    bottom, top = ax.get_ylim()
    ax_h = ax.get_window_extent(renderer).height          # axes height in pixels
    pad = pad_pts * fig.dpi / 72
    for a in labels:
        t = a.xy[1]
        hang = ax.transData.transform(a.xy)[1] - a.get_window_extent(renderer).y0 + pad
        f = hang / ax_h                                   # fraction of the axes the label needs
        bottom = min(bottom, (t - f * top) / (1 - f))
    ax.set_ylim(bottom, top)

FTIR = {
    "Ethanol Bottle_012": "Ethanol bottle", "Storage Bag": "Storage bag",
    "Salsa Container_010": "Salsa container", "nylon string": '"Nylon" string',
    "string_006": "String", "Zip Tie_014": "Zip tie", "Fishing Line_013": "Fishing line",
    "Water bottle_014": "Water bottle", "Heat Sealable Film_009": "Heat-sealable film",
    "Gloves": "Glove 1", "Glove 2": "Glove 2",
}

os.makedirs("IR_plots", exist_ok=True)

for name, label in FTIR.items():
    d = load_IR_csv(name)

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(d["wn"], d["T"], color="#2a78d6", lw=1.2)

    # IR peaks point down in %T, so look for peaks in -T
    # prominence = how deep a dip must be (10% of this spectrum's range)
    # distance = peaks closer than 50 cm-1 only keep the deeper one (stops labels overlapping)
    depth = 0.10 * (d["T"].max() - d["T"].min())
    peaks, _ = find_peaks(-d["T"], prominence=depth, distance=50)
    labels = []
    for p in peaks:
        wn, t = d["wn"][p], d["T"][p]
        labels.append(ax.annotate(f"{wn:.0f}", xy=(wn, t), xytext=(0, -3), textcoords="offset points",
                                  ha="center", va="top", rotation=90, fontsize=7))

    ax.set_xlim(4000, 400)                      # flipped x axis (high -> low wavenumber)
    ax.set_xlabel("Wavenumber (cm$^{-1}$)", fontsize=12)
    ax.set_ylabel("Transmittance (%)", fontsize=12)
    ax.set_title(f"{label} IR Spectrum", fontsize=14, pad=15)
    ax.spines[["top", "right"]].set_visible(False)

    plt.tight_layout()
    fit_labels(ax, labels)                      # make room under each dip for its label
    plt.savefig(f"IR_plots/{name}_IR.png", dpi=300)
    plt.close(fig)
    print(f"Saved {name}_IR.png")
