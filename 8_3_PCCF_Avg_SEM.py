# Note: Run this on the epi-paint kernel only.

# Import dependencies

from matplotlib.colors import LinearSegmentedColormap, to_hex
import matplotlib.pyplot as _plt
import pandas as _pd
import numpy as _np
import os as _os
from matplotlib.ticker import ScalarFormatter

def load_and_average_pccf(parent_folder, file_name="pccf_results.csv"):
    files_list = []
    for root, dirs, files in _os.walk(parent_folder):
        for file in files:
            if file == file_name:
                files_list.append(_os.path.join(root, file))

    if len(files_list) == 0:
        raise FileNotFoundError(f"No files named {file_name} found inside {parent_folder}")

    print(f"Found {len(files_list)} PCCF files.")

    # Load all CSV files
    dfs = [ _pd.read_csv(f).dropna(how="all").fillna(0) for f in files_list ]

    # Use r-grid from first file
    r_common = dfs[0]["r"].values

    # Interpolate onto common r-grid
    aligned = []
    for df in dfs:
        df2 = _pd.DataFrame({"r": r_common})
        for col in df.columns:
            if col == "r":
                continue
            df2[col] = _np.interp(r_common, df["r"], df[col])
        aligned.append(df2)

    # Stack all
    stacked = _pd.concat(aligned, axis=0)

    # Mean and SD
    pccf_mean = stacked.groupby("r").mean().reset_index()
    pccf_std  = stacked.groupby("r").std().reset_index()

    pccf_sem = pccf_std.copy()
    for col in pccf_std.columns:
        if col == "r":
            continue
        pccf_sem[col] = pccf_std[col] / _np.sqrt(len(files_list))

    return pccf_mean, pccf_sem, len(files_list)


folder = '/Users/mickyanand/Library/CloudStorage/OneDrive-IndianInstituteofScience/ActD/DoC_NB_Rev/ActD5H/PCCF' # <<< Set your folder path here

# Load PCCF averaged data (Mean + SEM)
pccf_mean, pccf_sem, num_reps = load_and_average_pccf(folder)
print(pccf_sem.describe())

# Colormaps
cmap_proteins_white = {
    'S2P': LinearSegmentedColormap.from_list('S2P', ['#FFFFFF', '#FF0000']),
    'S5P': LinearSegmentedColormap.from_list('S5P', ['#FFFFFF', '#FFAA00']),
    'SC35': LinearSegmentedColormap.from_list('SC35', ['#FFFFFF', '#AAFF00']),
    'H3K4me3': LinearSegmentedColormap.from_list('H3K4me3', ['#FFFFFF', '#00FF00']),
    'H3K27ac': LinearSegmentedColormap.from_list('H3K27ac', ['#FFFFFF', '#00FFAA']),
    'CTCF': LinearSegmentedColormap.from_list('CTCF', ['#FFFFFF', '#00AAFF']),
    'H3K27me3': LinearSegmentedColormap.from_list('H3K27me3', ['#FFFFFF', '#0000FF']),
    'H3K9me3': LinearSegmentedColormap.from_list('H3K9me3', ['#FFFFFF', '#AA00FF']),
    'LaminB1': LinearSegmentedColormap.from_list('LaminB1', ['#FFFFFF', '#FF00AA']),
}

pairs = [c for c in pccf_mean.columns if c != "r"]
targets = sorted(set([p.split("_vs_")[0] for p in pairs] +
                     [p.split("_vs_")[1] for p in pairs]))

fig, axes = _plt.subplots(3, 3, figsize=(12, 12), sharex=True, sharey=True)

for ax, target in zip(axes.flat, targets):

    relevant_cols = [p for p in pairs if target in p]

    for col in relevant_cols:

        t1, t2 = col.split("_vs_")
        other = t2 if t1 == target else t1
        if 'ActD' in folder:
            if target == 'H3K4me3' or target == 'H3K27ac' or other == 'H3K4me3' or other == 'H3K27ac':
                continue        
        color = to_hex(cmap_proteins_white[other](1.0))

        # Line (Mean)
        ax.plot(
            pccf_mean["r"], pccf_mean[col],
            label=other, color=color, linewidth=2
        )

        # Shaded region = SEM
        ax.fill_between(
            pccf_mean["r"],
            pccf_mean[col] - pccf_sem[col],
            pccf_mean[col] + pccf_sem[col],
            color=color,
            alpha=0.25,
            linewidth=0
        )

    ax.axhline(1.0, color="gray", linestyle="--", linewidth=1)
    ax.set_title(target, fontsize=12)
    ax.legend(fontsize=7, loc="upper right", frameon=False)

# Hide unused axes
for ax in axes.flat[len(targets):]:
    ax.axis("off")

fig.suptitle(f"Averaged PCCF (Mean ± SEM, N={num_reps})", fontsize=16)
_plt.xlabel("Distance (nm)", fontsize=14)
_plt.ylim(0.5,1.5)
_plt.tight_layout()
_plt.savefig(_os.path.join(folder, "PCCF_avg_plot_SEM.png"), dpi=300)
_plt.savefig(_os.path.join(folder, "PCCF_avg_plot_SEM.svg"), format="svg")
_plt.show()