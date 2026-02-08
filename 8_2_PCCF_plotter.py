# Note: Run this on the epi-paint kernel only.

# Import dependencies
 
from matplotlib.colors import LinearSegmentedColormap, to_hex
import matplotlib.pyplot as _plt
import pandas as _pd
import os as _os
from matplotlib.ticker import ScalarFormatter

folder = ''  # <<< Set your folder path here
folder = _os.path.join(folder, 'Analysis', 'PCCF')
file_extn = '.csv'
file_names = [f for f in _os.listdir(folder) if f.endswith(file_extn)]
pccf = _pd.read_csv(_os.path.join(folder, file_names[0]), delimiter=',')
print(pccf.head())

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

pairs = [c for c in pccf.columns if c != "r"]
targets = sorted(set([p.split("_vs_")[0] for p in pairs] + 
                     [p.split("_vs_")[1] for p in pairs]))

fig, axes = _plt.subplots(3, 3, figsize=(12, 12), sharex=True, sharey=True)

for ax, target in zip(axes.flat, targets):
    relevant_cols = [p for p in pairs if target in p]
    
    for col in relevant_cols:
        # Find the "other" protein in this pair
        t1, t2 = col.split("_vs_")
        other = t2 if t1 == target else t1
        
        # Pick a strong color from the colormap of the "other" protein
        color = to_hex(cmap_proteins_white[other](1.0))

        ax.plot(pccf["r"], pccf[col], label=other, color=color, linewidth=2)

    ax.axhline(1.0, color="gray", linestyle="--", linewidth=1)
    ax.set_title(target, fontsize=12)
    ax.legend(fontsize=7, loc="upper right", frameon=False)
    # ax.set_yscale('log')

    # # Use ScalarFormatter on log scale
    # formatter = ScalarFormatter()
    # formatter.set_scientific(False)   # disable scientific notation
    # formatter.set_useOffset(False)
    # formatter.set_powerlimits((0, 0))
    # ax.yaxis.set_major_formatter(formatter)


# Remove empty axes (if less than 9 targets)
for ax in axes.flat[len(targets):]:
    ax.axis("off")

fig.suptitle("Pair Cross-Correlation Functions per Target", fontsize=16)
_plt.xlabel("Distance (nm)", fontsize=14)
# _plt.yscale('log')
_plt.tight_layout()
_plt.savefig(_os.path.join(folder, 'PCCF_plot.png'), dpi=300)
_plt.savefig(_os.path.join(folder, 'PCCF_plot.svg'), format = 'svg')
_plt.show()