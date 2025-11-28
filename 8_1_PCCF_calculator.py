# Note: Run this on the epi-paint kernel only.

# Import dependencies

import os as _os
import os.path as _ospath
import numpy as _np
import h5py as _h5py
import yaml as _yaml
from shapely.geometry import Point, Polygon
from PyQt5.QtWidgets import QMessageBox as _QMessageBox
import matplotlib.pyplot as _plt
from scipy.spatial import cKDTree
import itertools
import time
from tqdm import tqdm

# Define the folder location and the file extension inside the folder.
folder = '/Users/abhinav/Library/CloudStorage/OneDrive-IndianInstituteofScience/AnalysisFolder/Epi/20251112/Mock/120825/cell5'  # <<< Set your folder path here
folder = _ospath.join(folder, 'Masked')
file_extn = '.hdf5'
pixel_size = 130
file_names = [f for f in _os.listdir(folder) if f.endswith(file_extn)]

# Define where the polygon file is saved.
polygon = _np.loadtxt(_os.path.join(folder, 'polygon_coords.csv'), delimiter = ',')
polygon = polygon * pixel_size
polygon = Polygon(polygon)
centroid = polygon.centroid

# Define parameters for the PCCF calculation.
r_max = 1000
dr = 10


# Define the output folder.
parent_folder, working_folder = _ospath.split(folder)
output_parent_folder = _ospath.join(parent_folder, 'Analysis', 'PCCF')
if not _ospath.exists(output_parent_folder):
    _os.makedirs(output_parent_folder)

# Define the list of proteins from the file names.
list_of_proteins = []

for file in file_names:
    protein = file.split('_')[0]
    if protein not in list_of_proteins:
        list_of_proteins.append(protein)

# Order the list of proteins
order = ['S2P', 'S5P', 'SC35', 'H3K4me3', 'H3K27ac', 'CTCF', 'H3K27me3', 'H3K9me3', 'LaminB1']
list_of_proteins = sorted(list_of_proteins, key = order.index)

def load_locs(path, qt_parent=None):
    with _h5py.File(path, "r") as locs_file:
        locs = locs_file["locs"][...]
    locs = _np.rec.array(
        locs, dtype=locs.dtype
    )  # Convert to rec array with fields as attributes
    info = load_info(path, qt_parent=qt_parent)
    return locs, info

class NoMetadataFileError(FileNotFoundError):
    pass

def load_info(path, qt_parent=None):
    path_base, path_extension = _ospath.splitext(path)
    filename = path_base + ".yaml"
    try:
        with open(filename, "r") as info_file:
            info = list(_yaml.load_all(info_file, Loader=_yaml.UnsafeLoader))
    except FileNotFoundError as e:
        print("\nAn error occured. Could not find metadata file:\n{}".format(filename))
        if qt_parent is not None:
            _QMessageBox.critical(
                qt_parent,
                "An error occured",
                "Could not find metadata file:\n{}".format(filename),
            )
        raise NoMetadataFileError(e)
    return info

def compute_pccf(points_A, points_B, r_max, dr, roi = None):
    points_A = _np.asarray(points_A)
    points_B = _np.asarray(points_B)

    # Build KDTree for fast distance queries
    tree_B = cKDTree(points_B)

    # Distance bins
    edges = _np.arange(0, r_max + dr, dr)
    r = 0.5 * (edges[:-1] + edges[1:])
    counts = _np.zeros_like(r)

    # Count actual pairs within each shell
    for p in tqdm(points_A, desc="Counting for points: "):
        idxs = tree_B.query_ball_point(p, r_max)
        if len(idxs) == 0:
            continue
        dists = _np.linalg.norm(points_B[idxs] - p, axis=1)
        hist, _ = _np.histogram(dists, bins=edges)
        counts += hist

    # Normalization
    area_total = polygon.area if polygon is not None else (
        (max(points_A[:,0].max(), points_B[:,0].max()) -
         min(points_A[:,0].min(), points_B[:,0].min())) *
        (max(points_A[:,1].max(), points_B[:,1].max()) -
         min(points_A[:,1].min(), points_B[:,1].min()))
    )

    density_B = len(points_B) / area_total
    norm = 2 * _np.pi * r * dr * density_B * len(points_A)
    g_AB = counts / norm
    return r, g_AB

def plot_doc_map(data, doc_values, protein_a, protein_b):
    _plt.figure(figsize=(8, 6))
    _plt.scatter(data[:, 0], data[:, 1], c=doc_values, cmap='coolwarm', marker='o', s=1)
    _plt.colorbar(label=f'DoC Score {protein_a} vs {protein_b}')
    _plt.title(f'DoC Map for {protein_a} vs {protein_b}')
    _plt.gca().invert_yaxis()
    _plt.show()

def compare_proteins(pair_args):
    (i, file_1), (j, file_2) = pair_args
    protein_1 = file_1.split('_')[0]
    protein_2 = file_2.split('_')[0]
    locs_1, info_1 = load_locs(_ospath.join(folder, file_1))
    locs_2, info_2 = load_locs(_ospath.join(folder, file_2))
    data_1 = _np.column_stack((locs_1.x, locs_1.y))
    data_1 = data_1 * pixel_size
    data_2 = _np.column_stack((locs_2.x, locs_2.y))
    data_2 = data_2 * pixel_size
    _, g_AB = compute_pccf(data_1, data_2, r_max=r_max, dr=dr, roi=polygon)
    _, g_BA = compute_pccf(data_2, data_1, r_max=r_max, dr=dr, roi=polygon)
    g_cross = 0.5 * (g_AB + g_BA)
    return g_cross


if __name__ == '__main__':
    _os.system('clear')
    start = time.time()
    pccf = {}
    pairs = list(itertools.combinations(enumerate(file_names), 2))
    # print(f'Total pairs to process: {len(pairs)}')
    pair_count = 1 
    for pair in pairs:
        print(f'Processing pair {pair_count}/{len(pairs)}: {pair[0][1].split("_")[0]} and {pair[1][1].split("_")[0]}')
        g_cross = compare_proteins(pair)
        pccf[pair] = g_cross
        print(f'Completed pair: {pair[0][1].split("_")[0]} and {pair[1][1].split("_")[0]}')
        pair_count += 1
        _os.system('clear')
    headers = ['r'] + [f'{pair[0][1].split("_")[0]}_vs_{pair[1][1].split("_")[0]}' for pair in pairs]
    edges = _np.arange(0, r_max + dr, dr)
    r = 0.5 * (edges[:-1] + edges[1:])
    # Stack results into one array
    out = [r]
    for pair in pairs:
        out.append(pccf[pair])
    out = _np.column_stack(out)

    # Save the results
    output_file = _ospath.join(output_parent_folder, 'pccf_results.csv')
    _np.savetxt(output_file, out, delimiter=',', header=','.join(headers), comments='')

    # with ProcessPoolExecutor() as executor:
    #     futures = [executor.submit(compare_proteins, pair) for pair in pairs]
    #     for f in as_completed(futures):
    #         try:
    #             result = f.result()  # will raise if worker errored
    #         except Exception as e:
    #             print("Worker failed:", e)

    end = time.time()
    print(f'Total time taken: {end - start:.2f} seconds')