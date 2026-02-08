# Note: Run this on the picasso kernel only.
# This script calculates the Degree of Colocalization (DoC) between different protein localizations and employs a parallel processing approach to speed up computations.

# Import dependencies

import os as _os
import os.path as _ospath
import numpy as _np
import pandas as _pd
import h5py as _h5py
import yaml as _yaml
from PyQt5.QtWidgets import QMessageBox as _QMessageBox
import matplotlib.pyplot as _plt
from scipy.stats import spearmanr
from scipy.spatial import cKDTree
import itertools
from tqdm.contrib.concurrent import process_map
from concurrent.futures import ProcessPoolExecutor, as_completed
import time
from numba import njit
import epi_paint_picasso_utilis as eppu

# Define the folder location and the file extension inside the folder.

folder = ''  # <<< Set your folder path here
min_radius = 100
step_size = 100
max_radius = 1000
folder = _ospath.join(folder, 'Sectored')
file_extn = '.hdf5'
pixel_size = 130
file_names = [f for f in _os.listdir(folder) if f.endswith(file_extn)]

# Define the output folder.

parent_folder, working_folder = _ospath.split(folder)
output_parent_folder = _ospath.join(parent_folder, 'Analysis', 'Correlations')
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

# Define radii for analysis in nanometers
radii = _np.arange(min_radius, max_radius + step_size, step_size)
# print(f'For every point, there will be {len(radii)} radii calculated!')

# def gradient_density(tree, spot, radii):
#     # distances, _ = tree.query(spot, k=len(tree.data), distance_upper_bound=max_radius)
#     distances = tree.query_ball_point(spot, max_radius)
#     distances = distances[distances < max_radius]  # Filter valid distances
#     densities = _np.array([_np.sum(distances <= r) / (_np.pi * r**2) for r in radii])
#     return densities

def gradient_density(tree, spot, radii):
    """
    Compute cumulative neighbor density within each radius in `radii`
    using KDTree.query_ball_point for fast neighbor lookup.

    Parameters
    ----------
    tree : scipy.spatial.cKDTree built on the point set
    spot : (2,) array-like, the query point
    radii : (R,) ndarray of radii (ascending)

    Returns
    -------
    densities : (R,) float64 ndarray
    """
    # indices of neighbors within max radius
    idxs = tree.query_ball_point(spot, max_radius)
    if len(idxs) == 0:
        return _np.zeros_like(radii, dtype=_np.float64)

    # distances to neighbors found
    # NOTE: we need the raw coordinates the tree was built on; pass them in via closure or store.
    # Here we access tree.data (cKDTree stores the points in .data)
    dists = _np.linalg.norm(tree.data[idxs] - spot, axis=1)
    dists.sort()  # ascending

    # cumulative counts at each radius
    counts = _np.searchsorted(dists, radii, side='right')

    areas = _np.pi * (radii ** 2)
    densities = counts / areas
    return densities

# @njit
# def spearman_correlation(a, b):
#     n = len(a)
#     rank_a = _np.argsort(_np.argsort(a))
#     rank_b = _np.argsort(_np.argsort(b))
#     mean_a = _np.mean(rank_a)
#     mean_b = _np.mean(rank_b)
#     num = _np.sum((rank_a - mean_a) * (rank_b - mean_b))
#     den = _np.sqrt(_np.sum((rank_a - mean_a)**2) * _np.sum((rank_b - mean_b)**2))
#     return num / den if den > 0 else 0.0

@njit
def _spearman_corr_numba(a, b):
    n = a.size
    # ranks via argsort(argsort(.))
    order_a = _np.argsort(a)
    order_b = _np.argsort(b)

    ranks_a = _np.empty(n, dtype=_np.int64)
    ranks_b = _np.empty(n, dtype=_np.int64)
    for r, idx in enumerate(order_a):
        ranks_a[idx] = r
    for r, idx in enumerate(order_b):
        ranks_b[idx] = r

    mean_a = ranks_a.mean()
    mean_b = ranks_b.mean()

    num = 0.0
    den_a = 0.0
    den_b = 0.0
    for i in range(n):
        da = ranks_a[i] - mean_a
        db = ranks_b[i] - mean_b
        num += da * db
        den_a += da * da
        den_b += db * db

    den = _np.sqrt(den_a * den_b)
    if den == 0.0:
        return 0.0
    return num / den

# def doc_calc(tree_a, tree_b, data_a, radii, protein_a, protein_b, output_parent_folder):
#     print(f'Calculating DoC for {protein_a} vs {protein_b} with PID {_os.getpid()}')
#     doc_value = []
#     for idx, spot in enumerate(data_a):
#         density_1 = gradient_density(tree_a, spot, radii)
#         density_2 = gradient_density(tree_b, spot, radii)
#         # corr = spearmanr(density_1, density_2) [0]
#         corr = spearman_correlation(density_1, density_2)

#         if corr == _np.nan:
#             corr = 0
#         # norm_density_1 = density_1 / density_1[-1] # This is not necessary as spearman's corr is scale invariant
#         # norm_density_2 = density_2 / density_2[-1] # This is not necessary as spearman's corr is scale invariant
#         # corr = spearmanr(norm_density_1, norm_density_2) [0]
#         doc_value.append(corr)

#         if idx % 10000 == 0:
#             print(f'Processed {idx}/{len(data_a)} points for {protein_a} vs {protein_b} with PID {_os.getpid()}')

#     doc_value = _np.array(doc_value)
#     output_folder = _ospath.join(output_parent_folder, str(min_radius) + '_' + str(step_size) + '_' + str(max_radius))
#     if not _ospath.exists(output_folder):
#         _os.makedirs(output_folder)
#     _np.savetxt(_ospath.join(output_folder, f'{protein_a}_vs_{protein_b}.csv'), doc_value, delimiter=',', fmt='%.4f')

#     print(f'Finished calculating DoC for {protein_a} vs {protein_b}. PID {_os.getpid()}')

#     return doc_value

def doc_calc(tree_a, tree_b, data_a, radii, protein_a, protein_b, output_parent_folder):
    print(f'Calculating DoC for {protein_a} vs {protein_b} with PID {_os.getpid()}')
    doc_value = _np.empty(len(data_a), dtype=_np.float64)

    for idx, spot in enumerate(data_a):
        density_1 = gradient_density(tree_a, spot, radii)
        density_2 = gradient_density(tree_b, spot, radii)

        # numba-jitted spearman
        corr = _spearman_corr_numba(density_1, density_2)
        # safety: sanitize any residual NaN/inf (shouldn't happen with the guard above)
        if not _np.isfinite(corr):
            corr = 0.0

        doc_value[idx] = corr

        if idx % 10000 == 0 and idx > 0:
            print(f'Processed {idx}/{len(data_a)} points for {protein_a} vs {protein_b} with PID {_os.getpid()}')

    output_folder = _ospath.join(output_parent_folder, str(min_radius) + '_' + str(step_size) + '_' + str(max_radius))
    if not _ospath.exists(output_folder):
        _os.makedirs(output_folder)
    _np.savetxt(_ospath.join(output_folder, f'{protein_a}_vs_{protein_b}.csv'),
                doc_value, delimiter=',', fmt='%.4f')

    print(f'Finished calculating DoC for {protein_a} vs {protein_b}. PID {_os.getpid()}')
    return doc_value

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
    print(f'Processing pair: {protein_1} and {protein_2} with PID {_os.getpid()}')
    # print(f'Comparing {protein_1} and {protein_2}. ({counter}/{int(len(file_names)*(len(file_names)-1))/2:.0f})')
    locs_1, info_1 = eppu.load_locs(_ospath.join(folder, file_1))
    locs_2, info_2 = eppu.load_locs(_ospath.join(folder, file_2))
    data_1 = _np.column_stack((locs_1.x, locs_1.y))
    data_1 = data_1 * pixel_size
    data_2 = _np.column_stack((locs_2.x, locs_2.y))
    data_2 = data_2 * pixel_size
    tree_1 = cKDTree(data_1)
    tree_2 = cKDTree(data_2)
    doc_1_vs_2 = doc_calc(tree_1, tree_2, data_1, radii, protein_1, protein_2, output_parent_folder)
    doc_2_vs_1 = doc_calc(tree_2, tree_1, data_2, radii, protein_2, protein_1, output_parent_folder)

if __name__ == '__main__':
    start = time.time()
    pairs = list(itertools.combinations(enumerate(file_names), 2)) 
    
    with ProcessPoolExecutor() as executor:
        futures = [executor.submit(compare_proteins, pair) for pair in pairs]
        for f in as_completed(futures):
            try:
                result = f.result()  # will raise if worker errored
            except Exception as e:
                print("Worker failed:", e)

    end = time.time()
    print(f'Total time taken: {end - start} seconds')