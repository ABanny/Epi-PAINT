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
import epi_paint_picasso_utilis as eppu
from concurrent.futures import ProcessPoolExecutor, as_completed

# Define the folder location and the file extension inside the folder.
folder = ''  # <<< Set your folder path here
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
file_names = sorted(file_names, key = lambda x: order.index(x.split('_')[0]))

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
    for p in points_A:
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

def compare_proteins_parallel(args):
    pair, data_1, data_2, r_max, dr, polygon_coords = args
    (i, file_1), (j, file_2) = pair
    protein_1 = file_1.split('_')[0]
    protein_2 = file_2.split('_')[0]
    
    print(f'Processing pair: {protein_1} and {protein_2} with PID {_os.getpid()}')
    
    polygon = Polygon(polygon_coords)
    _, g_AB = compute_pccf(data_1, data_2, r_max=r_max, dr=dr, roi=polygon)
    _, g_BA = compute_pccf(data_2, data_1, r_max=r_max, dr=dr, roi=polygon)
    
    print(f'Finished pair: {protein_1} and {protein_2} with PID {_os.getpid()}')
    
    return i, j, 0.5 * (g_AB + g_BA)


if __name__ == '__main__':
    start = time.time()
    _os.system('clear')
    print('Loading all localiztions...')
    all_data = {}
    for file in tqdm(file_names):
        locs, _ = eppu.load_locs(_ospath.join(folder, file))
        data = _np.column_stack((locs.x, locs.y)) * pixel_size
        all_data[file] = data
    print('All localizations loaded.')
    pairs = list(itertools.combinations(enumerate(file_names), 2))
    polygon_coords = list(polygon.exterior.coords)
    total_pairs = len(pairs)
    
    print(f'Total pairs to process PCCFs: {total_pairs}')
    
    task_args = [(pair, all_data[pair[0][1]], all_data[pair[1][1]], r_max, dr, polygon_coords) for pair in pairs]

    pccf = {}
    completed = 0
    with ProcessPoolExecutor() as executor:
        futures = {executor.submit(compare_proteins_parallel, arg): pair
                   for arg, pair in zip(task_args, pairs)}
        for future in as_completed(futures):
            pair = futures[future]
            try:
                i, j, g_cross = future.result()
                pccf[pair] = g_cross
                completed += 1
                print(f'Completed {completed}/{total_pairs}: {pair[0][1].split("_")[0]} vs {pair[1][1].split("_")[0]}')
            except Exception as e:
                print(f'Failed: {pair[0][1].split("_")[0]} vs {pair[1][1].split("_")[0]}: {e}')

    headers = ['r'] + [f'{pair[0][1].split("_")[0]}_vs_{pair[1][1].split("_")[0]}' for pair in pairs]
    edges = _np.arange(0, r_max + dr, dr)
    r = 0.5 * (edges[:-1] + edges[1:])
    out = _np.column_stack([r] + [pccf[pair] for pair in pairs])

    output_file = _ospath.join(output_parent_folder, 'pccf_results_parallel.csv')
    _np.savetxt(output_file, out, delimiter=',', header=','.join(headers), comments='')

    end = time.time()
    print(f'Total time taken to calculate PCCFs: {end - start:.2f} seconds')