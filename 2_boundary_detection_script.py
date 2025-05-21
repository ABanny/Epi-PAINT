# Note: Run this on the picasso kernel only.

import numpy as np
from scipy.ndimage import gaussian_filter
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from picasso import lib as _lib
from picasso import io as _io
import os as _os
import os.path as _ospath
import h5py as _h5py

def save_locs_withSuffix(path, locs, info, suffix=''):
    locs = _lib.ensure_sanity(locs, info)
    base, ext_locs = _ospath.splitext(path)
    output_locs_path = base + '_' + suffix + ext_locs    
    output_info_path = base + '_' + suffix + '.yaml'
    with _h5py.File(output_locs_path, "w") as locs_file:
        locs_file.create_dataset("locs", data=locs)
    _io.save_info(output_info_path, info, default_flow_style=False)


#Define folder
folder = '' # Insert the path to the Cell Folder. 
file_extn = '.hdf5'
file_names = [f for f in _os.listdir(folder) if f.endswith(file_extn)]
for file in file_names:
    if 'Lamin' in file:
        fpath = _ospath.join(folder, file)
        break

print(f'Loaded file {file}.')

locs, info = _io.load_locs(fpath)

# Example coordinates
x = locs['x']
y = locs['y']
xrange = locs['x'].max() - locs['x'].min()
yrange = locs['y'].max() - locs['y'].min()

# Create a 2D histogram (density map)
H, xedges, yedges = np.histogram2d(x, y, bins=(int(round(xrange)), int(round(yrange))))

# Smooth the density map
H_smooth = gaussian_filter(H, sigma=0.5)
plt.xlim(-25, int(round(xrange))+25)
plt.ylim(-25, int(round(yrange))+25)
plt.imshow(H_smooth.T, origin='lower', cmap='hot')
plt.colorbar()
plt.show()

# Thresholding
threshold = np.percentile(H_smooth, 95)
binary_mask = H_smooth > threshold


class MaskEditor:
    def __init__(self, binary_mask):
        self.binary_mask = binary_mask
        self.fig, self.ax = plt.subplots()
        self.im = self.ax.imshow(binary_mask, origin='lower', cmap='gray')
        self.rect = None
        self.start = None
        self.fig.canvas.mpl_connect('button_press_event', self.on_press)
        self.fig.canvas.mpl_connect('button_release_event', self.on_release)
        self.fig.canvas.mpl_connect('motion_notify_event', self.on_motion)
        plt.xlim(-25, int(round(xrange))+25)
        plt.ylim(-25, int(round(yrange))+25)
        plt.show()

    def on_press(self, event):
        if event.inaxes != self.ax:
            return
        self.start = (int(event.xdata), int(event.ydata))
        self.rect = Rectangle(self.start, 0, 0, linewidth=1, edgecolor='red', facecolor='none')
        self.ax.add_patch(self.rect)

    def on_motion(self, event):
        if self.start is None or event.inaxes != self.ax:
            return
        width = int(event.xdata) - self.start[0]
        height = int(event.ydata) - self.start[1]
        self.rect.set_width(width)
        self.rect.set_height(height)
        self.fig.canvas.draw()

    def on_release(self, event):
        if self.start is None or event.inaxes != self.ax:
            return
        x0, y0 = self.start
        x1, y1 = int(event.xdata), int(event.ydata)
        x_min, x_max = sorted((x0, x1))
        y_min, y_max = sorted((y0, y1))
        self.binary_mask[y_min:y_max, x_min:x_max] = False  # Deselect pixels
        self.im.set_data(self.binary_mask)
        self.start = None
        self.rect = None
        self.fig.canvas.draw()

# Run the Editor
editor = MaskEditor(binary_mask.T)

# Save the modified mask
cleaned_mask = editor.binary_mask

# Map coordinates to mask grid
x_indices = np.digitize(x, xedges) - 1
y_indices = np.digitize(y, yedges) - 1

# Clip indices to ensure they fall within the mask's bounds
x_indices = np.clip(x_indices, 0, cleaned_mask.shape[1] - 1)
y_indices = np.clip(y_indices, 0, cleaned_mask.shape[0] - 1)

# Extract retained localizations
mask_indices = cleaned_mask[y_indices, x_indices]
retained_locs = locs[mask_indices]

# Save retained localizations to a new file
save_locs_withSuffix(fpath, retained_locs, info, suffix='boundary')
print('Saved the boundary cleaned file')