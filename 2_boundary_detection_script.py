# Note: Run this on the picasso kernel only.

# This script allows manual boundary detection and removal of localizations.
# It creates a density map, applies thresholding, and provides an interactive
# interface to select and deselect regions. The cleaned localizations are saved
# to a new file with a '_boundary' suffix.

import numpy as np
from scipy.ndimage import gaussian_filter
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.widgets import Button, PolygonSelector
from matplotlib.path import Path
from picasso import io as _io
import os as _os
import os.path as _ospath
import epi_paint_picasso_utilis as eppu
import sys

folder = '' # <<< Set your folder path here
folder = _ospath.join(folder, 'Cleaned')
file_extn = '.hdf5'
file_names = [f for f in _os.listdir(folder) if f.endswith(file_extn)]
boundary_files = [f for f in file_names if 'boundary' in f]
if boundary_files:
    eppu.show_boundary_error(boundary_files)
    sys.exit()

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
        self.history = []

        self.fig = plt.figure(figsize=(10, 6))
        gs = self.fig.add_gridspec(
            1, 2,
            width_ratios=[4, 1],
            wspace=0.05
        )

        self.ax = self.fig.add_subplot(gs[0])
        self.ax_info = self.fig.add_subplot(gs[1])
        self.ax_info.axis('off')

        plt.subplots_adjust(bottom=0.2)
        self.im = self.ax.imshow(binary_mask, origin='lower', cmap='gray_r')

        instructions = (
        "Polygon:\n"
        "  • Click to add points.\n"
        "  • Connect back to starting vertex.\n\n"
        "Rectangle:\n"
        "  • Shift + drag to draw rectangles.\n\n"
        "Undo:\n"
        "  • If you made a mistake, UNDO it.\n\n"
        "Finish:\n"
        "  • Use the DONE button to finish.\n\n"
        )

        self.ax_info.text(
        0.0, 1.0, instructions,
        va='top', fontsize=8
        )

        # ---------- Buttons ----------
        ax_done = plt.axes([0.8, 0.05, 0.15, 0.07])
        ax_undo = plt.axes([0.6, 0.05, 0.15, 0.07])

        self.btn_done = Button(ax_done, 'Done')
        self.btn_undo = Button(ax_undo, 'Undo')

        self.btn_done.on_clicked(self.on_done)
        self.btn_undo.on_clicked(self.undo)

        # ---------- Polygon selector (default) ----------
        self.poly = PolygonSelector(
        self.ax,
        self.on_polygon_complete,
        useblit=True
        )

        # ---------- Rectangle state ----------
        self.rect_start = None
        self.rect_patch = None

        # ---------- Events ----------
        self.fig.canvas.mpl_connect('button_press_event', self.on_press)
        self.fig.canvas.mpl_connect('button_release_event', self.on_release)
        self.fig.canvas.mpl_connect('motion_notify_event', self.on_motion)
        self.fig.canvas.mpl_connect('key_press_event', self.on_key)
        self.fig.canvas.mpl_connect('close_event', self.on_close)

        self.finished = False

        plt.xlim(-25, int(round(xrange)) + 25)
        plt.ylim(-25, int(round(yrange)) + 25)
        plt.show()

    # ================= RECTANGLE (Shift + Drag) =================

    def on_press(self, event):
        if event.inaxes != self.ax or not event.key == 'shift':
            return

        self.rect_start = (int(event.xdata), int(event.ydata))
        self.rect_patch = Rectangle(
            self.rect_start, 0, 0,
            linewidth=1, edgecolor='red', facecolor='none'
        )
        self.ax.add_patch(self.rect_patch)

        # Disable polygon temporarily
        self.poly.set_active(False)

    def on_motion(self, event):
        if self.rect_start is None or event.inaxes != self.ax:
            return

        width = int(event.xdata) - self.rect_start[0]
        height = int(event.ydata) - self.rect_start[1]
        self.rect_patch.set_width(width)
        self.rect_patch.set_height(height)
        self.fig.canvas.draw_idle()

    def on_release(self, event):
        if self.rect_start is None or event.inaxes != self.ax:
            return

        self.history.append(self.binary_mask.copy())

        x0, y0 = self.rect_start
        x1, y1 = int(event.xdata), int(event.ydata)
        x_min, x_max = sorted((x0, x1))
        y_min, y_max = sorted((y0, y1))

        self.binary_mask[y_min:y_max, x_min:x_max] = False
        self.im.set_data(self.binary_mask)

        # Cleanup
        self.rect_start = None
        self.rect_patch.remove()
        self.rect_patch = None

        # Re-enable polygon
        self.poly.set_active(True)
        self.fig.canvas.draw_idle()

    # ================= POLYGON =================

    def on_polygon_complete(self, verts):
        if len(verts) < 3:
            return

        # Save history for undo
        self.history.append(self.binary_mask.copy())

        ny, nx = self.binary_mask.shape
        X, Y = np.meshgrid(np.arange(nx), np.arange(ny))
        points = np.vstack((X.ravel(), Y.ravel())).T

        path = Path(verts)
        mask = path.contains_points(points).reshape(ny, nx)

        # Remove selected region
        self.binary_mask[mask] = False
        self.im.set_data(self.binary_mask)
        self.fig.canvas.draw_idle()

        # 🔥 IMPORTANT PART 🔥
        # Reset selector so another polygon can be drawn
        self.poly.clear()
        self.poly.set_active(True)

    # ================= UNDO / DONE =================

    def undo(self, event):
        if not self.history:
            print("Nothing to undo")
            return
        self.binary_mask[:] = self.history.pop()
        self.im.set_data(self.binary_mask)
        self.fig.canvas.draw_idle()
        print("Undo")

    def on_key(self, event):
        if event.key == 'enter':
            self.on_done(None)

    def on_done(self, event):
        print("Selection finished.")
        self.finished = True
        plt.close(self.fig)

    def on_close(self, event):
        if not self.finished:
            print("Window closed without clicking Done. Changes discarded.")

# Run the Editor
editor = MaskEditor(binary_mask.T)

# Save the modified mask
if not editor.finished:
    eppu.show_boundary_not_saved_error()
    sys.exit()

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
eppu.save_locs_withSuffix(fpath, retained_locs, info, suffix='boundary')
print('Saved the boundary cleaned file')