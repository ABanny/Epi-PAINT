import numpy as _np
import os.path as _ospath
import h5py as _h5py
import yaml as _yaml
from PyQt5.QtWidgets import QApplication, QMessageBox
import sys
import pandas as _pd

### Input and Output Functions adapted from Picasso ###

def ensure_numpy_structured(locs):
    if isinstance(locs, _pd.DataFrame):
        return locs.to_records(index=False)
    return locs

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
            QMessageBox.critical(
                qt_parent,
                "An error occured",
                "Could not find metadata file:\n{}".format(filename),
            )
        raise NoMetadataFileError(e)
    return info

def save_info(path, info, default_flow_style=False):
    with open(path, "w") as file:
        _yaml.dump_all(info, file, default_flow_style=default_flow_style)

def ensure_sanity(locs, info):
    locs = ensure_numpy_structured(locs)
    names = locs.dtype.names

    # ---- finite check only on numeric fields ----
    numeric_fields = [
        name for name in names
        if _np.issubdtype(locs[name].dtype, _np.number)
    ]

    if numeric_fields:
        finite_mask = _np.all(
            _np.vstack([_np.isfinite(locs[name]) for name in numeric_fields]),
            axis=0
        )
        locs = locs[finite_mask]

    # ---- required fields ----
    required = ['x', 'y', 'lpx', 'lpy']
    for field in required:
        if field not in names:
            raise KeyError(
                f"ensure_sanity: required field '{field}' missing. "
                f"Available fields: {names}"
            )

    # ---- spatial sanity checks ----
    locs = locs[locs['x'] > 0]
    locs = locs[locs['y'] > 0]
    locs = locs[locs['x'] < info[0]["Width"]]
    locs = locs[locs['y'] < info[0]["Height"]]
    locs = locs[locs['lpx'] > 0]
    locs = locs[locs['lpy'] > 0]

    return locs

def save_locs_withSuffix(path, locs, info, suffix=''):

    locs = ensure_sanity(locs, info)

    # ensure DataFrame (required by ensure_sanity)
    if not isinstance(locs, _pd.DataFrame):
        locs = _pd.DataFrame.from_records(locs)

    # convert back to structured array for saving
    locs_np = locs.to_records(index=False)

    base, ext_locs = _ospath.splitext(path)
    output_locs_path = base + '_' + suffix + ext_locs    
    output_info_path = base + '_' + suffix + '.yaml'

    with _h5py.File(output_locs_path, "w") as locs_file:
        locs_file.create_dataset("locs", data=locs_np)

    save_info(output_info_path, info, default_flow_style=False)

### Error Message Functions ###

def show_boundary_error(boundary_files):
    app = QApplication.instance()
    if app is None:
        app = QApplication(sys.argv)

    msg = QMessageBox()
    msg.setIcon(QMessageBox.Critical)
    msg.setWindowTitle("Boundary Already Exists")
    msg.setText("Boundary file(s) detected.")
    msg.setInformativeText(
        "Please remove or move the existing boundary file(s) before proceeding:\n\n"
        + "\n".join(boundary_files)
    )
    msg.setStandardButtons(QMessageBox.Ok)

    msg.exec_()

def show_boundary_not_saved_error():
    app = QApplication.instance()
    if app is None:
        app = QApplication(sys.argv)

    msg = QMessageBox()
    msg.setIcon(QMessageBox.Critical)
    msg.setWindowTitle("Boundary Not Saved")
    msg.setText("Boundary changes were not saved.")
    msg.setInformativeText(
        "You closed the boundary editor without clicking Done. Please run the script again and click Done to save your changes."
    )
    msg.setStandardButtons(QMessageBox.Ok)

    msg.exec_()

### Code Separator ###

if __name__ == "__main__":
    pass