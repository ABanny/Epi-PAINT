# Epi-PAINT Analysis
These codes are utilized to analyze the multiplexed nuclear imaging data using novel DNA-PAINT speed sequences.

## Contents
- `1_file_cleaner.ipynb`: This notebook helps clean the data by discarding long binding events. 
- `2_boundary_detection_script.py`: This script helps manually mark regions that are not part of the Lamin to build a mask around the nuclei.
- `3_boundary_polygon_builder.ipynb`: This notebook builds the boundary polygon mask using the AlphaShape library. 
- `4_cell_sectioner.ipynb`: This notebook sections the cells into defined sections and provides each localization with a unique section ID for further analysis. 
- `5_feature_extractor.ipynb`: This notebook extracts the number of localizations within each sector for each species. This is the coarse grained analysis that has been used. 
- `6_correlation_plotting.ipynb`: This notebook plots the correlation confusion matrices. 
- `SpatialCorrelation`: Folder containing codes that were used for spatial correlation analysis.
    - `10_spatial_correlation_calculator.ipynb`: This notebook calculates the DoC (Degree of Correlation) for each localization in a given channel with localizations from all the other channels. 
    - `11_1_spatial_correlation_plotter_vs_single.ipynb`: This notebook builds scatter plots to visualize the spatial distribuition and the realtive DoC values for each localization. 
    - `11_2_spatial_correlation_plotter_all_vs_all.ipynb`: This notebook builds a whole cell scatter plot for all vs all channels. 
    - `11_3_spatial_correlation_box_plot.ipynb`: This notebook builds box plots to showcase the percentage of localization above a below a certain threshold for DoC values. 
    - `11_4_spatial_correlation_all_vs_all_KDE.ipynb`: This notebook builds KDE plots for all the localizations.
- `README.md`: This file.

## How to get started?
- Install [Conda](https://docs.conda.io/) to build and maintain virtual environments. 
- Install [Picasso](https://github.com/jungmannlab/picasso) into a virtual environment as described on the github page. This is necessary reconstruction of raw data, drift correction, alignment and visual rendering of final data as and when applicable. 
- Create a new virtual environment for most of the python notebooks used here:
```
conda env create -f environment.yml
```
- To use `picasso` tools, activate the environement using:
```
conda activate picasso
```
- Follow `readme.rst` on [Picasso](https://github.com/jungmannlab/picasso) for indepth guides.
- Note: Pick the right kernel while running each notebook and/or script. 

## How to organize the data? 
- Store all *.hdf5 and *.yaml files that were generated from `picasso localize` and `picasso render` post undrifting and channel alignment in a parent folder. 
- Run the `1_file_cleaner.ipynb` by providing the directory to the folder in the notebook. Set `max_bright_time` and `max_dark_time` on the `picasso` kernel. 
- This would make a folder named `Cleaned` under the parent folder with all the cleaned files. 
- Open all the cleaned channels in a single `picasso render` window and mark individual nuclei with a rectangular pick and export all channels and the pick region from `File > Save Picked Localizations` and `File > Save Picked Regions` for each nuclei separately into separate folders. You can name it `Cell_1`, `Cell_2`, etc.
- Provide the new `Cell_N` folder path for all subsequent analysis except for `6_correlation_plotting.ipynb` where the parent folder name needs to be provided.
- Ensure nomenclature of files to contain `ProteinName_*.hdf5` and `ProteinName_*.yaml` for detection of protein species within code.

## Points to note:
- The correct kernel is crucial for proper code running. 
- File names should match section in the code where the file name is read to determine the nature and categorization of data. 
- Code can be adapted to support more protein targets. 

## Credits and Acknowledgements
- Sections of code have been adapted from [Picasso](https://github.com/jungmannlab/picasso) for easier implementation. 
- Dr. Venkatareddy Dadireddy for code building support and guidance. 