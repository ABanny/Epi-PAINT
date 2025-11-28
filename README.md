# Epi-PAINT Analysis
These codes are utilized to analyze the multiplexed nuclear imaging data and origami data using novel DNA-PAINT speed sequences.

The pre-print is available [here](https://doi.org/10.1101/2024.12.21.629871)

## Contents
- `1_file_cleaner.ipynb`: This notebook helps clean the data by discarding long binding events. 
- `2_boundary_detection_script.py`: This script helps manually mark regions that are not part of the Lamin to build a mask around the nuclei.
- `3_boundary_polygon_builder.ipynb`: This notebook builds the boundary polygon mask using the AlphaShape library. 
- `4_cell_sectioner.ipynb`: This notebook sections the cells into defined sections and provides each localization with a unique section ID for further analysis. 
- `5_feature_extractor.ipynb`: This notebook extracts the number of localizations within each sector for each species. This is the coarse grained analysis that has been used. 
- `6_correlation_plotting.ipynb`: This notebook plots the correlation confusion matrices. 
- `7_1_SpatialCorrelation`: This script calculates the DoC (Degree of Correlation) for each localization in a given channel with localizations from all the other channels. This script has been parallelized for faster calculations.    
- `7_2_spatial_correlation_plotter_vs_single.ipynb`: This notebook builds scatter plots to visualize the spatial distribuition and the realtive DoC values for each localization. 
- `7_3_spatial_correlation_plotter_all_vs_all_KDE.ipynb`: This notebook builds KDE plots for all the localizations.
- `8_1_PCCF_calculator.py`: This script calculates PCCF for all localizations and stores the file in the appropriate directory. 
- `8_2_PCCF_plotter.py`: This script plots the PCCF for each cell mentioned.
- `8_3_PCCF_Avg_SEM.py`: This script averages the PCCF for all the cells, calculates the SEM for the average plots the same as line plots.
- `9_20nmOrigamiSpreadSDCalc.ipynb`: This notebook takes a 3D averaged origami file and uses multi-gaussian fitting to measure the spread of localizations on a predefined 20 nm DNA Origami template.
- `10_OrigamiFrameEvolution.ipynb`: This notebook plots the time evolution of target sampling of origamis.
- `11_1_crosstalk_analyzer.ipynb`: This notebook takes the *.hdf5 and *.yaml pick files and measures the crosstalk of imagers over different origmai species.
- `11_2_GenerateRandomPicks.ipynb`: This notebook takes many *.yaml finds the average number of picks and generate pick locations non-overlapping with any other picks.
- `11_3_crosstalk_strip_plotter.ipynb`: This notebook plots strip plots for the number of localization on each origami across all imaging channels. 
- `11_3_crosstalk_heatmap_plotter.ipynb`: This notebook plots the heatmap to showcase crosstalks.
- `12_kinetics_plotting.ipynb`: This notebook plots the kinetics of imager binding on individual origamis.
- `environment.yml`: This file contains information to build the epi-paint conda environment used extensively in these analysis codes. 
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