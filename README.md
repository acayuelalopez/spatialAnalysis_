# spatialAnalysis_
This repository contains a Python script designed to generate various types of plots from CSV files containing data related to T-cell and tumor distances. The script utilizes libraries such as `pandas`, `seaborn`, `matplotlib`, and `scipy` to create visualizations that aid in the analysis of spatial relationships between T-cells and tumors.

## Features

### Libraries Used
- **os**: For directory and file operations.
- **pandas**: For data manipulation and analysis.
- **seaborn**: For creating statistical plots.
- **matplotlib.pyplot**: For generating plots.
- **scipy.spatial.Voronoi**: For creating Voronoi diagrams.

### Functions
1. **generate_combined_bar_plots(input_directory, output_directory)**:
   - Reads CSV files from the input directory containing quantification data.
   - Combines data into a single DataFrame.
   - Generates bar plots for various parameters for each patient, categorized by infection status (`Inf` and `No_Inf`).

2. **generate_violin_box_plots(input_directory, output_directory)**:
   - Reads CSV files related to T-cell and tumor distances.
   - Combines data into a single DataFrame.
   - Generates violin and box plots for distance parameters for each patient, categorized by infection status.

3. **generate_violin_plots(input_directory, output_directory)**:
   - Similar to `generate_violin_box_plots`, but generates only violin plots.

4. **generate_scatter_plots(input_directory, tumor_border_directory, output_directory)**:
   - Reads CSV files related to T-cell and tumor distances.
   - Combines data into a single DataFrame.
   - Generates scatter plots of cell positions colored by distance to tumor and T-cell for each patient, categorized by infection status.

5. **generate_histograms(input_directory, output_directory)**:
   - Reads CSV files related to T-cell and tumor distances.
   - Combines data into a single DataFrame.
   - Generates histograms of minimum distances to tumor and T-cell for each patient, categorized by infection status.

6. **generate_voronoi_plots(input_directory, tumor_border_directory, output_directory)**:
   - Reads CSV files related to T-cell and tumor distances.
   - Combines data into a single DataFrame.
   - Generates Voronoi plots of cell positions colored by distance to tumor and T-cell for each patient, categorized by infection status.

7. **generate_heatmaps(input_directory, output_directory)**:
   - Reads CSV files related to T-cell and tumor distances.
   - Combines data into a single DataFrame.
   - Generates heatmaps of tumor distances for each patient, categorized by infection status.
   - Creates subplots for heatmaps of both distances to Tumor and TCell, considering the condition (`PRE` and `POST`). Only non-empty 'Classification' column values are plotted.

### Example Usage
The script includes example usage with specified directories for input, output, and tumor border data. The functions are commented out, indicating they can be called as needed.

### Key Parameters
- **PatientID**: Unique identifier for patients.
- **Status**: Infection status (`Inf` or `No_Inf`).
- **Condition**: Condition of the sample (`PRE` or `POST`).
- **Population**: Population type in the sample.
- **Distance Parameters**: Includes `Min Distance to TCell`, `Min Distance to Tumor`, etc.

### Plot Customizations
- **Figure Size**: Customizable figure sizes for plots.
- **Titles and Labels**: Titles and labels for plots are dynamically generated based on patient ID, parameter, and status.
- **Color Maps**: Various color maps are used for visualizing distances.

## Installation
To use this script, ensure you have the following Python libraries installed:
```bash
pip install pandas seaborn matplotlib scipy
