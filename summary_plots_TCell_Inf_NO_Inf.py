import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi, voronoi_plot_2d


def generate_combined_bar_plots(input_directory, output_directory):
    # Initialize an empty DataFrame to store all data
    combined_df = pd.DataFrame()

    # Iterate over each file in the input directory
    for filename in os.listdir(input_directory):
        if filename.endswith(".csv") and "quantification" in filename:
            # Read the CSV file
            file_path = os.path.join(input_directory, filename)
            print(f"Reading file: {file_path}")
            df = pd.read_csv(file_path)
            combined_df = pd.concat([combined_df, df], ignore_index=True)

    # Get unique patient IDs
    patient_ids = combined_df['PatientID'].unique()

    # Define the parameters to plot
    parameters = [
        'Count', 'TotalCells', 'Normalized by Tumor Anot Area µm^2', 'Normalized by CK+ Area µm^2',
        'Normalized by TotalCells', 'Average Distance to Tumor', 'Min Distance to Tumor',
        'Max Distance to Tumor', 'Average Distance to TCell', 'Min Distance to TCell',
        'Max Distance to TCell', 'Tumor Anot Area µm^2', 'CK+ Area µm^2'
    ]

    # Generate combined bar plots for each patient and parameter
    for patient_id in patient_ids:
        patient_df = combined_df[combined_df['PatientID'] == patient_id]

        # Create a directory for the plots if it doesn't exist
        patient_output_dir = os.path.join(output_directory, f'patient_{patient_id}' + '_bar_plots')
        os.makedirs(patient_output_dir, exist_ok=True)

        for parameter in parameters:
            for status in ['Inf', 'No_Inf']:
                status_df = patient_df[patient_df['Status'] == status]

                plt.figure(figsize=(36, 16))
                # Bar plot
                sns.barplot(x='Population', y=parameter, hue='Condition', data=status_df)
                plt.title(f'Bar Plot of {parameter} for Patient {patient_id} - Status {status}')
                plt.legend(title='Condition')
                plt.xticks(rotation=45, ha='right', fontsize=8)  # Rotate x-axis labels and adjust font size

                # Save the plot
                plot_filename = f'{patient_id}_{parameter}_{status}.png'
                plot_path = os.path.join(patient_output_dir, plot_filename)
                plt.savefig(plot_path)
                plt.close()


def generate_violin_box_plots(input_directory, output_directory):
    combined_df = pd.DataFrame()
    for filename in os.listdir(input_directory):
        if filename.endswith(".csv") and ("_TCell_distances.csv" in filename or "_tumor_distances.csv" in filename):
            parts = filename.split('_')
            condition = parts[0]
            patient_id = parts[1]
            population = parts[2]
            print(parts)
            print(filename)
            print(condition)
            print(patient_id)
            print(population)
            status = 'No_Inf' if 'No_Inf' in filename else 'Inf'
            print(status)
            file_path = os.path.join(input_directory, filename)
            df = pd.read_csv(file_path)
            df = df[df['Classification'].notna()]
            df['Condition'] = condition
            df['PatientID'] = patient_id
            df['Population'] = population
            df['Status'] = status
            combined_df = pd.concat([combined_df, df], ignore_index=True)

    patient_ids = combined_df['PatientID'].unique()
    distance_parameters = ['Min Distance to TCell', 'Min Distance to Tumor']

    for patient_id in patient_ids:
        patient_df = combined_df[combined_df['PatientID'] == patient_id]
        patient_output_dir = os.path.join(output_directory, f'patient_{patient_id}'+ '_violin_box_plots')
        os.makedirs(patient_output_dir, exist_ok=True)

        for parameter in distance_parameters:
            for status in ['Inf', 'No_Inf']:
                status_df = patient_df[patient_df['Status'] == status]

                plt.figure(figsize=(36, 16))
                plt.subplot(1, 2, 1)
                sns.violinplot(x='Population', y=parameter, hue='Condition', data=status_df, split=True)
                plt.title(f'Violin Plot of {parameter} for Patient {patient_id} - Status {status}')
                plt.legend(title='Condition')
                plt.xticks(rotation=45, ha='right', fontsize=8)

                plt.subplot(1, 2, 2)
                sns.boxplot(x='Population', y=parameter, hue='Condition', data=status_df)
                plt.title(f'Box Plot of {parameter} for Patient {patient_id} - Status {status}')
                plt.legend(title='Condition')
                plt.xticks(rotation=45, ha='right', fontsize=8)

                plot_filename = f'{patient_id}_{parameter}_{status}_violin_box.png'
                plot_path = os.path.join(patient_output_dir, plot_filename)
                plt.savefig(plot_path)
                plt.close()


def generate_violin_plots(input_directory, output_directory):
    combined_df = pd.DataFrame()
    for filename in os.listdir(input_directory):
        if filename.endswith(".csv") and ("_TCell_distances.csv" in filename or "_tumor_distances.csv" in filename):
            parts = filename.split('_')
            condition = parts[0]
            patient_id = parts[1]
            population = parts[2]
            print(parts)
            print(filename)
            print(condition)
            print(patient_id)
            print(population)
            status = 'No_Inf' if 'No_Inf' in filename else 'Inf'
            print(status)
            file_path = os.path.join(input_directory, filename)
            df = pd.read_csv(file_path)
            df = df[df['Classification'].notna()]
            df['Condition'] = condition
            df['PatientID'] = patient_id
            df['Population'] = population
            df['Status'] = status
            combined_df = pd.concat([combined_df, df], ignore_index=True)

    patient_ids = combined_df['PatientID'].unique()
    distance_parameters = ['Min Distance to TCell', 'Min Distance to Tumor']

    for patient_id in patient_ids:
        patient_df = combined_df[combined_df['PatientID'] == patient_id]
        patient_output_dir = os.path.join(output_directory, f'patient_{patient_id}'+ '_violin_plots')
        os.makedirs(patient_output_dir, exist_ok=True)

        for parameter in distance_parameters:
            for status in ['Inf', 'No_Inf']:
                status_df = patient_df[patient_df['Status'] == status]

                plt.figure(figsize=(26, 16))
                plt.subplot(1, 2, 1)
                sns.violinplot(x='Population', y=parameter, hue='Condition', data=status_df, split=True, inner="quart")
                plt.title(f'Violin Plot of {parameter} for Patient {patient_id} - Status {status}')
                plt.legend(title='Condition')
                plt.xticks(rotation=45, ha='right', fontsize=8)

                plot_filename = f'{patient_id}_{parameter}_{status}_violin_plot.png'
                plot_path = os.path.join(patient_output_dir, plot_filename)
                plt.savefig(plot_path)
                plt.close()

def generate_scatter_plots(input_directory, tumor_border_directory, output_directory):
    combined_df = pd.DataFrame()
    for filename in os.listdir(input_directory):
        if filename.endswith(".csv") and ("_TCell_distances.csv" in filename or "_tumor_distances.csv" in filename):
            parts = filename.split('_')
            condition = parts[0]
            patient_id = parts[1]
            population = parts[2]
            print(parts)
            print(filename)
            print(condition)
            print(patient_id)
            print(population)
            status = 'No_Inf' if 'No_Inf' in filename else 'Inf'
            print(status)
            file_path = os.path.join(input_directory, filename)
            df = pd.read_csv(file_path)
            df = df[df['Classification'].notna()]
            df['Condition'] = condition
            df['PatientID'] = patient_id
            df['Population'] = population
            df['Status'] = status
            combined_df = pd.concat([combined_df, df], ignore_index=True)

    patient_ids = combined_df['PatientID'].unique()

    for patient_id in patient_ids:
        patient_df = combined_df[combined_df['PatientID'] == patient_id]
        populations = patient_df['Population'].unique()

        for population in populations:
            population_df = patient_df[patient_df['Population'] == population]

            if not population_df.empty:
                patient_output_dir = os.path.join(output_directory, f'patient_{patient_id}'+'_scatter_plots')
                os.makedirs(patient_output_dir, exist_ok=True)

                for status in ['Inf', 'No_Inf']:
                    status_df = population_df[population_df['Status'] == status]

                    # Load tumor border data
                    tumor_border_file = next((tb_file for tb_file in os.listdir(tumor_border_directory) if
                                              tb_file.endswith(".csv") and patient_id in tb_file and condition in tb_file), None)
                    if tumor_border_file is None:
                        print(f"No corresponding tumor border file found for {filename}")
                        continue

                    tumor_border_path = os.path.join(tumor_border_directory, tumor_border_file)
                    tumor_border_df = pd.read_csv(tumor_border_path)
                    tumor_coords = tumor_border_df[['Centroid X (µm)', 'Centroid Y (µm)']].values

                    plt.figure(figsize=(36, 12))
                    distance_column_tumor = 'Min Distance to Tumor'
                    distance_column_tcell = 'Min Distance to TCell'

                    # Subfigure for PRE Tumor
                    plt.subplot(1, 2, 1)
                    pre_df = status_df[status_df['Condition'] == 'PRE']
                    if distance_column_tumor in pre_df.columns:
                        plt.scatter(tumor_coords[:, 0], tumor_coords[:, 1], c='gray', s=5, label='Tumor Border', alpha=0.2)
                        scatter_tumor_pre = plt.scatter(pre_df['Centroid X µm'],
                                                        pre_df['Centroid Y µm'],
                                                        c=pre_df[distance_column_tumor],
                                                        cmap='viridis', s=5)
                        plt.colorbar(scatter_tumor_pre).set_label(distance_column_tumor)
                        plt.title(f'Scatter Plot of Cell Positions Colored by {distance_column_tumor} for Patient {patient_id} (PRE) - Status {status}')
                        plt.legend(title='Condition')

                    # Subfigure for POST Tumor
                    plt.subplot(1, 2, 2)
                    post_df = status_df[status_df['Condition'] == 'POST']
                    if distance_column_tumor in post_df.columns:
                        plt.scatter(tumor_coords[:, 0], tumor_coords[:, 1], c='gray', s=5, label='Tumor Border', alpha=0.2)
                        scatter_tumor_post = plt.scatter(post_df['Centroid X µm'],
                                                         post_df['Centroid Y µm'],
                                                         c=post_df[distance_column_tumor],
                                                         cmap='viridis', s=5)
                        plt.colorbar(scatter_tumor_post).set_label(distance_column_tumor)
                        plt.title(f'Scatter Plot of Cell Positions Colored by {distance_column_tumor} for Patient {patient_id} (POST) - Status {status}')
                        plt.legend(title='Condition')

                    plot_filename_tumor = f'{patient_id}_{population}_{status}_scatter_plot_tumor.png'
                    plot_path_tumor = os.path.join(patient_output_dir, plot_filename_tumor)
                    plt.savefig(plot_path_tumor)
                    plt.close()

                    plt.figure(figsize=(36, 12))

                    # Subfigure for PRE TCell
                    plt.subplot(1, 2, 1)
                    pre_df = status_df[status_df['Condition'] == 'PRE']
                    if distance_column_tcell in pre_df.columns:
                        plt.scatter(tumor_coords[:, 0], tumor_coords[:, 1], c='gray', s=5, label='Tumor Border', alpha=0.2)
                        scatter_tcell_pre = plt.scatter(pre_df['Centroid X µm'],
                                                        pre_df['Centroid Y µm'],
                                                        c=pre_df[distance_column_tcell],
                                                        cmap='viridis', s=5)
                        plt.colorbar(scatter_tcell_pre).set_label(distance_column_tcell)
                        plt.title(f'Scatter Plot of Cell Positions Colored by {distance_column_tcell} for Patient {patient_id} (PRE) - Status {status}')
                        plt.legend(title='Condition')

                    # Subfigure for POST TCell
                    plt.subplot(1, 2, 2)
                    post_df = status_df[status_df['Condition'] == 'POST']
                    if distance_column_tcell in post_df.columns:
                        plt.scatter(tumor_coords[:, 0], tumor_coords[:, 1], c='gray', s=5, label='Tumor Border', alpha=0.2)
                        scatter_tcell_post = plt.scatter(post_df['Centroid X µm'],
                                                         post_df['Centroid Y µm'],
                                                         c=post_df[distance_column_tcell],
                                                         cmap='viridis', s=5)
                        plt.colorbar(scatter_tcell_post).set_label(distance_column_tcell)
                        plt.title(f'Scatter Plot of Cell Positions Colored by {distance_column_tcell} for Patient {patient_id} (POST) - Status {status}')
                        plt.legend(title='Condition')

                    plot_filename_tcell = f'{patient_id}_{population}_{status}_scatter_plot_tcell.png'
                    plot_path_tcell = os.path.join(patient_output_dir, plot_filename_tcell)
                    plt.savefig(plot_path_tcell)
                    plt.close()



def generate_histograms(input_directory, output_directory):
    # Initialize an empty DataFrame to store all data
    combined_df = pd.DataFrame()

    # Iterate over each file in the input directory
    for filename in os.listdir(input_directory):
        if filename.endswith(".csv") and ("_TCell_distances.csv" in filename or "_tumor_distances.csv" in filename):
            parts = filename.split('_')
            condition = parts[0]
            patient_id = parts[1]
            population = parts[2]
            status = 'No_Inf' if 'No_Inf' in filename else 'Inf'
            file_path = os.path.join(input_directory, filename)
            df = pd.read_csv(file_path)
            df = df[df['Classification'].notna()]
            df['Condition'] = condition
            df['PatientID'] = patient_id
            df['Population'] = population
            df['Status'] = status
            combined_df = pd.concat([combined_df, df], ignore_index=True)

    # Get unique patient IDs
    patient_ids = combined_df['PatientID'].unique()

    # Generate histograms for each patient and population
    for patient_id in patient_ids:
        patient_df = combined_df[combined_df['PatientID'] == patient_id]

        # Create a directory for the plots if it doesn't exist
        patient_output_dir = os.path.join(output_directory, f'patient_{patient_id}'+'_histogram_plots')
        os.makedirs(patient_output_dir, exist_ok=True)

        populations = patient_df['Population'].unique()

        for population in populations:
            population_df = patient_df[patient_df['Population'] == population]

            if population_df.empty:
                print(f"No data available for Patient {patient_id} - Population {population}")
                continue

            for status in ['Inf', 'No_Inf']:
                status_df = population_df[population_df['Status'] == status]

                if 'Min Distance to Tumor' not in status_df.columns or status_df['Min Distance to Tumor'].isna().all():
                    print(
                        f"No valid 'Min Distance to Tumor' data for Patient {patient_id} - Population {population} - Status {status}")
                    continue

                plt.figure(figsize=(12, 6))
                sns.histplot(data=status_df, x='Min Distance to Tumor', hue='Condition', element='step', stat='density',
                             common_norm=False)
                plt.xlabel('Min Distance to Tumor')
                plt.ylabel('Density')
                plt.title(
                    f'Histogram of Min Distance to Tumor for Patient {patient_id} - Population {population} - Status {status}')
                plt.legend(title='Condition')

                plot_filename = f'{patient_id}_{population}_{status}_histogram_tumor.png'
                plot_path = os.path.join(patient_output_dir, plot_filename)
                plt.savefig(plot_path)
                plt.close()

                if 'Min Distance to TCell' not in status_df.columns or status_df['Min Distance to TCell'].isna().all():
                    print(
                        f"No valid 'Min Distance to TCell' data for Patient {patient_id} - Population {population} - Status {status}")
                    continue

                plt.figure(figsize=(12, 6))
                sns.histplot(data=status_df, x='Min Distance to TCell', hue='Condition', element='step', stat='density',
                             common_norm=True, kde=True)
                plt.xlabel('Min Distance to TCell')
                plt.ylabel('Density')
                plt.title(
                    f'Histogram of Min Distance to TCell for Patient {patient_id} - Population {population} - Status {status}')
                plt.legend(title='Condition')

                plot_filename = f'{patient_id}_{population}_{status}_histogram_tcell.png'
                plot_path = os.path.join(patient_output_dir, plot_filename)
                plt.savefig(plot_path)
                plt.close()


def generate_voronoi_plots(input_directory, tumor_border_directory, output_directory):
    combined_df = pd.DataFrame()
    for filename in os.listdir(input_directory):
        if filename.endswith(".csv") and ("_TCell_distances.csv" in filename or "_tumor_distances.csv" in filename):
            parts = filename.split('_')
            condition = parts[0]
            patient_id = parts[1]
            population = parts[2]
            status = 'No_Inf' if 'No_Inf' in filename else 'Inf'
            file_path = os.path.join(input_directory, filename)
            df = pd.read_csv(file_path)
            df = df[df['Classification'].notna()]
            df['Condition'] = condition
            df['PatientID'] = patient_id
            df['Population'] = population
            df['Status'] = status
            combined_df = pd.concat([combined_df, df], ignore_index=True)
    patient_ids = combined_df['PatientID'].unique()
    for patient_id in patient_ids:
        patient_df = combined_df[combined_df['PatientID'] == patient_id]
        populations = patient_df['Population'].unique()
        for population in populations:
            population_df = patient_df[patient_df['Population'] == population]
            if not population_df.empty:
                patient_output_dir = os.path.join(output_directory, f'patient_{patient_id}' + '_voronoi_plots')
                os.makedirs(patient_output_dir, exist_ok=True)
                for status in ['Inf', 'No_Inf']:
                    status_df = population_df[population_df['Status'] == status]
                    tumor_border_file = next((tb_file for tb_file in os.listdir(tumor_border_directory) if
                                              tb_file.endswith(".csv") and patient_id in tb_file and condition in tb_file), None)
                    if tumor_border_file is None:
                        print(f"No corresponding tumor border file found for {filename}")
                        continue
                    tumor_border_path = os.path.join(tumor_border_directory, tumor_border_file)
                    tumor_border_df = pd.read_csv(tumor_border_path)
                    tumor_coords = tumor_border_df[['Centroid X (µm)', 'Centroid Y (µm)']].values

                    distance_column_tumor = 'Min Distance to Tumor'
                    distance_column_tcell = 'Min Distance to TCell'

                    plt.figure(figsize=(36, 12))

                    # Subfigure for PRE Tumor Voronoi plot
                    plt.subplot(1, 2, 1)
                    pre_df = status_df[status_df['Condition'] == 'PRE']
                    if not pre_df.empty:
                        points_pre = pre_df[['Centroid X µm', 'Centroid Y µm']].values
                        vor_pre = Voronoi(points_pre)
                        voronoi_plot_2d(vor_pre, ax=plt.gca(), show_vertices=False, line_colors='blue', line_width=2, line_alpha=0.6)
                        plt.scatter(tumor_coords[:, 0], tumor_coords[:, 1], c='gray', s=10, label='Tumor Border', alpha=0.5)
                        scatter_tumor_pre = plt.scatter(pre_df['Centroid X µm'],
                                                        pre_df['Centroid Y µm'],
                                                        c=pre_df[distance_column_tumor],
                                                        cmap='viridis', s=5)
                        plt.colorbar(scatter_tumor_pre).set_label(distance_column_tumor)
                        plt.title(f'Voronoi Plot of Cell Positions Colored by {distance_column_tumor} for Patient {patient_id} (PRE) - Status {status}')
                        plt.legend(title='Condition')

                    # Subfigure for POST Tumor Voronoi plot
                    plt.subplot(1, 2, 2)
                    post_df = status_df[status_df['Condition'] == 'POST']
                    if not post_df.empty:
                        points_post = post_df[['Centroid X µm', 'Centroid Y µm']].values
                        vor_post = Voronoi(points_post)
                        voronoi_plot_2d(vor_post, ax=plt.gca(), show_vertices=False, line_colors='blue', line_width=2, line_alpha=0.6)
                        plt.scatter(tumor_coords[:, 0], tumor_coords[:, 1], c='gray', s=10, label='Tumor Border', alpha=0.5)
                        scatter_tumor_post = plt.scatter(post_df['Centroid X µm'],
                                                         post_df['Centroid Y µm'],
                                                         c=post_df[distance_column_tumor],
                                                         cmap='viridis', s=5)
                        plt.colorbar(scatter_tumor_post).set_label(distance_column_tumor)
                        plt.title(f'Voronoi Plot of Cell Positions Colored by {distance_column_tumor} for Patient {patient_id} (POST) - Status {status}')
                        plt.legend(title='Condition')

                    plot_filename_tumor = f'{patient_id}_{population}_{status}_voronoi_plot_tumor.png'
                    plot_path_tumor = os.path.join(patient_output_dir, plot_filename_tumor)
                    plt.savefig(plot_path_tumor)
                    plt.close()

                    plt.figure(figsize=(36, 12))

                    # Subfigure for PRE TCell Voronoi plot
                    plt.subplot(1, 2, 1)
                    pre_df = status_df[status_df['Condition'] == 'PRE']
                    if not pre_df.empty:
                        points_pre = pre_df[['Centroid X µm', 'Centroid Y µm']].values
                        vor_pre = Voronoi(points_pre)
                        voronoi_plot_2d(vor_pre, ax=plt.gca(), show_vertices=False, line_colors='blue', line_width=2, line_alpha=0.6)
                        #plt.scatter(tumor_coords[:, 0], tumor_coords[:, 1], c='gray', s=10, label='Tumor Border', alpha=0.5)
                        scatter_tcell_pre = plt.scatter(pre_df['Centroid X µm'],
                                                        pre_df['Centroid Y µm'],
                                                        c=pre_df[distance_column_tcell],
                                                        cmap='viridis', s=5)
                        plt.colorbar(scatter_tcell_pre).set_label(distance_column_tcell)
                        plt.title(f'Voronoi Plot of Cell Positions Colored by {distance_column_tcell} for Patient {patient_id} (PRE) - Status {status}')
                        plt.legend(title='Condition')

                    # Subfigure for POST TCell Voronoi plot
                    plt.subplot(1, 2, 2)
                    post_df = status_df[status_df['Condition'] == 'POST']
                    if not post_df.empty:
                        points_post = post_df[['Centroid X µm', 'Centroid Y µm']].values
                        vor_post = Voronoi(points_post)
                        voronoi_plot_2d(vor_post, ax=plt.gca(), show_vertices=False, line_colors='blue', line_width=2, line_alpha=0.6)
                        #plt.scatter(tumor_coords[:, 0], tumor_coords[:, 1], c='gray', s=10, label='Tumor Border', alpha=0.5)
                        scatter_tcell_post = plt.scatter(post_df['Centroid X µm'],
                                                         post_df['Centroid Y µm'],
                                                         c=post_df[distance_column_tcell],
                                                         cmap='viridis', s=5)
                        plt.colorbar(scatter_tcell_post).set_label(distance_column_tcell)
                        plt.title(f'Voronoi Plot of Cell Positions Colored by {distance_column_tcell} for Patient {patient_id} (POST) - Status {status}')
                        plt.legend(title='Condition')

                    plot_filename_tcell = f'{patient_id}_{population}_{status}_voronoi_plot_tcell.png'
                    plot_path_tcell = os.path.join(patient_output_dir, plot_filename_tcell)
                    plt.savefig(plot_path_tcell)
                    plt.close()


def generate_heatmaps(input_directory, output_directory):
    # Initialize an empty DataFrame to store all data
    combined_df = pd.DataFrame()

    # Iterate over each file in the input directory
    for filename in os.listdir(input_directory):
        if filename.endswith(".csv") and ("_TCell_distances.csv" in filename or "_tumor_distances.csv" in filename):
            parts = filename.split('_')
            condition = parts[0]
            patient_id = parts[1]
            population = parts[2]
            status = 'No_Inf' if 'No_Inf' in filename else 'Inf'
            file_path = os.path.join(input_directory, filename)
            df = pd.read_csv(file_path)
            df = df[df['Classification'].notna()]
            df['Condition'] = condition
            df['PatientID'] = patient_id
            df['Population'] = population
            df['Status'] = status
            combined_df = pd.concat([combined_df, df], ignore_index=True)

    # Get unique patient IDs
    patient_ids = combined_df['PatientID'].unique()

    # Generate heatmaps for each patient and population
    for patient_id in patient_ids:
        patient_df = combined_df[combined_df['PatientID'] == patient_id]
        populations = patient_df['Population'].unique()
        for population in populations:
            population_df = patient_df[patient_df['Population'] == population]
            if not population_df.empty:
                patient_output_dir = os.path.join(output_directory, f'patient_{patient_id}' + '_heatmaps')
                os.makedirs(patient_output_dir, exist_ok=True)
                for status in ['Inf', 'No_Inf']:
                    status_df = population_df[population_df['Status'] == status]

                    plt.figure(figsize=(24, 12))

                    # Subplot for Min Distance to Tumor heatmap (PRE condition)
                    plt.subplot(2, 2, 1)
                    pre_df_tumor = status_df[status_df['Condition'] == 'PRE']
                    heatmap_data_tumor_pre = pre_df_tumor.pivot_table(index='Centroid Y µm', columns='Centroid X µm',
                                                                      values='Min Distance to Tumor', aggfunc='mean')
                    sns.heatmap(heatmap_data_tumor_pre, cmap='viridis')
                    plt.title(
                        f'Heatmap of Min Distance to Tumor (PRE) for Patient {patient_id} - Population {population} - Status {status}')

                    # Subplot for Min Distance to Tumor heatmap (POST condition)
                    plt.subplot(2, 2, 2)
                    post_df_tumor = status_df[status_df['Condition'] == 'POST']
                    heatmap_data_tumor_post = post_df_tumor.pivot_table(index='Centroid Y µm', columns='Centroid X µm',
                                                                        values='Min Distance to Tumor', aggfunc='mean')
                    sns.heatmap(heatmap_data_tumor_post, cmap='viridis')
                    plt.title(
                        f'Heatmap of Min Distance to Tumor (POST) for Patient {patient_id} - Population {population} - Status {status}')

                    # Subplot for Min Distance to TCell heatmap (PRE condition)
                    plt.subplot(2, 2, 3)
                    pre_df_tcell = status_df[status_df['Condition'] == 'PRE']
                    heatmap_data_tcell_pre = pre_df_tcell.pivot_table(index='Centroid Y µm', columns='Centroid X µm',
                                                                      values='Min Distance to TCell', aggfunc='mean')
                    sns.heatmap(heatmap_data_tcell_pre, cmap='viridis')
                    plt.title(
                        f'Heatmap of Min Distance to TCell (PRE) for Patient {patient_id} - Population {population} - Status {status}')

                    # Subplot for Min Distance to TCell heatmap (POST condition)
                    plt.subplot(2, 2, 4)
                    post_df_tcell = status_df[status_df['Condition'] == 'POST']
                    heatmap_data_tcell_post = post_df_tcell.pivot_table(index='Centroid Y µm', columns='Centroid X µm',
                                                                        values='Min Distance to TCell', aggfunc='mean')
                    sns.heatmap(heatmap_data_tcell_post, cmap='viridis')
                    plt.title(
                        f'Heatmap of Min Distance to TCell (POST) for Patient {patient_id} - Population {population} - Status {status}')

                    plot_filename = f'{patient_id}_{population}_{status}_heatmaps.png'
                    plot_path = os.path.join(patient_output_dir, plot_filename)
                    plt.savefig(plot_path)
                    plt.close()


# Example usage
input_directory = "/mnt/imgserver/CONFOCAL/IA/Projects/2025/2025_2_28_egonzalez/test/output_test/csv"
output_directory = "/mnt/imgserver/CONFOCAL/IA/Projects/2025/2025_2_28_egonzalez/test/output_test/figs"
tumor_border_directory = "/mnt/imgserver/CONFOCAL/IA/Projects/2024/2024_01_29_egonzalez/data/eva-gonzalez/spacio1/output/geojson"
#generate_violin_box_plots(input_directory, output_directory)
#generate_violin_plots(input_directory, output_directory)
#generate_scatter_plots(input_directory,tumor_border_directory, output_directory)
#generate_voronoi_plots(input_directory, tumor_border_directory, output_directory)
#generate_combined_bar_plots(input_directory, output_directory)
#generate_histograms(input_directory, output_directory)
generate_heatmaps(input_directory, output_directory)
