import os
import pandas as pd
import numpy as np
from scipy.spatial import cKDTree


def process_files(input_directory, tumor_border_directory, output_directory, data_file):
    # Read the data file
    data_df = pd.read_excel(data_file, engine='odf')

    # Initialize a dictionary to store quantification results for each population and patient
    quantification_results = {}

    # Iterate over each file in the input directory
    for filename in os.listdir(input_directory):
        if filename.endswith(".csv"):  # and "TCell" not in filename:
            print(filename)

            # Extract condition, patient ID, and population from the filename
            parts = filename.split('.ome.tif')
            condition_patient = parts[0]
            population_status = parts[1].replace('.csv', '')
            print(population_status)

            # Split population_status into population and status
            split_parts = population_status.split('_')
            if 'No' in split_parts:
                status_index = split_parts.index('No')
                status = 'No_Inf'
                population = '_'.join(split_parts[:status_index])
            elif 'Inf' in split_parts:
                status_index = split_parts.index('Inf')
                status = 'Inf'
                population = '_'.join(split_parts[:status_index])
            else:
                print(f"Unexpected format in filename: {filename}")
                continue

            # Now population should be 'Macrophage M1' and status should be 'Inf' or 'No_Inf'
            print(f"Population: {population}, Status: {status}")


            # Determine condition and patient ID based on the condition
            if condition_patient.startswith("POST"):
                condition = condition_patient[:4]
                patient_id = condition_patient[4:]
            else:  # Assuming the other condition is PRE
                condition = condition_patient[:3]
                patient_id = condition_patient[3:]

            # Read the CSV file
            file_path = os.path.join(input_directory, filename)
            df = pd.read_csv(file_path, sep='\t')

            # Print column names to debug
            print(f"Columns in {filename}: {df.columns.tolist()}")

            # Check if 'Classification' column exists
            if 'Classification' not in df.columns:
                print(f"'Classification' column not found in {filename}")
                continue

            # Count the total number of cells
            total_cells = len(df)

            # Count the number of cells for each classification
            classification_counts = df['Classification'].value_counts()

            # Get the tumor and CK+ areas from the data file
            tumor_area = data_df[(data_df['PatientID'] == int(patient_id)) & (data_df['Condition'] == condition)][
                'Tumor Anot Area µm^2'].values[0]
            ck_area = data_df[(data_df['PatientID'] == int(patient_id)) & (data_df['Condition'] == condition)][
                'CK+ Area µm^2'].values[0]

            # Find the corresponding tumor border file
            tumor_border_file = next((tb_file for tb_file in os.listdir(tumor_border_directory) if
                                      tb_file.endswith(".csv") and patient_id in tb_file), None)
            if tumor_border_file is None:
                print(f"No corresponding tumor border file found for {filename}")
                continue

            # Load tumor border data
            tumor_border_path = os.path.join(tumor_border_directory, tumor_border_file)
            tumor_border = pd.read_csv(tumor_border_path)
            print(f"Tumor border data loaded from {tumor_border_path}")

            # Calculate minimum distances from each cell centroid to the tumor border
            cell_coords = df[['Centroid X µm', 'Centroid Y µm']].values
            tumor_coords = tumor_border[['Centroid X (µm)', 'Centroid Y (µm)']].values
            tree = cKDTree(tumor_coords)
            distances, _ = tree.query(cell_coords)

            # Add distance columns to DataFrame
            df['Min Distance to Tumor'] = distances

            # Calculate average, min, and max distances to tumor border for each classification
            avg_distance_to_tumor = df.groupby('Classification')['Min Distance to Tumor'].mean()
            min_distance_to_tumor = df.groupby('Classification')['Min Distance to Tumor'].min()
            max_distance_to_tumor = df.groupby('Classification')['Min Distance to Tumor'].max()

            # Store the quantification results in the dictionary
            if (population, patient_id) not in quantification_results:
                quantification_results[(population, patient_id)] = []

            for classification, count in classification_counts.items():
                quantification_results[(population, patient_id)].append({
                    'Population': population,
                    'Condition': condition,
                    'PatientID': patient_id,
                    'Status': status,
                    'Classification': classification,
                    'Count': count,
                    'TotalCells': total_cells,
                    'Normalized by Tumor Anot Area µm^2': count / tumor_area,
                    'Normalized by CK+ Area µm^2': count / ck_area,
                    'Normalized by TotalCells': count / total_cells,
                    'Average Distance to Tumor': avg_distance_to_tumor[classification],
                    'Min Distance to Tumor': min_distance_to_tumor[classification],
                    'Max Distance to Tumor': max_distance_to_tumor[classification],
                    'Average Distance to TCell': np.nan,
                    'Min Distance to TCell': np.nan,
                    'Max Distance to TCell': np.nan,
                    'Tumor Anot Area µm^2': tumor_area,
                    'CK+ Area µm^2': ck_area
                })

            # Find the maximum number of elements separated by ':'
            max_elements_count = df['Classification'].dropna().apply(lambda x: len(x.split(':'))).max()

            # Remove values in Classification column that do not contain the maximum number of elements separated by ':'
            df.loc[df['Classification'].apply(lambda x: len(x.split(':')) != max_elements_count if isinstance(x,
                                                                                                              str) else False), 'Classification'] = ''

            # Save the output data with distances to a new CSV file
            output_csv_file = os.path.join(output_directory,
                                           f'{condition}_{patient_id}_{population}_{status}_tumor_distances.csv')
            df.to_csv(output_csv_file, index=False)
            print(f"Saved combined data with distances to {output_csv_file}")

    # Process TCell files separately and update quantification results with TCell distances
    for tcell_filename in os.listdir(input_directory):
        if "TCell" in tcell_filename and tcell_filename.endswith(".csv"):
            print(tcell_filename+"--------")
            # Extract condition, patient ID, and population from the filename
            parts_tcell = tcell_filename.split('.ome.tif')
            condition_patient_tcell = parts_tcell[0]
            population_status = parts_tcell[1].replace('.csv', '')

            # Split population_status into population and status
            split_parts = population_status.split('_')
            if 'No' in split_parts:
                status_index = split_parts.index('No')
                status = 'No_Inf'
                population = '_'.join(split_parts[:status_index])
            elif 'Inf' in split_parts:
                status_index = split_parts.index('Inf')
                status = 'Inf'
                population = '_'.join(split_parts[:status_index])
            else:
                print(f"Unexpected format in filename: {tcell_filename}")
                continue

            # Now population should be 'Macrophage M1' and status should be 'Inf' or 'No_Inf'
            print(f"Population: {population}, Status: {status}")


            # Determine condition and patient ID based on the condition
            if condition_patient_tcell.startswith("POST"):
                condition_tcell = condition_patient_tcell[:4]
                patient_id_tcell = condition_patient_tcell[4:]
            else:  # Assuming the other condition is PRE
                condition_tcell = condition_patient_tcell[:3]
                patient_id_tcell = condition_patient_tcell[3:]

            # Read the TCell CSV file
            tcell_file_path = os.path.join(input_directory, tcell_filename)
            df_tcell = pd.read_csv(tcell_file_path, sep='\t')

            # Print column names to debug
            print(f"Columns in {tcell_filename}: {df_tcell.columns.tolist()}")

            # Check if 'Classification' column exists and contains TCell label
            #if 'Classification' not in df_tcell.columns or "CD3+: CD16-: CD163-: CK-: HLA-DR-" not in df_tcell['Classification'].values:
            if 'Classification' not in df_tcell.columns or "CD3+: CD16-: CD206-: CD163-" not in df_tcell[
                'Classification'].values:
                print(f"'Classification' column not found or TCell label not found in {tcell_filename}")
                continue

            # Filter cells with TCell classification
            #tcell_cells = df_tcell[df_tcell['Classification'] == "CD3+: CD16-: CD163-: CK-: HLA-DR-"]
            tcell_cells = df_tcell[df_tcell['Classification'] == "CD3+: CD16-: CD206-: CD163-"]

            # Calculate minimum distances from each cell centroid to the TCell cells for each population file again
            for filename in os.listdir(input_directory):
                if filename.endswith(".csv"): #and "TCell" not in filename:
                    parts_pop = filename.split('.ome.tif')
                    condition_patient_pop = parts_pop[0]
                    population_status_pop = parts_pop[1].replace('.csv', '')

                    # Split population_status into population and status
                    split_parts = population_status_pop.split('_')
                    if 'No' in split_parts:
                        status_index = split_parts.index('No')
                        status_pop = 'No_Inf'
                        population_pop = '_'.join(split_parts[:status_index])
                    else:
                        print(f"Unexpected format in filename: {tcell_filename}")
                        continue

                    # Now population should be 'Macrophage M1' and status should be 'No_Inf'
                    print(f"Population: {population}, Status: {status_pop}")

                    if condition_patient_pop.startswith("POST"):
                        condition_pop = condition_patient_pop[:4]
                        patient_id_pop = condition_patient_pop[4:]
                    else:  # Assuming the other condition is PRE
                        condition_pop = condition_patient_pop[:3]
                        patient_id_pop = condition_patient_pop[3:]

                    # Check if the condition and patient ID match
                    if condition_pop == condition_tcell and patient_id_pop == patient_id_tcell:
                        pop_file_path = os.path.join(input_directory, filename)
                        df_pop = pd.read_csv(pop_file_path, sep='\t')
                        cell_coords_pop = df_pop[['Centroid X µm', 'Centroid Y µm']].values
                        tcell_coords = tcell_cells[['Centroid X µm', 'Centroid Y µm']].values
                        tree_tcell = cKDTree(tcell_coords)
                        distances_tcell, _ = tree_tcell.query(cell_coords_pop)

                        # Add distance columns to DataFrame for TCell cells
                        df_pop['Min Distance to TCell'] = distances_tcell

                        # Calculate average, min, and max distances to TCell cells for each classification
                        avg_distance_to_tcell = df_pop.groupby('Classification')['Min Distance to TCell'].mean()
                        min_distance_to_tcell = df_pop.groupby('Classification')['Min Distance to TCell'].min()
                        max_distance_to_tcell = df_pop.groupby('Classification')['Min Distance to TCell'].max()

                        # Find the maximum number of elements separated by ':'
                        max_elements_count = df_pop['Classification'].dropna().apply(lambda x: len(x.split(':'))).max()

                        # Remove values in Classification column that do not contain the maximum number of elements separated by ':'
                        df_pop.loc[df_pop['Classification'].apply(
                            lambda x: len(x.split(':')) != max_elements_count if isinstance(x,
                                                                                            str) else False), 'Classification'] = ''

                        # Save the output data with TCell distances to a new CSV file
                        output_csv_file_tcell = os.path.join(output_directory,
                                                             f'{condition_pop}_{patient_id_pop}_{population_pop}_{status_pop}_TCell_distances.csv')
                        df_pop.to_csv(output_csv_file_tcell, index=False)
                        print(f"Saved TCell data with distances to {output_csv_file_tcell}")

                        # Update quantification results with TCell distances for each classification count entry
                        if (population_pop, patient_id_pop) in quantification_results:
                            for result in quantification_results[(population_pop, patient_id_pop)]:
                                if result['Condition'] == condition_pop and result[
                                    'Classification'] in avg_distance_to_tcell.index:
                                    result.update({
                                        'Average Distance to TCell': avg_distance_to_tcell[result['Classification']],
                                        'Min Distance to TCell': min_distance_to_tcell[result['Classification']],
                                        'Max Distance to TCell': max_distance_to_tcell[result['Classification']]
                                    })
                        else:
                            print(f"Key ({population_pop}, {patient_id_pop}) not found in quantification_results")
    # Save the updated quantification results to CSV files
    for (population, patient_id), results in quantification_results.items():
        quantification_df = pd.DataFrame(results)
        # Ensure the population and status are separated correctly
        quantification_df['Population'] = quantification_df['Population'].apply(lambda x: x.split('_')[0])
        quantification_df['Status'] = quantification_df[
            'Status']  # .apply(lambda x: x.split('_')[1] if '_' in x else x)

        # Create rows for each combination of condition and status
        conditions = ['PRE', 'POST']
        statuses = ['Inf', 'No_Inf']
        combined_results = []
        for condition in conditions:
            for status in statuses:
                filtered_df = quantification_df[
                    (quantification_df['Condition'] == condition) & (quantification_df['Status'] == status)]
                combined_results.append(filtered_df)

        combined_df = pd.concat(combined_results)

        # Maintain rows with highest string length in 'Classification' column
        max_length = combined_df['Classification'].str.len().max()
        combined_df_filtered = combined_df[combined_df['Classification'].str.len() == max_length]

        # Save the filtered DataFrame to the CSV file using population_pop
        output_filename = f"{population}_quantification_{patient_id}.csv"
        output_path = os.path.join(output_directory, output_filename)
        combined_df_filtered.to_csv(output_path, index=False)

# Example usage
input_directory = "/mnt/imgserver/CONFOCAL/IA/Projects/2025/2025_2_28_egonzalez/test/input_test"
tumor_border_directory = "/mnt/imgserver/CONFOCAL/IA/Projects/2024/2024_01_29_egonzalez/data/eva-gonzalez/spacio1/output/geojson"
#output_directory = "/mnt/imgserver/CONFOCAL/IA/Projects/2024/2024_01_29_egonzalez/test/csv/csv"
output_directory = "/mnt/imgserver/CONFOCAL/IA/Projects/2025/2025_2_28_egonzalez/test/output_test"
data_file = "/mnt/imgserver/CONFOCAL/IA/Projects/2024/2024_01_29_egonzalez/count/csv/tumor_anot_area.ods"

process_files(input_directory, tumor_border_directory, output_directory, data_file)
