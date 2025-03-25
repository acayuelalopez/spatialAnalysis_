import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import Voronoi, voronoi_plot_2d
from sklearn.cluster import KMeans
from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

def generate_pca_plots(input_directory, output_directory):
    combined_df = pd.DataFrame()
    for filename in os.listdir(input_directory):
        if filename.endswith(".csv") and ("_TCell_distances.csv" in filename or "_tumor_distances.csv" in filename):
            parts = filename.split('_')
            condition = parts[0]
            patient_id = parts[1]
            population = '_'.join(parts[2:-2])
            status = 'Inf' if 'Inf' in filename else 'No_Inf'
            file_path = os.path.join(input_directory, filename)
            df = pd.read_csv(file_path)
            # Filter rows where the 'Classification' column is not empty
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
                patient_output_dir = os.path.join(output_directory, f'patient_{patient_id}' + '_pca_plots')
                os.makedirs(patient_output_dir, exist_ok=True)

                for status in ['Inf', 'No_Inf']:
                    status_df = population_df[population_df['Status'] == status]

                    # Perform PCA for Tumor distances
                    pca_tumor = PCA()
                    pca_tumor.fit(status_df[['Centroid X µm', 'Centroid Y µm', 'Min Distance to Tumor']])
                    explained_variance_tumor = np.cumsum(pca_tumor.explained_variance_ratio_)
                    plt.plot(range(1, len(explained_variance_tumor) + 1), explained_variance_tumor, marker='o')
                    plt.xlabel('Number of Components')
                    plt.ylabel('Cumulative Explained Variance')
                    plt.title('Explained Variance by PCA Components (Tumor)')
                    plt.show()

                    n_components_tumor = np.argmax(
                        explained_variance_tumor >= 0.95) + 1  # For example, 95% explained variance
                    print(f'Number of components to explain 95% variance (Tumor): {n_components_tumor}')

                    pca_tumor = PCA(n_components=n_components_tumor)
                    pca_result_tumor = pca_tumor.fit_transform(
                        status_df[['Centroid X µm', 'Centroid Y µm', 'Min Distance to Tumor']])
                    status_df['PCA1_Tumor'] = pca_result_tumor[:, 0]
                    status_df['PCA2_Tumor'] = pca_result_tumor[:, 1]

                    plt.figure(figsize=(12, 8))
                    sns.scatterplot(data=status_df, x='PCA1_Tumor', y='PCA2_Tumor', hue='Condition', palette='viridis')
                    plt.title(f'PCA Plot of Cell Positions (Tumor) for Patient {patient_id} - Status {status}')
                    plot_filename_tumor = f'{patient_id}_{population}_{status}_pca_plot_tumor.png'
                    plot_path_tumor = os.path.join(patient_output_dir, plot_filename_tumor)
                    plt.savefig(plot_path_tumor)
                    plt.close()

                    # Perform PCA for TCell distances
                    pca_tcell = PCA()
                    pca_tcell.fit(status_df[['Centroid X µm', 'Centroid Y µm', 'Min Distance to TCell']])
                    explained_variance_tcell = np.cumsum(pca_tcell.explained_variance_ratio_)
                    plt.plot(range(1, len(explained_variance_tcell) + 1), explained_variance_tcell, marker='o')
                    plt.xlabel('Number of Components')
                    plt.ylabel('Cumulative Explained Variance')
                    plt.title('Explained Variance by PCA Components (TCell)')
                    plt.show()

                    n_components_tcell = np.argmax(
                        explained_variance_tcell >= 0.95) + 1  # For example, 95% explained variance
                    print(f'Number of components to explain 95% variance (TCell): {n_components_tcell}')

                    pca_tcell = PCA(n_components=n_components_tcell)
                    pca_result_tcell = pca_tcell.fit_transform(
                        status_df[['Centroid X µm', 'Centroid Y µm', 'Min Distance to TCell']])
                    status_df['PCA1_TCell'] = pca_result_tcell[:, 0]
                    status_df['PCA2_TCell'] = pca_result_tcell[:, 1]

                    plt.figure(figsize=(12, 8))
                    sns.scatterplot(data=status_df, x='PCA1_TCell', y='PCA2_TCell', hue='Condition', palette='viridis')
                    plt.title(f'PCA Plot of Cell Positions (TCell) for Patient {patient_id} - Status {status}')
                    plot_filename_tcell = f'{patient_id}_{population}_{status}_pca_plot_tcell.png'
                    plot_path_tcell = os.path.join(patient_output_dir, plot_filename_tcell)
                    plt.savefig(plot_path_tcell)
                    plt.close()

def generate_tsne_plots(input_directory, output_directory):
    combined_df = pd.DataFrame()
    for filename in os.listdir(input_directory):
        if filename.endswith(".csv") and ("_TCell_distances.csv" in filename or "_tumor_distances.csv" in filename):
            parts = filename.split('_')
            condition = parts[0]
            patient_id = parts[1]
            population = '_'.join(parts[2:-2])
            status = 'Inf' if 'Inf' in filename else 'No_Inf'
            file_path = os.path.join(input_directory, filename)
            df = pd.read_csv(file_path)
            # Filter rows where the 'Classification' column is not empty
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
                patient_output_dir = os.path.join(output_directory, f'patient_{patient_id}'+'_tsne_plots')
                os.makedirs(patient_output_dir, exist_ok=True)

                for status in ['Inf', 'No_Inf']:
                    status_df = population_df[population_df['Status'] == status]

                    # Perform t-SNE for Tumor distances
                    tsne_tumor = TSNE(n_components=2, perplexity=30, n_iter=300)
                    tsne_result_tumor = tsne_tumor.fit_transform(status_df[['Centroid X µm', 'Centroid Y µm', 'Min Distance to Tumor']])
                    status_df['TSNE1_Tumor'] = tsne_result_tumor[:, 0]
                    status_df['TSNE2_Tumor'] = tsne_result_tumor[:, 1]

                    plt.figure(figsize=(12, 8))
                    sns.scatterplot(data=status_df, x='TSNE1_Tumor', y='TSNE2_Tumor', hue='Condition', palette='viridis')
                    plt.title(f't-SNE Plot of Cell Positions (Tumor) for Patient {patient_id} - Status {status}')
                    plot_filename_tumor = f'{patient_id}_{population}_{status}_tsne_plot_tumor.png'
                    plot_path_tumor = os.path.join(patient_output_dir, plot_filename_tumor)
                    plt.savefig(plot_path_tumor)
                    plt.close()

                    # Perform t-SNE for TCell distances
                    tsne_tcell = TSNE(n_components=2, perplexity=30, n_iter=300)
                    tsne_result_tcell = tsne_tcell.fit_transform(status_df[['Centroid X µm', 'Centroid Y µm', 'Min Distance to TCell']])
                    status_df['TSNE1_TCell'] = tsne_result_tcell[:, 0]
                    status_df['TSNE2_TCell'] = tsne_result_tcell[:, 1]

                    plt.figure(figsize=(12, 8))
                    sns.scatterplot(data=status_df, x='TSNE1_TCell', y='TSNE2_TCell', hue='Condition', palette='viridis')
                    plt.title(f't-SNE Plot of Cell Positions (TCell) for Patient {patient_id} - Status {status}')
                    plot_filename_tcell = f'{patient_id}_{population}_{status}_tsne_plot_tcell.png'
                    plot_path_tcell = os.path.join(patient_output_dir, plot_filename_tcell)
                    plt.savefig(plot_path_tcell)
                    plt.close()


input_directory = "/mnt/imgserver/CONFOCAL/IA/Projects/2025/2025_2_28_egonzalez/test/output_test/csv"
output_directory = "/mnt/imgserver/CONFOCAL/IA/Projects/2025/2025_2_28_egonzalez/test/output_test/figs"
tumor_border_directory = "/mnt/imgserver/CONFOCAL/IA/Projects/2024/2024_01_29_egonzalez/data/eva-gonzalez/spacio1/output/geojson"
generate_pca_plots(input_directory, output_directory)
#generate_tsne_plots(input_directory, output_directory)