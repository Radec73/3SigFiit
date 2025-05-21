import numpy as np
import pandas as pd
import math
from sigProfilerPlotting import plotActivity as plot_ac

def normalize_and_save(input_file, output_file):
    if input_file.endswith('.csv'):

        df = pd.read_csv(input_file, delimiter=",")
    elif input_file.endswith('.txt'):

        df = pd.read_csv(input_file, delimiter="\t")
    else:
        raise ValueError("Unsupported file format. Please provide a .csv or .txt file.")

    samples = df['Samples']

    data_values = df.drop('Samples', axis=1)

    min_values = data_values.min(axis=1)
    max_values = data_values.max(axis=1)

    normalized_data = data_values.apply(
        lambda row: ((row - min_values[row.name]) / (max_values[row.name] - min_values[row.name])).round(2), axis=1)

    normalized_data.insert(0, df.columns[0], samples)

    normalized_data.to_csv(output_file, sep=',', index=False)
    print(f"Processed data saved to {output_file}")


def frobenius_similarity(matrix1, matrix2):

    if matrix1.shape != matrix2.shape:
        raise ValueError("Matrices must have the same dimensions")

    matrix1 = np.round(matrix1, 2)
    matrix2 = np.round(matrix2, 2)

    frobenius_norm = np.linalg.norm(matrix2 - matrix1, 'fro')
    print(frobenius_norm)
    similarity = np.exp(-frobenius_norm)
    similarity_percent = similarity * 100
    print(f"Frobenius Similarity: {similarity_percent:.6f}%")
    return similarity


def convert_csv_to_txt(input_csv_path, output_txt_path):
    with open(input_csv_path, 'r') as csv_file:
        lines = csv_file.readlines()

    with open(output_txt_path, 'w') as txt_file:
        for line in lines:
            txt_file.write(line.replace(',', '\t'))


def format_sbs_csv(input_path, output_path, sw=''):
    if sw == 'SPA':

        df = pd.read_csv(input_path, sep='\t')

        df.to_csv(output_path, index=False)
    else:

        df = pd.read_csv(input_path, index_col=0)

        df = df.round(0).astype(int)

        df.reset_index(inplace=True)

        df.rename(columns={df.columns[0]: 'Samples'}, inplace=True)

        df['Samples'] = df['Samples'].str.replace('.bam', '', regex=False)

        df.to_csv(output_path, index=False)


def save_96_matrix_to_csv(matrices, output_path):

    matrix_96 = matrices.get('96')
    if matrix_96 is None or len(matrix_96) == 0:
        raise ValueError("'96' context not found or is empty in matrices.")


    if isinstance(matrix_96, dict):
        df = pd.DataFrame.from_dict(matrix_96, orient='index')
    else:
        df = matrix_96

    df = df.reindex(sorted(df.columns), axis=1)
    df.reset_index(inplace=True)
    df.rename(columns={'index': ''}, inplace=True)
    df.to_csv(output_path, index=False)


def aggregate_top_features(input_csv, output_csv):
    if input_csv.endswith('.csv'):

        df = pd.read_csv(input_csv)

    elif input_csv.endswith('.txt'):

        df = pd.read_csv(input_csv, delimiter="\t")

    else:
        raise ValueError("Unsupported file format. Please provide a .csv or .txt file.")

    feature_cols = [col for col in df.columns if col not in ['Samples', 'Total Count', 'Signature']]

    df[feature_cols] = df[feature_cols].apply(pd.to_numeric, errors='coerce')

    output_data = []

    for _, row in df.iterrows():
        sample_name = row['Samples']

        row_features = row[feature_cols].astype(float)
        top_features = row_features.nlargest(5).index.tolist()
        top_values = row_features.nlargest(5).values.tolist()

        output_data.append({
            'Sample': sample_name,
            'Top_Features': ', '.join(top_features),
            'Top_Values': ', '.join(map(str, top_values))
        })

    output_df = pd.DataFrame(output_data)
    output_df.to_csv(output_csv, index=False)
    print(f"Processed data saved to {output_csv}")


def aggregate_top_features_weighted(input_csv, sample_weights_csv, output_csv):

    if input_csv.endswith('.csv'):

        df = pd.read_csv(input_csv)
    elif input_csv.endswith('.txt'):

        df = pd.read_csv(input_csv, delimiter="\t")
    else:
        raise ValueError("Unsupported file format. Please provide a .csv or .txt file.")

    sample_weights = pd.read_csv(sample_weights_csv)

    metadata = pd.read_csv("Metadata.csv")

    df['Samples'] = df['Samples'].str.replace(r'\.bam$', '', regex=True)

    df['Patient'] = df['Samples'].str.extract(r'(CRUK\d+)')

    feature_cols = [col for col in df.columns if col not in ['Samples', 'Total Count', 'Signature', 'Patient']]

    df[feature_cols] = df[feature_cols].apply(pd.to_numeric, errors='coerce')

    df = df.merge(metadata[['TRACERxID', 'Smoking status']], left_on='Patient', right_on='TRACERxID', how='left')

    df = df.merge(sample_weights[['File', 'Weight', 'LineCount']], left_on='Samples', right_on='File', how='left')

    if df['Weight'].isna().sum() > 0:
        print("Warning: Some samples do not have a weight assigned.")

    df.drop(columns=['File'], inplace=True)

    output_data = []

    for patient, group in df.groupby('Patient'):
        weighted_features = {col: 0 for col in feature_cols}
        total_linecount = group['LineCount'].sum()
        for _, sample_row in group.iterrows():
            sample_weight = sample_row['Weight']

            for feature in feature_cols:
                weighted_features[feature] += sample_row[feature] * sample_weight

        top_features = sorted(weighted_features, key=weighted_features.get, reverse=True)[:5]
        #top_values = [weighted_features[feature] for feature in top_features]
        top_values = [math.floor(weighted_features[feature] * 100) / 100 for feature in top_features]

        smoking_status = group['Smoking status'].iloc[0] if not group['Smoking status'].isnull().all() else "Unknown"

        output_data.append({
            'Patient': patient,
            'Smoking_Status': smoking_status,
            'Total_LineCount': total_linecount,
            'Top_Features': ', '.join(top_features),
            'Top_Values': ', '.join(map(str, top_values))
        })

    output_df = pd.DataFrame(output_data)
    output_df.to_csv(output_csv, index=False)
    print(f"Processed data saved to {output_csv}")


def refitting_analysis_pipeline(
        input_txt_path,
        formatted_csv_path,
        normalized_csv_path,
        sorted_csv_path,
        weights_csv_path,
        weighted_csv_path,
        output_txt_path,
        output_plot_path,
        mutation_counts_path,
        mutation_counts_dest_path, sw=''):
    format_sbs_csv(input_txt_path, formatted_csv_path, sw)
    format_sbs_csv(mutation_counts_path, mutation_counts_dest_path)
    normalize_and_save(formatted_csv_path, normalized_csv_path)
    aggregate_top_features(normalized_csv_path, sorted_csv_path)
    aggregate_top_features_weighted(normalized_csv_path, weights_csv_path, weighted_csv_path)
    convert_csv_to_txt(formatted_csv_path, output_txt_path)
    plot_ac.plotActivity(output_txt_path, output_plot_path, bin_size=50)
