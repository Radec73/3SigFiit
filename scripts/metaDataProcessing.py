import os
import re
import pandas as pd


def extract_patient_id(filename):
    match = re.match(r"^(CRUK\d+)_.*$", filename)
    return match.group(1) if match else None


def count_lines(file_path):
    with open(file_path, "r") as f:
        return sum(1 for line in f if line.strip() and not line.startswith("#"))


def calculate_sample_weights(folder_path, output_csv):
    vcf_files = [f for f in os.listdir(folder_path) if f.endswith('.vcf')]

    data = []

    for file in vcf_files:
        file_path = os.path.join(folder_path, file)
        line_count = count_lines(file_path)
        patient_id = extract_patient_id(file)

        if patient_id:
            data.append({"File": file[:17], "PatientID": patient_id, "LineCount": line_count})

    df = pd.DataFrame(data)

    if not df.empty:

        df["Weight"] = df.groupby("PatientID")["LineCount"].transform(lambda x: x / x.sum())
        outer_folder = os.path.dirname(folder_path)
        output_path = os.path.join(outer_folder, output_csv)
        df.to_csv(output_path, index=False)
        print(f"Sample weights saved to {output_path}")
    else:
        print("No valid patient data found. No file was saved.")
