import os
import shutil

def delete_files_and_folders(folder):
    if os.path.exists(folder):
        for item in os.listdir(folder):
            item_path = os.path.join(folder, item)
            if os.path.isfile(item_path) or os.path.islink(item_path):
                os.remove(item_path)
            elif os.path.isdir(item_path):
                shutil.rmtree(item_path)
        print(f"All contents of '{folder}' have been deleted.")
    else:
        print(f"Folder '{folder}' does not exist.")


def SigProfilerVcfInput(input_folder, output_folder):
    if not os.path.exists(input_folder):
        os.makedirs(input_folder)

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    else:
        delete_files_and_folders(output_folder)

    for file_name in os.listdir(input_folder):
        if file_name.endswith('.vcf'):
            input_file_path = os.path.join(input_folder, file_name)
            output_file_path = os.path.join(output_folder, file_name[:17]) + '.vcf'
            with open(input_file_path, 'r') as vcf_file, open(output_file_path, 'w') as output_file:
                for line in vcf_file:
                    if not line.startswith('##') and not line.startswith('#'):
                        columns = line.strip().split('\t')
                        selected_columns = [columns[i] for i in [0, 1, 2, 3, 4]]
                        selected_columns[0] = columns[0].replace('chr', '')
                        selected_columns[2] = columns[2].replace('.', file_name[:17])
                        output_file.write('\t'.join(selected_columns) + '\n')


def SigMinerVcfInput(input_folder, output_folder):
    if not os.path.exists(input_folder):
        os.makedirs(input_folder)

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    else:
        delete_files_and_folders(output_folder)

    for file_name in os.listdir(input_folder):
        if file_name.endswith('.vcf'):
            input_file_path = os.path.join(input_folder, file_name)
            output_file_path = os.path.join(output_folder, file_name[:17] + '.vcf')

            with open(input_file_path, 'r') as infile, open(output_file_path, 'w') as outfile:
                outfile.write(infile.read())


def MutSignaturesVcfInput(input_folder, output_folder):
    if not os.path.exists(input_folder):
        os.makedirs(input_folder)

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    else:
        delete_files_and_folders(output_folder)

    for file_name in os.listdir(input_folder):
        if file_name.endswith('.vcf'):
            input_file_path = os.path.join(input_folder, file_name)
            output_file_path = os.path.join(output_folder, file_name[:17]) + '.vcf'

            with open(input_file_path, 'r') as f:
                header_lines = []
                data_lines = []
                for line in f:
                    if line.startswith("##"):
                        header_lines.append(line)
                    else:
                        data_lines.append(line)

            header = data_lines[0].strip().split('\t')
            sample_name = header[-1]
            header[-1] = "XTR1"
            header.append("SAMPLEID")

            data = []
            for line in data_lines[1:]:
                parts = line.strip().split('\t')
                parts[1] = str(parts[1])
                parts.append(sample_name)
                data.append(parts)

            with open(output_file_path, 'w') as f:
                # f.writelines(header_lines)
                f.write('\t'.join(header).replace("#", "") + '\n')
                for row in data:
                    f.write('\t'.join(row) + '\n')
