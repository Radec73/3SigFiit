import os
import subprocess
from scripts import inputDataProcessing as dp, \
    resultsProcessing as rp, \
    parametersProcessing as pp, \
    metaDataProcessing as mp, \
    run

from SigProfilerMatrixGenerator.scripts import (SigProfilerMatrixGeneratorFunc as matrix)

input_raw_data_folder = 'InputFiles/InputData'
output_processed_data_folder = 'InputFiles/ProcessedData'

def option_1():
    print("*** SigProfiler Tools ***")
    dp.SigProfilerVcfInput(input_raw_data_folder, output_processed_data_folder)
    print("\nDATA SUCCESSFULLY PROCESSED AND PREPARED FOR ANALYSIS!\n")

    while True:
        print("Choose your option:")
        print("h - Help and instructions for running SigProfiler modules")
        print("d - Run De Novo Signature Extraction with SigProfilerExtractor module")
        print("r - Run Signature Exposures Refitting with SigProfilerAssignment module")
        print("b - Back to menu")
        choice = input()
        if choice == 'h':
            with open("InputFiles/InstructionsFiles/sigprofiler.txt", "r", encoding="utf-8") as file:
                help = file.read()
                print(help)
        elif choice == 'r':
            dp.delete_files_and_folders("OutputFiles/SigProfilerOutput/Refitting")
            print("Running *** SigProfilerAssignment *** Signature Exposures Refitting with these parameters:\n")
            sigProfilerAssignmentParameters = pp.load_config("InputFiles/Settings/SigProfilerRefitting.txt","Assignment")
            pp.running_parameters(sigProfilerAssignmentParameters)
            try:
                run.run_SigProfiler_Assignment(sigProfilerAssignmentParameters,output_processed_data_folder,'OutputFiles/SigProfilerOutput/Refitting')
                print("Run was successful. You can find results in Refitting folder.\n")

                counts = matrix.SigProfilerMatrixGeneratorFunc(
                    project="",
                    reference_genome="GRCh38",
                    path_to_input_files=output_processed_data_folder,
                    exome='genome',
                    bed_file=None,
                    chrom_based=False,
                    plot=False,
                    gs=False,
                )
                rp.save_96_matrix_to_csv(counts, "OutputFiles/SigProfilerOutput/Mutation_Counts.csv")
                rp.refitting_analysis_pipeline(
                    input_txt_path='OutputFiles/SigProfilerOutput/Refitting/Assignment_Solution/Activities/Assignment_Solution_Activities.txt',
                    formatted_csv_path='OutputFiles/Refitting Analysis/SigProfiler/Assignment_Solution_Activities.csv',
                    normalized_csv_path='OutputFiles/Refitting Analysis/SigProfiler/Normalized_Assignment_Solution_Activities.csv',
                    sorted_csv_path='OutputFiles/Refitting Analysis/SigProfiler/Sorted_Assignment_Solution_Activities.csv',
                    weights_csv_path='OutputFiles/Refitting Analysis/sample_weights.csv',
                    weighted_csv_path='OutputFiles/Refitting Analysis/SigProfiler/Weighted_Assignment_Solution_Activities.csv',
                    output_txt_path='OutputFiles/Refitting Analysis/SigProfiler/Assignment_Solution_Activities.txt',
                    output_plot_path='OutputFiles/Refitting Analysis/SigProfiler/Activities_in_samples.pdf',
                    mutation_counts_path='OutputFiles/SigProfilerOutput/Mutation_Counts.csv',
                    mutation_counts_dest_path='OutputFiles/Refitting Analysis/SigProfiler/Mutation_Counts.csv',
                    sw="SPA")
                print("Refitting Analysis files are stored in Refitting Analysis/SigProfiler\n")

            except Exception as e:
                print(f"An error occurred: {e}")
                return

        elif choice == 'd':
            dp.delete_files_and_folders("OutputFiles/SigProfilerOutput/Denovo")
            print("Running *** SigProfilerExtractor *** De Novo Signature Extraction with these parameters:\n")
            sigProfilerExtractorParameters = pp.load_config("InputFiles/Settings/SigProfilerDeNovo.txt", "Extractor")
            pp.running_parameters(sigProfilerExtractorParameters)
            try:
                run.run_SigProfiler_Extractor(sigProfilerExtractorParameters,output_processed_data_folder,'OutputFiles/SigProfilerOutput/Denovo')
                print("Run was successful. You can find results in Denovo folder.\n")
            except Exception as e:
                print(f"An error occurred: {e}")
                return

        elif choice == 'b':
            print("Thank you for using this tool.\n")
            break
        else:
            print("Invalid choice. Please choose from available choices.\n")


def option_2():
    r_script_path = "Rtools/mutsignatures.R"
    print("*** MutSignatures ***")
    dp.MutSignaturesVcfInput(input_raw_data_folder, output_processed_data_folder)
    print("\nDATA SUCCESSFULLY PROCESSED AND PREPARED FOR ANALYSIS!\n")

    while True:
        print("Choose your option:")
        print("h - Help and instructions for running MutSignatures")
        print("d - Run De Novo Signature Extraction with MutSignatures framework")
        print("r - Run Signature Exposures Refitting with MutSignatures framework")
        print("b - Back to menu")
        choice = input()
        if choice == 'h':
            with open("InputFiles/InstructionsFiles/mutsignatures.txt", "r", encoding="utf-8") as file:
                help = file.read()
                print(help)

        elif choice == 'd':
            dp.delete_files_and_folders("OutputFiles/MutSignaturesOutput/Denovo")
            print("\nRunning *** MutSignatures *** De Novo Extraction of Signatures  *** with these parameters:\n")
            mutSignatures_parameters = pp.load_config("InputFiles/Settings/MutSignaturesDeNovo.txt", "MutSignatures")
            pp.running_parameters(mutSignatures_parameters)
            try:
                subprocess.run(["Rscript", r_script_path,
                                str(mutSignatures_parameters["num_processesToExtract"]),
                                str(mutSignatures_parameters["num_totIterations"]),
                                str(mutSignatures_parameters["thresh_removeWeakMutTypes"]),
                                str(mutSignatures_parameters["thresh_removeLastPercent"]),
                                str(mutSignatures_parameters["distanceFunction"]),
                                str(mutSignatures_parameters["num_totReplicates"]),
                                str(mutSignatures_parameters["eps"]),
                                str(mutSignatures_parameters["stopconv"]),
                                str(mutSignatures_parameters["niter"]),
                                str(mutSignatures_parameters["approach"]),
                                str(mutSignatures_parameters["stopRule"]),
                                str(mutSignatures_parameters["algorithm"]),
                                str(mutSignatures_parameters["logIterations"]),
                                str(mutSignatures_parameters["seed"]),
                                str(mutSignatures_parameters["sig_db"]),
                                'denovo'], check=True)
                print("\nR script execution completed successfully.\n")
                print("You can find results in MutSignaturesOutput.\n")
            except subprocess.CalledProcessError as e:
                print(f"\nError: R script execution failed with return code {e.returncode}.\n")
        elif choice == 'r':
            dp.delete_files_and_folders("OutputFiles/MutSignaturesOutput/Refitting")
            print("\nRunning *** MutSignatures *** Signature Exposures Refitting\n")
            mutSignatures_parameters = pp.load_config("InputFiles/Settings/MutSignaturesRefitting.txt", "MutSignatures")
            pp.running_parameters(mutSignatures_parameters)
            try:
                subprocess.run(["Rscript", r_script_path,
                                '', '', '', '', '', '', '', '', '', '', '', '', '', '',
                                str(mutSignatures_parameters["sig_db"]),
                                'refitting'], check=True)
                print("\nR script execution completed successfully.\n")
                print("You can find results in MutSignaturesOutput.\n")

                rp.refitting_analysis_pipeline(
                    input_txt_path='OutputFiles/MutSignaturesOutput/Refitting/Refitting_absolute_exposure.csv',
                    formatted_csv_path='OutputFiles/Refitting Analysis/MutSignatures/Assignment_Solution_Activities.csv',
                    normalized_csv_path='OutputFiles/Refitting Analysis/MutSignatures/Normalized_Assignment_Solution_Activities.csv',
                    sorted_csv_path='OutputFiles/Refitting Analysis/MutSignatures/Sorted_Assignment_Solution_Activities.csv',
                    weights_csv_path='OutputFiles/Refitting Analysis/sample_weights.csv',
                    weighted_csv_path='OutputFiles/Refitting Analysis/MutSignatures/Weighted_Assignment_Solution_Activities.csv',
                    output_txt_path='OutputFiles/Refitting Analysis/MutSignatures/Assignment_Solution_Activities.txt',
                    output_plot_path='OutputFiles/Refitting Analysis/MutSignatures/Activities_in_samples.pdf',
                    mutation_counts_path='OutputFiles/MutSignaturesOutput/Mutation_Counts.csv',
                    mutation_counts_dest_path='OutputFiles/Refitting Analysis/MutSignatures/Mutation_Counts.csv')

                print("Refitting Analysis files are stored in Refitting Analysis/MutSignatures\n")
            except subprocess.CalledProcessError as e:
                print(f"\nError: R script execution failed with return code {e.returncode}.\n")

        elif choice == 'b':
            print("Thank you for using this tool.\n")
            break
        else:
            print("Invalid choice. Please choose from available choices.\n")


def option_3():
    r_script_path = "Rtools/sigminer.R"
    print("*** SigMiner ***")
    dp.SigMinerVcfInput(input_raw_data_folder, output_processed_data_folder)
    print("\nDATA SUCCESSFULLY PROCESSED AND PREPARED FOR ANALYSIS!\n")

    while True:
        print("Choose your option:")
        print("h - Help and instructions for running Sigminer")
        print("d - Run De Novo Signature Extraction with Sigminer")
        print("r - Run Signature Exposures Refitting with Sigminer")
        print("b - Back to menu")
        choice = input()
        if choice == 'h':
            with open("InputFiles/InstructionsFiles/sigminer.txt", "r", encoding="utf-8") as file:
                help = file.read()
                print(help)

        elif choice == 'd':
            dp.delete_files_and_folders("OutputFiles/SigMinerOutput/Denovo")
            print(
                "\nRunning *** SigMiner *** De Novo Extraction of Signatures using method sig_extract() with these parameters:\n")
            sigMiner_parameters = pp.load_config("InputFiles/Settings/SigMinerDenovo.txt", "SigMiner")
            pp.running_parameters(sigMiner_parameters)

            try:
                subprocess.run(["Rscript", r_script_path,
                                str(sigMiner_parameters["nruns"]),
                                str(sigMiner_parameters["nsigs"]),
                                '',
                                str(sigMiner_parameters["method"]),
                                str(sigMiner_parameters["mode"]),
                                str(sigMiner_parameters["genome_build"]),
                                '',
                                '',
                                str(sigMiner_parameters["mut_type"]),
                                str(sigMiner_parameters["seed"]),
                                'denovo'], check=True)
                print("\nR script execution completed successfully.\n")
                print("You can find results in SigMinerOutput/Denovo.\n")
            except subprocess.CalledProcessError as e:
                print(f"\nError: R script execution failed with return code {e.returncode}.\n")
                return

        elif choice == 'r':
            dp.delete_files_and_folders("OutputFiles/SigMinerOutput/Refitting")
            print(
                "Running *** SigMiner *** Signature Exposure Refitting using method sig_fit() with these parameters:\n")
            sigMiner_parameters = pp.load_config("InputFiles/Settings/SigMinerRefitting.txt", "SigMiner")
            pp.running_parameters(sigMiner_parameters)

            try:
                subprocess.run(["Rscript", r_script_path,
                                '',
                                '',
                                str(sigMiner_parameters["threshold"]),
                                str(sigMiner_parameters["method"]),
                                str(sigMiner_parameters["mode"]),
                                str(sigMiner_parameters["genome_build"]),
                                str(sigMiner_parameters["sig_db"]),
                                str(sigMiner_parameters["db_type"]),
                                str(sigMiner_parameters["mut_type"]),
                                '',
                                'refitting'], check=True)
                print("\nR script execution completed successfully.\n")
                print("You can find results in SigMinerOutput/Refitting.\n")
                rp.refitting_analysis_pipeline(
                    input_txt_path='OutputFiles/SigMinerOutput/Refitting/SBS_fitting_absolute_exposure.csv',
                    formatted_csv_path='OutputFiles/Refitting Analysis/Sigminer/Assignment_Solution_Activities.csv',
                    normalized_csv_path='OutputFiles/Refitting Analysis/SigMiner/Normalized_Assignment_Solution_Activities.csv',
                    sorted_csv_path='OutputFiles/Refitting Analysis/Sigminer/Sorted_Assignment_Solution_Activities.csv',
                    weights_csv_path='OutputFiles/Refitting Analysis/sample_weights.csv',
                    weighted_csv_path='OutputFiles/Refitting Analysis/Sigminer/Weighted_Assignment_Solution_Activities.csv',
                    output_txt_path='OutputFiles/Refitting Analysis/Sigminer/Assignment_Solution_Activities.txt',
                    output_plot_path='OutputFiles/Refitting Analysis/Sigminer/Activities_in_samples.pdf',
                    mutation_counts_path='OutputFiles/SigMinerOutput/Mutation_Counts.csv',
                    mutation_counts_dest_path='OutputFiles/Refitting Analysis/Sigminer/Mutation_Counts.csv')
                print("Refitting Analysis files are stored in Refitting Analysis/Sigminer\n")
            except subprocess.CalledProcessError as e:
                print(f"\nError: R script execution failed with return code {e.returncode}.\n")
                return

        elif choice == 'b':
            print("Thank you for using this tool.\n")
            break
        else:
            print("Invalid choice. Please choose from available choices.\n")


def main():
    while True:
        print("\nMutational Signatures Analysis Tools\n")
        print("1 - SigProfiler")
        print("2 - MutSignatures")
        print("3 - SigMiner")
        print("E - Exit")

        choice = input("\nChoose software for analyzing samples by picking number from 1 to 3:\n")

        if choice == '1':
            option_1()
        elif choice == '2':
            option_2()
        elif choice == '3':
            option_3()
        elif choice == 'E':
            print("Exiting program. Thank you for using this tool.")
            break
        else:
            print("Invalid choice. Please enter a number between 1 and 3 or E for exiting.")


if __name__ == "__main__":
    print("\n" + "=" * 100)
    print("  WELCOME TO THE".center(95))
    print("   GENOMIC DATA ANALYZER TOOL   ".center(95))
    print("=" * 100)
    print("""
       Powered by:
            Python implementation of SigProfiler tools
            R packages: Sigminer & MutSignatures

        Designed for:
            De novo signature extraction
            Refitting to known COSMIC signatures

        Reminder:
          Make sure your input files are:
            - Correctly formatted in vcf format
            - Located in the InputData folder
    """)
    print("=" * 100 + "\n")
    if not os.path.exists("OutputFiles/Refitting Analysis/sample_weights.csv"):
        print("Calculating sample weights... ")
        mp.calculate_sample_weights("InputFiles/InputData", "sample_weights.csv")

        print("File with sample weights was saved to Refitting Analysis/sample_weights.csv")

    else:
        print("File with sample weights is stored in Refitting Analysis/sample_weights.csv")
    main()
