keys_extractor = ['input_type','reference_genome', 'opportunity_genome', 'cosmic_version',
                  'context_type', 'exome', 'minimum_signatures', 'maximum_signatures','nmf_replicates','resample',
                  'batch_size','cpu','gpu','nmf_init','precision','matrix_normalization','seeds','min_nmf_iterations',
                  'max_nmf_iterations','nmf_test_conv','nmf_tolerance','nnls_add_penalty','nnls_remove_penalty',
                  'initial_remove_penalty','collapse_to_SBS96','clustering_distance','export_probabilities',
                  'make_decomposition_plots','stability','min_stability','combined_stability','allow_stability_drop',
                  'get_all_signature_matrices']
keys_assignment=['input_type','context_type','cosmic_version',
                'exome','genome_build','signature_database','exclude_signature_subgroups','export_probabilities',
                'export_probabilities_per_mutation','make_plots','sample_reconstruction_plots','verbose',
                'collapse_to_SBS96','nnls_add_penalty','nnls_remove_penalty','initial_remove_penalty']
keys_sigminer=['nruns','nsigs','data_folder','method','threshold','mode',
               'genome_build','sig_db','db_type','mut_type','seed']
keys_mutsignatures=['num_processesToExtract','num_totIterations','thresh_removeWeakMutTypes','thresh_removeLastPercent',
                    'distanceFunction','num_totReplicates','eps','stopconv','niter','approach','stopRule','algorithm','logIterations','sig_db','seed']

def running_parameters(config):
    for key, value in config.items():
        print(f"{key}: {repr(value)}")


def parse_value(value):
    if value.lower() == "true":
        return True
    elif value.lower() == "false":
        return False
    elif value.isdigit():
        return int(value)
    try:
        float_value = float(value)
        return float_value
    except ValueError:
        return value if value.lower() != "none" else None


def load_config(filename,sw):
    config = {}
    if sw == "SigMiner":
        keys = keys_sigminer
    elif sw == "Extractor":
        keys = keys_extractor
    elif sw == "Assignment":
        keys = keys_assignment
    elif sw == 'MutSignatures':
        keys = keys_mutsignatures
    else:
        keys = []
    with open(filename, "r") as file:
        for line in file:
            line = line.strip().rstrip(",")
            if line and "=" in line:
                key, value = line.split("=", 1)
                key = key.strip()
                if key not in keys:
                    raise ValueError(f"Error: Unexpected key '{key}' found in config file. Please change it according to instruction file.")
                config[key.strip()] = parse_value(value.strip())
    return config

