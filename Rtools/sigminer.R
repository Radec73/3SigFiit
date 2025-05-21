suppressPackageStartupMessages({
library(sigminer)
library(readr)
library(NMF)
library(BSgenome.Hsapiens.UCSC.hg38)
library(nnls)
})
args <- commandArgs(trailingOnly = TRUE)

nmfruns <- as.integer(args[1])
nsigs <- as.integer(args[2])
threshold <- as.numeric(args[3])
method <- args[4]
mode <- args[5]
genome_build <- args[6]
sig_db <- args[7]
db_type <- args[8]
mut_type <- args[9]
seed <- as.integer(args[10])
analysis <- args[11]


# if (sig_db == 'latest_SBS_GRCh38'){
#     exx <- read.table("COSMIC_v3.4_SBS_GRCh38_exome.txt", header = TRUE, sep = "\t", stringsAsFactors = TRUE)
#     rownames(exx) <- exx[, 1]
#     sig_matrix <- exx[, -1]
#     cosmic_catalog$db <- as.mutation.signatures(sig_matrix)
# } else {
#     cosmic_catalog <- get_sig_db(sig_db = sig_db)
# }

ref_genome_lib <- "BSgenome.Hsapiens.UCSC.hg38"

# Set base paths
wd <- getwd()
data_dir <- file.path(wd, "InputFiles", "ProcessedData")
output_base <- file.path(wd, "OutputFiles", "SigMinerOutput")
counts_output <- file.path(output_base)
denovo_tally_output <- file.path(output_base, "Denovo", "TallyResults")
denovo_sig_output <- file.path(output_base, "Denovo", "Signatures")
refitting_output <- file.path(output_base, "Refitting")

# Create necessary folders
dir.create(counts_output, recursive = TRUE, showWarnings = FALSE)
dir.create(denovo_tally_output, recursive = TRUE, showWarnings = FALSE)
dir.create(denovo_sig_output, recursive = TRUE, showWarnings = FALSE)
dir.create(refitting_output, recursive = TRUE, showWarnings = FALSE)

# Load data
vcf_files <- list.files(path = data_dir, pattern = "\\.vcf$", full.names = TRUE)
maf <- read_vcf(vcf_files, genome_build = genome_build, keep_only_pass = FALSE)

# Tally mutations
mt_tally <- sig_tally(maf, mode = mode, ref_genome = ref_genome_lib, genome_build = genome_build, useSyn = TRUE)

if (analysis == "denovo") {

  # De novo extraction
  mt_sig <- sig_extract(
    mt_tally$nmf_matrix,
    n_sig = nsigs,
    nrun = nmfruns,
    seed = seed,
    method = method,
    cores = 4
  )

  # Save de novo matrices
  for (i in seq_along(mt_sig)) {
    matrix_name <- names(mt_sig)[i]
    sig_matrix <- basis(mt_sig[[i]])  # EXTRACT basis matrix

    write.table(
      sig_matrix,
      file = file.path(denovo_tally_output, paste0(matrix_name, ".txt")),
      sep = "\t",
      quote = FALSE
    )
  }

  # Save tallied mutations
  output_tally(mt_tally, denovo_tally_output, mut_type = mut_type)

  # Save extracted signatures
  output_sig(mt_sig, denovo_sig_output, mut_type = mut_type)

} else if (analysis == "refitting") {

  mat_obj <- mt_tally$nmf_matrix
  mat_obj_t <- t(mat_obj)

  # Signature fitting
  fit_res <- sig_fit(
    catalogue_matrix = mat_obj_t,
    sig_db = sig_db,
    sig_index = "ALL",
    db_type = db_type,
    show_index = TRUE,
    return_class = "data.table",
    method = "NNLS",
    auto_reduce = FALSE,
    return_error = TRUE,
    rel_threshold = threshold,
    mode = mode
  )

  # Save mutation counts
  write.csv(
    t(as.data.frame(mat_obj)),
    file = file.path(counts_output, "Mutation_Counts.csv"),
    row.names = TRUE,
    quote = FALSE
  )

  # Save refitting results
  output_fit(fit_res, refitting_output, mut_type = mut_type, sig_db = sig_db)
}

quit(status = 0)
