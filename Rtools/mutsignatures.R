suppressPackageStartupMessages({
  library(readr)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(mutSignatures)
  library(sigminer)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
num_processesToExtract <- as.integer(args[1])
num_totIterations <- as.integer(args[2])
thresh_removeWeakMutTypes <- as.numeric(args[3])
thresh_removeLastPercent <- as.numeric(args[4])
distanceFunction <- args[5]
num_totReplicates <- as.integer(args[6])
eps <- as.numeric(args[7])
stopconv <- as.integer(args[8])
niter <- as.integer(args[9])
approach <- args[10]
stopRule <- args[11]
algorithm <- args[12]
logIterations <- args[13]
seed <- as.integer(args[14])
sig_db <- args[15]
analysis <- args[16]

hg38 <- BSgenome.Hsapiens.UCSC.hg38
cosmic_catalog <- get_sig_db(sig_db = sig_db)

# Set paths
wd <- getwd()
data_dir <- file.path(wd, "InputFiles", "ProcessedData")
output_dir <- file.path(wd, "OutputFiles", "MutSignaturesOutput")
denovo_outputs <- file.path(output_dir, "Denovo")
fit_outputs <- file.path(output_dir, "Refitting")

# Create output directories if missing
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(denovo_outputs, recursive = TRUE, showWarnings = FALSE)
dir.create(fit_outputs, recursive = TRUE, showWarnings = FALSE)

# Read VCF files
vcf_files <- list.files(path = data_dir, pattern = "\\.vcf$", full.names = TRUE)
vcf_data <- lapply(vcf_files, read_delim, delim = "\t", comment = "##", col_names = TRUE)
vcf <- do.call(rbind, vcf_data)

# Step 1: Attach context
context <- attachContext(
  mutData = vcf,
  chr_colName = "CHROM",
  start_colName = "POS",
  end_colName = "POS",
  nucl_contextN = 3,
  BSGenomeDb = hg38,
  context_colName = "CONTEXT"
)

# Step 2: Remove mismatches
context_clean <- removeMismatchMut(
  mutData = context,
  refMut_colName = "REF",
  context_colName = "CONTEXT",
  refMut_format = "N"
)

# Step 3: Compute mutation types
mut_type <- attachMutType(
  mutData = context_clean,
  ref_colName = "REF",
  var_colName = "ALT",
  context_colName = "CONTEXT",
  format = 1,
  mutType_dict = "alexa",
  mutType_colName = "MUT_TYPE"
)

# Step 4: Count mutation types per sample
allVCF_counts <- countMutTypes(
  mutTable = mut_type,
  mutType_colName = "MUT_TYPE",
  sample_colName = "SAMPLEID"
)

# Save mutation counts
write.csv(x = as.data.frame(allVCF_counts), file = file.path(output_dir, "Mutation_Counts.csv"), row.names = TRUE, quote = FALSE)

if (analysis == "refitting") {
  analysis_catalog <- as.mutation.signatures(as.data.frame(cosmic_catalog$db))

  refitting_activities_analysis <- resolveMutSignatures(
    mutCountData = allVCF_counts,
    signFreqData = analysis_catalog
  )

  count_res <- refitting_activities_analysis$Results$count.result
  freq_res <- refitting_activities_analysis$Results$freq.result

  write.csv(x = t(as.data.frame(count_res)), file = file.path(fit_outputs, "Refitting_absolute_exposure.csv"), row.names = TRUE, quote = FALSE)
  write.csv(x = t(as.data.frame(freq_res)), file = file.path(fit_outputs, "Refitting_relative_freq_exposure.csv"), row.names = TRUE, quote = FALSE)

} else if (analysis == "denovo") {
  # Step 6: Perform De Novo Signature Extraction
  params <- mutSignatures::setMutClusterParams(
    num_processesToExtract = num_processesToExtract,
    num_totIterations = num_totIterations,
    num_parallelCores = 4,
    thresh_removeWeakMutTypes = thresh_removeWeakMutTypes,
    thresh_removeLastPercent = thresh_removeLastPercent,
    distanceFunction = distanceFunction,
    num_totReplicates = num_totReplicates,
    eps = eps,
    stopconv = stopconv,
    niter = niter,
    guided = TRUE,
    debug = FALSE,
    approach = approach,
    stopRule = stopRule,
    algorithm = algorithm,
    logIterations = logIterations,
    seed = seed
  )

  denovo <- decipherMutationalProcesses(input = allVCF_counts, params = params)

  xprt_sig <- as.data.frame(denovo$Results$signatures)
  xprt_exp <- as.data.frame(denovo$Results$exposures)

  write.csv(xprt_sig, file = file.path(denovo_outputs, "Signatures.csv"), row.names = TRUE, quote = FALSE)
  write.csv(xprt_exp, file = file.path(denovo_outputs, "Exposures.csv"), row.names = TRUE, quote = FALSE)

  analysis_catalog <- as.mutation.signatures(as.data.frame(cosmic_catalog$db))

  msig1 <- matchSignatures(
    mutSign = denovo$Results$signatures,
    reference = analysis_catalog,
    threshold = 0.5,
    method = "cosine",
    plot = TRUE
  )

  heatmap_plot <- msig1$plot + theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7)
  )
  ggplot2::ggsave(
    filename = file.path(denovo_outputs, "Cosine_Similarity_Heatmap.png"),
    plot = heatmap_plot,
    width = 10,
    height = 5,
    dpi = 300
  )


  exposure_plot <- msigPlot(denovo$Results$exposures, top = 15) +
    ggplot2::scale_fill_manual(values = c('#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854'))

  ggplot2::ggsave(
    filename = file.path(denovo_outputs, "Exposures_Barplot.png"),
    plot = exposure_plot,
    width = 6.5,
    height = 4.5,
    dpi = 300
  )
  # Retrieve signatures (results)
  signatures <- denovo$Results$signatures

  # Define output folder and image size
  for (i in 1:5) {
    png(filename = file.path(denovo_outputs, paste0("Denovo_Signature", i, "_Barplot.png")),
        width = 800, height = 600, res = 120)

    # Generate and draw the plot
    msigPlot(signatures, signature = i, ylim = c(0, 0.16),
             main = paste("De Novo Signature", i))

    dev.off()
  }


  write.csv(x = as.data.frame(msig1$distanceDataFrame), file = file.path(denovo_outputs, "Cosine_Distance_Dataframe.csv"), row.names = TRUE, quote = FALSE)
  write.csv(x = as.data.frame(msig1$distanceMatrix), file = file.path(denovo_outputs, "Cosine_Distance_Matrix.csv"), row.names = TRUE, quote = FALSE)


}

quit(status = 0)
  #
  # signature_plot1 <- msigPlot(denovo$Results$signatures, signature = 1, ylim = c(0, 0.16),
  #          main = 'De Novo Signature 1')
  # ggplot2::ggsave(
  #   filename = file.path(denovo_outputs, "Denovo_Signature1_Barplot.png"),
  #   plot = signature_plot1,
  #   width = 6.5,
  #   height = 4.5,
  #   dpi = 300
  # )
  #
  # signature_plot2 <- msigPlot(denovo$Results$signatures, signature = 2, ylim = c(0, 0.16),
  #          main = 'De Novo Signature 2') + ggplot2::scale_fill_manual(values = c('#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854'))
  # ggplot2::ggsave(
  #   filename = file.path(denovo_outputs, "Denovo_Signature2_Barplot.png"),
  #   plot = signature_plot2,
  #   width = 6.5,
  #   height = 4.5,
  #   dpi = 300
  # )
  #
  # signature_plot3 <- msigPlot(denovo$Results$signatures, signature = 3, ylim = c(0, 0.16),
  #          main = 'De Novo Signature 3') + ggplot2::scale_fill_manual(values = c('#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854'))
  # ggplot2::ggsave(
  #   filename = file.path(denovo_outputs, "Denovo_Signature3_Barplot.png"),
  #   plot = signature_plot3,
  #   width = 6.5,
  #   height = 4.5,
  #   dpi = 300
  # )
  #
  # signature_plot4 <- msigPlot(denovo$Results$signatures, signature = 4, ylim = c(0, 0.16),
  #          main = 'De Novo Signature 4') + ggplot2::scale_fill_manual(values = c('#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854'))
  # ggplot2::ggsave(
  #   filename = file.path(denovo_outputs, "Denovo_Signature4_Barplot.png"),
  #   plot = signature_plot4,
  #   width = 6.5,
  #   height = 4.5,
  #   dpi = 300
  # )
  #
  # signature_plot5 <- msigPlot(denovo$Results$signatures, signature = 5, ylim = c(0, 0.16),
  #          main = 'De Novo Signature 5') + ggplot2::scale_fill_manual(values = c('#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854'))
  #
  # ggplot2::ggsave(
  #   filename = file.path(denovo_outputs, "Denovo_Signature5_Barplot.png"),
  #   plot = signature_plot5,
  #   width = 6.5,
  #   height = 4.5,
  #   dpi = 300
  # )