suppressPackageStartupMessages({
  library(readr)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(mutSignatures)
  library(sigminer)
  library(ggplot2)
})
getwd()
# Step 1: Read your CSV file
df <- read.csv("OutputFiles/MutSignaturesOutput/Denovo/Signatures.csv", row.names = 1, check.names = FALSE)

# Step 2: Convert to mutationSignatures object
signatures <- as.mutation.signatures(df)

cosmic_catalog <- get_sig_db(sig_db = "latest_SBS_GRCh38")
analysis_catalog <- as.mutation.signatures(as.data.frame(cosmic_catalog$db))

  msig1 <- matchSignatures(
    mutSign = signatures,
    reference = analysis_catalog,
    threshold = 0.5,
    method = "cosine",
    plot = TRUE
  )

  heatmap_plot <- msig1$plot + theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7)
  )
  ggplot2::ggsave(
    filename = file.path("Cosine_Similarity_Heatmap.png"),
    plot = heatmap_plot,
    width = 10,
    height = 5,
    dpi = 300
  )