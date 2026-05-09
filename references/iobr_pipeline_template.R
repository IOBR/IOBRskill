###############################################################################
# IOBR Full TME Deconvolution and Signature Scoring Pipeline
# Reference template adapted from the author's production workflow
#
# IOBR Citation:
# Zeng DQ, Fang YR, … , Yu GC, Liao WJ.
# Enhancing Immuno-Oncology Investigations Through Multidimensional Decoding of
# Tumour Microenvironment with IOBR 2.0, Cell Reports Methods, 2024
# https://doi.org/10.1016/j.crmeth.2024.100910
###############################################################################

library(IOBR)
library(tidyverse)

############################
# User-configurable parameters
ProjectID  <- "PROJECT001"
tumor_type <- "stad"         # Cancer type abbreviation for TIMER
index      <- 1
############################


# ===========================================================================
# Step 0: Data preparation — log2 transform if needed
# ===========================================================================
# eset should be a gene × sample matrix (rows = genes, cols = samples)
eset[1:5, 1:5]

# Ensure log2 scale (auto-detects and transforms if necessary)
eset <- log2eset(eset)


# ===========================================================================
# Step 1: TME deconvolution — run all 8 methods
# ===========================================================================
# Each method may take 1-5 minutes depending on sample size.
# If TIMER throws an error for your cancer type, comment it out.

cibersort <- deconvo_tme(eset = eset, method = "cibersort",   arrays = FALSE, perm = 1000)
epic      <- deconvo_tme(eset = eset, method = "epic",         arrays = FALSE)
mcp       <- deconvo_tme(eset = eset, method = "mcpcounter")
xcell     <- deconvo_tme(eset = eset, method = "xcell",        arrays = FALSE)
estimate  <- deconvo_tme(eset = eset, method = "estimate")
timer     <- deconvo_tme(eset = eset, method = "timer",
                          group_list = rep(tumor_type, ncol(eset)))
quantiseq <- deconvo_tme(eset = eset, method = "quantiseq",
                          tumor = TRUE, arrays = FALSE, scale_mrna = TRUE)
ips       <- deconvo_tme(eset = eset, method = "ips", plot = FALSE)

# Merge all deconvolution results by sample ID
tme_combine <- cibersort %>%
  inner_join(., mcp,       by = "ID") %>%
  inner_join(., xcell,     by = "ID") %>%
  inner_join(., epic,      by = "ID") %>%
  inner_join(., estimate,  by = "ID") %>%
  inner_join(., quantiseq, by = "ID") %>%
  inner_join(., timer,     by = "ID") %>%
  inner_join(., ips,       by = "ID")

save(tme_combine,
     file = paste0("03-tme/", index, "-1-", ProjectID, "-TME-Cell-fraction.RData"))


# ===========================================================================
# Step 2: Signature scoring — IOBR collection (PCA method)
# ===========================================================================
sig_res <- calculate_sig_score(
  pdata          = NULL,
  eset           = eset,
  adjust_eset    = TRUE,
  signature      = signature_collection,   # 323 curated signatures
  method         = "PCA",
  mini_gene_count = 2
)

save(sig_res,
     file = paste0("03-tme/", index, "-2-", ProjectID, "-Signature-score-mycollection.RData"))


# ===========================================================================
# Step 3: Pathway scoring — Hallmark + GO + KEGG + Reactome (ssGSEA)
# ===========================================================================
# load_data() downloads and caches signature files from GitHub on first use
hallmark  <- IOBR::load_data("hallmark")    # 50 pathways
go_bp     <- IOBR::load_data("go_bp")       # 7,658 gene sets
go_cc     <- IOBR::load_data("go_cc")
go_mf     <- IOBR::load_data("go_mf")
kegg      <- IOBR::load_data("kegg")        # 186 pathways
reactome  <- IOBR::load_data("reactome")    # 1,615 pathways

sig_go_kegg <- calculate_sig_score(
  pdata          = NULL,
  eset           = eset,
  signature      = c(hallmark, go_bp, go_cc, go_mf, kegg, reactome),
  method         = "ssgsea",
  mini_gene_count = 3,
  parallel.size  = 1
)

save(sig_go_kegg,
     file = paste0("03-tme/", index, "-3-", ProjectID,
                   "-Signature-score-Hallmark-GO-KEGG.RData"))


# ===========================================================================
# Step 4: Merge all results into one master table
# ===========================================================================
tme_sig_combine <- tme_combine %>%
  inner_join(., sig_res,     by = "ID") %>%
  inner_join(., sig_go_kegg, by = "ID")

save(tme_sig_combine,
     file = paste0("03-tme/", index, "-0-", ProjectID, "-Merge-all-signature.RData"))

cat(paste0(">>>>> Pipeline complete: ", ProjectID, "\n"))
cat(paste0(">>>>> TME methods: 8\n"))
cat(paste0(">>>>> Signature scores: ", ncol(sig_res) + ncol(sig_go_kegg) - 1, "\n"))
cat(paste0(">>>>> Combined columns: ", ncol(tme_sig_combine), "\n"))
