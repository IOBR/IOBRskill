# IOBR Built-in Data Reference

This document catalogs all built-in data shipped with the IOBR package (v2.2.0).
These datasets are available without any external download — use them directly.

---

## 1. Annotation / Probe Data

### 1.1 `extdata/genes.txt` — Gene Length & ID Mapping
- **Purpose**: Used by `count2tpm()` for TPM normalization
- **Columns**: `HUGO.symbols`, `Cell.population`, `ENTREZID`, `ENSEMBL.ID`
- **Rows**: 111 genes (reference gene set for cell populations)
- **Usage**: Automatically loaded by `count2tpm()`, no manual action needed

### 1.2 `extdata/probesets.txt` — Probe Annotations for Microarray
- **Purpose**: Probe-to-cell-type mapping for microarray deconvolution
- **Columns**: Probe ID, Cell population
- **Rows**: 176 probes
- **Usage**: Used internally by MCPcounter and related functions

### 1.3 `IOBR:::mcp_probesets` — MCPcounter Probe Sets
- **Purpose**: Pre-compiled probe sets for MCPcounter deconvolution
- **Type**: data.frame
- **Access**: `IOBR:::mcp_probesets` (internal)

### 1.4 `IOBR:::Top_probe` — Top Probes for Cell Types
- **Purpose**: Top discriminative probes for cell type identification
- **Type**: data.frame with probe IDs and cell type labels
- **Access**: `IOBR:::Top_probe` (internal)

---

## 2. Signature Gene Sets

### 2.1 `IOBR:::signature_tme` — TME Signatures (186 sets)
- **Purpose**: Core TME-related gene signature collection
- **Access**: `IOBR:::signature_tme` (triple colon — internal object)
- **Contents**: 186 gene signatures covering:
  - T cell signatures (CD8 effector, Tregs, exhaustion)
  - B cell signatures
  - Macrophage polarization (M1/M2)
  - NK cell signatures
  - EMT signatures (EMT1, EMT2, EMT3)
  - Immune checkpoint genes
  - DDR, cell cycle, mismatch repair
  - TME scores (TMEscoreA/B_CIR)
  - Signatures from Rooney et al, Bindea et al, Li et al, Peng et al
  - Single-cell derived signatures (Chu et al, Xue et al, Kieffer et al)
  - CAF subtypes (myCAF, iCAF)

### 2.2 `signature_collection` — Comprehensive Collection (323 sets)
- **Purpose**: Full curated gene signature collection
- **Access**: `data("signature_collection")` then use `signature_collection`
- **Contains all signatures from `signature_tme` plus**:
  - Hallmark gene sets (50 pathways) — groups: `hallmark`, `hallmark1`, `hallmark2`, `hallmark3`
  - Metabolism signatures (113 sets) — fatty acid, hypoxia, cholesterol, etc.
  - Additional cell-type-specific signatures
  - CAF single-cell signatures

### 2.3 `sig_group` — Signature Groupings
- **Purpose**: Categorization of signatures by biological theme
- **Access**: `data("sig_group")` then use `sig_group`
- **Categories**:

| Category | Count | Description |
|----------|-------|-------------|
| `tumor_signature` | 12 | Tumor biology (DDR, cell cycle, EMT) |
| `EMT` | 5 | Epithelial-mesenchymal transition |
| `io_biomarkers` | 11 | Immuno-oncology biomarkers |
| `immu_microenvironment` | 10 | Immune microenvironment scores |
| `immu_suppression` | 8 | Immunosuppressive signatures |
| `immu_exclusion` | 12 | Immune exclusion mechanisms |
| `immu_exhaustion` | 7 | T cell exhaustion |
| `TCR_BCR` | 5 | TCR/BCR signaling |
| `tme_signatures1` | 41 | TME signatures subset 1 |
| `tme_signatures2` | 94 | TME signatures subset 2 |
| `Bcells` | 13 | B cell signatures across methods |
| `Tcells` | 64 | T cell signatures across methods |
| `DCs` | 14 | Dendritic cell signatures |
| `Macrophages` | 22 | Macrophage signatures (M0/M1/M2) |
| `Neutrophils` | 18 | Neutrophil signatures |
| `Monocytes` | 11 | Monocyte/MDSC signatures |
| `CAFs` | 12 | Cancer-associated fibroblast signatures |
| `NK` | 10 | NK cell signatures |
| `tme_cell_types` | 164 | All cell type signatures combined |
| `CIBERSORT` | 22 | CIBERSORT LM22 cell types |
| `MCPcounter` | 10 | MCPcounter populations |
| `EPIC` | 8 | EPIC cell types |
| `xCell` | 67 | xCell cell types |
| `quanTIseq` | 11 | quanTIseq cell types |
| `ESTIMATE` | 4 | ESTIMATE scores |
| `IPS` | 6 | Immunophenoscore axes |
| `TIMER` | 6 | TIMER cell types |
| `hallmark` | 50 | MSigDB Hallmark (50 pathways) |
| `Metabolism` | 113 | Metabolism pathways |
| `Rooney_et_al` | 15 | Rooney et al immune signatures |
| `Bindea_et_al` | 24 | Bindea et al immune signatures |
| `Li_et_al` | 17 | Li et al immune pathways |
| `Peng_et_al` | 7 | Peng et al ICB signatures |
| `CAF_singleCell` | 9 | CAF single-cell subtypes |
| `Chu_et_al` | 26 | Chu et al single-cell CD4/CD8 |
| `Xue_et_al` | 33 | Xue et al macrophage single-cell |
| `Kieffer_et_al` | 8 | Kieffer et al CAF subtypes |

---

## 3. Example / Toy Datasets

| Dataset | Access | Description |
|---------|--------|-------------|
| `imvigor210_pdata` | `data("imvigor210_pdata")` | IMvigor210 bladder cancer immunotherapy cohort (phenotype data) |
| `pdata_stad` | `data("pdata_stad")` | TCGA-STAD gastric cancer toy phenotype data |
| `tcga_stad_pdata` | `data("tcga_stad_pdata")` | TCGA-STAD clinical and molecular annotation |
| `stad_group` | `data("stad_group")` | TCGA-STAD cancer type grouping vector |
| `subgroup_data` | `data("subgroup_data")` | Example dataset for subgroup survival analysis |
| `null_models` | `data("null_models")` | NULL model coefficients for MCPcounter |

---

## 4. Downloaded Reference Data (auto-cached)

IOBR automatically downloads and caches these reference files on first use:

| File | Function | Description |
|------|----------|-------------|
| `lm22.rda` | `deconvo_tme(method="cibersort")` | CIBERSORT LM22 signature matrix (22 immune cell types) |
| `common_genes.rda` | `deconvo_tme(method="estimate")` | Common gene set for ESTIMATE scoring |
| `SI_geneset.rda` | `deconvo_tme(method="estimate")` | Stromal/Immune gene signatures for ESTIMATE |

Cached location: `~/.cache/iobr_data/` or temp directory.

---

## 4b. On-demand Pathway Data via `load_data()`

IOBR provides `IOBR::load_data()` to download and cache pathway gene sets from GitHub on first use. These are essential for comprehensive pathway scoring and GSEA.

```r
hallmark  <- IOBR::load_data("hallmark")   # 50 pathways
go_bp     <- IOBR::load_data("go_bp")      # 7,658 gene sets
go_cc     <- IOBR::load_data("go_cc")      # Cellular component
go_mf     <- IOBR::load_data("go_mf")      # Molecular function
kegg      <- IOBR::load_data("kegg")       # 186 pathways
reactome  <- IOBR::load_data("reactome")   # 1,615 pathways
```

| Dataset | Gene Sets | Description |
|---------|-----------|-------------|
| `"hallmark"` | 50 | MSigDB Hallmark gene sets — well-curated biological pathways |
| `"go_bp"` | 7,658 | Gene Ontology Biological Process |
| `"go_cc"` | ~1,000 | Gene Ontology Cellular Component |
| `"go_mf"` | ~1,700 | Gene Ontology Molecular Function |
| `"kegg"` | 186 | KEGG pathway gene sets |
| `"reactome"` | 1,615 | Reactome pathway gene sets |

**Usage pattern with `calculate_sig_score()`:**
```r
hallmark <- IOBR::load_data("hallmark")
kegg     <- IOBR::load_data("kegg")

sig_res <- calculate_sig_score(
  pdata = NULL, eset = eset,
  signature = c(hallmark, kegg),
  method = "ssgsea",
  mini_gene_count = 3,
  adjust_eset = TRUE
)
```

**Note**: First call to `load_data()` for each dataset requires internet access to download from GitHub. Subsequent calls use the cached version.

---

## 5. Quick Reference: How to Access Signatures

```r
library(IOBR)

# TME signatures (internal — use :::)
sig_tme <- IOBR:::signature_tme       # 186 gene sets
length(sig_tme)                        # 186
names(sig_tme)[1:5]                    # CD_8_T_effector, DDR, APM, ...

# Full collection (lazy-loaded with data())
data("signature_collection")           # 323 gene sets
data("sig_group")                      # Category groupings

# Get signatures by category
sig_group$CIBERSORT                    # 22 CIBERSORT cell type column names
sig_group$Macrophages                  # 22 macrophage-related signatures
sig_group$hallmark                     # 50 Hallmark pathway names

# Select specific signatures for scoring
my_sigs <- signature_collection[sig_group$io_biomarkers]  # IO biomarker signatures
```
