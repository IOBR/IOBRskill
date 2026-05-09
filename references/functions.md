# IOBR Core Function Reference

## Data Preprocessing Functions

### `count2tpm()`
Convert raw count matrix to TPM.
```r
count2tpm(countMat, idType = "Ensembl", org = "hsa", source = "local")
```
- `countMat`: Gene × sample count matrix. Row names = gene IDs.
- `idType`: `"Ensembl"`, `"Entrez"`, or `"Symbol"`
- `org`: `"hsa"` (human) or `"mmu"` (mouse)
- `source`: `"local"` (cached gene length) or `"online"` (download)

### `log2eset()`
Auto-detect and apply log2 transformation if data is not already log-scaled.
```r
log2eset(eset)
```

### `anno_eset()`
Annotate expression matrix with gene symbols.
```r
anno_eset(eset, annotation, symbol = "gene_symbol", probe = "probe_id", method = "mean")
```
- `eset`: Expression matrix or ExpressionSet
- `annotation`: Data frame with probe-to-gene mapping
- `method`: `"mean"` or `"sd"` for handling duplicate probes

### `remove_batcheffect()`
Remove batch effects across datasets.
```r
remove_batcheffect(eset_list, pdata_list, idType = "ensembl", dataType = "count")
```
- `eset_list`: List of expression matrices
- `dataType`: `"Array"`, `"Count"`, or `"Tpm"`

### `find_outlier_samples()`
Detect outlier samples using WGCNA connectivity.
```r
find_outlier_samples(eset, plot = TRUE)
```

### `iobr_pca()`
PCA visualization for batch effect assessment.
```r
iobr_pca(pdata, eset, group, scale = TRUE)
```

### `mouse2human_eset()`
Convert mouse gene symbols to human orthologs.
```r
mouse2human_eset(eset)
```

---

## TME Deconvolution Functions

### `deconvo_tme()`
Main deconvolution wrapper. Supports 10 methods.
```r
deconvo_tme(eset, method, arrays = FALSE, perm = 1000,
            absolute.mode = FALSE, reference = NULL,
            tumor = TRUE, scale_mrna = TRUE,
            group_list = NULL, platform = "affymetrix")
```

Method-specific parameters:

| Method | Key Parameters | Notes |
|--------|---------------|-------|
| `"cibersort"` | `arrays`, `perm` (100–10000), `absolute` | `perm=1000` recommended |
| `"mcpcounter"` | None extra | Simplest to run |
| `"epic"` | `tumor`, `scale_mrna` | Estimates cancer cell fraction |
| `"xcell"` | `arrays` | 64 cell types, can be slow |
| `"estimate"` | `platform` | Returns stromal/immune/tumor purity scores |
| `"timer"` | `group_list` (required) | Must provide cancer type per sample |
| `"quantiseq"` | `arrays`, `tumor` | M1/M2 macrophage resolution |
| `"ips"` | None extra | 4-axis immunophenoscore |
| `"integration"` | All applicable | Runs multiple methods together |

---

## Signature Analysis Functions

### `calculate_sig_score()`
Calculate gene signature scores per sample.
```r
calculate_sig_score(pdata = NULL, eset, signature, method = "ssgsea",
                    mini_gene_count = 3, adjust_eset = TRUE)
```

### `format_signatures()`
Convert external gene sets to IOBR signature format.
```r
format_signatures(file, sheet = 1)
```
- Input: CSV or Excel with columns: signature_name, gene

### `format_msigdb()`
Convert MSigDB GMT files to IOBR signature format for `calculate_sig_score()`.
```r
format_msigdb(gmt_file)
```

### `get_sig_sc()`
Extract top markers from single-cell DE results for `calculate_sig_score()`.
```r
get_sig_sc(sc_markers_matrix, ...)
```

### `tme_cluster()`
Cluster samples based on TME cell composition.
```r
tme_cluster(input, features, id = "sample_id", scale = TRUE,
            method = "kmeans", max.nc = 6)
```
- `method`: `"kmeans"`, `"hclust"`, `"consensus"`

---

## Batch Statistical Analysis Functions

### `batch_surv()`
Batch survival analysis — Cox regression for multiple variables.
```r
batch_surv(pdata, variable, time = "time", status = "status", best_cutoff = FALSE)
```
- `variable`: Vector of **column names** in pdata to test. Must merge pdata + TME first.
- Returns data frame with columns: `P`, `HR`, `CI_low_0.95`, `CI_up_0.95`, `ID`
- Returns hazard ratios, confidence intervals, p-values for all variables.

### `subgroup_survival()`
Subgroup Cox analysis — extract HR/CI within subgroups.
```r
subgroup_survival(pdata, variable, group, time, status)
```

### `batch_cor()`
Batch correlation analysis — correlates one target variable against multiple features.
```r
batch_cor(data, target, feature, method = c("spearman", "pearson", "kendall"))
```
- `target`: Column name (string) of the target variable to correlate against
- `feature`: Vector of column names to correlate with target
- `method`: `"spearman"` (default), `"pearson"`, or `"kendall"`
- Returns tibble: `sig_names`, `p.value`, `statistic`, `p.adj`, `log10pvalue`, `stars`
- **NOT** a pairwise matrix — use `cor()` or `get_cor_matrix()` for full pairwise

### `batch_pcc()`
Batch partial correlation coefficient — control for a third variable.
```r
batch_pcc(data, target, variable, control)
```

### `batch_wilcoxon()`
Batch Wilcoxon rank-sum test — compare feature distributions between two groups.
```r
batch_wilcoxon(data, target = "group", feature = NULL, feature_manipulation = FALSE)
```
- `target`: Column name (string) of the grouping variable. Must have exactly 2 levels.
- `feature`: Character vector of feature column names to test. If NULL, interactive selection.
- Returns tibble: `sig_names`, `p.value`, `statistic`, `p.adj`, `log10pvalue`, `stars`

### `iobr_deg()`
Differential expression analysis (limma or DESeq2), with built-in volcano and heatmap.
```r
iobr_deg(eset, pdata, group_id, array = TRUE,
         padj_cutoff = 0.05, logfc_cutoff = 1)
```

### `find_markers_in_bulk()`
Multi-group marker finding via Seurat methods.
```r
find_markers_in_bulk(pdata, eset, group, nfeatures = 50,
                     top_n = 100, thresh.use = 0.25)
```

### `sig_gsea()`
Gene Set Enrichment Analysis via fgsea.
```r
sig_gsea(genesets, gene_symbol, logfc, org = "hsa",
         msigdb = TRUE, category = "H", subcategory = NULL)
```

### `LR_cal()`
Ligand-receptor interaction analysis.
```r
LR_cal(eset, pdata, group, species = "human")
```

### `find_mutations()`
Identify mutations associated with TME signatures.
```r
find_mutations(mutation_matrix, signature_matrix, signature = NULL,
               method = "wilcoxon", min_mut_freq = 0.02)
```

### `feature_manipulation()`
Clean features: remove missing values, outliers, zero-variance variables.
```r
feature_manipulation(data, ...)
```

---

## Visualization Functions

### TME Composition

### `cell_bar_plot()`
Stacked bar plot for cell fractions (CIBERSORT/EPIC/quanTIseq).
```r
cell_bar_plot(input, id = "ID", title = "Cell Fraction", features = NULL,
              pattern = NULL, legend.position = "bottom", coord_flip = TRUE,
              palette = 3, show_col = FALSE, cols = NULL)
```
- Use `features` to specify cell columns; filter out QC columns (P.value/Correlation/RMSE).

### `percent_bar_plot()`
Stacked percentage bar plot.
```r
percent_bar_plot(input, x, y, subset.x = NULL, color = NULL, palette = NULL,
                 title = NULL, axis_angle = 0, coord_flip = FALSE,
                 add_Freq = TRUE, add_sum = TRUE, round.num = 2)
```

### Boxplots

### `sig_box()`
Boxplot with statistical comparisons.
```r
sig_box(data, signature, variable, palette = "nrc", cols = NULL,
        jitter = FALSE, point_size = 5, show_pairwise_p = TRUE)
```
- No `method` parameter — test is auto-selected. Use `palette` (not `palette_group`).

### Heatmaps

### `sig_heatmap()`
Comprehensive heatmap with grouping, scaling, annotation. Requires `tidyHeatmap` package.
```r
sig_heatmap(input, id = "ID", features, group, condition = NULL,
            scale = FALSE, palette = 2, palette_group = "jama",
            size_row = 8, path = NULL, index = NULL)
```
- `features`: Character vector of feature column names to plot (REQUIRED)
- `group`: **Column name** (string) of grouping variable in `input`, NOT a vector. E.g., `"TME_subtype"`.
- `path`: Directory to save output (NULL = display only)
- `index`: Filename prefix for saved output
- Requires: `install.packages("tidyHeatmap")`

### Correlation

### `get_cor()`
Scatter plot with correlation and regression line for a single pair.
```r
get_cor(data, variable1, variable2, method = "spearman")
```

### `get_cor_matrix()`
Correlation matrix heatmap between two variable sets.
```r
get_cor_matrix(data, variable1, variable2, method = "spearman", palette = 2)
```

### `iobr_cor_plot()`
Batch correlation visualization — signature vs signature/phenotype.
```r
iobr_cor_plot(data, signature, target, method = "spearman")
```

### Survival

### `sig_surv_plot()`
Single Kaplan-Meier curve for one signature/gene.
```r
sig_surv_plot(input_pdata, signature, project = "KM", ID = "ID",
              time = "time", status = "status", time_type = "month",
              break_month = "auto", cols = NULL, palette = "jama",
              mini_sig = "score", fig.type = "png", save_path = NULL, index = 1)
```
- `mini_sig`: `"score"` (continuous, auto-split by median) or `"category"` (already grouped)
- `time_type`: `"day"`, `"month"`, or `"year"`

### `batch_sig_surv_plot()`
Batch KM curves for multiple signatures, optionally across multiple projects.
```r
batch_sig_surv_plot(input_pdata, signature, id = "ID",
                     column_of_project = "ProjectID", project = NULL,
                     time = "time", status = "status", time_type = "day",
                     break_month = "auto", palette = "jama", cols = NULL,
                     mini_sig = "score", save_path, show_col = TRUE, fig_type = "pdf")
```

### `surv_group()`
KM curve comparing high/low groups with risk table.
```r
surv_group(input_pdata, target_group, ID = "ID", levels = c("High", "Low"),
           reference_group = NULL, project = NULL, time = "time", status = "status",
           time_type = "month", break_month = "auto", cols = NULL, palette = "jama",
           mini_sig = "score", save_path = NULL, fig.type = "pdf", index = 1,
           width = 6, height = 6.5, font.size.table = 3)
```

### `sig_forest()`
Forest plot for `batch_surv()` results.
```r
sig_forest(data, signature, pvalue = "P", HR = "HR",
           CI_low_0.95 = "CI_low_0.95", CI_up_0.95 = "CI_up_0.95",
           n = 10, max_character = 25, discrete_width = 35,
           color_option = 1, cols = NULL, text.size = 13)
```
- `data`: Output from `batch_surv()` — must contain columns matching `signature`, `pvalue`, `HR`, etc.
- `signature`: **Column name** in data containing variable names (e.g., `"ID"`)

### `roc_time()`
Time-dependent ROC curves with AUC annotation at multiple time points.
```r
roc_time(input, vars, time = "time", status = "status", time_point = 12,
         time_type = "month", palette = "jama", cols = "normal", seed = 1234,
         show_col = FALSE, path = NULL, main = "PFS", index = 1,
         fig.type = "pdf", width = 5, height = 5.2)
```

### `sig_roc()`
Multiple ROC curves for binary response prediction.
```r
sig_roc(data, response, variables, fig.path = NULL, main = NULL,
        file.name = NULL, palette = "jama", cols = NULL, alpha = 1,
        compare = FALSE, smooth = TRUE, compare_method = "bootstrap", boot.n = 100)
```

### Dimensionality Reduction

### `iobr_pca()` (visualization mode)
PCA scatter plot.
```r
iobr_pca(pdata, eset, group, scale = TRUE)
```

### DEG / Volcano

### `iobr_deg()` (with visualization)
DEG with built-in volcano plot and heatmap output.
```r
iobr_deg(eset, pdata, group_id, array = TRUE,
         padj_cutoff = 0.05, logfc_cutoff = 1)
```

---

## Modeling Functions

### `PrognosticModel()`
Build survival prediction model with Lasso/Ridge.
```r
PrognosticModel(x, y, scale = TRUE, seed = 42,
                train_ratio = 0.7, nfold = 5)
```

### `BinomialModel()`
Binary classification model.
```r
BinomialModel(x, y, scale = TRUE, seed = 42,
              train_ratio = 0.7, nfold = 5)
```

---

## Custom Helper Functions

### `find_top_feature_per_subtype()`
For each TME subtype, find the cell type with highest mean scaled z-score. Used in Fig09 to select the most representative feature per subtype for boxplot visualization.

```r
find_top_feature_per_subtype <- function(data, cell_cols, subtype_col = "TME_subtype") {
  subtypes <- sort(unique(data[[subtype_col]]))
  result <- character(0)
  for (st in subtypes) {
    st_data <- data[data[[subtype_col]] == st, cell_cols, drop = FALSE]
    cell_means <- colMeans(st_data, na.rm = TRUE)
    result[st] <- names(sort(cell_means, decreasing = TRUE))[1]
  }
  return(result)
}
```
- `data`: scaled data frame with subtype column and cell columns
- `cell_cols`: character vector of cell column names to search
- `subtype_col`: name of the subtype grouping column (default `"TME_subtype"`)
- Returns: named character vector — `c(C1 = "top_cell_col", C2 = "top_cell_col", ...)`
