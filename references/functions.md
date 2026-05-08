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
calculate_sig_score(pdata, eset, signature, method = "ssgsea",
                    mini_gene_count = 3, adjust_eset = TRUE)
```

### `format_signatures()`
Convert external gene sets to IOBR signature format.
```r
format_signatures(file, sheet = 1)
```
- Input: CSV or Excel with columns: signature_name, gene

### `tme_cluster()`
Cluster samples based on TME cell composition.
```r
tme_cluster(input, features, id = "sample_id", scale = TRUE,
            method = "kmeans", max.nc = 6)
```
- `method`: `"kmeans"`, `"hclust"`, `"consensus"`

### `LR_cal()`
Ligand-receptor interaction analysis.
```r
LR_cal(eset, pdata, group, species = "human")
```

---

## Statistical Functions

### `find_markers_in_bulk()`
Differential expression between groups.
```r
find_markers_in_bulk(pdata, eset, group, nfeatures = 50,
                     top_n = 100, thresh.use = 0.25)
```

### `iobr_deg()`
Differential expression analysis (limma or DESeq2).
```r
iobr_deg(eset, pdata, group_id, array = TRUE,
         padj_cutoff = 0.05, logfc_cutoff = 1)
```

### `sig_gsea()`
Gene Set Enrichment Analysis.
```r
sig_gsea(genesets, gene_symbol, logfc, org = "hsa",
         msigdb = TRUE, category = "H", subcategory = NULL)
```

### `batch_surv()`
Batch survival analysis across multiple features.
```r
batch_surv(pdata, tme_data, time, status, method = "cox")
```

### `batch_cor()`
Batch correlation analysis.
```r
batch_cor(x, y, method = "spearman")
```

### `find_mutations()`
Identify mutations associated with TME signatures.
```r
find_mutations(mutation_matrix, signature_matrix, signature = NULL,
               method = "wilcoxon", min_mut_freq = 0.02)
```

---

## Visualization Functions

### `sig_box()`
Boxplots with statistical testing.
```r
sig_box(pdata, signature, variable, palette_group = "jama",
        show_pvalue = TRUE, method = "wilcox.test")
```

### `sig_heatmap()`
Comprehensive heatmap with annotations.
```r
sig_heatmap(pdata, signature, scale = TRUE, palette = 2,
            palette_group = "jama", size_row = 8,
            custom_cols = NULL, custom_heatmap_cols = NULL)
```

### `percent_bar_plot()`
Stacked percentage bar plot.
```r
percent_bar_plot(data, palette_group = "jama")
```

### `cell_bar_plot()`
Cell fraction bar plot.
```r
cell_bar_plot(data, palette_group = "jama")
```

### `get_cor()`
Scatter plot with correlation and regression.
```r
get_cor(x, y, method = "spearman", palette = "jama")
```

### `get_cor_matrix()`
Correlation matrix heatmap.
```r
get_cor_matrix(data, method = "spearman", palette = 2)
```

### `sig_forest()`
Forest plot for Cox regression.
```r
sig_forest(pdata, signature, time, status)
```

### `sig_surv_plot()`
Kaplan-Meier survival curves.
```r
sig_surv_plot(pdata, signature, time, status, group = NULL,
              palette_group = "jama")
```

### `roc_time()`
Time-dependent ROC curves.
```r
roc_time(pdata, signature, time, status,
         time_point = c(12, 36, 60))
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
