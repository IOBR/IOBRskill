# IOBRskill

A Claude Code skill for end-to-end tumor microenvironment (TME) analysis of bulk transcriptome data using the [IOBR](https://github.com/IOBR/IOBR) R package.

IOBRskill automates the complete TME analysis pipeline — from data preprocessing and annotation through immune cell deconvolution, signature scoring, statistical analysis, and publication-quality visualization — all through natural language interaction with Claude Code.

![IOBRskill Pipeline Overview](scripts/IOBRskill.png)

## What It Does

IOBRskill guides Claude Code through a standardized 6-phase TME analysis workflow:

| Phase | Step | Description |
|-------|------|-------------|
| 0 | Plan | Generate `IOBR-pipeline.md` — ASCII tree of all scripts, inputs/outputs, and expected results |
| 1 | `01-data_preprocessing.R` | Data loading, normalization, probe annotation, QC, batch correction, sample matching, pdata summary |
| 2 | `02-tme_deconvolution.R` | Immune cell deconvolution (8 methods: CIBERSORT, CIBERSORT-ABS, MCPcounter, EPIC, xCell, ESTIMATE, TIMER, quanTIseq, IPS) |
| 3 | `03-signature_analysis.R` | Gene signature scoring (ssGSEA/PCA/Z-score), Hallmark/KEGG/GO/Reactome pathway scoring, TME subtype clustering |
| 4 | `04-statistical_analysis.R` | Merge → scale → cluster → Wilcoxon/Kruskal tests → Cox survival → correlation matrix. All results saved to `04-figs/data/` |
| 5 | `05-visualization.R` | Nature-style publication figures: Fig01–Fig09 (barplot, heatmap, forest, correlation, boxplot, KM) — PNG 300dpi + PDF |
| 6 | Wrap-Up Note | Generate `05-note/IOBR-analysis-README.md` with actual output tree, figure count, and key findings summary |

### Decision Points

IOBRskill pauses at key steps to let the user choose:

1. **Phase 1**: Data type (Count/TPM/Array) and species (human/mouse)
2. **Phase 2**: Deconvolution method(s) — CIBERSORT + MCPcounter + ESTIMATE recommended
3. **Phase 3**: Scoring method (ssGSEA default) and signature collection
4. **Phase 5**: Color palette preference

## Requirements

- [Claude Code](https://docs.anthropic.com/en/docs/claude-code) CLI (latest version recommended)
- R >= 4.2.0
- [IOBR](https://github.com/IOBR/IOBR) R package >= 2.2.0

## Installation

### Step 1 -- Install the Skill

Choose **one** of the following methods:

<details>
<summary><b>Method A: Clone to global skills directory (Recommended)</b></summary>

This makes IOBRskill available in **all** your Claude Code sessions:

```bash
git clone https://github.com/IOBR/IOBRskill.git ~/.claude/skills/IOBRskill
```

That's it. Claude Code automatically discovers skills under `~/.claude/skills/`.
</details>

<details>
<summary><b>Method B: Install as a project-level skill</b></summary>

Install into a specific project so the skill only activates in that workspace:

```bash
cd /path/to/your/project
git clone https://github.com/IOBR/IOBRskill.git
```

Or add it as a git submodule:

```bash
cd /path/to/your/project
git submodule add https://github.com/IOBR/IOBRskill.git IOBRskill
```
</details>

<details>
<summary><b>Method C: Manual copy</b></summary>

Download the repository and copy the `IOBRskill/` folder to your preferred location:

```bash
# Download
wget https://github.com/IOBR/IOBRskill/archive/refs/heads/main.zip
unzip main.zip

# Copy to global skills
cp -r IOBRskill-main/ ~/.claude/skills/IOBRskill

# Or copy to a specific project
cp -r IOBRskill-main/ /path/to/your/project/IOBRskill
```
</details>

### Step 2 -- Install IOBR R Package

IOBRskill depends on the IOBR R package. Install it in your R environment:

```r
# Install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install IOBR from GitHub (latest version)
BiocManager::install("IOBR/IOBR")

# Or install from CRAN (stable release)
install.packages("IOBR")

# Verify
packageVersion("IOBR")
# Should print: [1] '2.2.0' or higher
```

If you plan to use CIBERSORT deconvolution, also install:

```r
BiocManager::install("preprocessCore")
```

### Step 3 -- Verify

Open Claude Code and test the skill:

```
/IOBRskill
```

You should see the IOBR analysis pipeline launched. If you get an IOBR not found error, run:

```bash
Rscript -e 'library(IOBR); cat("IOBR", as.character(packageVersion("IOBR")), "OK\n")'
```

## Trigger Conditions

IOBRskill activates when you:

- Type `/IOBRskill` directly
- Mention any of these keywords in your prompt:
  - **English**: tumor microenvironment analysis, TME deconvolution, immune infiltration, ligand-receptor analysis, CIBERSORT, MCPcounter, ESTIMATE, immune cell deconvolution
  - **Chinese**: 肿瘤微环境分析, 肿瘤微环境解析, 免疫浸润分析, 受体配体分析, 通路分析

## How to Use

Simply describe what you want to analyze. IOBRskill will interactively guide you through the pipeline:

### Example prompts

```
# Start a full TME analysis
/IOBRskill

# With your own data
"I have a raw count matrix at ~/data/expr.csv from 50 human tumor samples.
I want to do a comprehensive TME analysis."

# Specific method
"Run CIBERSORT and MCPcounter deconvolution on my expression matrix"

# Focus on specific analysis
"I want to analyze ligand-receptor interactions in my TME data"
```

## Output Structure

IOBRskill creates a standardized directory structure. Here is an example from IMvigor210 (bladder cancer immunotherapy cohort, 348 samples, RNA-seq):

```
IMvigor210/
├── 01-script/
│   ├── 01-data_preprocessing.R       # Count → TPM, gene annotation, log2, QC
│   ├── 02-tme_deconvolution.R        # 9-method deconvolution + merge
│   ├── 03-signature_analysis.R       # ssGSEA: TME signatures + Hallmark + KEGG
│   ├── 04-statistical_analysis.R     # Merge, scale, cluster, Wilcoxon, survival
│   └── 05-visualization.R           # Fig01–Fig09 (all figures)
├── 02-input/
│   ├── annotated_eset.csv            # Log2-normalized expression (genes x samples)
│   ├── pdata.csv                     # Sample phenotype data
│   └── pdata_summary.csv             # Variable type/level summary
├── 03-tme/
│   ├── tme_cibersort.csv             # 22 immune cell fractions per sample
│   ├── tme_mcpcounter.csv            # 8 cell population scores per sample
│   ├── tme_estimate.csv              # Stromal/Immune/Tumor purity scores
│   ├── tme_epic.csv                  # 8 cell types incl. cancer fraction
│   ├── tme_xcell.csv                 # 64 cell type scores
│   ├── tme_quantiseq.csv             # 10 cell type fractions
│   ├── tme_ips.csv                   # 4 immunotherapy response axes
│   ├── sig_score_tme.csv             # TME signature scores (ssGSEA)
│   ├── sig_score_pathway.csv         # Hallmark + GO + KEGG + Reactome scores
│   └── tme_sig_combine.csv           # All results merged by sample ID
├── 04-figs/
│   ├── Fig01-barplot_cibersort.png/pdf        # CIBERSORT cell composition barplot
│   ├── Fig02-barplot_epic.png/pdf             # EPIC cell composition barplot
│   ├── Fig03a-heatmap_cibersort.png/pdf       # CIBERSORT heatmap (z-score)
│   ├── Fig03b-heatmap_mcpcounter.png/pdf      # MCPcounter heatmap (z-score)
│   ├── Fig03c-heatmap_xcell.png/pdf           # xCell heatmap (z-score)
│   ├── Fig04a-forest_cibersort.png/pdf        # CIBERSORT survival forest (if survival)
│   ├── Fig04b-forest_other_tme.png/pdf        # MCPcounter+EPIC+xCell forest (if survival)
│   ├── Fig04c-forest_tme_signature.png/pdf    # TME signatures forest (if survival)
│   ├── Fig04d-forest_go_kegg.png/pdf          # GO/KEGG/Hallmark forest (if survival)
│   ├── Fig05-cor_matrix.png/pdf               # TME cell correlation matrix heatmap
│   ├── Fig06a-heatmap_cibersort_subtype.png/pdf  # CIBERSORT heatmap by TME subtype
│   ├── Fig06b-heatmap_mcpcounter_subtype.png/pdf # MCPcounter heatmap by TME subtype
│   ├── Fig06c-heatmap_xcell_subtype.png/pdf   # xCell heatmap by TME subtype
│   ├── Fig06d-km_subtype.png/pdf              # KM plot by TME subtype (if survival)
│   ├── Fig07a-boxplot_cibersort_subtype.png/pdf   # CIBERSORT subtype boxplots
│   ├── Fig07b-boxplot_mcpcounter_subtype.png/pdf  # MCPcounter subtype boxplots
│   ├── Fig07c-boxplot_xcell_subtype.png/pdf   # xCell subtype boxplots
│   ├── Fig07d-boxplot_tme_signature_subtype.png/pdf # TME signatures subtype boxplots
│   ├── Fig07e-boxplot_pathway_subtype.png/pdf  # GO/KEGG/Hallmark subtype boxplots
│   ├── Fig08a-top10_<variable>.png/pdf        # Top 10 differential boxplots (if categorical)
│   ├── Fig08b-top10_<variable>.png/pdf        # Top 10 differential boxplots (if categorical)
│   ├── kmplot/                                # KM plots (most/least significant variables)
│   └── data/                                  # Statistical result tables
│       ├── 01-tme_pdata_merged.csv
│       ├── 02-tme_scaled.csv
│       ├── 03-tme_subtype.csv
│       ├── 04-diff_<variable>.csv
│       ├── 05-surv_cibersort.csv
│       ├── 05-surv_other_tme.csv
│       ├── 05-surv_tme_signatures.csv
│       ├── 05-surv_pathways.csv
│       └── 06-cor_matrix.csv
├── 05-note/
│   ├── IOBR-pipeline.md               # Analysis plan (generated in Phase 0)
│   └── IOBR-analysis-README.md        # Output tree + key findings summary (Phase 6)
└── 06-log/
    ├── 01-data_preprocessing.log      # Execution logs with IOBR citation header
    ├── 02-tme_deconvolution.log
    ├── 03-signature_analysis.log
    ├── 04-statistical_analysis.log
    └── 05-visualization.log
```

### Figure Logic

| Condition | Figure | Content |
|-----------|--------|---------|
| Always | Fig01–02 | CIBERSORT + EPIC fraction barplots |
| Always | Fig03a–c | Heatmaps without subtype (CIBERSORT, MCPcounter, xCell) |
| IF survival data | Fig04a–d | Forest plots (CIBERSORT, other TME, signatures, pathways) |
| Always | Fig05 | TME cell correlation matrix |
| Always | Fig06a–c | Heatmaps by TME subtype |
| IF survival data | Fig06d | KM plot by TME subtype |
| Always | Fig07a–e | Most positive cell per subtype boxplots |
| IF categorical vars | Fig08a–b... | Top 10 differential boxplots (Wilcoxon/Kruskal) |
| IF survival data | Fig09/kmplot/ | 8 KM plots (extreme HR from each group) |

### Key Rules

- **sig_forest**: Must filter `abs(HR) > 10` before passing data -- extreme HR causes blank plots
- **ALL boxplots**: Use `sig_box()` with fixed params: `palette = "paired1"`, `show_pvalue = TRUE`, `angle_x_text = 60`, `hjust = 1`, `size_of_pvalue = 3`, `size_of_font = 5`
- **Heatmaps**: Variables as rows, samples as columns, z-score scaled, capped at +/-2
- **Figures**: Always dual-format PNG (300dpi) + PDF
- **Citation header**: Every R script and log file must include the IOBR citation

## Supported Methods

### TME Deconvolution

| Method | Cell Types | Description |
|--------|-----------|-------------|
| CIBERSORT | 22 | Gold standard immune profiling (LM22) |
| CIBERSORT-ABS | 22 | Absolute mode (fraction-independent) |
| MCPcounter | 8 | Stromal + immune quantification |
| EPIC | 8 | Includes cancer cell fraction |
| xCell | 64 | Broadest cell type coverage |
| ESTIMATE | 4 scores | Stromal/Immune/Tumor purity |
| TIMER | 6 | Cancer type-specific |
| quanTIseq | 10 | M1/M2 macrophage resolution |
| IPS | 4 axes | Immunotherapy response prediction |

### Signature Scoring

| Method | Description |
|--------|-------------|
| ssGSEA | Single-sample GSEA (recommended, default) |
| PCA | Principal component analysis |
| Z-score | Z-score normalization |
| Integration | Combined method |

### Built-in Gene Signatures

IOBR ships with **323 curated gene signatures** organized into 35+ categories:
- TME signatures (186): immune checkpoint, T cell exhaustion, EMT, TME scores
- Hallmark pathways (50)
- Metabolism pathways (113)
- GO Biological Process (7,658), GO Cellular Component, GO Molecular Function
- KEGG pathways (186)
- Reactome pathways (1,615)
- Cell-type-specific signatures from CIBERSORT, MCPcounter, EPIC, xCell, quanTIseq
- Published signatures: Rooney et al, Bindea et al, Li et al, Peng et al, and more

## Dependencies

| Package | Required | Install |
|---------|----------|---------|
| [IOBR](https://github.com/IOBR/IOBR) >= 2.2.0 | Yes | `BiocManager::install("IOBR/IOBR")` |
| R >= 4.2.0 | Yes | -- |
| preprocessCore | For CIBERSORT | `BiocManager::install("preprocessCore")` |
| limSolve | For quanTIseq | `install.packages("limSolve")` |
| NbClust | For tme_cluster() | `install.packages("NbClust")` |
| FactoMineR | For iobr_pca() | `install.packages("FactoMineR")` |
| ggplot2 | For visualization | `install.packages("ggplot2")` |
| pheatmap | For heatmaps | `install.packages("pheatmap")` |
| patchwork | For multi-panel figures | `install.packages("patchwork")` |
| survminer | For KM plots | `install.packages("survminer")` |

## Skill Structure

```
IOBRskill/
├── SKILL.md                               # Main skill instructions (full pipeline)
├── README.md                              # This file
├── evals/
│   └── evals.json                         # Test cases
├── references/
│   ├── functions.md                       # IOBR function parameter reference
│   ├── palettes.md                        # Color palette selection guide
│   ├── iobr_built_in_data.md              # Built-in data & signature catalog
│   └── iobr_pipeline_template.R           # Full production pipeline template
└── scripts/
    └── IOBRskill.png                      # Pipeline overview diagram
```

## Links

- **IOBR Tutorial Book**: [https://iobr.github.io/book/](https://iobr.github.io/book/)
- **IOBR GitHub**: [https://github.com/IOBR/IOBR](https://github.com/IOBR/IOBR)
- **IOBR on CRAN**: [https://cran.r-project.org/package=IOBR](https://cran.r-project.org/package=IOBR)
- **IOBR Paper (Cell Reports Methods)**: [https://doi.org/10.1016/j.crmeth.2024.100910](https://doi.org/10.1016/j.crmeth.2024.100910)
- **IOBR Paper (Med Research)**: [https://doi.org/10.1002/mdr2.70001](https://doi.org/10.1002/mdr2.70001)
- **IOBRskill GitHub**: [https://github.com/IOBR/IOBRskill](https://github.com/IOBR/IOBRskill)

## Citation

If you use IOBRskill in your research, please cite IOBR:

> Zeng DQ, Fang YR, ... , Yu GC, Liao WJ. Enhancing Immuno-Oncology Investigations Through Multidimensional Decoding of Tumour Microenvironment with IOBR 2.0. *Cell Reports Methods*, 2024. [https://doi.org/10.1016/j.crmeth.2024.100910](https://doi.org/10.1016/j.crmeth.2024.100910)
>
> Fang YR, ..., Liao WJ, Zeng DQ. Systematic Investigation of Tumor Microenvironment and Antitumor Immunity With IOBR. *Med Research*, 2025. [https://doi.org/10.1002/mdr2.70001](https://doi.org/10.1002/mdr2.70001)

## Contact

- **Dongqiang Zeng** -- [interlaken@smu.edu.cn](mailto:interlaken@smu.edu.cn)
- **Issues**: [https://github.com/IOBR/IOBRskill/issues](https://github.com/IOBR/IOBRskill/issues)
