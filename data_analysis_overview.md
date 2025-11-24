---
output:
  word_document: default
  html_document: default
---
# Methods: Protein-Protein Interaction Network Analysis of Bladder Cancer Signaling Pathways

## Overview

This analysis integrates protein-protein interaction (PPI) network topology, bulk RNA-seq gene expression, single-cell RNA-seq classification, and functional enrichment to characterize genes in three receptor tyrosine kinase signaling pathways (EGFR, FGFR3, and ERBB2) in bladder cancer.

---

## 1. Pathway Gene Retrieval and Network Expansion

### 1.1 Seed Gene Selection
**Source**: Reactome Pathway Database (https://reactome.org/)

- Downloaded Ensembl2Reactome mapping file containing human pathway annotations
- Filtered for pathway-specific genes:
  - **EGFR pathway**: Genes annotated to EGFR-related pathways
  - **FGFR3 pathway**: Genes annotated to FGFR3-related pathways
  - **ERBB2 pathway**: Genes annotated to ERBB2-related pathways
- Converted Ensembl IDs to HGNC gene symbols using biomaRt (Ensembl v11.5)

### 1.2 Protein-Protein Interaction Network Construction
**Source**: STRING Database v11.5 (https://string-db.org/)

**Parameters**:
- Species: Homo sapiens (taxonomy ID: 9606)
- Confidence score threshold: 500 (medium-high confidence)
- Network expansion: Up to 100 neighbor genes per pathway

**Algorithm** (`expand_signature_network.R`):
1. Map seed genes to STRING protein IDs
2. Retrieve all interacting neighbors from STRING database
3. Rank neighbors by degree (number of connections to seed genes)
4. Select top N neighbors (max_added_genes = 100) with highest connectivity
5. Build expanded network including seed genes + top neighbors
6. Filter to genes present in expression dataset
7. Extract network as igraph object for downstream analysis

### 1.3 Network Topology Metrics
**Tool**: igraph R package

Calculated centrality metrics for each gene node:
- **Degree**: Number of direct protein-protein interactions
- **Betweenness**: Number of shortest paths passing through the node (bridge importance)
- **Closeness**: Average distance to all other nodes (centrality)
- **Eigenvector centrality**: Influence based on connections to other highly-connected nodes
- **Hub score**: HITS algorithm hub score (regulatory potential)
- **Authority score**: HITS algorithm authority score (downstream effector potential)
- **Community**: Module assignment using Louvain clustering algorithm

---

## 2. Gene Expression Analysis

### 2.1 Bulk RNA-Seq Dataset
**Source**: Uroscanseq bladder cancer cohort, 482 samples (UC high quality) subset to the following subtypes; UroA, UroB, UroC, Gu and BaSq (omit Mes and ScNE).

**Samples**: Bladder cancer tumor samples classified into molecular subtypes:
- **UroA**: Urothelial-A subtype
- **UroB**: Urothelial-B subtype
- **UroC**: Urothelial-C subtype
- **GU**: Genomically Unstable subtype
- **BaSq**: Basal-Squamous subtype

**Expression Processing**:
- Gene expression data normalized and Z-score transformed
- Calculated mean expression per gene for each subtype
- Calculated global mean across all samples

**Expression Threshold**:
- Genes classified as "well-expressed" if global mean expression exceeds dataset-specific threshold (< -0.081, 25th percentile) 
- Used for filtering low/unreliable expression genes

### 2.2 Differential Expression Analysis
**Method**: One-way ANOVA with effect size calculations (`run_anova.R`)

**Statistical Framework**:
- **Null hypothesis**: Gene expression does not differ across bladder cancer subtypes
- **Test**: One-way ANOVA comparing expression across 5 subtypes
- **Degrees of freedom**: 
  - Between groups: k-1 = 4 (5 subtypes)
  - Within groups: N-k (samples - subtypes)

**Effect Size Metrics**:
1. **Eta-squared (η²)**: Proportion of total variance explained by subtype
   - Formula: SS_between / SS_total
   - Interpretation: 0.01 = small, 0.06 = medium, 0.14 = large

2. **Partial eta-squared (η²p)**: Variance explained excluding other factors
   - Formula: SS_between / (SS_between + SS_within)

3. **Omega-squared (ω²)**: Unbiased estimate of effect size
   - Formula: (SS_between - df_between × MS_within) / (SS_total + MS_within)
   - Less biased than η² for small sample sizes

4. **Cohen's f**: Standardized measure of group separation
   - Formula: √(η² / (1 - η²))
   - Interpretation: <0.1 = small, 0.25 = medium, 0.4 = large

5. **Expression Range**: Maximum - minimum subtype mean expression

6. **Effect Size Ratio**: Expression range / overall standard deviation
   - Represents signal-to-noise ratio

**Multiple Testing Correction**:
- **Benjamini-Hochberg (FDR)**: Controls false discovery rate (α = 0.05)
- **Bonferroni**: Conservative family-wise error rate control

**Significance Categories**:
- `significant_05`: Adjusted p-value < 0.05
- `significant_01`: Adjusted p-value < 0.01
- `sig_and_meaningful`: Significant + η² ≥ 0.01 (at least small effect)
- `sig_and_large_effect`: Significant + η² ≥ 0.06 (medium/large effect)

### 2.3 Subtype Expression Ranking
For each gene, ranked the 5 subtypes by mean expression level:
- Rank 1 = highest expression subtype
- Rank 5 = lowest expression subtype
- Identifies subtype-specific activation patterns

---

## 3. Single-Cell RNA-Seq Analysis

### 3.1 Datasets
**Normal Bladder**: Single-cell RNA-seq from healthy bladder tissue
**Bladder Cancer**: Single-cell RNA-seq from bladder cancer tumors

**Processing**: Seurat R package
- Data normalized and scaled
- Unsupervised clustering to identify cell populations
- Urothelial cell cluster identified using canonical markers:
  - GATA3, FOXA1, KRT7, KRT8, PPARG, ELF3, FGFR3, CCND1, TP63, KRT5, KRT14, EPCAM, CDH1

### 3.2 Gene Classification in Single Cells
**Method**: Binary classification based on expression thresholds (`run_single_cell_classification.R`)

**Thresholds**:
- **Detection rate**: ≥5% of urothelial cells must express the gene
- **Mean expression**: ≥0.1 in urothelial cells

**Classification**:
- **"Urothelial"**: Meets both thresholds (expressed in bladder epithelium)
- **"Non-urothelial"**: Below thresholds (not typically expressed in urothelium)

**Metrics Calculated**:
- Mean expression in urothelial cells
- Detection rate (proportion of cells with non-zero expression)
- Maximum expression (peak expression in any cell)
- Number of cells expressing
- Total urothelial cells in dataset

**Analysis performed separately for**:
- Normal bladder tissue (`sc_norm_*` columns)
- Cancer tissue (`sc_cancer_*` columns)

---

## 4. Functional Enrichment Analysis

### 4.1 Gene Ontology Enrichment
**Tool**: clusterProfiler R package with org.Hs.eg.db annotation database

**Method**:
1. Convert gene symbols to Entrez IDs using biomaRt
2. Perform GO enrichment analysis for Biological Process (BP) ontology
3. Multiple testing correction: Benjamini-Hochberg (FDR)
4. Significance thresholds:
   - p-value < 0.05
   - q-value < 0.2

### 4.2 Functional Annotation Metrics

**Per-gene metrics**:
- **min_p_adjust**: Minimum adjusted p-value across all GO terms for this gene
  - Lower = stronger functional annotation evidence

- **functional_diversity**: Number of distinct GO Biological Process terms
  - Measures pleiotropy (multi-functionality)

- **functional_significance**: -log10(min_p_adjust)
  - Higher = stronger confidence in functional role
  - Values >20 indicate very strong annotation

- **functional_specificity_ratio**: Gene's functional_significance / median
  - >1 indicates above-average functional characterization

### 4.3 Functional Categorization

**Functional Categories** (based on GO term pattern matching):
- **Translation**: Ribosomal and translation-related processes
- **MAPK_Signaling**: ERK/MAPK cascade and signaling
- **Growth_Factor_Response**: ERBB, peptide hormone, and growth factor signaling
- **Phosphorylation**: Tyrosine phosphorylation and kinase activity
- **Transferase_Regulation**: Transferase and kinase activity regulation
- **Other**: All other or unannotated functions

### 4.4 Pathway Role Assignment

**Integrated network-function classification**:

Rules based on network topology + functional category:

1. **Functional_Hub**: High degree (≥80th percentile) + functional annotation
2. **Functional_Bridge**: High betweenness (≥80th percentile) + functional annotation
3. **Core_Process**: Translation or MAPK_Signaling categories
4. **Signal_Receptor**: Growth_Factor_Response category
5. **Signal_Modifier**: Phosphorylation category
6. **Peripheral**: All other genes

### 4.5 Functional-Structural Importance Score

**Composite metric** integrating multiple dimensions:

Formula: Z(functional_significance) + Z(degree) + Z(betweenness)

- Combines standardized Z-scores of:
  - Functional annotation strength
  - Network connectivity (degree)
  - Network bridging capacity (betweenness)
  
- Higher scores = genes important in both structure AND function
- Identifies key driver genes in the pathway

---

## 5. Data Integration and Consolidation

### 5.1 Merging Data Layers
**Script**: `consoldiate_node_meta.R`

**Integration steps**:

1. **Load individual pathway network metrics**:
   - `egfr_node_metrics.Rdata`
   - `fgfr3_node_metrics.Rdata`
   - `erbb2_node_metrics.Rdata`

2. **Add pathway identifiers**:
   ```r
   egfr_node_metrics$pathway = "EGFR"
   fgfr3_node_metrics$pathway = "FGFR3"
   erbb2_node_metrics$pathway = "ERBB2"
   ```

3. **Combine network metrics** from all three pathways:
   ```r
   node_metrics_combined = rbind(egfr_node_metrics, fgfr3_node_metrics, erbb2_node_metrics)
   ```

4. **Add expression threshold classification** (`is_well_expressed`):
   - Load filtered pathway genes that passed expression threshold
   - Mark genes as TRUE/FALSE based on presence in filtered gene list
   - TRUE = global mean expression > threshold (25th percentile = -0.081)

5. **Add global mean expression**:
   - Calculate mean expression across all samples (all subtypes)
   - Join to node metrics by gene name

6. **Join ANOVA statistical results and effect sizes**:
   - Left join `anova_with_effects` by gene name
   - Adds: F-statistic, p-value, degrees of freedom, effect sizes (η², ω², Cohen's f), significance flags

7. **Calculate expression rankings across subtypes**:
   - Rank genes within each subtype by mean expression (1 = highest)
   - Creates: `rank_UroA`, `rank_UroB`, `rank_UroC`, `rank_GU`, `rank_BaSq`
   - Uses descending rank (highest expression = rank 1)

8. **Split by pathway for subsequent annotations**:
   ```r
   egfr_metrics = node_metrics_combined %>% filter(pathway == "EGFR")
   fgfr3_metrics = node_metrics_combined %>% filter(pathway == "FGFR3")
   erbb2_metrics = node_metrics_combined %>% filter(pathway == "ERBB2")
   ```

9. **Merge single-cell classifications** (normal and cancer):
   - Left join single-cell data for each pathway separately
   - Adds: `sc_norm_*` columns (normal bladder single-cell data)
   - Adds: `sc_cancer_*` columns (cancer bladder single-cell data)

10. **Combine pathways after single-cell merge**:
    ```r
    node_metrics_combined = rbind(egfr_metrics, fgfr3_metrics, erbb2_metrics)
    ```

11. **Join functional enrichment annotations**:
    - Left join functional annotations for each pathway
    - Adds: GO enrichment metrics, functional categories, pathway roles

12. **Merge Human Protein Atlas (HPA) IHC data**:
    - Left join HPA annotations for each pathway (see Section 7)
    - Adds: Normal tissue IHC data (`hpa_normal_*` columns)
    - Adds: Cancer tissue IHC data (`hpa_cancer_*` columns)
    - Adds: Derived classifications (expression patterns, quality flags, filter flags)
    - **Note**: NO filtering applied; all genes with HPA data retained

13. **Final consolidation**:
    - Combine all three pathways into final integrated dataset
    - Rename for clarity:
      - `egfr_meta_consolidated`
      - `fgfr3_meta_consolidated`
      - `erbb2_meta_consolidated`
      - `node_meta_consolidated` (all pathways combined)

### 5.2 Final Consolidated Dataset

**Output files**:
- `results_out/egfr_meta_consolidated.Rdata`: EGFR pathway with all annotations
- `results_out/fgfr3_meta_consolidated.Rdata`: FGFR3 pathway with all annotations
- `results_out/erbb2_meta_consolidated.Rdata`: ERBB2 pathway with all annotations
- `results_out/node_meta_consolidated.Rdata`: All three pathways combined

**Data structure**: Each gene includes:
- **Network topology** (9 columns): degree, betweenness, closeness, eigenvector, hub score, authority score, community, is_seed, pathway
- **Expression profiles** (11 columns): mean expression per subtype, global mean, expression threshold flag
- **Statistical testing** (18 columns): ANOVA results, effect sizes, multiple testing corrections, significance flags
- **Expression rankings** (5 columns): rank within each subtype (UroA, UroB, UroC, GU, BaSq)
- **Single-cell data - normal** (6 columns): mean expression, detection rate, max expression, cells expressing, total cells, classification
- **Single-cell data - cancer** (6 columns): mean expression, detection rate, max expression, cells expressing, total cells, classification
- **Functional annotations** (6 columns): GO enrichment metrics, functional category, pathway role, functional-structural importance
- **HPA normal tissue** (4 columns): tissue, cell type, expression level, antibody reliability
- **HPA cancer tissue** (5 columns): cancer type, high/medium/low/not detected patient counts
- **HPA derived classifications** (6 columns): normal/cancer expressed, reliability, confidence, expression pattern, data quality
- **HPA filter flags** (4 columns): has data, normal quality, cancer quality, passes filter
- **HPA biological flags** (3 columns): high confidence, therapeutic target, potential suppressor

**Total columns**: 76+ (may vary by pathway)

---

## 6. Software and Packages

### R Version
R version 4.x.x

### Key R Packages
- **Network analysis**: igraph, STRINGdb
- **Expression analysis**: dplyr, tidyr
- **Statistical testing**: base R stats
- **Single-cell analysis**: Seurat
- **Functional enrichment**: clusterProfiler, org.Hs.eg.db
- **Gene annotation**: biomaRt

### External Databases
- **STRING v11.5**: Protein-protein interaction networks
- **Reactome**: Pathway annotations
- **Gene Ontology**: Functional annotations via clusterProfiler
- **Ensembl**: Gene ID mapping via biomaRt
- **Human Protein Atlas**: Protein expression validation (IHC data)

---

## 7. Protein Expression Validation Using Human Protein Atlas

### 7.1 Data Acquisition
**Source**: Human Protein Atlas (https://www.proteinatlas.org/)
**Script**: `get_hpa.R`
**Functions**: `download_hpa_data()` and `add_hpa_annotations()`

Downloaded immunohistochemistry (IHC) data to validate gene expression at the protein level:

**Normal Tissue Data**:
- URL: `https://www.proteinatlas.org/download/tsv/normal_ihc_data.tsv.zip`
- Function: `download_hpa_data(normal_ihc_url, "Normal Tissue IHC Data")`
- Processing:
  - Downloads and extracts TSV file to temporary directory
  - Reads data using `read_tsv()`
  - Filters for: `Tissue == "Urinary bladder"`
- Columns extracted and renamed:
  - `Gene name` → `name`: Gene symbol
  - `Tissue` → `hpa_normal_tissue`: Tissue type (Urinary bladder)
  - `Cell type` → `hpa_normal_cell_type`: Cell type where protein detected (e.g., "urothelial cells")
  - `Level` → `hpa_normal_level`: Expression level (High, Medium, Low, Not detected)
  - `Reliability` → `hpa_normal_reliability`: Antibody validation status (Enhanced, Supported, Approved, Uncertain)
- Removed columns: `Gene` (Ensembl ID), `IHC tissue name` (redundant)
- Saved as: `results_data/bladder_normal_ihc.Rdata`

**Cancer Tissue Data**:
- URL: `https://www.proteinatlas.org/download/tsv/cancer_data.tsv.zip`
- Function: `download_hpa_data(cancer_url, "Cancer Data")`
- Processing:
  - Downloads and extracts TSV file to temporary directory
  - Reads data using `read_tsv()`
  - Filters for: `Cancer == "urothelial cancer"`
- Columns extracted and renamed:
  - `Gene name` → `name`: Gene symbol
  - `Cancer` → `hpa_cancer_type`: Cancer type (urothelial cancer)
  - `High` → `hpa_cancer_high`: Number of patients with high protein expression
  - `Medium` → `hpa_cancer_medium`: Number of patients with medium expression
  - `Low` → `hpa_cancer_low`: Number of patients with low expression
  - `Not detected` → `hpa_cancer_not_detected`: Number of patients with no detectable protein
- Removed columns: `Gene` (Ensembl ID)
- Saved as: `results_data/bladder_cancer_ihc.Rdata`

**Note**: Protein Atlas typically includes ~12 patient samples per cancer type

### 7.2 Data Integration and Annotation
**Function**: `add_hpa_annotations(pathway_metrics, bladder_normal_ihc, bladder_cancer_ihc)`

Applied HPA annotations to each pathway's network metrics separately:
- `egfr_hpa <- add_hpa_annotations(egfr_metrics, bladder_normal_ihc, bladder_cancer_ihc)`
- `fgfr3_hpa <- add_hpa_annotations(fgfr3_metrics, bladder_normal_ihc, bladder_cancer_ihc)`
- `erbb2_hpa <- add_hpa_annotations(erbb2_metrics, bladder_normal_ihc, bladder_cancer_ihc)`

**Process**:
1. Left joins normal IHC data by gene name
2. Left joins cancer IHC data by gene name
3. Creates 15 derived classification columns (see below)
4. Returns annotated data frame with all HPA columns

Saved outputs:
- `results_data/egfr_hpa.Rdata`
- `results_data/fgfr3_hpa.Rdata`
- `results_data/erbb2_hpa.Rdata`

### 7.3 Derived Classifications and Metrics

The `add_hpa_annotations()` function creates the following classifications:

**Normal Expression Classification**:
- **hpa_normal_expressed** (Boolean):
  - TRUE: Level = "High" or "Medium"
  - FALSE: Level = "Low" or "Not detected"
  - NA: No data available

**Normal Reliability Classification**:
- **hpa_normal_reliable** (Categorical):
  - "High": Antibody reliability = Enhanced or Supported (most trustworthy)
  - "Medium": Antibody reliability = Approved
  - "Low": Antibody reliability = Uncertain
  - NA: No data available

**Cancer Expression Classification**:
- **hpa_cancer_expressed** (Boolean):
  - TRUE: (hpa_cancer_high + hpa_cancer_medium) > 0 patients
  - FALSE: Only Low expression OR ≥10 patients with Not detected
  - NA: No data available

**Cancer Expression Confidence**:
- **hpa_cancer_confidence** (Categorical):
  - "High": (High + Medium) ≥ 6 patients (>50% of cohort)
  - "Medium": (High + Medium) 3-5 patients (25-50%)
  - "Low": (High + Medium) 1-2 patients (<25%)
  - "None": No High or Medium expression

### 7.4 Biological Expression Patterns

Combined normal and cancer classifications to identify biological patterns:

**hpa_expression_pattern** (Categorical):
- **"Constitutive"**: Expressed in normal urothelium AND cancer
  - Interpretation: Core urothelial genes maintained in malignancy
  
- **"Lost_in_cancer"**: Expressed in normal but absent/reduced in cancer
  - Interpretation: Potential tumor suppressors; silenced during transformation
  
- **"Cancer_upregulated"**: Not expressed in normal but detected in cancer
  - Interpretation: Cancer-specific activation; potential therapeutic targets
  
- **"Not_expressed"**: Absent in both normal and cancer
  - Interpretation: Not relevant to urothelial biology
  
- **"Insufficient_data"**: Missing normal or cancer data
  - Interpretation: Cannot classify expression pattern

- **"Unknown"**: Edge cases not fitting other categories

### 7.5 Data Quality Assessment

**hpa_data_quality** assigned based on completeness and reliability:

- **"High"**: 
  - Has both normal and cancer IHC data
  - Antibody reliability = Enhanced, Supported, or Approved
  
- **"Medium"**: 
  - Has both normal and cancer data
  - Any antibody reliability level (including Uncertain)
  
- **"Low"**: 
  - Missing either normal or cancer dataset

- **"None"**:
  - No data available for either normal or cancer

### 7.6 Quality Control Filter Flags

The `add_hpa_annotations()` function creates four filter flags for quality control:

**hpa_filter_has_data** (Boolean):
- TRUE: `hpa_data_quality` is "High" or "Medium"
- FALSE: Missing normal or cancer data
- Purpose: Ensures both datasets are present

**hpa_filter_normal_quality** (Boolean):
- TRUE: `hpa_normal_expressed` = TRUE AND `hpa_normal_reliable` = "High" or "Medium"
- FALSE: Otherwise
- Purpose: Identifies genes reliably detected in normal bladder with validated antibodies

**hpa_filter_cancer_quality** (Boolean):
- TRUE: `hpa_cancer_confidence` = "High" or "Medium" (≥3 patients)
- FALSE: Otherwise
- Purpose: Identifies genes with reliable cancer expression (≥25% of patients)

**hpa_passes_filter** (Boolean):
- TRUE: `hpa_filter_has_data` AND (`hpa_filter_normal_quality` OR `hpa_filter_cancer_quality`)
- FALSE: Otherwise
- Purpose: Combined quality criteria - has complete data AND reliable expression in either normal or cancer
- **Note**: This flag is created but NOT used for filtering in `get_hpa.R`; all annotations are retained

### 7.7 Biological Interpretation Flags

Three additional flags for identifying genes of biological interest:

**hpa_high_confidence** (Boolean):
- TRUE: `hpa_data_quality` = "High" AND `hpa_cancer_confidence` = "High" AND `hpa_normal_reliable` = "High"
- FALSE: Otherwise
- Purpose: Gold standard protein validation with maximum confidence

**hpa_therapeutic_target** (Boolean):
- TRUE: `hpa_expression_pattern` = "Cancer_upregulated" AND `hpa_cancer_confidence` = "High" or "Medium"
- FALSE: Otherwise
- Purpose: Flags proteins upregulated in cancer with good patient evidence; candidates for targeted therapy

**hpa_potential_suppressor** (Boolean):
- TRUE: `hpa_expression_pattern` = "Lost_in_cancer" AND `hpa_normal_expressed` = TRUE
- FALSE: Otherwise
- Purpose: Flags proteins lost during cancer transformation; candidate tumor suppressors

### 7.8 Integration into Consolidated Metrics

**Script**: `consoldiate_node_meta.R`

HPA annotations added to consolidated pathway metrics:
1. Load HPA-annotated data for each pathway
2. Left join with pathway metrics by gene name:
   ```r
   egfr_metrics <- egfr_metrics %>% left_join(egfr_hpa, by = "name")
   fgfr3_metrics <- fgfr3_metrics %>% left_join(fgfr3_hpa, by = "name")
   erbb2_metrics <- erbb2_metrics %>% left_join(erbb2_hpa, by = "name")
   ```
3. Combine all pathways into final consolidated dataset

**Output files**:
- `results_out/egfr_meta_consolidated.Rdata`: EGFR pathway with all annotations including HPA
- `results_out/fgfr3_meta_consolidated.Rdata`: FGFR3 pathway with all annotations including HPA
- `results_out/erbb2_meta_consolidated.Rdata`: ERBB2 pathway with all annotations including HPA
- `results_out/node_meta_consolidated.Rdata`: All three pathways combined

**Important Note**: The `get_hpa.R` script creates filter flags but does NOT filter the data. All genes with any level of HPA data are retained in the output files. Filtering decisions can be made downstream based on the `hpa_passes_filter` flag or custom criteria.

---

## 8. Quality Control and Filtering

### Network Construction
- Only included genes present in expression dataset
- Required minimum STRING confidence score of 500
- Removed isolated nodes (degree < 1)

### Expression Analysis
- Filtered low-expressed genes based on global mean threshold
- Required successful ANOVA convergence for statistical metrics
- Applied multiple testing correction for all p-values

### Single-Cell Analysis
- Identified urothelial cells using 13 canonical marker genes
- Required minimum detection rate (5%) and expression (0.1) for classification
- Analyzed only urothelial cells (excluded other cell types)

### Functional Enrichment
- Required FDR-adjusted p-value < 0.05 for GO term inclusion
- Handled unannotated genes with default values (0 or "Unannotated")

---

## 8. Interpretation Guidelines

### Network Importance
- **High degree + high betweenness** = Critical hub and bridge gene
- **High eigenvector centrality** = Connected to other important genes
- **Community assignment** = Functional module membership

### Expression Patterns
- **High eta-squared** = Subtype-specific expression (biologically important)
- **High Cohen's f** = Strong separation between subtypes
- **Rank = 1 in subtype** = Highest expression in that subtype

### Single-Cell Context
- **Urothelial in normal, Non-urothelial in cancer** = Downregulated in cancer
- **Non-urothelial in normal, Urothelial in cancer** = Upregulated/activated in cancer
- **High detection rate** = Widely expressed across cell population

### Functional Role
- **High functional_structural_importance** = Key driver gene
- **Functional_Hub or Functional_Bridge** = Central to pathway function
- **High functional_diversity** = Multi-functional/pleiotropic gene

### Protein Atlas Expression Patterns
- **Constitutive (HPA)** = Core urothelial genes; maintained from normal to cancer
- **Cancer_upregulated (HPA)** = Cancer-specific activation; potential therapeutic targets
- **Lost_in_cancer (HPA)** = Tumor suppressor candidates; silenced in malignancy
- **hpa_data_quality = High** = Both normal and cancer IHC data with reliable antibodies
- **hpa_cancer_confidence = High/Medium** = Detected in ≥25% of patient samples (≥3/12)

---

## References

- STRING Database: https://string-db.org/
- Reactome Pathway Database: https://reactome.org/
- Human Protein Atlas: https://www.proteinatlas.org/
  - Uhlen et al., Science 2015
- clusterProfiler: Yu et al., OMICS 2012
- Seurat: Hao et al., Cell 2021
- biomaRt: Durinck et al., Bioinformatics 2009
