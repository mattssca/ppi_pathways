---
output:
  word_document: default
  html_document: default
---
# Column Documentation for Pathway Metrics Data

This data frame integrates **network topology metrics**, **gene expression data**, **statistical analysis**, **single-cell RNA-seq classifications**, **functional enrichment annotations**, and **protein-level validation** for genes in three receptor tyrosine kinase signaling pathways (EGFR, FGFR3, and ERBB2) in bladder cancer.

**Total columns**: 76+

**Data layers**:
- Network topology (9 columns)
- Expression profiles (11 columns)
- Statistical testing & effect sizes (18 columns)
- Expression rankings (5 columns)
- Single-cell RNA-seq - Normal (6 columns)
- Single-cell RNA-seq - Cancer (6 columns)
- Functional enrichment (6 columns)
- Human Protein Atlas - Normal tissue (4 columns)
- Human Protein Atlas - Cancer tissue (5 columns)
- HPA derived classifications (6 columns)
- HPA filter flags (4 columns)
- HPA biological interpretation flags (3 columns)

---

## **1. Basic Gene & Network Information**

- **`name`** 
  - **What**: Gene symbol (HGNC name)
  - **Source**: Original seed genes from Reactome pathway database + expanded network neighbors from STRING database
  - **Interpretation**: Unique identifier for each gene in the protein-protein interaction network

- **`degree`** 
  - **What**: Number of direct connections (edges) a gene has in the network
  - **Source**: Calculated from STRING protein-protein interaction network using igraph
  - **Interpretation**: Higher values indicate hub genes with many interaction partners; important for network connectivity

- **`is_seed`** 
  - **What**: Boolean indicating if gene was an original pathway member
  - **Source**: Reactome pathway database (ERBB2 pathway genes)
  - **Interpretation**: TRUE = original pathway gene; FALSE = neighbor gene added during network expansion

---

## **2. Network Centrality Metrics**

- **`betweenness`** 
  - **What**: Number of shortest paths between other nodes that pass through this gene
  - **Source**: Calculated from STRING PPI network using igraph
  - **Interpretation**: Higher values indicate "bridge" genes connecting different network modules; important for information flow

- **`closeness`** 
  - **What**: Average shortest path distance from this gene to all other genes
  - **Source**: Calculated from STRING PPI network using igraph
  - **Interpretation**: Higher values indicate genes centrally positioned in the network; can quickly influence other genes

- **`eigenvector`** 
  - **What**: Influence score based on connections to other well-connected genes
  - **Source**: Calculated from STRING PPI network using igraph
  - **Interpretation**: High values indicate genes connected to other important hub genes; reflects network prestige

- **`hub_score`** 
  - **What**: HITS algorithm hub score (authority-hub analysis)
  - **Source**: Calculated from STRING PPI network using igraph
  - **Interpretation**: High scores indicate genes that point to many authoritative genes; regulatory potential

- **`authority_score`** 
  - **What**: HITS algorithm authority score
  - **Source**: Calculated from STRING PPI network using igraph
  - **Interpretation**: High scores indicate genes pointed to by many hubs; downstream effectors

- **`community`** 
  - **What**: Network module/cluster assignment
  - **Source**: Community detection algorithm (likely Louvain) on STRING PPI network
  - **Interpretation**: Genes with same number belong to same functional module; groups genes with similar roles

---

## **3. Bulk RNA-Seq Expression Data (UroScan cohort)**

- **`mean_expr_UroB`** 
  - **What**: Mean Z-scored expression in Urothelial-B subtype samples
  - **Source**: UroScan bladder cancer bulk RNA-seq dataset
  - **Interpretation**: Positive = higher than cohort average; negative = lower; reflects subtype-specific activity

- **`mean_expr_GU`** 
  - **What**: Mean Z-scored expression in Genomically Unstable subtype
  - **Source**: UroScan bladder cancer bulk RNA-seq dataset
  - **Interpretation**: Same as above for GU subtype

- **`mean_expr_UroA`** 
  - **What**: Mean Z-scored expression in Urothelial-A subtype
  - **Source**: UroScan bladder cancer bulk RNA-seq dataset
  - **Interpretation**: Same as above for UroA subtype

- **`mean_expr_BaSq`** 
  - **What**: Mean Z-scored expression in Basal-Squamous subtype
  - **Source**: UroScan bladder cancer bulk RNA-seq dataset
  - **Interpretation**: Same as above for BaSq subtype

- **`mean_expr_UroC`** 
  - **What**: Mean Z-scored expression in Urothelial-C subtype
  - **Source**: UroScan bladder cancer bulk RNA-seq dataset
  - **Interpretation**: Same as above for UroC subtype

---

## **4. Pathway & Expression Threshold**

- **`pathway`** 
  - **What**: Which signaling pathway this gene belongs to
  - **Source**: Manually assigned during consolidation
  - **Interpretation**: "ERBB2", "EGFR", or "FGFR3" - indicates which pathway network this gene came from

- **`is_well_expressed`** 
  - **What**: Boolean indicating if gene passes expression threshold
  - **Source**: Filtering based on mean expression across all samples in UroScan dataset
  - **Interpretation**: TRUE = gene has sufficient expression for reliable analysis; FALSE = low/unreliable expression

- **`global_mean`** 
  - **What**: Mean expression across ALL samples (all subtypes combined)
  - **Source**: Calculated from UroScan bulk RNA-seq data
  - **Interpretation**: Overall expression level; positive = above average; used for filtering and normalization

---

## **5. ANOVA Statistical Results**

- **`f_statistic`** 
  - **What**: F-statistic from one-way ANOVA testing expression differences across 5 subtypes
  - **Source**: ANOVA calculated on UroScan expression data
  - **Interpretation**: Higher values indicate larger differences between subtypes; measures overall variance ratio

- **`p_value`** 
  - **What**: Raw p-value from ANOVA test
  - **Source**: ANOVA calculated on UroScan expression data
  - **Interpretation**: Probability of observing these expression differences by chance; <0.05 suggests real differences

- **`df_between`** 
  - **What**: Degrees of freedom between groups (subtypes - 1 = 4)
  - **Source**: ANOVA calculation
  - **Interpretation**: Number of subtypes minus 1; used in F-test calculation

- **`df_within`** 
  - **What**: Degrees of freedom within groups (total samples - subtypes)
  - **Source**: ANOVA calculation  
  - **Interpretation**: Sample size minus number of groups; reflects statistical power

---

## **6. Effect Size Metrics**

- **`eta_squared`** 
  - **What**: Proportion of total variance explained by subtype
  - **Source**: Calculated from ANOVA results
  - **Interpretation**: 0-1 scale; 0.01=small, 0.06=medium, 0.14=large effect; indicates biological importance

- **`partial_eta_squared`** 
  - **What**: Proportion of variance explained excluding other factors
  - **Source**: Calculated from ANOVA results
  - **Interpretation**: Similar to eta_squared but adjusts for other variables (same as eta_squared in one-way ANOVA)

- **`omega_squared`** 
  - **What**: Unbiased estimate of effect size (less biased than eta_squared)
  - **Source**: Calculated from ANOVA results
  - **Interpretation**: More conservative than eta_squared; better for small samples; same interpretation scale

- **`cohens_f`** 
  - **What**: Cohen's f effect size metric
  - **Source**: Calculated from eta_squared (sqrt(eta_squared / (1 - eta_squared)))
  - **Interpretation**: <0.1=small, 0.25=medium, 0.4=large effect; standardized measure of group separation

- **`expression_range`** 
  - **What**: Difference between highest and lowest subtype mean expression
  - **Source**: Calculated from subtype-specific means
  - **Interpretation**: Larger values = more dramatic expression differences across subtypes

- **`overall_sd`** 
  - **What**: Standard deviation of expression across all samples
  - **Source**: Calculated from UroScan expression data
  - **Interpretation**: Overall expression variability; higher = more variable gene

- **`effect_size_ratio`** 
  - **What**: Ratio of expression range to overall standard deviation
  - **Source**: expression_range / overall_sd
  - **Interpretation**: >1 indicates subtype differences exceed overall variability; signal-to-noise ratio

---

## **7. Multiple Testing Corrections & Significance**

- **`p_adjusted_bh`** 
  - **What**: Benjamini-Hochberg corrected p-value (FDR control)
  - **Source**: p.adjust() on ANOVA p-values
  - **Interpretation**: Controls false discovery rate; <0.05 indicates significant after accounting for multiple testing

- **`p_adjusted_bonf`** 
  - **What**: Bonferroni corrected p-value (FWER control)
  - **Source**: p.adjust() on ANOVA p-values
  - **Interpretation**: Very conservative; <0.05 indicates high confidence in significance

- **`significant_05`** 
  - **What**: Boolean for statistical significance at α=0.05
  - **Source**: p_adjusted_bh < 0.05
  - **Interpretation**: TRUE = statistically significant differential expression across subtypes

- **`significant_01`** 
  - **What**: Boolean for statistical significance at α=0.01
  - **Source**: p_adjusted_bh < 0.01
  - **Interpretation**: TRUE = highly significant differential expression

- **`effect_size_category`** 
  - **What**: Categorical effect size classification
  - **Source**: Based on eta_squared thresholds
  - **Interpretation**: "Negligible"/<1%, "Small"/1-6%, "Medium"/6-14%, "Large"/≥14% variance explained

- **`cohens_f_category`** 
  - **What**: Categorical Cohen's f classification
  - **Source**: Based on cohens_f thresholds
  - **Interpretation**: "Small"/<0.1, "Medium"/0.1-0.25, "Large"/0.25-0.4, "Very Large"/≥0.4

- **`sig_and_meaningful`** 
  - **What**: Boolean for significance + at least small effect
  - **Source**: significant_05 & eta_squared ≥ 0.01
  - **Interpretation**: TRUE = both statistically significant AND biologically meaningful

- **`sig_and_large_effect`** 
  - **What**: Boolean for significance + medium/large effect
  - **Source**: significant_05 & eta_squared ≥ 0.06
  - **Interpretation**: TRUE = significant with substantial biological importance

---

## **8. Expression Ranking Across Subtypes**

- **`rank_UroA`**, **`rank_UroB`**, **`rank_UroC`**, **`rank_GU`**, **`rank_BaSq`** 
  - **What**: Rank of expression level for this gene in each subtype (1=highest)
  - **Source**: Ranked transformation of mean_expr columns
  - **Interpretation**: 1 = this subtype has highest expression; 5 = lowest; identifies subtype preference

---

## **9. Single-Cell RNA-Seq - Normal Bladder Tissue**

- **`sc_norm_mean_expression`** 
  - **What**: Mean expression in urothelial cells from normal bladder
  - **Source**: Single-cell RNA-seq analysis (bladder_normal_seurat_processed)
  - **Interpretation**: Expression level in healthy urothelial cells; baseline activity

- **`sc_norm_detection_rate`** 
  - **What**: Proportion of normal urothelial cells expressing this gene
  - **Source**: Single-cell RNA-seq analysis
  - **Interpretation**: 0-1 scale; cellular prevalence; >0.05 suggests widespread expression

- **`sc_norm_max_expression`** 
  - **What**: Maximum expression value observed in any normal urothelial cell
  - **Source**: Single-cell RNA-seq analysis
  - **Interpretation**: Peak expression capacity; identifies highly active cells

- **`sc_norm_cells_expressing`** 
  - **What**: Absolute number of normal urothelial cells expressing the gene
  - **Source**: Single-cell RNA-seq analysis
  - **Interpretation**: Raw count; depends on total cells analyzed

- **`sc_norm_total_urothelial_cells`** 
  - **What**: Total number of urothelial cells in normal dataset
  - **Source**: Single-cell RNA-seq dataset
  - **Interpretation**: Denominator for calculating detection rate; sample size

- **`sc_norm_classification`** 
  - **What**: Binary classification of gene expression pattern in normal tissue
  - **Source**: Based on detection_rate ≥5% & mean_expression ≥0.1 thresholds
  - **Interpretation**: "Urothelial" = expressed in normal bladder cells; "Non-urothelial" = not typically expressed

---

## **10. Single-Cell RNA-Seq - Bladder Cancer Tissue**

- **`sc_cancer_mean_expression`** 
  - **What**: Mean expression in urothelial cells from bladder cancer
  - **Source**: Single-cell RNA-seq analysis (bladder_cancer_seurat_processed)
  - **Interpretation**: Expression in cancer cells; compare to normal to see changes

- **`sc_cancer_detection_rate`** 
  - **What**: Proportion of cancer urothelial cells expressing this gene
  - **Source**: Single-cell RNA-seq analysis
  - **Interpretation**: Cellular prevalence in tumors; upregulation vs normal indicates activation

- **`sc_cancer_max_expression`** 
  - **What**: Maximum expression in any cancer urothelial cell
  - **Source**: Single-cell RNA-seq analysis
  - **Interpretation**: Peak cancer cell expression; identifies highly expressing subpopulations

- **`sc_cancer_cells_expressing`** 
  - **What**: Absolute number of cancer cells expressing the gene
  - **Source**: Single-cell RNA-seq analysis
  - **Interpretation**: Raw count of expressing cells

- **`sc_cancer_total_urothelial_cells`** 
  - **What**: Total urothelial cells in cancer dataset
  - **Source**: Single-cell RNA-seq dataset
  - **Interpretation**: Sample size for cancer analysis

- **`sc_cancer_classification`** 
  - **What**: Binary classification in cancer tissue
  - **Source**: Based on detection_rate ≥5% & mean_expression ≥0.1
  - **Interpretation**: "Urothelial" = expressed in cancer; "Non-urothelial" = not expressed; compare to normal classification

---

## **11. Functional Enrichment Analysis**

- **`min_p_adjust`** 
  - **What**: Minimum adjusted p-value across all GO terms this gene belongs to
  - **Source**: clusterProfiler GO enrichment analysis (Biological Process)
  - **Interpretation**: Strength of functional annotation; lower = stronger evidence for functional role

- **`functional_diversity`** 
  - **What**: Number of distinct GO Biological Process terms this gene is annotated to
  - **Source**: Count of unique GO terms from enrichment analysis
  - **Interpretation**: Higher values = multi-functional gene; pleiotropy indicator

- **`functional_significance`** 
  - **What**: -log10(min_p_adjust)
  - **Source**: Transformation of GO enrichment p-value
  - **Interpretation**: Higher = stronger functional annotation confidence; >20 is very strong

- **`functional_category`** 
  - **What**: Broad functional classification based on GO terms
  - **Source**: Pattern matching on GO term descriptions
  - **Interpretation**: Categories include: "Translation", "MAPK_Signaling", "Growth_Factor_Response", "Phosphorylation", "Transferase_Regulation", "Other"

- **`functional_specificity_ratio`** 
  - **What**: Gene's functional_significance divided by median across all genes
  - **Source**: Normalized functional_significance
  - **Interpretation**: >1 = above-average functional annotation; indicates well-characterized genes

---

## **12. Network-Function Integration**

- **`pathway_role`** 
  - **What**: Combined role based on network topology + functional category
  - **Source**: Conditional assignment using degree, betweenness, and functional_category
  - **Interpretation**: 
    - "Functional_Hub" = high degree + functional annotation
    - "Functional_Bridge" = high betweenness + functional annotation
    - "Core_Process" = Translation/MAPK signaling genes
    - "Signal_Receptor" = Growth factor response genes
    - "Signal_Modifier" = Phosphorylation genes
    - "Peripheral" = other genes

- **`functional_structural_importance`** 
  - **What**: Composite Z-score combining functional + network metrics
  - **Source**: Sum of Z-scores: functional_significance + degree + betweenness
  - **Interpretation**: Higher = more important in both network structure AND biological function; identifies key driver genes

---

## **13. Human Protein Atlas - Normal Tissue IHC**

- **`hpa_normal_tissue`** 
  - **What**: Tissue type analyzed (should be "Urinary bladder" for all rows)
  - **Source**: Human Protein Atlas normal tissue IHC database (normal_ihc_data.tsv)
  - **Script**: Downloaded and filtered by `get_hpa.R` using `download_hpa_data()` function
  - **Interpretation**: Organ/tissue where protein expression was measured

- **`hpa_normal_cell_type`** 
  - **What**: Specific cell type(s) where protein expression was detected
  - **Source**: Pathologist annotation of IHC staining patterns from Human Protein Atlas
  - **Script**: Renamed from "Cell type" column in `get_hpa.R`
  - **Interpretation**: 
    - "urothelial cells" = Expression in bladder epithelial lining (primary interest)
    - May also include: "smooth muscle cells", "endothelial cells", etc.
    - Multiple cell types listed if protein is expressed in multiple compartments

- **`hpa_normal_level`** 
  - **What**: Intensity/level of protein expression from antibody staining
  - **Source**: Semi-quantitative pathologist scoring of staining intensity (Human Protein Atlas)
  - **Script**: Renamed from "Level" column in `get_hpa.R`
  - **Interpretation**:
    - "High" = Strong staining; abundant protein
    - "Medium" = Moderate staining
    - "Low" = Weak staining; protein present but at low levels
    - "Not detected" = No visible staining; protein absent or below detection threshold

- **`hpa_normal_reliability`** 
  - **What**: Quality/validation status of the antibody used for IHC
  - **Source**: Human Protein Atlas antibody validation pipeline
  - **Script**: Renamed from "Reliability" column in `get_hpa.R`
  - **Interpretation**:
    - "Enhanced" = Highest reliability; validated by protein/gene characterization + independent antibodies
    - "Supported" = High reliability; consistent with RNA expression data
    - "Approved" = Medium reliability; consistent with literature
    - "Uncertain" = Lower reliability; unclear or conflicting data

---

## **14. Human Protein Atlas - Cancer Tissue IHC**

- **`hpa_cancer_type`** 
  - **What**: Type of cancer analyzed (should be "urothelial cancer" for all rows)
  - **Source**: Human Protein Atlas cancer database (cancer_data.tsv)
  - **Script**: Downloaded and filtered by `get_hpa.R` using `download_hpa_data()` function
  - **Interpretation**: Cancer type; urothelial cancer = bladder cancer

- **`hpa_cancer_high`** 
  - **What**: Number of patient samples showing HIGH protein expression
  - **Source**: Count from HPA patient cohort (~12 samples per cancer type)
  - **Script**: Renamed from "High" column in `get_hpa.R`
  - **Interpretation**: Higher numbers = more patients with strong protein expression

- **`hpa_cancer_medium`** 
  - **What**: Number of patient samples showing MEDIUM protein expression
  - **Source**: Count from HPA patient cohort
  - **Script**: Renamed from "Medium" column in `get_hpa.R`
  - **Interpretation**: Moderate expression frequency across patients

- **`hpa_cancer_low`** 
  - **What**: Number of patient samples showing LOW protein expression
  - **Source**: Count from HPA patient cohort
  - **Script**: Renamed from "Low" column in `get_hpa.R`
  - **Interpretation**: Weak but detectable expression

- **`hpa_cancer_not_detected`** 
  - **What**: Number of patient samples with NO detectable protein
  - **Source**: Count from HPA patient cohort
  - **Script**: Renamed from "Not detected" column in `get_hpa.R`
  - **Interpretation**: Protein absent in these patients; compare to normal to identify cancer-specific changes

---

## **15. Human Protein Atlas - Derived Classifications**

All derived classifications created by `add_hpa_annotations()` function in `get_hpa.R`:

- **`hpa_normal_expressed`** 
  - **What**: Boolean classification of normal tissue expression
  - **Source**: Derived from hpa_normal_level using `add_hpa_annotations()` function
  - **Script**: `get_hpa.R` (applied via `add_hpa_annotations()`)
  - **Interpretation**:
    - TRUE = Level is "High" or "Medium" (reliably expressed in normal urothelium)
    - FALSE = Level is "Low" or "Not detected" (not expressed or unreliable)
    - NA = No data available

- **`hpa_normal_reliable`** 
  - **What**: Categorical antibody reliability classification
  - **Source**: Derived from hpa_normal_reliability using `add_hpa_annotations()` function
  - **Script**: `get_hpa.R` (applied via `add_hpa_annotations()`)
  - **Interpretation**:
    - "High" = Enhanced or Supported antibodies (most trustworthy)
    - "Medium" = Approved antibodies
    - "Low" = Uncertain antibodies (interpret cautiously)

- **`hpa_cancer_expressed`** 
  - **What**: Boolean classification of cancer tissue expression
  - **Source**: Derived from hpa_cancer_high + hpa_cancer_medium counts
  - **Script**: `get_hpa.R` (applied via `add_hpa_annotations()`)
  - **Interpretation**:
    - TRUE = At least 1 patient with High or Medium expression
    - FALSE = Only Low expression OR ≥10 patients with Not detected
    - NA = No data available

- **`hpa_cancer_confidence`** 
  - **What**: Confidence level of cancer expression based on patient frequency
  - **Source**: Calculated from (hpa_cancer_high + hpa_cancer_medium) counts
  - **Script**: `get_hpa.R` (applied via `add_hpa_annotations()`)
  - **Interpretation**:
    - "High" = ≥6 patients with High/Medium (>50% of cohort)
    - "Medium" = 3-5 patients with High/Medium (25-50% of cohort)
    - "Low" = 1-2 patients with High/Medium (<25% of cohort)
    - "None" = No patients with High/Medium expression

- **`hpa_expression_pattern`** 
  - **What**: Biological expression pattern across normal and cancer tissues
  - **Source**: Combined classification from hpa_normal_expressed + hpa_cancer_expressed
  - **Script**: `get_hpa.R` (applied via `add_hpa_annotations()`)
  - **Interpretation**:
    - "Constitutive" = Expressed in both normal and cancer (core urothelial genes)
    - "Lost_in_cancer" = Expressed in normal but absent in cancer (potential tumor suppressors)
    - "Cancer_upregulated" = Not in normal but detected in cancer (potential therapeutic targets)
    - "Not_expressed" = Absent in both (not relevant to urothelial biology)
    - "Insufficient_data" = Missing normal or cancer data
    - "Unknown" = Edge cases not fitting other categories

- **`hpa_data_quality`** 
  - **What**: Overall quality assessment of HPA data for this gene
  - **Source**: Based on data completeness and antibody reliability
  - **Script**: `get_hpa.R` (applied via `add_hpa_annotations()`)
  - **Interpretation**:
    - "High" = Both normal and cancer data available with Enhanced/Supported/Approved antibodies
    - "Medium" = Both datasets available but lower reliability
    - "Low" = Missing either normal or cancer dataset
    - "None" = No data available for either normal or cancer

---

## **16. Human Protein Atlas - Filter Flags**

Quality control flags created by `add_hpa_annotations()` function:

- **`hpa_filter_has_data`** 
  - **What**: Boolean indicating presence of usable HPA data
  - **Source**: TRUE if hpa_data_quality is "High" or "Medium"
  - **Script**: `get_hpa.R` (applied via `add_hpa_annotations()`)
  - **Interpretation**: TRUE = has both normal AND cancer data; FALSE = incomplete data

- **`hpa_filter_normal_quality`** 
  - **What**: Boolean for high-quality normal tissue protein expression
  - **Source**: TRUE if hpa_normal_expressed = TRUE AND hpa_normal_reliable = "High" or "Medium"
  - **Script**: `get_hpa.R` (applied via `add_hpa_annotations()`)
  - **Interpretation**: TRUE = reliably detected in normal bladder with validated antibody

- **`hpa_filter_cancer_quality`** 
  - **What**: Boolean for reliable cancer tissue protein expression
  - **Source**: TRUE if hpa_cancer_confidence = "High" or "Medium" (≥3 patients)
  - **Script**: `get_hpa.R` (applied via `add_hpa_annotations()`)
  - **Interpretation**: TRUE = detected in ≥25% of cancer patient samples

- **`hpa_passes_filter`** 
  - **What**: Combined quality control flag for HPA inclusion
  - **Source**: TRUE if has data AND (normal quality OR cancer quality)
  - **Script**: `get_hpa.R` (applied via `add_hpa_annotations()`)
  - **Interpretation**: TRUE = meets minimum quality criteria for HPA-based analysis

---

## **17. Human Protein Atlas - Biological Interpretation Flags**

Additional annotation flags for biological insights:

- **`hpa_high_confidence`** 
  - **What**: Boolean for highest quality protein expression data
  - **Source**: TRUE if hpa_data_quality = "High" AND hpa_cancer_confidence = "High" AND hpa_normal_reliable = "High"
  - **Script**: `get_hpa.R` (applied via `add_hpa_annotations()`)
  - **Interpretation**: TRUE = gold standard protein validation; maximum confidence

- **`hpa_therapeutic_target`** 
  - **What**: Boolean flag for potential therapeutic target genes
  - **Source**: TRUE if hpa_expression_pattern = "Cancer_upregulated" AND hpa_cancer_confidence = "High" or "Medium"
  - **Script**: `get_hpa.R` (applied via `add_hpa_annotations()`)
  - **Interpretation**: TRUE = protein upregulated in cancer with good patient evidence; candidate for targeted therapy

- **`hpa_potential_suppressor`** 
  - **What**: Boolean flag for potential tumor suppressor genes
  - **Source**: TRUE if hpa_expression_pattern = "Lost_in_cancer" AND hpa_normal_expressed = TRUE
  - **Script**: `get_hpa.R` (applied via `add_hpa_annotations()`)
  - **Interpretation**: TRUE = protein lost during cancer transformation; candidate tumor suppressor

---

## **Summary of Data Sources:**

1. **STRING Database** → Network topology metrics (degree, betweenness, etc.)
2. **Reactome Pathways** → Seed genes (is_seed)
3. **UroScan Bulk RNA-seq** → Expression values, ANOVA statistics, subtypes
4. **Single-cell RNA-seq (Normal)** → sc_norm_* columns
5. **Single-cell RNA-seq (Cancer)** → sc_cancer_* columns
6. **Gene Ontology (via clusterProfiler)** → Functional annotations
7. **Human Protein Atlas** → Protein-level validation (hpa_* columns)
   - Normal tissue IHC: `hpa_normal_*` columns
   - Cancer tissue IHC: `hpa_cancer_*` columns
   - Derived classifications: Expression patterns, quality metrics, filter flags
   - Processing: Downloaded via `download_hpa_data()`, annotated via `add_hpa_annotations()`

---

## **Data Processing Scripts:**

1. **`get_pathway_genes.R`** → Retrieves seed genes from Reactome
2. **`expand_network.R`** → Expands PPI network using STRING database
3. **`create_seurat_objects.R`** → Processes single-cell RNA-seq data
4. **`run_single_cell_annotations.R`** → Classifies genes in single-cell data
5. **`anova_sig_genes.R`** → Calculates ANOVA statistics and effect sizes
6. **`functional_enrichment_analysis.R`** → GO enrichment and functional roles
7. **`get_hpa.R`** → Downloads HPA data and adds protein annotations
   - Functions: `download_hpa_data()`, `add_hpa_annotations()`
   - Creates: Normal/cancer IHC data, derived classifications, filter flags
8. **`consoldiate_node_meta.R`** → Integrates all data layers into final dataset
   - Combines: Network, expression, statistics, single-cell, functional, HPA
   - Outputs: Pathway-specific and combined consolidated data files

---

This integrated dataset allows comprehensive characterization of genes in the receptor tyrosine kinase pathway networks from multiple biological perspectives: network importance, expression patterns across cancer subtypes, single-cell behavior, molecular function, and protein-level validation in normal and cancer tissues.

