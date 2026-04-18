# 1. Executive Summary

This report presents a comprehensive bulk RNA-seq analysis of MCF-7
breast cancer cells cultured under control conditions and in coculture
with adipose-derived stem cells (ADSCs).

The analysis identifies transcriptional changes associated with tumor
microenvironment interactions, including differential gene expression,
functional enrichment, and pathway-level reprogramming.

Key findings indicate strong modulation of extracellular matrix
organization, hypoxia response, metabolic reprogramming, and pro-growth
signaling pathways.

------------------------------------------------------------------------

# 2. Experimental Design

The dataset consists of 6 RNA-seq samples:

<table>
<thead>
<tr>
<th>Condition</th>
<th>Replicates</th>
<th>Description</th>
</tr>
</thead>
<tbody>
<tr>
<td>Control</td>
<td>3</td>
<td>MCF-7 cells cultured alone</td>
</tr>
<tr>
<td>ADSC</td>
<td>3</td>
<td>MCF-7 cells cocultured with ADSCs</td>
</tr>
</tbody>
</table>

No batch effects were reported unless otherwise specified.

------------------------------------------------------------------------

# 3. Data Processing Workflow

RNA-seq preprocessing was performed using a standard pipeline:

-   Quality control: FastQC / MultiQC
-   Adapter trimming: Trimmomatic
-   Alignment: HISAT2
-   Gene quantification: featureCounts

Raw counts were used for downstream statistical analysis.

Analysis was performed following established Bioconductor workflows for
bulk RNA-seq differential expression analysis, ensuring reproducibility,
analytical transparency, and adherence to widely accepted computational
standards.

## 3.1 Statistical framework

Differential expression analysis was performed using edgeR, which models
count data using a negative binomial distribution. Library normalization
was performed using the TMM method to account for compositional
differences between samples.

Statistical significance was defined using the false discovery rate
(FDR) framework to control for multiple hypothesis testing.

------------------------------------------------------------------------

# 4. Exploratory Data Analysis (EDA)

## 4.1 Principal Component Analysis (PCA)

PCA was used to assess global transcriptional differences between
conditions.

<img src="../downstream/plots/ExactTest/EDA/04-PCA_log2CPM.png" alt="" width="1200" style="display: block; margin: auto;" />

PCA revealed a clear separation between control and ADSC-treated
samples, indicating strong global transcriptional reprogramming induced
by tumor microenvironment interactions.

## 4.2 Sample Clustering

Sample-distance heatmap was performed to assess sample similarity.

<img src="../downstream/plots/ExactTest/EDA/05-sample_distance_heatmap.png" alt="" width="1200" style="display: block; margin: auto;" />

------------------------------------------------------------------------

# 5. Differential Expression Analysis (DEGs)

Differential expression analysis was performed using edgeR with TMM
normalization.

**Criteria for significance:**

    ## - FDR < 0.05

    ## - |log2FC| > 1

## 5.1 Summary of DEGs

    ## Total DEGs identified: 1013

    ## Upregulated genes: 887

    ## Downregulated genes: 126

## 5.2 Volcano Plot of DEGs

<img src="../downstream/plots/ExactTest/Differential_Expression/09-volcano_plot_DEGs.png" alt="" width="1200" style="display: block; margin: auto;" />

## 5.3 Heatmap of DEGs

<img src="../downstream/plots/ExactTest/Differential_Expression/10-DEGs_heatmap.png" alt="" width="1200" style="display: block; margin: auto;" />

------------------------------------------------------------------------

# 6. Functional Enrichment Analysis

Functional enrichment analysis was performed using GO, KEGG, Reactome,
and GSEA.

## 6.1 Gene Ontology terms

### **Biological Process**

-   Response to hypoxia
-   Extracellular matrix organization
-   Cellular response to oxygen levels

### **Molecular Functions**

-   Extracellular matrix structural constituent
-   Monosaccharide binding
-   Oxidoreductase activity

### **Cellular Component**

-   Vesicle lumen
-   Basement membrane
-   Extracellular matrix

## 🧠 Biological interpretation (GO enrichment)

The GO enrichment analysis highlights a coordinated transcriptional
response involving hypoxia signaling, extracellular matrix remodeling,
and metabolic adaptation. The enrichment of response to hypoxia and
cellular response to oxygen levels suggests activation of oxygen-sensing
pathways consistent with tumor microenvironment stress conditions.

In parallel, extracellular matrix organization and related structural
components indicate active remodeling of the tumor microenvironment, a
key feature associated with increased cell adhesion, migration, and
stromal interaction in breast cancer progression.

At the molecular level, enrichment of extracellular matrix structural
components and oxidoreductase activity reflects extracellular remodeling
processes coupled with metabolic and redox adaptation. Finally, the
presence of extracellular and vesicle-associated compartments supports
enhanced intercellular communication within the tumor microenvironment.

Overall, these results are consistent with early pro-tumorigenic
transcriptional reprogramming driven by stromal interaction, as
described in ADSC-cocultured MCF-7 models.

## 6.2 KEGG Pathways

Significant pathways include:

-   HIF-1 signaling pathway
-   Glycolysis / gluconeogenesis
-   Complement and coagulation cascades

## 6.3 Reactome Pathways

-   Interleukin-4 and Interleukin-13 signaling
-   Fibronectin matrix formation
-   Assembly of collagen fibrils and other multimeric structures

## 🧠 Biological interpretation (KEGG and Reactome pathway)

The KEGG and Reactome enrichment analyses consistently highlight major
biological processes associated with tumor microenvironment remodeling
and metabolic adaptation. The enrichment of the HIF-1 signaling pathway
together with glycolysis/gluconeogenesis suggests a metabolic shift
toward hypoxia-driven glycolytic reprogramming, a well-established
feature of early tumor progression and cellular adaptation to low oxygen
conditions.

In parallel, pathways such as ECM-receptor interaction, fibronectin
matrix formation, and collagen organization indicate strong
extracellular matrix remodeling and reorganization of cell–matrix
interactions, supporting increased cell adhesion, migration, and stromal
communication.

Additionally, the activation of immune-related pathways, including
complement and coagulation cascades and interleukin-4/13 signaling,
suggests a pro-inflammatory microenvironment potentially driven by
stromal–tumor crosstalk.

Overall, these pathways converge on a biological program characterized
by hypoxia adaptation, metabolic reprogramming, and extracellular matrix
remodeling, consistent with ADSC-induced pro-tumorigenic transcriptional
changes in MCF-7 cells.

------------------------------------------------------------------------

# 7. Network Analysis

Protein-protein interaction (PPI) analysis was performed using STRING
database.

Hub genes were identified using network topology analysis (cytoHubba MCC
method)

<img src="../downstream/plots/PPI/Hub_genes.png" alt="" width="1837" style="display: block; margin: auto;" />

------------------------------------------------------------------------

# 8. Integrated Biological Interpretation

The transcriptional profile of ADSC-treated MCF-7 cells reveals
coordinated changes across multiple biological systems.

A strong enrichment of extracellular matrix organization and integrin
signaling pathways suggests active remodeling of the tumor
microenvironment, potentially enhancing cell adhesion and migratory
capacity.

In parallel, hypoxia-related and metabolic pathways indicate a shift
toward glycolytic metabolism and cellular stress adaptation, consistent
with early tumor progression mechanisms.

Additionally, activation of cytokine signaling and JAK-STAT pathways
highlights crosstalk between tumor cells and stromal components,
suggesting a pro-inflammatory microenvironment.

Finally, enrichment of growth factor signaling and FOXO-mediated
transcription suggests enhanced proliferative and survival signaling.

Overall, these results indicate that ADSC coculture promotes a
pro-tumorigenic transcriptional program in MCF-7 cells.

------------------------------------------------------------------------

# 9. Conclusions

1.  ADSC coculture induces strong transcriptional reprogramming in MCF-7
    cells
2.  ECM remodeling and hypoxia response are key biological features
3.  Metabolic and inflammatory pathways are significantly altered
4.  Results are consistent with a pro-tumorigenic tumor microenvironment
    effect

------------------------------------------------------------------------

# 10. Reproducibility

-   R version: v.4.5.1
-   Main Packages: edgeR, clusterProfiler, ReactomePA, fgsea, WGCNA
-   Session information available upon request
-   Recommended: use renv for environment reconstruction

------------------------------------------------------------------------

# 11. Study limitations and considerations

This analysis is based on a limited sample size (n = 3 per group), which
is common in RNA-seq experimental designs but may limit statistical
power for detecting low-effect genes.

Differential expression analysis was performed using edgeR, which
assumes negative binomial distribution and is robust for small sample
sizes, but results should be interpreted in the context of biological
variability.

No batch effects were detected or reported; however, batch correction
methods (e.g., surrogate variable analysis) may be required in more
complex experimental designs.

Finally, functional enrichment results depend on existing pathway
annotations and should be interpreted as hypothesis-generating rather
than definitive biological conclusions

------------------------------------------------------------------------

# 12. Deliverables (client-ready output)

This analysis provides structured outputs designed for downstream
interpretation and decision-making, including:

-   Differential expression results with statistical confidence (log2FC,
    FDR, CPM)
-   Ranked gene lists for biomarker prioritization
-   Functional pathway summaries highlighting key dysregulated
    biological processes
-   Network-based identification of hub genes potentially involved in
    regulatory control
-   Visual reports (PCA, volcano plots, heatmaps) suitable for
    publication or presentations
-   Biological interpretation summary to support hypothesis generation

------------------------------------------------------------------------

# 13. Biological and translational relevance

The observed transcriptional changes reflect key features of tumor
microenvironment-driven cancer progression, including extracellular
matrix remodeling, hypoxia adaptation, and metabolic reprogramming.

These processes are consistent with features associated with increased
tumor aggressiveness and may provide insights into early molecular
events linked to breast cancer progression in stromal-rich environments.

Rather than isolated pathway effects, the results indicate a coordinated
tumor microenvironment-driven transcriptional program organized into
three functional biological modules:

1.  **Microenvironment remodeling module** (ECM organization, integrin
    signaling, collagen formation)  
2.  **Metabolic adaptation module** (HIF-1 signaling, glycolysis,
    oxidative stress response)  
3.  **Immune-stromal communication module** (cytokine signaling,
    complement activation, JAK-STAT pathway)

These findings may provide actionable insights into stromal–tumor
interactions and could help prioritize candidate pathways for further
experimental validation or therapeutic targeting.

------------------------------------------------------------------------

# 14. Potential applications

This type of analysis can support:

-   Biomarker discovery for cancer progression studies
-   Identification of therapeutic targets in tumor microenvironment
    interactions
-   Functional characterization of treatment-induced transcriptional
    changes
-   Preclinical hypothesis generation for breast cancer research
-   Integration with proteomics or clinical datasets for multi-omics
    analysis
