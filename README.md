# Identification of Expression Quantitative Super Enhancers (EqSEs) using Single-cell Multiomic Data

## Author: Guoshuai Cai

![Image_of_workflow](https://github.com/GuoshuaiCai/EqSE/blob/main/workflow.png|width=100)

### Tool dependencies:

reachtools, https://github.com/cxzhu/Paired-Tag

macs, https://github.com/macs3-project/MACS, macs2 was used in this pipeline.

ROSE, https://github.com/stjude/ROSE

Seurat, https://satijalab.org/seurat/, Seurat v4.0.2 was used in this pipeline.

Others: Tool dependencies of all above tools.

### STEP 1 (sequencing data processing)
1. Process same-cell multi-omic data of each sub-libarary, following the protocol in https://github.com/cxzhu/Paired-Tag.
2. Merge bam files.
3. (Optional) Extracting H3K27ac data from pooled data of different histones.

### STEP 2 (cell type identification)
1. Constrcut Seurat object with single-cell transcriptome data.
2. Identify cell types and generate corresponding barcode lists.

In the mouse brain study, the RNA sequencing read counts were normalized using the “LogNormalize” method and further centered and scaled by standard deviations. 2000 highly variable feature were selected, and dimension of their space was reduced using the top 20 principal components. Further, cell clusters were identified using the Louvain optimization-based clustering method, and were visualized in reduced dimensions of Uniform Manifold Approximation and Projection (UMAP). Based on the expression of markers that were used in the original study of GSE152020, 18 cell types were identified, including 7 FC cell types (FC L2/3, FC L4, FC L5, FC PT, FC NP, FC CT, FC L6), 3 HC cell types (HC CA1, HC CA2/3, HC DG), 3 Interneurons (InNeu-CGE, InNeu-Sst, InNeu-Pvalb), Oligodendrocyte precursor cells (OPC), Oligodendrocytes (Oligo), Astrocytes, Endothelial cells and Ependymal cells (Fig. S1).

### STEP 3 (SE identification and quantification)
1. Extract data for the specific cell type.
2. Call peaks of H3K27ac modification.
3. Identify SE based on the rank of the H3K27ac signal.
4. Pool SEs identified in all cell types.
5. Quantify histone modification at all SEs in each single cell.

### STEP4 (data conjugation)
Conjugate histone modification data into above seurat object.

### STEP5 (EqSE identification)
1. Perform standard correlation analysis to identify EqSEs.
2. Integrative analysis and visualization of EqSEs.
