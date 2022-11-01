# Identification of Expression Quantitative Super Enhancers (EqSEs) using Single-cell Multiomic Data

#### Tool dependencies:

reachtools, https://github.com/cxzhu/Paired-Tag

macs, https://github.com/macs3-project/MACS, macs2 was used in this pipeline.

ROSE, https://github.com/stjude/ROSE

Seurat, https://satijalab.org/seurat/, Seurat v4.0.2 was used in this pipeline.

Others: Tool dependencies of all above tools.

#### STEP 1 (sequencing data processing)
1. Process same-cell multi-omic data of each sub-libarary, following the protocol in https://github.com/cxzhu/Paired-Tag.
2. Merge bam files.
3. (Optional) Extracting H3K27ac data from pooled data of different histones.

#### STEP 2 (cell type identification)
1. Constrcut Seurat object with single-cell transcriptome data.
2. Identify cell types and generate corresponding barcode lists.

#### STEP 3 (SE identification and quantification)
1. Extract data for the specific cell type.
2. Call peaks of H3K27ac modification.
3. Identify SE based on the rank of the H3K27ac signal.
4. Pool SEs identified in all cell types.
5. Quantify histone modification at all SEs in each single cell.

#### STEP4 (data conjugation)
Conjugate histone modification data into above seurat object.

#### STEP5 (EqSE identification)
#Perform standard correlation analysis to identify EqSEs (Fig. )
#Integrative analysis and visualization of EqSEs (Fig. )
