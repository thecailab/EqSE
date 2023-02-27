# Identification of Expression Quantitative Super Enhancers (eSEs) using Single-cell Multiomic Data

## Author: Guoshuai Cai

<img src="https://github.com/GuoshuaiCai/KEEP/blob/main/workflow.png" width="800" height="800">

### Platform/data Compatibility
With necessary modifications, the proposed workflow is versatile and is logically adaptable to current established single-cell multi-omic measures of transcriptome and histone modification, including Paired-Tag,  CoTECH, SET-seq, EpiDamID and others. The below pipeline is constructed for Paired-Tag data processing. For data from other platforms, specific tools (e.g, reachtools for Paired-Tag) are likely needed for processing the data with their specific formats and settings such as barcodes and demultiplexing. Users can check with their orginal papers for appropriate tools and protocols and accordingly modify the STEPs 1-2 of this pipeline.

### Tool dependencies:

reachtools, https://github.com/cxzhu/Paired-Tag

macs, https://github.com/macs3-project/MACS, macs2 was used in this pipeline

ROSE, https://github.com/stjude/ROSE

Seurat, https://satijalab.org/seurat/, Seurat v4.0.2 was used in this pipeline

zinbwave, https://github.com/drisso/zinbwave

wcorr, https://cran.r-project.org/web/packages/wCorr/index.html

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

### STEP5 (eSE identification)
1. Perform standard correlation analysis to identify eSEs.
2. Integrative analysis and visualization of eSEs.


## Example of correlation calculation
```{r}
load("wcorr_func.R")
load("data_demo.RData")
```

### weight calculation
```{r}
w.exp<-weight.calc(count.exp, celltype=celltype)
w.histone<-weight.calc(count.histone, celltype=celltype)
```

### normalization
```{r}
dat.exp<-log2(t(t(count.exp)/colSums(count.exp)*median(colSums(count.exp)))+1)
dat.histone<-log2(t(t(count.histone)/colSums(count.histone)*median(colSums(count.histone)))+1)
```

### weighted corrlation calculation
```{r}
wcorr<-wcorr.calc(dat.exp, dat.histone,w.exp,w.histone,celltype=celltype)
```

### plot
```{r}
wcorr.plot(dat.exp, dat.histone,w.exp,w.histone,celltype=celltype, cell="FcNeu", peak="chr2:73584986-73596983_Chn1",color="red")
```
<img src="https://github.com/GuoshuaiCai/KEEP/blob/main/plot_demo.png" width="800" height="800">
