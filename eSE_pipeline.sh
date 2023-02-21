#Tool dependencies:

#reachtools, https://github.com/cxzhu/Paired-Tag

#macs, https://github.com/macs3-project/MACS, macs2 was used in this pipeline.

#ROSE, https://github.com/stjude/ROSE

#Seurat, https://satijalab.org/seurat/, Seurat v4.0.2 was used in this pipeline.

#Others: Tool dependencies of all above tools.


#STEP 1
#processing same-cell multi-omic data, following the protocol in https://github.com/cxzhu/Paired-Tag

module load samtools/gcc/1.5
samtools merge -b bamlist_dna.txt -o merged_dna_mm10_sorted_rmdup.bam

#split
perl splitbam.pl merged_dna_mm10_sorted_rmdup.bam

#STEP 2
#Seurat analysis to identify cell populations with lists of barcodes.

#STEP3: Identify SEs in each cell type. Using FcNeu as an example,

perl filterbam.pl H3K27ac_merged_dna_mm10_sorted_rmdup.bam FcNeu_barcodes.txt

samtools sort -o FcNeu_H3K27ac_merged_dna_mm10_sorted_rmdup_sorted.bam FcNeu_H3K27ac_merged_dna_mm10_sorted_rmdup.bam
samtools index FcNeu_H3K27ac_merged_dna_mm10_sorted_rmdup_sorted.bam

module load python3/anaconda/2021.11
macs2 callpeak -t FcNeu_H3K27ac_merged_dna_mm10_sorted_rmdup_sorted.bam -f BAM -g mm -n FcNeu_H3K27ac_merged_dna_macs2 -B -q 0.01

python ROSE_main.py -g MM10 -i FcNeu_H3K27ac_merged_dna_macs2_summits.bed -r FcNeu_H3K27ac_merged_dna_mm10_sorted_rmdup_sorted.bam -o ROSE


#pool super enhancers
cat ROSE/FcNeu_H3K27ac_merged_dna_macs2_summits_SuperStitched.bed ROSE/HcNeu_H3K27ac_merged_dna_macs2_summits_SuperStitched.bed ROSE/InNeu_H3K27ac_merged_dna_macs2_summits_SuperStitched.bed ROSE/Other_H3K27ac_merged_dna_macs2_summits_SuperStitched.bed > ROSE/All_H3K27ac_merged_dna_macs2_summits_SuperStitched.bed
cat ROSE/FcNeu_H3K27ac_merged_dna_macs2_summits_SuperStitched_GENE_TO_REGION.txt ROSE/HcNeu_H3K27ac_merged_dna_macs2_summits_SuperStitched_GENE_TO_REGION.txt ROSE/InNeu_H3K27ac_merged_dna_macs2_summits_SuperStitched_GENE_TO_REGION.txt ROSE/Other_H3K27ac_merged_dna_macs2_summits_SuperStitched_GENE_TO_REGION.txt > ROSE/All_H3K27ac_merged_dna_macs2_summits_SuperStitched_GENE_TO_REGION.txt
cat ROSE/FcNeu_H3K27ac_merged_dna_macs2_summits_SuperStitched_REGION_TO_GENE.txt ROSE/HcNeu_H3K27ac_merged_dna_macs2_summits_SuperStitched_REGION_TO_GENE.txt ROSE/InNeu_H3K27ac_merged_dna_macs2_summits_SuperStitched_REGION_TO_GENE.txt ROSE/Other_H3K27ac_merged_dna_macs2_summits_SuperStitched_REGION_TO_GENE.txt > ROSE/All_H3K27ac_merged_dna_macs2_summits_SuperStitched_REGION_TO_GENE.txt

reachtools bam2Mtx2 H3K27ac_merged_dna_mm10_sorted_rmdup_sorted.bam ROSE/All_H3K27ac_merged_dna_macs2_summits_SuperStitched.bed
#generating H3K27ac_merged_dna_mm10_sorted_rmdup_sorted_mtx2

#STEP5
#conjugating histone modification data into above seurat object
#perform standard correlation analysis to identify eSEs (Fig. )
#integrative analysis and visualization of eSEs (Fig. )

