# Metastasis-scRNA-seq-
Run this in R-studios after downloading Seurat library and the cell sample. Change the link from line 6 to where the downloaded cell info is on the computer.  

Seurat can be dowloaded through the link: https://satijalab.org/seurat/ 

Also download the ggplot2, dplyr, and patchwork packages

Sample sample can be downloaded by the link: https://www.10xgenomics.com/resources/datasets/20-k-mixture-of-nsclc-dt-cs-from-7-donors-3-v-3-1-3-1-standard-6-1-0
Click 'donor 2', scroll to the bottom and click "gene expression - feature / cell matrix (per sample)"
the dowload zipped file will appear as "20k_NSCLC_DTC_3p_nextgem_donor_2_count_sample_feature_bc_matrix.tar.gz" in your files. 
Unzip it and you should have a file folder named "sample_feature_bc_matrix" with three files inside: barcodes.tsv.gz , features.tsv.gz , and matrix.mtx.gz
Acquire the pathname of the folder through double click + OPTION + Copy "sample_feature_bc_matrix" as Pathname.
Paste the pathname to line 6 of the code
