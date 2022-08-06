

library(dplyr)
library(Seurat)

data_dir <- '/Users/peterwang/Downloads/sample_feature_bc_matrix'
data <- Read10X(data.dir = data_dir)
cells = CreateSeuratObject(counts = data$`Gene Expression`)

cells[["percent.mt"]] <- PercentageFeatureSet(cells, pattern = "^MT-")
VlnPlot(cells, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

cells <- subset(cells, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 25)
VlnPlot(cells, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

cells <- NormalizeData(cells, normalization.method = "LogNormalize", scale.factor = 10000) 
cells <- FindVariableFeatures(cells, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(cells), 10) 
all.genes <- rownames(cells)
plot1 <- VariableFeaturePlot(cells)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
cells <- ScaleData(cells, features = all.genes)
cells <- RunPCA(cells, features = VariableFeatures(object = cells))

DimHeatmap(cells, dims = 1, cells = 500, balanced = TRUE)

cells <- JackStraw(cells, num.replicate = 200)
cells <- ScoreJackStraw(cells, dims = 1:20)
JackStrawPlot(cells, dims = 1:20)

cells <- FindNeighbors(cells, dims = 1:15)
cells <- FindClusters(cells, resolution = 0.5)
cells <- RunUMAP(cells, dims = 1:15)

DimPlot(cells, reduction = 'umap')

markers <- FindAllMarkers(cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
genes <- markers %>% group_by(cluster) %>% slice_max(n = 10, order_by = avg_log2FC)

write.csv(genes, '/Users/peterwang/Downloads/Gene_Clusters_info2.csv')

FeaturePlot(cells, features = c('CD3E', 'CD8A', 'MS4A1', 'IGKC', 'CD14')) # immune cells (macrophages, B-cells, plasma Cells, T-cells)

new.cluster.ids <- c('0', 'T Cells', 'T Cells', '3', 'Monocytes', '5', '6', 'B Cells', '8', 'Plasma Cells', '10', '11', '12') # keep the variables that we have set
names(new.cluster.ids) <- levels(cells)
cells <- RenameIdents(cells, new.cluster.ids)
cells <- RunUMAP(cells, dims = 1:15)
DimPlot(cells, reduction = 'umap', label = TRUE, pt.size = 0.5) + NoLegend()



FeaturePlot(cells, features = c('ADAM12', 'CDH11', 'CDH2', 'COL1A1', 'COL3A1', 'COL5A1', 'COL6A1', 'COL6A2', 'CTGF')) # upregulated in cell adhesion and migration
FeaturePlot(cells, features = c('CYP1B1', 'CYP1B1', 'CYP1B1', 'DLC1', 'FBLN1', 'FBLN5', 'FGFR1', 'FN1', 'HAS2', 'LUM'))
FeaturePlot(cells,features = c('MMP2', 'MYL9', 'NID2', 'NR2F1', 'NRP1', 'PLAT', 'PPAP2B', 'PRKCA', 'RECK', 'SERPINE1'))
FeaturePlot(cells,features = c('SERPINE2', 'SPOCK1', 'TGM2', 'TNFAIP6', 'TPM1', 'VCAN', 'WNT5A'))

FeaturePlot(cells, features = c('CD24', 'CDH1', 'CXADR', 'CXCL16', 'DSG3', 'DSG3', 'ELF3')) # down regulated in cell adhesion and migration
FeaturePlot(Cells, features = c( 'EPCAM', 'EPHA', 'JUP', 'MPZL2', 'OVOL2', 'PLXNB1', 'S100P', 'SLC7A5', 'SYK'))

FeaturePlot(cells, features = c('CDKN2C','EMP3','FBN1', 'IGFBP3', 'IL1R1', 'LTBP1', 'MME', 'PMP22' )) # upregulated in development, cell differentiatioin, and proliferation
FeaturePlot(cells, features = c('PTGER2', 'PTX3', 'SRGN', 'SULF1', 'SYNE1', 'TAGLN', 'TUBA1A', 'TUBA1A', 'VIM', 'ZEB1'))

FeaturePlot(cells, features = c('ABLIM1', 'ADRB2', 'ALDH1A3', 'ANK3', 'BIK CA2', 'CTSL2', 'FGFR2', 'FGFR3', 'FST'))  # down regulated in development, cell differentiation, and proliferation
FeaturePlot(cells,features = c('GJB3', 'IFI30', 'IL18', 'KLK7', 'KRT15', 'KRT17', 'LSR'))
FeaturePlot(cells, features = c('MAP7', 'MBP', 'OCLN', 'PKP2', 'PPL', 'PRSS8', 'RAPGEF5', 'SPINT1'))

FeaturePlot(cells, features = c('DCN','LOX','TEP1')) # upregulated in angiogenesis and wound healing

FeaturePlot(cells, features = c('ABCA1','GALNT10','SLC22A4')) # upregulated in metabolism

FeaturePlot(cells, features = c('GPX3','SLC27A2','SMPDL3B', 'SORL1', 'ST6GALNAC2')) # down regulated in metabolism

FeaturePlot(cells, features = c('C5ORF13','CDK14','EML1', 'FSTL1', 'LTBP2', 'MAP1B', 'RGS4', 'SYT11', 'TMEM158')) # upregulated in others or unclassified

FeaturePlot(cells, features = c('AGR2', 'C10ORF10', 'CDS1', 'FAM169A', 'FXYD3', 'KLK10', 'LAD1', 'MTUS1', 'PLS1')) # down regulated in others or unclassififed
FeaturePlot(cells, features = c('PRRG4', 'RHOD', 'SERPINB1', 'SLPI', 'TMEM30B', 'TPD52L1', 'TSPAN1', 'ZHX2', 'ZNF165'))

FeaturePlot(cells, features = c('AGR2','C10ORF10','CDS1', 'FAM169A', 'FXYD3', 'KLK10', 'LAD1', 'MTUS1', 'PLS1', 'PRRG4', 'RHOD', 'SERPINB1', 'SLP1', 'TMEM30B', 'TPD52L1', 'TSPAN1', 'ZHX2', 'ZNF165')) # downregulated in others or unclassified
