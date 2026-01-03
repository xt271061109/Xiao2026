## Loading packages
library(Seurat)
library(Signac)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(ComplexHeatmap)
library(EnhancedVolcano)
library(colorRamp2)
library(clusterProfiler)
library(enrichplot)
library(SingleCellExperiment)
library(slingshot)
library(scater)
library(AUCell)
library(msigdbr)
library(biomaRt)
library(httr)
library(org.Mm.eg.db)
library(tidyr)
library(viridis)
library(GenomeInfoDb)
library(GenomicRanges)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(chromVAR)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(biovizBase)
library(RPresto)
library(glmGamPoi)
library(GenomicRanges)
library(EnrichedHeatmap)
library(RColorBrewer)



#read WT
## Change working directory
WT_dir <-
  '/fs/ess/PAS2205/Tong_xiao/Single-cell_multiome/processed/WT_filtered_feature_bc_matrix'
list.files(WT_dir)
## Read data
WT <- Read10X(data.dir = WT_dir)
## extract RNA and ATAC data
WT_rna_counts <- WT$`Gene Expression`
WT_atac_counts <- WT$Peaks
## Create Seurat object for WT
WT <- CreateSeuratObject(counts = WT_rna_counts, project = "WT")
WT[["percent.mt"]] <- PercentageFeatureSet(WT, pattern = "^mt-")
## Now add in the ATAC-seq data
## we'll only use peaks in standard chromosomes
grange.counts <-
  StringToGRanges(rownames(WT_atac_counts), sep = c(":", "-"))
grange.use <-
  seqnames(grange.counts) %in% standardChromosomes(grange.counts)
WT_atac_counts <- WT_atac_counts[as.vector(grange.use),]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"


frag.file <-
  "/fs/ess/PAS2205/Tong_xiao/Single-cell_multiome/processed/WT_atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
  counts = WT_atac_counts,
  sep = c(":", "-"),
  genome = 'mm10',
  fragments = frag.file,
  min.cells = 10,
  annotation = annotations
)
WT[["ATAC"]] <- chrom_assay

#read KO
## Change working directory
KO_dir <-
  '/fs/ess/PAS2205/Tong_xiao/Single-cell_multiome/processed/KO_filtered_feature_bc_matrix'
list.files(KO_dir)
## Read data
KO <- Read10X(data.dir = KO_dir)
## extract RNA and ATAC data
KO_rna_counts <- KO$`Gene Expression`
KO_atac_counts <- KO$Peaks
## Create Seurat object for KO
KO <- CreateSeuratObject(counts = KO_rna_counts, project = "KO")
KO[["percent.mt"]] <- PercentageFeatureSet(KO, pattern = "^mt-")
## Now add in the ATAC-seq data
## we'll only use peaks in standard chromosomes
grange.counts <-
  StringToGRanges(rownames(KO_atac_counts), sep = c(":", "-"))
grange.use <-
  seqnames(grange.counts) %in% standardChromosomes(grange.counts)
KO_atac_counts <- KO_atac_counts[as.vector(grange.use),]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

frag.file <-
  "/fs/ess/PAS2205/Tong_xiao/Single-cell_multiome/processed/KO_atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
  counts = KO_atac_counts,
  sep = c(":", "-"),
  genome = 'mm10',
  fragments = frag.file,
  min.cells = 10,
  annotation = annotations
)
KO[["ATAC"]] <- chrom_assay

#combine WT and KO into one seurat object
combined <- merge(
  WT,
  y = KO,
  add.cell.ids = c("WT", "KO"),
  project = "LCMV"
)
table(combined$orig.ident)
# Join layers (assays) from both objects, for example, RNA and atac assays
combined <-
  JoinLayers(combined,
             layers = list("RNA", "ATAC"),
             merge.data = TRUE)

#QC
VlnPlot(
  combined,
  features = c("nCount_ATAC", "nCount_RNA", "percent.mt"),
  ncol = 3,
  log = TRUE,
  pt.size = 0
) + NoLegend()
combined <- subset(
  x = combined,
  subset = nCount_ATAC < 7e4 &
    nCount_ATAC > 5e3 &
    nCount_RNA < 25000 &
    nCount_RNA > 1000 &
    percent.mt < 20
)
combined
#preprocessing
# RNA analysis
DefaultAssay(combined) <- "RNA"
combined <-
  SCTransform(combined, verbose = FALSE) %>% RunPCA() %>% RunUMAP(
    dims = 1:50,
    reduction.name = 'umap.rna',
    reduction.key = 'rnaUMAP_'
  )

head(combined$SCT@SCTModel.list)

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(combined) <- "ATAC"
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 'q0')
combined <- RunSVD(combined)
combined <- RunUMAP(
  combined,
  reduction = 'lsi',
  dims = 2:50,
  reduction.name = "umap.atac",
  reduction.key = "atacUMAP_"
)

#calculate a WNN graph
combined <- FindMultiModalNeighbors(
  combined,
  reduction.list = list("pca", "lsi"),
  dims.list = list(1:50, 2:50)
)
combined <- RunUMAP(
  combined,
  nn.name = "weighted.nn",
  reduction.name = "wnn.umap",
  reduction.key = "wnnUMAP_"
)
combined <- FindClusters(
  combined,
  graph.name = "wsnn",
  algorithm = 3,
  verbose = FALSE,
  resolution = 0.1
)


#--------------------------------------------------------------------------------
#Single-cell ATAC-seq analysis
DefaultAssay(combined) <- "ATAC"

# first compute the GC content for each peak
combined <-
  RegionStats(combined, genome = BSgenome.Mmusculus.UCSC.mm10)
# link peaks to genes
combined <- LinkPeaks(object = combined,
                      peak.assay = "ATAC",
                      expression.assay = "SCT",
)
# compute gene activities
#gene.activities <- GeneActivity(combined)
gene.activities <- GeneActivity(combined, extend.upstream = 5000,
                                extend.downstream = 0,)

# add the gene activity matrix to the Seurat object as a new assay
combined[['geneactivity']] <-
  CreateAssayObject(counts = gene.activities)
combined <- NormalizeData(
  object = combined,
  assay = 'geneactivity',
  normalization.method = 'LogNormalize',
  scale.factor = median(combined$nCount_RNA)
)
