library('scRNAseq')
library('BiocFileCache')
bfc <- BiocFileCache()
raw.path <-
  bfcrpath(
    bfc,
    file.path(
      "http://cf.10xgenomics.com/samples",
      "cell-exp/2.1.0/pbmc4k/pbmc4k_raw_gene_bc_matrices.tar.gz"
    )
  )
untar(raw.path, exdir = file.path(tempdir(), "pbmc4k"))

library('DropletUtils')
library('Matrix')
fname <- file.path(tempdir(), "pbmc4k/raw_gene_bc_matrices/GRCh38")
sce.pbmc <- read10xCounts(fname, col.names = TRUE)
library('scater')
rownames(sce.pbmc) <- uniquifyFeatureNames(rowData(sce.pbmc)$ID, rowData(sce.pbmc)$Symbol)
library('EnsDb.Hsapiens.v86')
location <- mapIds(
  EnsDb.Hsapiens.v86,
  keys = rowData(sce.pbmc)$ID,
  column = "SEQNAME",
  keytype = "GENEID"
)
# cell-detection
set.seed(100)
e.out <- emptyDrops(counts(sce.pbmc))
sce.pbmc <- sce.pbmc[, which(e.out$FDR <= 0.001)]
#The function emptyDrops was applied to determine which droplets are empty
#based on the counts for each gene in each cell. A threshold of .001 was set
#to make sure that at most .001 of the droplets which show a significant 
#deviation from the null model, which assumes that transcripts are passed from the
#environment into the droplets in a way that can be modelled by multinomial sampling, 
#are false positives. Given that the algorithm behind this function depends on 
#simulations, a fixed seed is needed to make the code reproducible. 

# quality-control
stats <- perCellQCMetrics(sce.pbmc, subsets = list(Mito = which(location == "MT")))

high.mito <- isOutlier(stats$subsets_Mito_percent,
                       type = "higher")
sce.pbmc <- sce.pbmc[, !high.mito]
#From the droplets assumed to be not empty given the FDR threshold of .001,
# we calculated the number of mitochondrial transcripts in each of these 
# as a quality measure. Outliers were detected based on an interval given by the median 
# of these transcript numbers and 3 median deviations around this value, everything
# outside this interval is classified as an outlier. Finally, we kept those cells which 
# satisfy that they are not an outlier based on this procedure.

# normalization
library('scran')
set.seed(1000)
clusters <- quickCluster(sce.pbmc)
sce.pbmc <- computeSumFactors(sce.pbmc, cluster = clusters)
sce.pbmc <- logNormCounts(sce.pbmc)
# Just like in the emptyDrops procedure, a seed is fixed for reproducibility purposes.
# Clusters are obtained using the quickCluster function. Given the similarity of the
# clusters, a size factor value for each of them is obtained by taking the median value
# of calculations computed on multiple pools taken from the cluster,using the function 
# computeSumFactors. Finally, we can divide the number of counts between the library 
#size factor of the gene in the cluster to  normalize, and use the log function to 
#observe significant differences in the number of counts.
