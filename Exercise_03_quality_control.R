## Data
library('scRNAseq')
sce.416b <- LunSpikeInData(which = "416b")
sce.416b$block <- factor(sce.416b$block)

# Download the relevant Ensembl annotation database
# using AnnotationHub resources
library('AnnotationHub')
ah <- AnnotationHub()
query(ah, c("Mus musculus", "Ensembl", "v97"))
# Annotate each gene with its chromosome location
ens.mm.v97 <- ah[["AH73905"]]
location <- mapIds(
  ens.mm.v97,
  keys = rownames(sce.416b),
  keytype = "GENEID",
  column = "SEQNAME"
)
# Identify the mitochondrial genes
is.mito <- which(location == "MT")

library('scater')
sce.416b <- addPerCellQC(sce.416b,
                         subsets = list(Mito = is.mito))

###Exercises to descard outliers using thresholds

qc.lib <- sce.416b$sum < 100000
qc.nexprs <- sce.416b$detected < 5000
qc.spike <- sce.416b$altexps_ERCC_percent > 10
qc.mito <- sce.416b$subsets_Mito_percent > 10
discard <- qc.lib | qc.nexprs | qc.spike | qc.mito

qc.lib2 <- isOutlier(sce.416b$sum, log = TRUE, type = "lower")
qc.nexprs2 <- isOutlier(sce.416b$detected, log = TRUE,
                        type = "lower")
qc.spike2 <- isOutlier(sce.416b$altexps_ERCC_percent,
                       type = "higher")
qc.mito2 <- isOutlier(sce.416b$subsets_Mito_percent,
                      type = "higher")
discard2 <- qc.lib2 | qc.nexprs2 | qc.spike2 | qc.mito2

batch <- paste0(sce.416b$phenotype, "-", sce.416b$block)
qc.lib3 <- isOutlier(sce.416b$sum,
                     log = TRUE,
                     type = "lower",
                     batch = batch)
qc.nexprs3 <- isOutlier(sce.416b$detected,
                        log = TRUE,
                        type = "lower",
                        batch = batch)
qc.spike3 <- isOutlier(sce.416b$altexps_ERCC_percent,
                       type = "higher",
                       batch = batch)
qc.mito3 <- isOutlier(sce.416b$subsets_Mito_percent,
                      type = "higher",
                      batch = batch)
discard3 <- qc.lib3 | qc.nexprs3 | qc.spike3 | qc.mito3

# Summarize the number of cells removed for each reason
DataFrame(
  LibSize = sum(qc.lib),
  NExprs = sum(qc.nexprs),
  SpikeProp = sum(qc.spike),
  MitoProp = sum(qc.mito),
  Total = sum(discard)
)
DataFrame(
  LibSize = sum(qc.lib2),
  NExprs = sum(qc.nexprs2),
  SpikeProp = sum(qc.spike2),
  MitoProp = sum(qc.mito2),
  Total = sum(discard2)
)

DataFrame(
  LibSize = sum(qc.lib3),
  NExprs = sum(qc.nexprs3),
  SpikeProp = sum(qc.spike3),
  MitoProp = sum(qc.mito3),
  Total = sum(discard3)
)

# Was qc.lib necessary for creating discord?
#Option 1
y=seq(1,length(qc.lib))
y=y[(qc.lib=="TRUE")&(qc.spike=="FALSE")&(qc.mito=="FALSE")&(qc.nexprs=="FALSE")]
#Option 2
which(qc.lib & (!qc.spike & !qc.mito & !qc.nexprs))
#Since 2 cells where discarded due to qc.lib and not to the other options, it was important

#Which filter discarded more cells? discard or discard2?
length(which(discard=="TRUE"))
length(which(discard2=="TRUE"))
#Since the true values in discard a more than in discard2, the first filter discarded more cells.

#By considering the sample batch, did we discard more cells using automatic threshold detection?
length(which(discard3=="TRUE"))

###Exercises with EmptyDrops
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

set.seed(100)
e.out <- emptyDrops(counts(sce.pbmc))

set.seed(100)
limit <- 100
all.out <-
  emptyDrops(counts(sce.pbmc), lower = limit, test.ambient = TRUE)

bcrank <- barcodeRanks(counts(sce.pbmc))

#Why does emptyDrops() return NA values?
#Using help 
?emptyDrops
##NA values are assigned to empty droplets, which are determined by those libraries below a limit of gene counts or no gene counts.

#Are the p-values the same for e.out and all.out?
identical(all.out$PValue,e.out$PValue)
#FALSE, due to the presence of NA values

#What if you subset to the non-NA entries?
identical(all.out$PValue[which(!is.na(e.out$FDR))],e.out$PValue[which(!is.na(e.out$FDR))])
#TRUE