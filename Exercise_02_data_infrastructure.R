library('scRNAseq')
sce.416b <- LunSpikeInData(which = "416b")

# Load the SingleCellExperiment package
library('SingleCellExperiment')
# Extract the count matrix from the 416b dataset
counts.416b <- counts(sce.416b)
# Construct a new SCE from the counts matrix
sce <- SingleCellExperiment(assays = list(counts = counts.416b))

# Inspect the object we just created
sce

## How big is it?
pryr::object_size(sce)

# Access the counts matrix from the assays slot
# WARNING: This will flood RStudio with output!

# 1. The general method
assay(sce, "counts")[1:6, 1:3]
# 2. The special method for the assay named "counts"
counts(sce)[1:6, 1:3]

sce <- scater::logNormCounts(sce)
# Inspect the object we just updated
sce

## How big is it?
pryr::object_size(sce)

# 1. The general method
assay(sce, "logcounts")[1:6, 1:3]
# 2. The special method for the assay named "logcounts"
logcounts(sce)[1:6, 1:3]

# assign a new entry to assays slot
assay(sce, "counts_100") <- assay(sce, "counts") + 100
# List the assays in the object
assays(sce)
assayNames(sce)

## How big is it?
pryr::object_size(sce)

# Extract the sample metadata from the 416b dataset
colData.416b <- colData(sce.416b)
# Add some of the sample metadata to our SCE
colData(sce) <- colData.416b[, c("phenotype", "block")]
# Inspect the object we just updated
sce
# Access the sample metadata from our SCE
colData(sce)
# Access a specific column of sample metadata from our SCE
table(sce$block)

# Example of function that adds extra fields to colData
sce <- scater::addPerCellQC(sce.416b)
# Access the sample metadata from our updated SCE
colData(sce)

# Inspect the object we just updated
sce

## How big is it?
pryr::object_size(sce)

## Add the lognorm counts again
sce <- scater::logNormCounts(sce)

## How big is it?
pryr::object_size(sce)

# E.g., subset data to just wild type cells
# Remember, cells are columns of the SCE
x <- sce$phenotype == "wild type phenotype"
class(x)
table(x)
table(sce$phenotype)

sce[, sce$phenotype == "wild type phenotype"]

# Access the feature metadata from our SCE
# It's currently empty!
rowData(sce)

# Example of function that adds extra fields to rowData
sce <- scater::addPerFeatureQC(sce)
# Access the feature metadata from our updated SCE
rowData(sce)

## How big is it?
pryr::object_size(sce)


# Download the relevant Ensembl annotation database
# using AnnotationHub resources
library('AnnotationHub')
ah <- AnnotationHub()
query(ah, c("Mus musculus", "Ensembl", "v97"))

# Annotate each gene with its chromosome location
#download_data <- function(ah, id) { access_data(sh, id) }
ensdb <- ah[["AH73905"]]
chromosome <- mapIds(ensdb,
                     keys = rownames(sce),
                     keytype = "GENEID",
                     column = "SEQNAME")

class(chromosome)
rowData(sce)$chromosome <- chromosome

# Access the feature metadata from our updated SCE
rowData(sce)

## How big is it?
pryr::object_size(sce)

# E.g., subset data to just genes on chromosome 3
# NOTE: which() needed to cope with NA chromosome names
x <- rowData(sce)$chromosome
table(is.na(x))
y <- x == '3'
table(is.na(y))
table(is.na(which(y)))
head(which(y))
sce[which(rowData(sce)$chromosome == "3"), ]

# Access the metadata from our SCE
# It's currently empty!
metadata(sce)

# The metadata slot is Vegas - anything goes
metadata(sce) <- list(favourite_genes = c("Shh", "Nck1", "Diablo"),
                      analyst = c("Pete"))

# Access the metadata from our updated SCE
metadata(sce)

# E.g., add the PCA of logcounts
# NOTE: We'll learn more about PCA later
sce <- scater::runPCA(sce)
# Inspect the object we just updated
sce
# Access the PCA matrix from the reducedDims slot
reducedDim(sce, "PCA")[1:6, 1:3]

# E.g., add a t-SNE representation of logcounts
# NOTE: We'll learn more about t-SNE later
sce <- scater::runTSNE(sce)
# Inspect the object we just updated
sce
# Access the t-SNE matrix from the reducedDims slot
head(reducedDim(sce, "TSNE"))

# E.g., add a 'manual' UMAP representation of logcounts
# NOTE: We'll learn more about UMAP later and a
# 		  simpler way to compute it.
u <- uwot::umap(t(logcounts(sce)), n_components = 2)
# Add the UMAP matrix to the reducedDims slot
# Access the UMAP matrix from the reducedDims slot
reducedDim(sce, "UMAP") <- u

# List the dimensionality reduction results stored in # the object
reducedDims(sce)

# Extract the ERCC SCE from the 416b dataset
ercc.sce.416b <- altExp(sce.416b, "ERCC")
# Inspect the ERCC SCE
ercc.sce.416b

# Add the ERCC SCE as an alternative experiment to our SCE
altExp(sce, "ERCC") <- ercc.sce.416b
# Inspect the object we just updated
sce

## How big is it?
pryr::object_size(sce)

# List the alternative experiments stored in the object
altExps(sce)

# Subsetting the SCE by sample also subsets the
# alternative experiments
sce.subset <- sce[, 1:10]
ncol(sce.subset)
ncol(altExp(sce.subset))

## How big is it?
pryr::object_size(sce.subset)

# Extract existing size factors (these were added
# when we ran scater::logNormCounts(sce))
head(sizeFactors(sce))

# 'Automatically' replace size factors
sce <- scran::computeSumFactors(sce)
head(sizeFactors(sce))

# 'Manually' replace size factors
sizeFactors(sce) <- scater::librarySizeFactors(sce)
head(sizeFactors(sce))


##Exercise with sce objetc
#Which function defines the sce class?
#SingleCellExperiment::SingleCellExperiment
#SingleCellExperiment funnction (after ::), from SingleCellExperiment library

#What are the minimum type of tables an sce object contains?
#'Assays' for primary data containing the counts for each gene (row) in each cell (column)
#'rowData' for gene metadata, e.g. ID, chromosome, symbol, mean expresion, ...
#'colData' for cell metadata, e.g. Barcode, donor, treatment, ...
  
#Where are the colnames(sce) used?
colnames(sce)
#To identify cells when searching in the Assays table (columns)

#Similarly, where are the rownames(sce) used?
rownames(sce)
#To identify genes when searching in the Assays table (rows)

#How many principal components did we compute?
dim(reducedDim(sce))
#50

#Which three chromosomes have the highest mean gene expression?
sort(tapply(rowData(sce)$mean,rowData(sce)$chromosome,mean), decreasing=T)
# With 'tapply' we will apply the function mean (last argument) to subsets of the rowData(sce)$mean vector 
# (containing information about mean expression for each gene), in a way that these subsets are definided 
# by the information in the rowData(sce)$chromosome (number of chromosome). Sort allows to order the subsets
# mean in a decreasing level of mean gene expression, so that it is easier to identify the chromosomes with
# higher mean gene expression.

##Exercises with ERCC
## Read the data from the web
ercc_info <-
  read.delim(
    'https://tools.thermofisher.com/content/sfs/manuals/cms_095046.txt',
    as.is = TRUE,
    row.names = 2,
    check.names = FALSE
  )
#Use ERCC ID to align this table with sce object (ERCC alt experiment)
y=match(rownames(altExp(sce, 'ERCC')),(rownames(ercc_info)))
table(is.na(y))
ercc_info=ercc_info[y,]
identical(rownames(altExp(sce, "ERCC")), rownames(ercc_info))
#We downloaded from the internet a table containing information about multiple ERCC genes, but we're only
#interested in the info from the subset of them that were used in our alt experiment. In order to do this, 
#we'll need the row numbers from the internet table which contain the info from the genes used in the alt experiment.
#We can use the function match, which will return the row numbers from the internet table where the ERCC genes 
# we know we used in the alt experiment coincide with the ERCC genes in the general table. After that, we check
# that every ERCC gene in the alt experiment matched with at least one gene in the general internet table by 
#seeing if there was a NA, but there wasn't. Then we can save the information from the rows in the general
#internet table we were interested, and check with identical if they efectively correspond to the ERCC genes used
#in the alt experiment.

#Use plot() to plot concentration in Mix 1 (attomoles/ul) vs ERCC counts of the sce object (in alt exp) 
pdf("/home/igarnica/igarnica/counts_vs_mix1_conc.pdf")
for(i in 1:dim(ercc_info[2])){
  plot(log(ercc_info[,3]),log(assay(altExp(sce,'ERCC'),'counts')[,i]), main=paste("",rownames(altExp(sce, "ERCC"))[i],""),xlab="log concentration", ylab="log counts")    
}
dev.off()
#We saved the multiple plots containing the information of the concentration in Mix 1 (attomoles/ul) 
# for the ERCC gene counts in each cell of the alt experiment (saved in the columns of the Assay table), 
#in a pdf document.