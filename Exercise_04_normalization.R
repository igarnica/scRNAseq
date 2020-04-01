#Data
library('scRNAseq')
sce.zeisel <- ZeiselBrainData(ensembl = TRUE)
#Data quality
library('scater')
is.mito <- which(rowData(sce.zeisel)$featureType == "mito")
stats <- perCellQCMetrics(sce.zeisel, subsets = list(Mt = is.mito))
qc <-
  quickPerCellQC(stats,
                 percent_subsets = c("altexps_ERCC_percent", "subsets_Mt_percent"))
sce.zeisel <- sce.zeisel[, !qc$discard]

#Calculate library size and library size factor
lib.sf.zeisel <- librarySizeFactors(sce.zeisel)
ls.zeisel <- colSums(counts(sce.zeisel))

##Exercises about library size and library size factor
#Are ls.zeisel and lib.sf.zeisel identical?
identical(lib.sf.zeisel,ls.zeisel)
#False

#Are they proportional?
unique(lib.sf.zeisel/ls.zeisel)
#We can check if they are proportional by dividing each corresponding cell in both vectors
#to make sure they differ by a proportionality constant. Two responses are given by the command
#but they're pretty similar, so they are probably due to rounding errors during division.

identical(lib.sf.zeisel,ls.zeisel/mean(ls.zeisel))
lib.sf.zeisel=ls.zeisel/mean(ls.zeisel)
#Since we expect the mean library size factor, representing the normalized expression value, to be 1
# mean(k*ls.zeisel)=1  => k*mean(ls.zeisel)=1 => k= 1/mean(ls.zeisel)
#From this last derivation, we can see that the expected proportionality constant is 1/mean(ls.zeisel), so that
#the normalized expression value is 1. It is sufficient to multiply by this proportionality constant in each
#cell of the vector containing the counts for all genes in a cell. By running the identical command in line 28,
#we can see that our function gives the same result that calculating the library size factor the tradional way.

##Exercises with quickCluster

library('scran')
# Pre-clustering
set.seed(100)
# Every time we run a simulation, where the answer is non-deterministic, it is important to set the seed number 
# to a fixed number, in order to always use the same random numbers and get the same answer (even though the 
# answer can vary using other random numbers) for reproducibility purposes.

#How many quick clusters did we get?
clust.zeisel <- quickCluster(sce.zeisel)
summary(clust.zeisel)

# 12 clusters.

#How many cells per quick cluster did we get?
# 243,281,252,324,224,236,440,325,231,259
  
#How many quick clusters will we get if we set the minimum size to 200? Use 100 as the seed.
set.seed(100)
clust.zeisel_2 <- quickCluster(sce.zeisel, min.size=200)
summary(clust.zeisel_2)

#10 clusters.

deconv.sf.zeisel <-
  calculateSumFactors(sce.zeisel, clusters = clust.zeisel, min.mean = 0.1)
summary(deconv.sf.zeisel)
plot(
  ls.zeisel,
  deconv.sf.zeisel,
  log = "xy",
  xlab = "Library size",
  ylab = "Size factor"
)


#How many lines do you see?
#From the plot we can see aproximately 2 lines, which could be due to multiple factors like samples 
#coming from different individuals or groups of individuals.
