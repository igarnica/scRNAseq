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