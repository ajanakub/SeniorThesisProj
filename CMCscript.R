## "CMCscript.R" Bukola Ajanaku (September 22, 2021)
## This is my 399 project script.
## GOAL: Investigat whetger glial co-accessible chromatin regions associated
## with Schizophrenia are spatially located with neuronal chromatin rather
## than being isolated glial chromatin.
## Linux screen: screen 4, "runCMCScript.R"
## ml R/4.0.2 geos/3.5.0 udunits/2.2.26 gdal/2.4.1 proj/6.3.1 pandoc
## R

.libPaths(c("/sc/arion/projects/psychgen/scratch/temp_from_kiran/RLib",.libPaths()))

suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(BiocParallel))
suppressPackageStartupMessages(library(gtools))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(qvalue))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(decorate))

options(xtable.type="html")

knitr::opts_chunk$set(
  echo=TRUE,
  warning=FALSE,
  message=TRUE,
  error = FALSE,
  tidy = FALSE,
  cache = TRUE,
  cache.lazy = FALSE,
  dev = c("png", "pdf"),
  fig.width=7, fig.height=7)

options(markdown.HTML.stylesheet = 'css/custom.css')

register(SnowParam(4, "SOCK", progressbar=TRUE))

# loading data and naming needed values
  load("/sc/arion/projects/epigenAD/Bukola/CMC_ATAC/count_matrix_glia_imputed_olig_residualized_biology_kept_PFC.Rdata")
    gliaCountsMatrix = counts_matrix
    gliaMetadata = METADATA
    gliaMetadataPeaks = metadata_peaks

  load("/sc/arion/projects/epigenAD/Bukola/CMC_ATAC/count_matrix_neuron_imputed_glu_residualized_biology_kept_PFC.Rdata")
    neuronCountsMatrix = counts_matrix
    neuronMetadata = METADATA
    neuronMetadataPeaks = metadata_peaks

# getting rid of these variables so things don't get confusing
    rm(counts_matrix)
    rm(METADATA)
    rm(metadata_peaks)

# merging both metadata by individual id (finding shared ids)
  collapser = merge(neuronMetadata, gliaMetadata, by = "Individual_ID")

# Subset cell type matrices by their IDs
  neuronCountsMatrix = neuronCountsMatrix[, unlist(collapser["ID.x"])]
  neuronMetadata = neuronMetadata[unlist(collapser["ID.x"]),]
  gliaCountsMatrix = gliaCountsMatrix[, unlist(collapser["ID.y"])]
  gliaMetadata = gliaMetadata[unlist(collapser["ID.y"]),]

# Making metadata peaks into grange objects
  nMETAPEAKS = as(neuronMetadataPeaks, "GRanges")
  gMETAPEAKS = as(gliaMetadataPeaks, "GRanges")

# rbinding glia and neuronal cell types to create mother matrix
  combinedDF = rbind(neuronCountsMatrix, gliaCountsMatrix)

## Disease as factor
  gliaMetadata$Dx2 = factor(gliaMetadata$Dx2, c("Control", "Case"))
  neuronMetadata$Dx2 = factor(neuronMetadata$Dx2, c("Control", "Case"))

## Notes:
## 1) "Learn local correlation structure" from tutorial (call CRDS using runclusteredgenome function)
##          residValues is just combinedDF, 651851 from dim of combinedDF
## 2) Create function that permutes the peaks in one chromosome (one main and 10 permutations)
## 3) use

# grabbing function from Gabe Hoffman's package
runOrderedClusteringGenome(
  X,
  gr,
  method = c("adjclust", "hclustgeo"),
  quiet = FALSE,
  alpha = 0.5,
  adjacentCount = 500,
  setNANtoZero = FALSE,
  method.corr = c("pearson", "spearman")
)

# Consolidated code from Kiran to make calling CRDs easier.
call_clusters_collapse=function(input_mat,peaksGR,meanCluster_Size,permutation){
  if (permutation > 0) {
    all_flt_peaks_split=split(peaksGR,seqnames(peaksGR))
    peaks_shuffle=as.character(unlist(sapply(lapply(all_flt_peaks_split,names),sample)))
    input_data_shuffle=input_data[peaks_shuffle,]
    rownames(input_data_shuffle)=peaks_shuffle
    colnames(input_data_shuffle)=colnames(input_mat)
    input_mat=input_data_shuffle
    }

    else {
    input_mat=input_mat
	 }

	  treeList = runOrderedClusteringGenome(input_mat,peaksGR)
    names(peaksGR)=peaksGR$names
    treeListClusters = createClusters(treeList, method = "meanClusterSize", meanClusterSize=meanCluster_Size)
    clstScore =scoreClusters(treeList, treeListClusters)
    return(list(treeList,treeListClusters,clstScore))
}

# Original calling of CRDs. No shuffling (so permuatation 0 of course).
# Must save as R.data

  call_clusters_collapse(combinedDF, nMETAPEAKS, c(10, 25, 50, 100), 0)
  save(c(treeList,treeListClusters,clstScore), file = "/sc/arion/projects/epigenAD/Bukola/CMC_ATAC/permutations/Original.RDATA")

# now doing the for loop of permutations (shuffled treelists)
for(val %in% 1:10){
  peaks_shuffle=as.character(unlist(sapply(lapply(all_flt_peaks_split,names),sample)))

}
