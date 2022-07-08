## Bukola Ajanaku
## March 1, 2022
## simplify499.R : Taking information from for499.R document, cleaning, and
## simplifying the code.
## Screen 2 as 408431.simplify499
## R

options(warn=-1)
.libPaths(c("/sc/arion/projects/epigenAD/Bukola/RLib",.libPaths()))

suppressPackageStartupMessages(library(decorate))
suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(BiocParallel))
suppressPackageStartupMessages(library(glasso))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))


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

load("Processed.RDATA") # getting genetic information from temp_from_kiran

# START FROM DECORATE ##########################################################

# filtering metaData for only Scz and Control
keep = metaData$Dx %in% c('SCZ', 'Control')
metaData = metaData[keep,]
chipCounts = chipCounts[,keep]
metaData$Dx = factor(metaData$Dx, c("Control", "SCZ"))

## Processing data to find data that is being expressed: *Not sure of line*
isexpr = rowSums(cpm(chipCounts)>1) >= 0.2*ncol(chipCounts)
peakLocs2 = peakLocs[which(isexpr)]

## Functional use of limma/voom *Not sure of line*
countObj = DGEList( chipCounts[isexpr,] )
countObj = calcNormFactors( countObj )
design = model.matrix( ~ as.character(`ATACSeq_report:Sequencing_Order_ID`) +
`ATACSeq_report:Mean_GC_content`+
`ATACSeq_report:Mapped_Reads` +
`Age of Death` +
`PMI (in hours)` + Sex , metaData)
vobj = voom( countObj, design, plot= FALSE)

## Other variables used to residualize data.
dcmp = svd(vobj$E, nv=5, nu=0)
frac = dcmp$d^2 / sum(dcmp$d^2) * 100
xlab = paste0('PC1: ', round(frac[1], 1), '%')
ylab = paste0('PC2: ', round(frac[2], 1), '%')

## Residualizing the data.
dsgn = model.matrix( ~ dcmp$v[,1:2] + as.character(`ATACSeq_report:Sequencing_Order_ID`) +
`ATACSeq_report:Mean_GC_content`+
`ATACSeq_report:Mapped_Reads` +
`Age of Death` +
`PMI (in hours)` + Sex , metaData)
fitPC = lmFit(vobj, dsgn)
quantResid = residuals( fitPC, vobj )

vobj2 = voom( countObj, dsgn, plot=FALSE)

fitPC2 = lmFit(vobj2, dsgn)
quantResid2 = residuals( fitPC2, vobj2 )

dsgn = model.matrix( ~ dcmp$v[,1:2] + as.character(`ATACSeq_report:Sequencing_Order_ID`) +
`ATACSeq_report:Mean_GC_content`+
`ATACSeq_report:Mapped_Reads` +
`Age of Death` +
`PMI (in hours)` + Sex + Dx, metaData)

fitDE = lmFit(vobj2, dsgn)

fitDE = eBayes(fitDE)

# save(list = ls(), file = "/sc/arion/projects/epigenAD/Bukola/copyCoder.RDATA")
# END FROM DECORATE ############################################################

# In this example, I will be using a cluster size of 25: -
treeList = runOrderedClusteringGenome( quantResid2, peakLocs2)
treeListClusters = createClusters( treeList, method='meanClusterSize',
  meanClusterSize=c(10, 25, 50, 75, 80, 100)
n_clusters = countClusters( treeListClusters )
clstScore = scoreClusters(treeList, treeListClusters, BPPARAM=SerialParam() )

  # Now, dropping all the N = 1s
originalClstScore = clstScore
cleanclst10 = clstScore[["10"]][clstScore[["10"]]$N > 2, ]
cleanclst25 = clstScore[["25"]][clstScore[["25"]]$N > 2, ]
cleanclst50 = clstScore[["50"]][clstScore[["50"]]$N > 2, ]
cleanclst75 = clstScore[["75"]][clstScore[["75"]]$N > 2, ]
cleanclst80 = clstScore[["80"]][clstScore[["80"]]$N > 2, ]
cleanclst100 = clstScore[["100"]][clstScore[["100"]]$N > 2, ]

clstScore = list(cleanclst10, cleanclst25, cleanclst50, cleanclst75, cleanclst80,
  cleanclst100)
names(clstScore) = c("10", "25", "50", "75", "80", "100")

    # Retaining clusters based on strength and getting these retained clusters. Then collapsing similar clusters.
clustInclude = retainClusters(clstScore, "LEF", 0.1 )
treeListClusters_filter = filterClusters( treeListClusters, clustInclude )
treeListClusters_collapse = collapseClusters( treeListClusters_filter, peakLocs2 )

# FROM KIRAN ###################################################################
newCRD = list()

    # ADDING THE METHYLATION NAMES PER CLST SIZE (FOR BOTH DISEASE AND CONTROL: ALL)
for (m in 1:length(treeListClusters_collapse)){
      newCRD_temp = list()
  for (mm in 1:length(treeListClusters_collapse[[m]])){
      df = as.data.frame(treeListClusters_collapse[[m]][[mm]])
      df$ID=rownames(df)
      colnames(df) <- c("cluster", "ID")
      peakLocs2_df = as.data.frame(peakLocs2)
      peakLocs2_df$ID = rownames(peakLocs2_df)
      peakLocs2_cluster = merge(peakLocs2_df, df, by = "ID")
      peakLocs2_cluster$uniqueID = paste0(names(treeListClusters_collapse)[m],"_",
        names(treeListClusters_collapse[[m]])[mm],"_",peakLocs2_cluster$cluster) # new
      newCRD_temp[[mm]] = peakLocs2_cluster
}

  names(newCRD_temp) = names(treeListClusters_collapse[[m]])
  newCRD_temp_df = do.call(rbind, newCRD_temp)
  newCRD[[m]] = newCRD_temp_df
}

names(newCRD) = names(treeListClusters_collapse)

# Saving variables so far:
# save(list = ls(all.names = TRUE), file =  "newVarSimp499.RDATA")

## Create grange object combining peakLocs and clstScore ##################################################
  for(i in 1:length(names(newCRD))) {
    clusterset = newCRD[[i]]
    names(clusterset)[names(clusterset) == "uniqueID"] = "forComb"

    clstScore[[i]]$forComb = paste0(clstScore[[i]]$id, "_", clstScore[[i]]$chrom, "_", clstScore[[i]]$cluster)

    firstMerge = merge(clusterset, clstScore[[i]], by = "forComb")

    mergePeakLocs = peakLocs2[peakLocs2$ID %in% firstMerge$ID,]
    secondMerge = merge(firstMerge, mergePeakLocs, by = "ID")

    colnames(secondMerge) = gsub(".x","",colnames(secondMerge))

    secondMerge = dplyr::select(secondMerge, -c("seqnames.y", "chrom", "cluster.y", "start.y",
      "end.y", "width.y", "strand.y"))

    assign(paste0("finalObject","_",sub( "_.*$", "", clusterset$forComb[1])), as(secondMerge, "GRanges"))
  }

# Just getting distance values for each cluster in each clustering scheme. CONTINUE
for(i in 1:length(ls(patt="finalObject_")){
  finalObject = allFinalObjects[[p]]
  clusterobj = unlist(range(split(finalObject_25,~forComb)))
}

# RUNNING GLASSO ON CLUSTER:
  counter = seq(0.1, 1.0, length.out = 10)
  ct = 1
  finalResults = list()

for (i in counter){

    distance_parameter = i

    for (i in counter){
        s = i

    m = "80_chr7_195"

    cluster_peaks = finalObject[finalObject$forComb == m ,]$ID
    clusterMat = input_mat[rownames(input_mat) %in% cluster_peaks,]

    peakIDS = rownames(clusterMat) # peakID of the specific cluster
    peakLocs2$ID = names(peakLocs2)
    newCRD_subset = peakLocs2[peakLocs2$ID %in% peakIDS]
    clusterMat_subset = clusterMat[rownames(clusterMat) %in% peakIDS,]

    cov_mat = cov(t(clusterMat_subset))
    diag(cov_mat) = diag(cov_mat) + 1e-4

    distance_df = matrix(0, ncol=length(peakIDS), nrow = length(peakIDS))
    rownames(distance_df) = peakIDS
    colnames(distance_df) = peakIDS
    distance_df_melt = melt(distance_df)
    dfDist_distance = distance(newCRD_subset[match(distance_df_melt$Var1,
      newCRD_subset$ID)], newCRD_subset[match(distance_df_melt$Var2, newCRD_subset$ID)])
    distance_df_melt$distance = dfDist_distance
    distance_df_melt = distance_df_melt[, c(1,2,4)]
    distance_df = xtabs(distance~., distance_df_melt)
    dist_matrix = distance_df

    xmin <- 10e3 # minimum pairwise peak distance

    out <- (1-(xmin/dist_matrix)^s) * distance_parameter

    out[!is.finite(out)] <- 0
    out[out < 0] <- 0

    GL <- glasso::glasso(cov_mat, out)
    rownames(GL$wi) = rownames(cov_mat)
    colnames(GL$wi) = rownames(cov_mat)
    orginal_GL = melt(GL$wi)

    GL$wi[lower.tri(GL$wi)] <- NA
    GL_data = melt(GL$wi)
    GL_data=GL_data[!is.na(GL_data$value),]

    distance_df_melt_temp = distance_df_melt
    distance_df_melt_temp$ID = paste0(distance_df_melt_temp$Var1,"_",distance_df_melt_temp$Var2)
    cor_mat = cor(t(clusterMat_subset))
    cor_mat_melt = melt(cor_mat)
    cor_mat_melt$ID = paste0(cor_mat_melt$Var1,"_",cor_mat_melt$Var2)
    cor_mat_melt_merge = merge(cor_mat_melt, distance_df_melt_temp, by="ID")


    GL_data$ID = paste0(GL_data$Var1, "_", GL_data$Var2)
    GL_results = merge(GL_data, cor_mat_melt_merge, by = "ID")
    GL_results = GL_results[GL_results$Var1.x != GL_results$Var2.x,]
    GL_results$signif = 0
    GL_results$signif[GL_results$value.x != 0] = 1
    GL_results$sval = s
    GL_results$distval = distance_parameter

    finalResults[[ct]] = GL_results
    ct = ct + 1
    }
}

finalResults1=do.call(rbind,finalResults)

################################################################################
#### Plotting correlation matrix: (all other s in screen 4)

highlight_df = GL_results[GL_results$signif == 0,]

pdf(file = "/sc/arion/projects/epigenAD/Bukola/s0.1corrmat.pdf")
ggplot(cor_mat_melt_merge, aes(distance,value))+geom_point()+xlab("distance")+ylab("correlation")+theme_classic()
dev.off()

# scp -r ajanab02@minerva.hpc.mssm.edu:/sc/arion/projects/epigenAD/Bukola/s0.1corrmat.pdf /Users/bukola/Documents/plotsforKiran

################################################################################

for(i in 1:length(names(allFinalObjects))) {
    assign(paste0("clustersIn",names(allFinalObjects)[i]),
      do.call(rbind, lapply(ls(patt= paste0("holder.", names(allFinalObjects)[i], "_")), get)))
}

clustersIn10 = as.data.frame(clustersIn10)
clustersIn25 = as.data.frame(clustersIn25)
clustersIn50 = as.data.frame(clustersIn50)
clustersIn75 = as.data.frame(clustersIn75)
clustersIn80 = as.data.frame(clustersIn80)
clustersIn100 = as.data.frame(clustersIn100)
colnames(clustersIn10) = c("clstUniqueID", "numberOfPeaks", "pairsOfPeaks", "signifPairs")
colnames(clustersIn25) = c("clstUniqueID", "numberOfPeaks", "pairsOfPeaks", "signifPairs")
colnames(clustersIn50) = c("clstUniqueID", "numberOfPeaks", "pairsOfPeaks", "signifPairs")
colnames(clustersIn75) = c("clstUniqueID", "numberOfPeaks", "pairsOfPeaks", "signifPairs")
colnames(clustersIn80) = c("clstUniqueID", "numberOfPeaks", "pairsOfPeaks", "signifPairs")
colnames(clustersIn100) = c("clstUniqueID", "numberOfPeaks", "pairsOfPeaks", "signifPairs")
clustersIn10$percSig = as.numeric(clustersIn10$signifPairs) / as.numeric(clustersIn10$pairsOfPeaks)
clustersIn25$percSig = as.numeric(clustersIn25$signifPairs) / as.numeric(clustersIn25$pairsOfPeaks)
clustersIn50$percSig = as.numeric(clustersIn50$signifPairs) / as.numeric(clustersIn50$pairsOfPeaks)
clustersIn75$percSig = as.numeric(clustersIn75$signifPairs) / as.numeric(clustersIn75$pairsOfPeaks)
clustersIn80$percSig = as.numeric(clustersIn80$signifPairs) / as.numeric(clustersIn80$pairsOfPeaks)
clustersIn100$percSig = as.numeric(clustersIn100$signifPairs) / as.numeric(clustersIn100$pairsOfPeaks)

# HW:
# try and plot the above clustersin25 through directory from KIRAN
# Gives start and end of each cluster by splitting against forComb: clusterobj = unlist(range(split(finalObject_25,~forComb)))
# clusterobj$names = names(clusterobj)
# width(clusterobj[1:10])

# PLOT 1: Significant Peak Pairs vs Total Peak Pairs

pdf(file = "/sc/arion/projects/epigenAD/Bukola/sigPairsOverTotalPairs.pdf")
  plot(clustersIn10$pairsOfPeaks, clustersIn10$signifPairs)
  plot(clustersIn25$pairsOfPeaks, clustersIn25$signifPairs)
  plot(clustersIn50$pairsOfPeaks, clustersIn50$signifPairs)
  plot(clustersIn75$pairsOfPeaks, clustersIn75$signifPairs)
  plot(clustersIn80$pairsOfPeaks, clustersIn80$signifPairs)
  plot(clustersIn100$pairsOfPeaks, clustersIn100$signifPairs)
dev.off()

pdf(file = "/sc/arion/projects/epigenAD/Bukola/sigPercOverTotalPairs.pdf")
  plot(clustersIn10$percSig, clustersIn10$pairsofPeaks)
  plot(clustersIn25$percSig, clustersIn25$pairsofPeaks)
  plot(clustersIn50$percSig, clustersIn50$pairsofPeaks)
  plot(clustersIn75$percSig, clustersIn75$pairsofPeaks)
  plot(clustersIn80$percSig, clustersIn80$pairsofPeaks)
  plot(clustersIn100$percSig, clustersIn100$pairsofPeaks)
dev.off()


pdf(file = "/sc/arion/projects/epigenAD/Bukola/hist.pdf")
  hist(as.numeric(clustersIn10[clustersIn10$percSig == 1, "numberOfPeaks"]))
dev.off()

# PLOT 2: Significant Peak Pairs in CRD vs Width of CRD
clusterobj10 = unlist(range(split(finalObject_10,~forComb)))
  clusterobj10$names = names(clusterobj10)
  clusterobj10$width = width(clusterobj10)

clusterobj25 = unlist(range(split(finalObject_25,~forComb)))
  clusterobj25$names = names(clusterobj25)
  clusterobj25$width = width(clusterobj25)

clusterobj50 = unlist(range(split(finalObject_50,~forComb)))
  clusterobj50$names = names(clusterobj50)
  clusterobj50$width = width(clusterobj50)

clusterobj75 = unlist(range(split(finalObject_75,~forComb)))
  clusterobj75$names = names(clusterobj75)
  clusterobj75$width = width(clusterobj75)

clusterobj80 = unlist(range(split(finalObject_80,~forComb)))
  clusterobj80$names = names(clusterobj80)
  clusterobj80$width = width(clusterobj80)

clusterobj100 = unlist(range(split(finalObject_100,~forComb)))
  clusterobj100$names = names(clusterobj100)
  clusterobj100$width = width(clusterobj100)

colnames(clustersIn10)[which(names(clustersIn10) == "clstUniqueID")] = "forComb"
colnames(clustersIn25)[which(names(clustersIn25) == "clstUniqueID")] = "forComb"
colnames(clustersIn50)[which(names(clustersIn50) == "clstUniqueID")] = "forComb"
colnames(clustersIn75)[which(names(clustersIn75) == "clstUniqueID")] = "forComb"
colnames(clustersIn80)[which(names(clustersIn80) == "clstUniqueID")] = "forComb"
colnames(clustersIn100)[which(names(clustersIn100) == "clstUniqueID")] = "forComb"

dfclusterobj10 = as.data.frame(clusterobj10)
dfclusterobj25 = as.data.frame(clusterobj25)
dfclusterobj50 = as.data.frame(clusterobj50)
dfclusterobj75 = as.data.frame(clusterobj75)
dfclusterobj80 = as.data.frame(clusterobj80)
dfclusterobj100 = as.data.frame(clusterobj100)

colnames(dfclusterobj10)[which(names(dfclusterobj10) == "names")] = "forComb"
colnames(dfclusterobj25)[which(names(dfclusterobj25) == "names")] = "forComb"
colnames(dfclusterobj50)[which(names(dfclusterobj50) == "names")] = "forComb"
colnames(dfclusterobj75)[which(names(dfclusterobj75) == "names")] = "forComb"
colnames(dfclusterobj80)[which(names(dfclusterobj80) == "names")] = "forComb"
colnames(dfclusterobj100)[which(names(dfclusterobj100) == "names")] = "forComb"

# creating the merged clusters for distplots
  clustersIn10 = as.data.frame(clustersIn10)
  dfclusterobj10 = as.data.frame(dfclusterobj10)
  rownames(dfclusterobj10) = NULL
  dfclusterobj10$width.1 = NULL
  clustersIn10$forComb = as.character(clustersIn10$forComb)
  mergedClusters10 = merge(clustersIn10, dfclusterobj10, by = "forComb")

  clustersIn25 = as.data.frame(clustersIn25)
  dfclusterobj25 = as.data.frame(dfclusterobj25)
  rownames(dfclusterobj25) = NULL
  dfclusterobj25$width.1 = NULL
  clustersIn25$forComb = as.character(clustersIn25$forComb)
  mergedClusters25 = merge(clustersIn25, dfclusterobj25, by = "forComb")

  clustersIn50 = as.data.frame(clustersIn50)
  dfclusterobj50 = as.data.frame(dfclusterobj50)
  rownames(dfclusterobj50) = NULL
  dfclusterobj50$width.1 = NULL
  clustersIn50$forComb = as.character(clustersIn50$forComb)
  mergedClusters50 = merge(clustersIn50, dfclusterobj50, by = "forComb")

  clustersIn75 = as.data.frame(clustersIn75)
  dfclusterobj75 = as.data.frame(dfclusterobj75)
  rownames(dfclusterobj75) = NULL
  dfclusterobj75$width.1 = NULL
  clustersIn75$forComb = as.character(clustersIn75$forComb)
  mergedClusters75 = merge(clustersIn75, dfclusterobj75, by = "forComb")

  clustersIn80 = as.data.frame(clustersIn80)
  dfclusterobj80 = as.data.frame(dfclusterobj80)
  rownames(dfclusterobj80) = NULL
  dfclusterobj80$width.1 = NULL
  clustersIn80$forComb = as.character(clustersIn80$forComb)
  mergedClusters80 = merge(clustersIn80, dfclusterobj80, by = "forComb")


  clustersIn100 = as.data.frame(clustersIn100)
  dfclusterobj100 = as.data.frame(dfclusterobj100)
  rownames(dfclusterobj100) = NULL
  dfclusterobj100$width.1 = NULL
  clustersIn100$forComb = as.character(clustersIn100$forComb)
  mergedClusters100 = merge(clustersIn100, dfclusterobj100, by = "forComb")

pdf(file = "/sc/arion/projects/epigenAD/Bukola/alldistplots.pdf")
  plot(mergedClusters10$width, mergedClusters10$signifPairs)
  plot(mergedClusters25$width, mergedClusters25$signifPairs)
  plot(mergedClusters50$width, mergedClusters50$signifPairs)
  plot(mergedClusters75$width, mergedClusters75$signifPairs)
  plot(mergedClusters80$width, mergedClusters80$signifPairs)
  plot(mergedClusters100$width, mergedClusters100$signifPairs)
dev.off()

pdf(file = "/sc/arion/projects/epigenAD/Bukola/callHiswhenPerc1.pdf")
  median10 = median(as.numeric(clustersIn10[clustersIn10$percSig == 1, "numberOfPeaks"]))
  median25 = median(as.numeric(clustersIn25[clustersIn25$percSig == 1, "numberOfPeaks"]))
  median50 = median(as.numeric(clustersIn50[clustersIn50$percSig == 1, "numberOfPeaks"]))
  median75 = median(as.numeric(clustersIn75[clustersIn75$percSig == 1, "numberOfPeaks"]))
  median80 = median(as.numeric(clustersIn80[clustersIn80$percSig > 0.9, "numberOfPeaks"]))
  median100 = median(as.numeric(clustersIn100[clustersIn100$percSig == 1, "numberOfPeaks"]))

  hist(as.numeric(clustersIn10[clustersIn10$percSig == 1, "numberOfPeaks"]))
      abline(v = as.numeric(median10), col = 'blue', lwd = 3)
      text(x = as.numeric(median10) + 2 , y = 2000, paste("Median =", median10), col = "blue", cex = 2)
  hist(as.numeric(clustersIn25[clustersIn25$percSig == 1, "numberOfPeaks"]))
      abline(v = as.numeric(median25), col = 'blue', lwd = 3)
      text(x = as.numeric(median25) + 1, y = 20, paste("Median =", median25), col = "blue", cex = 2)
  hist(as.numeric(clustersIn50[clustersIn50$percSig == 1, "numberOfPeaks"]))
      abline(v = as.numeric(median50), col = 'blue', lwd = 3)
      text(x = as.numeric(median50) + 1 , y = 6, paste("Median =", median50), col = "blue", cex = 2)
  hist(as.numeric(clustersIn75[clustersIn75$percSig == 1, "numberOfPeaks"]))
      abline(v = as.numeric(median75), col = 'blue', lwd = 3)
      text(x = as.numeric(median75) - 1 , y = 0.8, paste("Median =", median75), col = "blue", cex = 2)
  hist(as.numeric(clustersIn80[clustersIn80$percSig > 0.9, "numberOfPeaks"])) # for 1, 0. for 0.9, 1. for 0.8, 2
      abline(v = as.numeric(median80), col = 'blue', lwd = 3)
      text(x = as.numeric(median80) - 1 , y = 0.8, paste("Median =", median80), col = "blue", cex = 2)
  hist(as.numeric(clustersIn100[clustersIn100$percSig == 1, "numberOfPeaks"]))
      abline(v = as.numeric(median100), col = 'blue', lwd = 3)
      text(x = as.numeric(median100) - 1, y = 0.8, paste("Median =", median100), col = "blue", cex = 2)
dev.off()

pdf(file = "/sc/arion/projects/epigenAD/Bukola/testHis.pdf")
  ggplot(data.frame(clustersIn10[clustersIn10$percSig == 1,]),
    aes(clustersIn10[clustersIn10$percSig == 1, "numberOfPeaks"])) + geom_histogram(bins = 10)
dev.off()

# scp -r ajanab02@minerva.hpc.mssm.edu:/sc/arion/projects/epigenAD/Bukola/alldistplots.pdf /Users/bukola/Documents/plotsforKiran

# treelist plot for sig and large cluster in cluster 10
  mergedClusters10[mergedClusters10$percSig == 1, "forComb"] # "10_chr1_1104"
  query = range(finalObject_10[finalObject_10$forComb == "10_chr1_1104"])

  library(EnsDb.Hsapiens.v86) # is this right???
  ensdb = EnsDb.Hsapiens.v86

pdf(file = "/sc/arion/projects/epigenAD/Bukola/decorateSignifClust10.pdf")
  plotDecorate( ensdb, treeList, treeListClusters_collapse, peakLocs2, query, showGenes = TRUE)
dev.off()

To-Do:
1) make histograms of percentages
2) identify medians per graph (
  make hist of distribution of number of peaks when percentage = 1
3)
