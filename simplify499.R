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

# save(list = ls(), file = "tester.RDATA")
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
      df = as.data.frame(treeListClusters_collapse[[1]][[1]])
      df$ID = rownames(df)
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
# save(list = ls(all.names = TRUE), file =  "/sc/arion/projects/epigenAD/Bukola/newVarSimp499.RDATA")

## Create grange object combining peakLocs and clstScore ##################################################
  for(i in 1:length(names(newCRD))) {
    clusterset = newCRD[[i]]
    names(clusterset)[names(clusterset) == "uniqueID"] = "forComb"

    clstScore[[i]]$forComb = paste0(clstScore[[i]]$id, "_", clstScore[[i]]$chrom, "_", clstScore[[i]]$cluster)

    firstMerge = merge(clusterset, clstScore[[i]], by = "forComb")

#    peakLocs2_df


#    mergePeakLocs = peakLocs2[peakLocs2$ID %in% firstMerge$ID,]
#    secondMerge = merge(firstMerge, mergePeakLocs, by = "ID")
#    colnames(secondMerge) = gsub(".x","",colnames(secondMerge))
#    secondMerge = dplyr::select(secondMerge, -c("seqnames.y", "chrom", "cluster.y", "start.y",
#      "end.y", "width.y", "strand.y"))

    firstMerge$chrom = NULL
    firstMerge$cluster.x = NULL

    assign(paste0("finalObject","_",sub( "_.*$", "", clusterset$forComb[1])), as(firstMerge, "GRanges"))
  }

# Just getting distance values for each cluster in each clustering scheme. CONTINUE
for(i in 1:length(ls(patt="finalObject_"))){
  finalObject = allFinalObjects[[i]]
  clusterobj = unlist(range(split(finalObject,~forComb)))
}

# RUNNING GLASSO ON CLUSTER:
input_mat = quantResid2 # our peaks matrix
input_mat = peakLocs2

allFinalObjects = do.call(list, lapply(ls(patt="finalObject_"), get))

names(allFinalObjects) = c(sub( "_.*$", "", allFinalObjects[[1]]$forComb)[1],
  sub( "_.*$", "", allFinalObjects[[2]]$forComb)[1], sub( "_.*$", "", allFinalObjects[[3]]$forComb)[1],
  sub( "_.*$", "", allFinalObjects[[4]]$forComb)[1], sub( "_.*$", "", allFinalObjects[[5]]$forComb)[1],
  sub( "_.*$", "", allFinalObjects[[6]]$forComb)[1])

finalObject_10$GL_result = NA
finalObject_25$GL_result = NA
finalObject_50$GL_result = NA
finalObject_75$GL_result = NA
finalObject_80$GL_result = NA
finalObject_100$GL_result = NA

GL_results_allsave = list()

for (m in 1:length(ls(patt="finalObject_"))) {

  finalObject = finalObject_100
    GLs_for_100 = list()

    for(mm in 1:length(finalObject$forComb)) {

      clst = unique(finalObject$forComb)[mm]

      cluster_peaks = finalObject[finalObject$forComb == clst,]$ID
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

        distance_parameter = 0.5 # look into
        s = 0.2 # power law, fast vs slow decay of peak correlations
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

      GL_results$cluster = unique(finalObject$forComb)[mm]

      GLs_for_100[[mm]] = GL_results
      }
#      GL_for_25_df = do.call(rbind, GLs_for_25)
#     GL_for_10_df = do.call(rbind, GLs_for_10)
#     GL_for_50_df = do.call(rbind, GLs_for_50)
#     GL_for_80_df = do.call(rbind, GLs_for_80)
#     GL_for_100_df = do.call(rbind, GLs_for_100)

}

#      uniqueIDdf = m

      #### CHECK MATHHHHH:
#      pairPeaks = dim(GL_results)[1]
#      sigPairPeaks = sum(GL_results$signif == 1)
#      totalPeaks = length(cluster_peaks)

#      holder = list(uniqueIDdf, totalPeaks, pairPeaks, sigPairPeaks)

#      assign(paste0("holder", ".", uniqueIDdf), holder)

#    }
#  }

#      checking10 = GL_results
# GL100_chr11_71 = GL_results
# GL10_chr1_1117 = GL_results
# save(checking10, checking100, file= "/sc/arion/projects/epigenAD/Bukola/newChecks.RDATA")
# save(GL_for_75_df, file= "/sc/arion/projects/epigenAD/Bukola/for75GL.RDATA")


# setwd("/sc/arion/projects/epigenAD/Bukola")
# save(list = ls(), file =
# load("/sc/arion/projects/epigenAD/Bukola/afterHoldr.RDATA")

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

pdf(file = "/sc/arion/projects/epigenAD/Bukola/forrdistplots.pdf")
  plot(mergedClusters10$width, mergedClusters10$signifPairs, main= "10: Signifant Pairs vs Distance", xlab = "Distance", ylab = "Significant Pairs", col = ifelse(clustersIn10$percSig > 0.5,'green','red'))
    legend("topright", legend = c("Significance > 0.5", "Significance < 0.5"), fill = c("green", "red"))
  plot(mergedClusters25$width, mergedClusters25$signifPairs, main= "25: Signifant Pairs vs Distance", xlab = "Distance", ylab = "Significant Pairs", col = ifelse(clustersIn25$percSig > 0.5,'green','red'))
    legend("topright", legend = c("Significance > 0.5", "Significance < 0.5"), fill = c("green", "red"))
  plot(mergedClusters50$width, mergedClusters50$signifPairs, main= "50: Signifant Pairs vs Distance", xlab = "Distance", ylab = "Significant Pairs", col = ifelse(clustersIn50$percSig > 0.5,'green','red'))
    legend("topright", legend = c("Significance > 0.5", "Significance < 0.5"), fill = c("green", "red"))
  plot(mergedClusters75$width, mergedClusters75$signifPairs, main= "75: Signifant Pairs vs Distance", xlab = "Distance", ylab = "Significant Pairs", col = ifelse(clustersIn75$percSig > 0.5,'green','red'))
    legend("topright", legend = c("Significance > 0.5", "Significance < 0.5"), fill = c("green", "red"))
  plot(mergedClusters80$width, mergedClusters80$signifPairs, main= "80: Signifant Pairs vs Distance", xlab = "Distance", ylab = "Significant Pairs", col = ifelse(clustersIn80$percSig > 0.5,'green','red'))
    legend("topright", legend = c("Significance > 0.5", "Significance < 0.5"), fill = c("green", "red"))
  plot(mergedClusters100$width, mergedClusters100$signifPairs, main= "100: Signifant Pairs vs Distance", xlab = "Distance", ylab = "Significant Pairs", col = ifelse(clustersIn100$percSig > 0.5,'green','red'))
    legend("topright", legend = c("Significance > 0.5", "Significance < 0.5"), fill = c("green", "red"))
dev.off()
# scp -r ajanab02@minerva.hpc.mssm.edu:/sc/arion/projects/epigenAD/Bukola/forrdistplots.pdf /Users/bukola/Documents/finalplots


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

################################################################################
# LAST RESULTS:

# 10: percSig = 1
checkClust10 = as(mergedClusters10[mergedClusters10$percSig == 1,], "GRanges")
mainClustersfor10 = dim(mergedClusters10[mergedClusters10$percSig == 1,])[1]
overlayClustersfor10 = GenomicRanges::intersect(checkClust10, HIC_loops)
confirmedMainClustersfor10 = length(ranges(overlayClustersfor10))
percBiologicallValidated10 = confirmedMainClustersfor10 / mainClustersfor10
### of the 3923 clusters of maintained significance, only 1547 could be biologically validated using

# 25: percSig = 1
checkClust25 = as(mergedClusters25[mergedClusters25$percSig == 1,], "GRanges")
mainClustersfor25 = dim(mergedClusters25[mergedClusters25$percSig == 1,])[1]
overlayClustersfor25 = GenomicRanges::intersect(checkClust25, HIC_loops)
confirmedMainClustersfor25 = length(ranges(overlayClustersfor25))
percBiologicallValidated25 = confirmedMainClustersfor25 / mainClustersfor25

# 50: percSig = 1
checkClust50 = as(mergedClusters50[mergedClusters50$percSig == 1,], "GRanges")
mainClustersfor50 = dim(mergedClusters50[mergedClusters50$percSig == 1,])[1]
overlayClustersfor50 = GenomicRanges::intersect(checkClust50, HIC_loops)
confirmedMainClustersfor50 = length(ranges(overlayClustersfor50))
percBiologicallValidated50 = confirmedMainClustersfor50 / mainClustersfor50

# 75: percSig = 1
checkClust75 = as(mergedClusters75[mergedClusters75$percSig == 1,], "GRanges")
mainClustersfor75 = dim(mergedClusters75[mergedClusters75$percSig == 1,])[1]
overlayClustersfor75 = GenomicRanges::intersect(checkClust75, HIC_loops)
confirmedMainClustersfor75 = length(ranges(overlayClustersfor75))
percBiologicallValidated75 = confirmedMainClustersfor75 / mainClustersfor75

# 80: percSig = 1
checkClust80 = as(mergedClusters80[mergedClusters80$percSig == 1,], "GRanges")
mainClustersfor80 = dim(mergedClusters80[mergedClusters80$percSig == 1,])[1]
overlayClustersfor80 = GenomicRanges::intersect(checkClust80, HIC_loops)
confirmedMainClustersfor80 = length(ranges(overlayClustersfor80))
percBiologicallValidated80 = confirmedMainClustersfor80 / mainClustersfor80

# 100: percSig = 1
checkClust100 = as(mergedClusters100[mergedClusters100$percSig == 1,], "GRanges")
mainClustersfor100 = dim(mergedClusters100[mergedClusters100$percSig == 1,])[1]
overlayClustersfor100 = GenomicRanges::intersect(checkClust100, HIC_loops)
confirmedMainClustersfor100 = length(ranges(overlayClustersfor100))
percBiologicallValidated100 = confirmedMainClustersfor100 / mainClustersfor100

################################################################################
Before_glasso = allFinalObjects[["10"]] # before glasso

dropPeaks1 = as.character(highlight_df$Var1)
dropPeaks2 = as.character(highlight_df$Var2)

`%!in%` <- Negate(`%in%`)

After_glasso1 = as.data.frame(Before_glasso)
After_glasso2 = After_glasso1[After_glasso1$ID %!in% dropPeaks1,]
After_glasso = After_glasso2[After_glasso2$ID %!in% dropPeaks2,]

After_glasso = as(After_glasso, "GRanges")

load("/sc/arion/projects/epigenAD/Bukola/loop_CRD_results.RDATA")

HIC_loops = Loop_GR$neuron_loop_anchors

all_peaks = peakLocs2
# input = HIC_loops
# BG = all_peaks

A_vs_B_fisher_test=function(A,B,input,BG){
C1 = input
mcols(C1)=NULL
names(C1)=NULL
C2=BG
mcols(C2)=NULL
names(C2)=NULL
mcols(A)=NULL
names(A)=NULL
mcols(B)=NULL
names(B)=NULL
strand(B)="*"
strand(A)="*"
gr1=GenomicRanges::intersect(C1,A)
gr2=GenomicRanges::intersect(C1,B)
gr3=GenomicRanges::intersect(C2,A)
gr4=GenomicRanges::intersect(C2,B)
c1_in_A_in_B=sum(as.numeric(GenomicRanges::width(gr1)))
c1_not_A_in_B=sum(as.numeric(GenomicRanges::width(gr2)))
c1_in_A_not_B=sum(as.numeric(GenomicRanges::width(gr3)))
c1_not_A_not_B=sum(as.numeric(GenomicRanges::width(gr4)))
matA=matrix(c(c1_in_A_in_B,c1_not_A_in_B,c1_in_A_not_B,c1_not_A_not_B),nrow=2)
colnames(matA)=c("A","not_A")
rownames(matA)=c("B","not_B")
#print(matA[1,2])
#print(matA[2,2])
#if (matA[2,2]==0 | matA[1,1]==0 ){
#compt_score=data.frame("OR"=0,"Pvalue"=1,LL=0,UL=0)
#} else if (matA[1,2]==0 | matA[2,1]==0) {
#compt_score=data.frame("OR"="Inf","Pvalue"=.01,LL=0,UL=0)
#print(compt_score)
#} else {
dt=ifelse(nchar(min(matA))>3,10^(nchar(min(matA))-4),nchar(min(matA)))
outputA=fisher.test((matA)/dt)
compt_score=data.frame("OR"=outputA$estimate[1],"Pvalue"=format(outputA$p.value[1],digits=6),"LL"=outputA$conf.int[1],"UL"=outputA$conf.int[2])
#}
#print(compt_score)
compt_score
}

################################################################################
compartment_enrichment=function(dysCRD,BG,compartment_GR){
c1_in_compt_in_diff=sum(as.numeric(GenomicRanges::width(GenomicRanges::intersect(dysCRD,compartment_GR))))
c1_not_compt_in_diff=sum(as.numeric(GenomicRanges::width(GenomicRanges::setdiff(dysCRD,compartment_GR))))
non_dys_CRD=GenomicRanges::setdiff(BG,dysCRD)
c1_in_compt_not_diff=sum(as.numeric(GenomicRanges::width(GenomicRanges::intersect(non_dys_CRD,compartment_GR))))
c1_not_compt_not_diff=sum(as.numeric(GenomicRanges::width(GenomicRanges::setdiff(non_dys_CRD,compartment_GR))))
matA=matrix(c(c1_in_compt_in_diff,c1_not_compt_in_diff,c1_in_compt_not_diff,c1_not_compt_not_diff),nrow=2)
colnames(matA)=c("diff","not_diff")
rownames(matA)=c("compt","not_compt")
#if (matA[1,2]==0 | matA[1,1]==0| matA[2,1]==0 | matA[2,2]==0  ){
#compt_score=data.frame("OR"=0,"Pvalue"=1,LL=0,UL=0)
#} else {
#dt=ifelse(nchar(min(matA))>3,10^(nchar(min(matA))-4),nchar(min(matA)))
#outputA=oddsratio((matA)/dt)
#compt_score=data.frame("OR"=outputA$measure[2],"Pvalue"=format(outputA$p.value[2],digits=6),"LL"=outputA$measure[2,2],"UL"=outputA$measure[2,3])
#}
dt=ifelse(nchar(min(matA))>3,10^(nchar(min(matA))-4),nchar(min(matA)))
outputA=fisher.test((matA)/dt)
compt_score=data.frame("OR"=outputA$estimate[1],"Pvalue"=format(outputA$p.value[1],digits=6),"LL"=outputA$conf.int[1],"UL"=outputA$conf.int[2])
compt_score
}

################################## >0.5 vs 0.5>
# 10
Gr1in10 = as(mergedClusters10[mergedClusters10$percSig > 0.5,], "GRanges")
fromGr1in10 = compartment_enrichment(Gr1in10, peakLocs2, HIC_loops)
Gr2in10 = as(mergedClusters10[mergedClusters10$percSig <= 0.5,], "GRanges")
fromGr2in10 = compartment_enrichment(Gr2in10, peakLocs2, HIC_loops)

   # for KIRAN
   Gr1in10 = as(mergedClusters10[mergedClusters10$percSig > 0.75,], "GRanges")
   Gr2in10 = as(mergedClusters10[mergedClusters10$percSig <= 0.10,], "GRanges")
   fromGr1in10 = compartment_enrichment(Gr1in10, peakLocs2, HIC_loops)
   fromGr2in10 = compartment_enrichment(Gr2in10, peakLocs2, HIC_loops)

# for percSig = 1 of cluster sizes 50
checkClust50 = as(mergedClusters50[mergedClusters50$percSig == 1,], "GRanges")
forClust50compartment = compartment_enrichment(checkClust50, peakLocs2, HIC_loops)

  # for percSig of 0.5
  Gr1in50 = as(mergedClusters50[mergedClusters50$percSig > 0.50,], "GRanges")
  Gr2in50 = as(mergedClusters50[mergedClusters50$percSig <= 0.50,], "GRanges")
  fromGr1in50 = compartment_enrichment(Gr1in50, peakLocs2, HIC_loops)
  fromGr2in50 = compartment_enrichment(Gr2in50, peakLocs2, HIC_loops)
  fromGr1in50
  fromGr2in50

# for percSig = 1 of cluster sizes 100
checkClust100 = as(mergedClusters100[mergedClusters100$percSig == 1,], "GRanges")
forClust100compartment = compartment_enrichment(checkClust100, peakLocs2, HIC_loops)

# save(mergedClusters100, mergedClusters10, mergedClusters25, mergedClusters50, mergedClusters75, mergedClusters80, file= "/sc/arion/projects/epigenAD/Bukola/mergedClustersforKiran.RDATA")


######## Results:
# 1
compartment_enrichment(finalObject_75,peakLocs2,HIC_loops)

# 2 Crd Comparison

# 3 Kiran's Version of CRD Comparison

# save(GL_results, file= "/sc/arion/projects/epigenAD/Bukola/origGL_results.RDATA")

# save(list = ls(), file = "/sc/arion/projects/epigenAD/Bukola/fixedVariabless.RDATA")

# save(GL_for_100_df, GL_for_80_df, file = "/sc/arion/projects/epigenAD/Bukola/GLsfor100and80.RDATA")

# save(GL_for_10_df, file = "/sc/arion/projects/epigenAD/Bukola/GLsfor10.RDATA")

######## confirming HIC_loops
.libPaths(c("/sc/arion/projects/psychgen/scratch/temp_from_kiran/RLib",.libPaths()))
source("/sc/arion/work/girdhk01/scripts/chipseq_files.R")
#source(‘/sc/arion/projects/roussp01a/pengfei/hicchip/scripts/hic_helper.R’)
#source(“/sc/arion/work/girdhk01/scripts/myscripts/CRD/CMC_SV_help.R”)
source("/sc/arion/work/girdhk01/scripts/myscripts/CRD/TRH_help_functions.R")
source("/sc/arion/work/girdhk01/scripts/myscripts/CRD/HiC_pengfei.R")

source("/sc/arion/work/girdhk01/scripts/myscripts/CRD/TRH_help_functions.R")
# "/sc/arion/projects/epigenAD/Bukola/GLs25and50.RDATA")
# "/sc/arion/projects/epigenAD/Bukola/GLsfor100and80.RDATA")
