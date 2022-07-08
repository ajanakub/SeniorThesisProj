## Bukola Ajanaku
## February 11, 2022
## for499.R : This is my code for using regression methods on the chipseq data
## that was provided by Kiran. Need information regarding data type again.
## on screen 1 as for499
## only run if decorate cannot be loaded: ml R/4.0.2 geos/3.5.0 udunits/2.2.26 gdal/2.4.1 proj/6.3.1 pandoc
## R
## based on my other code: finalizing.R
## screen 4 has the limited clstScore (only 10, 25, 80, 100)

.libPaths(c("/sc/arion/projects/psychgen/scratch/temp_from_kiran/RLib",.libPaths()))

suppressPackageStartupMessages(library(decorate))
suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(BiocParallel))
suppressPackageStartupMessages(library(glasso))


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

load("Processed.RDATA")

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

    # Now, for clustering the correlated ATAC sequences
treeList = runOrderedClusteringGenome( quantResid2, peakLocs2)
treeListClusters = createClusters( treeList, method='meanClusterSize', meanClusterSize=c(10, 25, 50, 75, 80, 100))
forc80 = countClusters( treeListClusters )
clstScore = scoreClusters(treeList, treeListClusters, BPPARAM=SerialParam() )

    # Now, dropping all the N = 1s
originalClstScore = clstScore
    # clearing out all the clusters with only 1 feature
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

newCRD = list()

    # ADDING THE METHYLATION NAMES PER CLST SIZE (FOR BOTH DISEASE AND CONTROL: ALL)
for (m in 1:length(treeListClusters_collapse)){
      newCRD_temp = list()
  for (mm in 1:length(treeListClusters_collapse[[m]])){
      df = as.data.frame(treeListClusters_collapse[[m]][[mm]])
      df$ID=rownames(df)
      colnames(df) <- c("cluster", "ID")
      peakLocs2_df = as.data.frame(peakLocs2)
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

# save(list = ls(all.names = TRUE), file =  "originalVariables.RDATA") from screen 4
# save(list = ls(all.names = TRUE), file =  "newVariables.RDATA") messed up

input_mat = quantResid2 #our peaks matrix
cluster_peaks = newCRD[[1]][newCRD[[1]]$uniqueID == "10_chr1_5", "ID"]
clusterMat = input_mat[rownames(input_mat) %in% cluster_peaks,]


# Glasso: But starting with Cicro example
peakIDS = rownames(clusterMat) # peakID of the specific cluster
newCRD_subset = peakLocs2[peakLocs2$ID %in% peakIDS] # GRanges object like peaklocs 2 input_CRD
clusterMat_subset = clusterMat[rownames(clusterMat) %in% peakIDS,]
​
cov_mat <- cov(t(clusterMat_subset))
diag(cov_mat) <- diag(cov_mat) + 1e-4
​
distance_df = matrix(0,ncol=length(peakIDS),nrow=length(peakIDS))
rownames(distance_df) = peakIDS
colnames(distance_df) = peakIDS
distance_df_melt = melt(distance_df)
dfDist_distance = distance(newCRD_subset[match(distance_df_melt$Var1,newCRD_subset$ID)], newCRD_subset[match(distance_df_melt$Var2,newCRD_subset$ID)])
distance_df_melt$distance = dfDist_distance
distance_df_melt = distance_df_melt[,c(1,2,4)]
distance_df = xtabs(distance~.,distance_df_melt)

distance_df = xtabs(distance~., data = distance_df_melt)​

dist_matrix = distance_df
​
distance_parameter = 0.5
s = 0.1
xmin <- 10e3
out <- (1-(xmin/dist_matrix)^s) * distance_parameter
  # temp=melt(out)
  # plot(temp$value)
out[!is.finite(out)] <- 0
out[out < 0] <- 0
out
​
rho_mat = out
GL <- glasso::glasso(cov_mat, rho_mat)
# heatmap(GL$wi)
rownames(GL$wi) = rownames(cov_mat)
colnames(GL$wi) = rownames(cov_mat)
GL_data = melt(GL$wi)
# GL_data$value[GL_data$value!=0]=1

highlight_df = GL_data[GL_data$value == 0.000000,]

library(dplyr)
library(ggplot2)

## yes
ggplot(GL_data, aes(Var1, Var2, fill= value)) +
    geom_tile()+ scale_fill_gradient(low="white", high="green") +
    theme_classic() + geom_tile(data = highlight_df, fill = "red")

plot(sort(distance_df_melt$distance),type= "b")

distance_df_melt_temp = distance_df_melt
distance_df_melt_temp$ID = paste0(distance_df_melt_temp$Var1,"_",distance_df_melt_temp$Var2)
cor_mat = cor(t(clusterMat_subset))
cor_mat_melt = melt(cor_mat)
cor_mat_melt$ID = paste0(cor_mat_melt$Var1,"_",cor_mat_melt$Var2)
cor_mat_melt_merge = merge(cor_mat_melt,distance_df_melt_temp,by="ID")
ggplot(cor_mat_melt_merge,aes(distance,value))+geom_point()+xlab("distance")+ylab("correlation")+theme_classic()

GL_data = melt(GL$wi)
GL_data$ID = paste0(GL_data$Var1, "_", GL_data$Var2)
GL_results = merge(GL_data, cor_mat_melt_merge, by = "ID")
GL_results$signif = 0
GL_results$signif[GL_results$value.x != 0] = 1

ggplot(GL_results, aes(distance,value.y, color = factor(signif)))+geom_point()+xlab("Distance (bp)")+ylab("Correlation")+theme_classic() +
  ggtitle("Comparing Cor & Signif over Distance ") + labs(color = "Significance")

  highlight_df = GL_results[GL_results$signif == 0,]

  (ggplot(GL_results, aes(Var1, Var2, fill= value.y)) +
      geom_tile()+ scale_fill_gradient(low="white", high="green") +
      theme_classic() + geom_tile(data = highlight_df, fill = "red")
)
      ggplot(GL_results, aes(Var1, Var2, fill= value.y)) +
          geom_tile()+ scale_fill_gradient(low="white", high="green") +
          theme_classic() + geom_tile(data = highlight_df, fill = "red")
---

distance_df_melt_temp = distance_df_melt
distance_df_melt_temp$ID = paste0(distance_df_melt_temp$Var1,"_",distance_df_melt_temp$Var2)
cor_mat=cor(t(res_mat_subset))
cor_mat_melt=melt(cor_mat)
cor_mat_melt$ID=paste0(cor_mat_melt$Var1,"_",cor_mat_melt$Var2)
cor_mat_melt_merge=merge(cor_mat_melt,distance_df_melt_temp,by="ID")

ggplot(cor_mat_melt_merge,aes(distance,value))+geom_point()+xlab("distance")+ylab("correlation")+theme_classic()

ggplot(GL_data, aes(Var1, Var2, fill= value)) +
    geom_tile()+ scale_fill_gradient(low="white", high="green") +
    theme_classic() + geom_tile(data = highlight_df, fill = "red")
