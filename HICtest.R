## HICtest.R - checking simplify499 data with HIC_loops
## (neuronal anchorloops gather from HIC data).
## May 2, 2021 - Bukola Ajanaku
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

load("/sc/arion/projects/epigenAD/Bukola/fixedVariabless.RDATA")

load("/sc/arion/projects/epigenAD/Bukola/GL_results_HIC_loops_75_15546_16840.RDATA")
GL_75_pt1 = Glasso_results_list
load("/sc/arion/projects/epigenAD/Bukola/GL_results_HIC_loops_75_14250_15545.RDATA")
GL_75_pt2 = Glasso_results_list
load("/sc/arion/projects/epigenAD/Bukola/GL_results_HIC_loops_75_22022_23317.RDATA")
GL_75_pt3 = Glasso_results_list
load("/sc/arion/projects/epigenAD/Bukola/GL_results_HIC_loops_75_16841_18135.RDATA")
GL_75_pt4 = Glasso_results_list
load("/sc/arion/projects/epigenAD/Bukola/GL_results_HIC_loops_75_19432_20726.RDATA")
GL_75_pt5 = Glasso_results_list
load("/sc/arion/projects/epigenAD/Bukola/GL_results_HIC_loops_75_20727_22021.RDATA")
GL_75_pt6 = Glasso_results_list
load("/sc/arion/projects/epigenAD/Bukola/GL_results_HIC_loops_75_10364_11658.RDATA")
GL_75_pt7 = Glasso_results_list
load("/sc/arion/projects/epigenAD/Bukola/GL_results_HIC_loops_75_6478_7772.RDATA")
GL_75_pt8 = Glasso_results_list
load("/sc/arion/projects/epigenAD/Bukola/GL_results_HIC_loops_75_33681_34975.RDATA")
GL_75_pt9 = Glasso_results_list
load("/sc/arion/projects/epigenAD/Bukola/GL_results_HIC_loops_75_11659_12954.RDATA")
GL_75_pt10 = Glasso_results_list
load("/sc/arion/projects/epigenAD/Bukola/GL_results_HIC_loops_75_5183_6477.RDATA")
GL_75_pt11 = Glasso_results_list
load("/sc/arion/projects/epigenAD/Bukola/GL_results_HIC_loops_75_25909_27203.RDATA")
GL_75_pt12 = Glasso_results_list
load("/sc/arion/projects/epigenAD/Bukola/GL_results_HIC_loops_75_29795_31089.RDATA")
GL_75_pt13 = Glasso_results_list
load("/sc/arion/projects/epigenAD/Bukola/GL_results_HIC_loops_75_3887_5182.RDATA")
GL_75_pt14 = Glasso_results_list
load("/sc/arion/projects/epigenAD/Bukola/GL_results_HIC_loops_75_12955_14249.RDATA")
GL_75_pt15 = Glasso_results_list
load("/sc/arion/projects/epigenAD/Bukola/GL_results_HIC_loops_75_18136_19431.RDATA")
GL_75_pt16 = Glasso_results_list
load("/sc/arion/projects/epigenAD/Bukola/GL_results_HIC_loops_75_32385_33680.RDATA")
GL_75_pt17 = Glasso_results_list
load("/sc/arion/projects/epigenAD/Bukola/GL_results_HIC_loops_75_9069_10363.RDATA")
GL_75_pt18 = Glasso_results_list
load("/sc/arion/projects/epigenAD/Bukola/GL_results_HIC_loops_75_1296_2591.RDATA")
GL_75_pt19 = Glasso_results_list
load("/sc/arion/projects/epigenAD/Bukola/GL_results_HIC_loops_75_27204_28498.RDATA")
GL_75_pt20 = Glasso_results_list
load("/sc/arion/projects/epigenAD/Bukola/GL_results_HIC_loops_75_1_1295.RDATA")
GL_75_pt21 = Glasso_results_list
load("/sc/arion/projects/epigenAD/Bukola/GL_results_HIC_loops_75_31090_32384.RDATA")
GL_75_pt22 = Glasso_results_list
load("/sc/arion/projects/epigenAD/Bukola/GL_results_HIC_loops_75_28499_29794.RDATA")
GL_75_pt23 = Glasso_results_list
load("/sc/arion/projects/epigenAD/Bukola/GL_results_HIC_loops_75_7773_9068.RDATA")
GL_75_pt24 = Glasso_results_list



GL_results_75 = rbind(GL_75_pt1, GL_75_pt2, GL_75_pt3, GL_75_pt4, GL_75_pt5, GL_75_pt6,
  GL_75_pt7, GL_75_pt8, GL_75_pt9, GL_75_pt10, GL_75_pt11, GL_75_pt12,GL_75_pt13,
  GL_75_pt14, GL_75_pt15, GL_75_pt16, GL_75_pt17, GL_75_pt18, GL_75_pt19,
  GL_75_pt19, GL_75_pt20, GL_75_pt21, GL_75_pt21, GL_75_pt22, GL_75_pt23,
  GL_75_pt24)

save(GL_75_pt1, GL_75_pt2, GL_75_pt3, GL_75_pt4, GL_75_pt5, GL_75_pt6,
    GL_75_pt7, GL_75_pt8, GL_75_pt9, GL_75_pt10, GL_75_pt11, GL_75_pt12,GL_75_pt13,
    GL_75_pt14, GL_75_pt15, GL_75_pt16, GL_75_pt17, GL_75_pt18, GL_75_pt19,
    GL_75_pt19, GL_75_pt20, GL_75_pt21, GL_75_pt21, GL_75_pt22, GL_75_pt23,
    GL_75_pt24, GL_results_75, file = "/sc/arion/projects/epigenAD/Bukola/finished75GL.pdf")

##### Kiran Equations (2) ######################################################
compartment_enrichment = function(dysCRD,BG,compartment_GR){
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

A_vs_B_fisher_test = function(A,B,input,BG){
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
Glasso_results_list = GL_results_75
# Enriched: After GLasso:
table(Glasso_results_list$OR>1 & Glasso_results_list$Pvalue<.05 & Glasso_results_list$OR!= 0 & Glasso_results_list$OR!="Inf")

# Depleted: Before GLasso:
table(Glasso_results_list$OR<1 & Glasso_results_list$Pvalue<.05 & Glasso_results_list$OR!= 0 & Glasso_results_list$OR!="Inf")
################################################################################
# Plot 1:

ce10 = compartment_enrichment(finalObject_10, peakLocs2, HIC_loops)
ce25 = compartment_enrichment(finalObject_25, peakLocs2, HIC_loops)
ce50 = compartment_enrichment(finalObject_50, peakLocs2, HIC_loops)
ce75 = compartment_enrichment(finalObject_75, peakLocs2, HIC_loops)
ce80 = compartment_enrichment(finalObject_80, peakLocs2, HIC_loops)
ce100 = compartment_enrichment(finalObject_100, peakLocs2, HIC_loops)

rownames(ce10) = "10"
rownames(ce25) = "25"
rownames(ce50) = "50"
rownames(ce75) = "75"
rownames(ce80) = "80"
rownames(ce100) = "100"

df_ce = rbind(ce10, ce25, ce50, ce75, ce80, ce100)
df_ce$cluster_size = rownames(df_ce)
df_ce$cluster_size <- factor(df_ce$cluster_size, levels = c("10", "25", "50", "75", "80", "100"))


zz = ggplot(df_ce, aes(OR, cluster_size)) + geom_point() + ylab("Clustering Size") +
  xlab("Odds Ratio") + theme_classic(base_size = 12)+ geom_errorbarh(aes(xmin = LL, xmax = UL),size=1) +
  geom_vline(xintercept = 1, color = "green", size = 0.7, linetype = "dashed") + ggtitle("Odds Ratio by Clustering Scheme") +
  theme(axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 13),
  axis.title.x = element_text(size = 15), axis.title.y = element_text(size= 15), plot.title=element_text(hjust=0.5))


pdf("/sc/arion/projects/epigenAD/Bukola/save_ce.pdf")
  zz
dev.off()
# scp -r ajanab02@minerva.hpc.mssm.edu:/sc/arion/projects/epigenAD/Bukola/save_ce.pdf /Users/bukola/Documents/finalplots

# Plot 2:
load("/sc/arion/projects/epigenAD/Bukola/GLsfor75.RDATA")
GL_75 = GL_df
dropper = as.character(unique(GL_75[GL_75$signif == 0, "Var1"]))
dropper2 = as.character(unique(GL_75[GL_75$signif == 0, "Var2"]))
finalObject_75_df = as.data.frame(finalObject_75)
`%!in%` <- Negate(`%in%`)
enrichedFinObj75df = finalObject_75_df[finalObject_75_df$ID %!in% dropper,]
enrichedFinObj75df = enrichedFinObj75df[enrichedFinObj75df$ID %!in% dropper2,]
enrichedFinObj75 = as(enrichedFinObj75df, "GRanges")
compartment_enrichment(enrichedFinObj75, peakLocs2, HIC_loops)

#enrclusterobj75 = unlist(range(split(enrichedFinObj75,~forComb)))
#  enrclusterobj75$names = names(enrclusterobj75)
#  enrclusterobj75$width = width(enrclusterobj75)
#clusterobj75 = enrclusterobj75

clusterobj75 = unlist(range(split(finalObject_75,~forComb)))
  clusterobj75$names = names(clusterobj75)
  clusterobj75$width = width(clusterobj75)

colnames(clustersIn75)[which(names(clustersIn75) == "clstUniqueID")] = "forComb"
  dfclusterobj75 = as.data.frame(clusterobj75)
  colnames(dfclusterobj75)[which(names(dfclusterobj75) == "names")] = "forComb"

  clustersIn75 = as.data.frame(clustersIn75)
  dfclusterobj75 = as.data.frame(dfclusterobj75)
  rownames(dfclusterobj75) = NULL
  dfclusterobj75$width.1 = NULL
  clustersIn75$forComb = as.character(clustersIn75$forComb)
  mergedClusters75 = merge(clustersIn75, dfclusterobj75, by = "forComb")

pdf(file = "/sc/arion/projects/epigenAD/Bukola/origDist75.pdf")
  plot(mergedClusters75$width, mergedClusters75$signifPairs)
dev.off()

# scp -r ajanab02@minerva.hpc.mssm.edu:/sc/arion/projects/epigenAD/Bukola/origDist75.pdf /Users/bukola/Documents/finalplots


pair_GR=list()
HIC_pair=list()
ct = 1
GL_results = GL_75[GL_75$cluster == “75_chr17_23”]

for (i in (1:nrow(GL_results))){
pair_GR[[i]]=c(peakLocs2[peakLocs2$ID %in% GL_results$Var1[i]],peakLocs2[peakLocs2$ID %in% GL_results$Var2[i]])
pair_GR[[i]]$anchors=paste0("loop_",i)
#pair_GR[[i]]$HIC_loops_id=HIC_loops$loop_id[subjectHits(findOverlaps(pair_GR[[i]],HIC_loops))]
pair_GR[[i]]$GL_status=GL_results$signif[i]
len=uniqueN(subjectHits(findOverlaps(pair_GR[[i]],HIC_loops)))
if( len > 0 ){
HIC_pair[[ct]]=HIC_loops[unique(subjectHits(findOverlaps(pair_GR[[i]],HIC_loops)))]
HIC_pair[[ct]]$GL_status=GL_results$signif[i]
ct = ct + 1
}
}
pair_GR_list=lapply(pair_GR,as.data.frame)
pair_GR_list=do.call(rbind,pair_GR_list)
pair_GR_list_notinGL=as(pair_GR_list[pair_GR_list$GL_status==0,],"GRanges")
pair_GR_list_inGL=as(pair_GR_list[pair_GR_list$GL_status==1,],"GRanges")
HIC_GR=HIC_loops[HIC_loops$loop_id %in% unique(HIC_loops$loop_id[subjectHits(findOverlaps(pair_GR_list_inGL,HIC_loops))])]
HIC_GR_list=split(HIC_GR,~loop_id)
GL_results_list=split(pair_GR_list_inGL,~anchors)
df_results=list()
ct = 1
for (i in (1:length(GL_results_list))){
for (ii in (1:length(HIC_GR_list))){
gr1=split(pair_GR_list_inGL,~anchors)[[i]]
gr2=split(pair_GR_list_notinGL,~anchors)[[i]]
df=A_vs_B_fisher_test(gr1,gr2,HIC_GR_list[[ii]],peakLocs2)
df$loop_id=names(HIC_GR_list)[ii]
df$GL_loop_id=names(GL_results_list)[i]
df_results[[ct]]=df
ct = ct + 1
}
}
Glasso_results_list=do.call(rbind,df_results)
Glasso_results_list$Pvalue=as.numeric(as.character(Glasso_results_list$Pvalue))

finalObject_10$width = as.numeric(width(finalObject_10))
finalObject_25$width = as.numeric(width(finalObject_25))
finalObject_50$width = as.numeric(width(finalObject_50))
finalObject_75$width = as.numeric(width(finalObject_75))
finalObject_80$width = as.numeric(width(finalObject_80))
finalObject_100$width = as.numeric(width(finalObject_100))

dfFinalObject_10 = as.data.frame(finalObject_10)
dfFinalObject_25 = as.data.frame(finalObject_25)
dfFinalObject_50 = as.data.frame(finalObject_50)
dfFinalObject_75 = as.data.frame(finalObject_75)
dfFinalObject_80 = as.data.frame(finalObject_80)
dfFinalObject_100 = as.data.frame(finalObject_100)

dfFinalObject_10$width = as.numeric(dfFinalObject_10$width)
dfFinalObject_25$width = as.numeric(dfFinalObject_25$width)
dfFinalObject_50$width = as.numeric(dfFinalObject_50$width)
dfFinalObject_75$width = as.numeric(dfFinalObject_75$width)
dfFinalObject_80$width = as.numeric(dfFinalObject_80$width)
dfFinalObject_100$width = as.numeric(dfFinalObject_100$width)

library(gridExtra)

pdf(file = "/sc/arion/projects/epigenAD/Bukola/forrdistplots.pdf")

  plot10 = ggplot(dfFinalObject_10, aes(width, mean_abs_corr)) + geom_point(colour = "lightgray") +
  ggtitle("10: Correlation vs Distance") + theme_classic() +
  labs(x = "Distance (bp)", y = "Peak Correlation") + geom_hline(yintercept=0.75, color = 'green', size = 0.7, linetype = "dashed") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

  plot25 = ggplot(dfFinalObject_25, aes(width, mean_abs_corr)) + geom_point(colour = "lightgray") +
  ggtitle("25: Correlation vs Distance") + theme_classic() +
  labs(x = "Distance (bp)", y = "Peak Correlation") + geom_hline(yintercept=0.75, color = 'green', size = 0.7, linetype = "dashed") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

  plot50 = ggplot(dfFinalObject_50, aes(width, mean_abs_corr)) + geom_point(colour = "lightgray") +
  ggtitle("50: Correlation vs Distance") + theme_classic() +
  labs(x = "Distance (bp)", y = "Peak Correlation") + geom_hline(yintercept=0.75, color = 'green', size = 0.7, linetype = "dashed") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

  plot75 = ggplot(dfFinalObject_75, aes(width, mean_abs_corr)) + geom_point(colour = "lightgray") +
  ggtitle("75: Correlation vs Distance") + theme_classic() +
  labs(x = "Distance (bp)", y = "Peak Correlation") + geom_hline(yintercept=0.75, color = 'green', size = 0.7, linetype = "dashed") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

  plot80 = ggplot(dfFinalObject_80, aes(width, mean_abs_corr)) + geom_point(colour = "lightgray") +
  ggtitle("80: Correlation vs Distance") + theme_classic() +
  labs(x = "Distance (bp)", y = "Peak Correlation") + geom_hline(yintercept=0.75, color = 'green', size = 0.7, linetype = "dashed") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

  plot100 = ggplot(dfFinalObject_100, aes(width, mean_abs_corr)) + geom_point(colour = "lightgray") +
  ggtitle("100: Correlation vs Distance") + theme_classic() +
  labs(x = "Distance (bp)", y = "Peak Correlation") + geom_hline(yintercept=0.75, color = 'green', size = 0.7, linetype = "dashed") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

  grid.arrange(plot10, plot25, plot50, plot75, plot80, plot100, nrow = 3, ncol = 2)

  dev.off()

scp -r ajanab02@minerva.hpc.mssm.edu:/sc/arion/projects/epigenAD/Bukola/forrdistplots.pdf /Users/bukola/Documents/finalplots


# Plot 3: histograms:

merging75$forComb = merging75$cluster

finalObject_75[startsWith(as.character(finalObject_75$cluster), "75_chr1_"), ]

mchr1 = merging75[startsWith(as.character(merging75$cluster), "75_chr1_"), ]
mchr2 = merging75[startsWith(as.character(merging75$cluster), "75_chr2_"), ]
mchr3 = merging75[startsWith(as.character(merging75$cluster), "75_chr3_"), ]
mchr4  = merging75[startsWith(as.character(merging75$cluster), "75_chr4_"), ]
mchr5 = merging75[startsWith(as.character(merging75$cluster), "75_chr5_"), ]
mchr6 = merging75[startsWith(as.character(merging75$cluster), "75_chr6_"), ]
mchr7 = merging75[startsWith(as.character(merging75$cluster), "75_chr7_"), ]
mchr8 = merging75[startsWith(as.character(merging75$cluster), "75_chr8_"), ]
mchr9 = merging75[startsWith(as.character(merging75$cluster), "75_chr9_"), ]
mchr10 = merging75[startsWith(as.character(merging75$cluster), "75_chr10_"), ]
mchr11 = merging75[startsWith(as.character(merging75$cluster), "75_chr11_"), ]
mchr12 = merging75[startsWith(as.character(merging75$cluster), "75_chr12_"), ]
mchr13 = merging75[startsWith(as.character(merging75$cluster), "75_chr13_"), ]
mchr14 = merging75[startsWith(as.character(merging75$cluster), "75_chr14_"), ]
mchr15 = merging75[startsWith(as.character(merging75$cluster), "75_chr15_"), ]
mchr16 = merging75[startsWith(as.character(merging75$cluster), "75_chr16_"), ]
mchr17 = merging75[startsWith(as.character(merging75$cluster), "75_chr17_"), ]
mchr18 = merging75[startsWith(as.character(merging75$cluster), "75_chr18_"), ]
mchr19 = merging75[startsWith(as.character(merging75$cluster), "75_chr19_"), ]
mchr20 = merging75[startsWith(as.character(merging75$cluster), "75_chr20_"), ]
mchr21 = merging75[startsWith(as.character(merging75$cluster), "75_chr21_"), ]
mchr22 = merging75[startsWith(as.character(merging75$cluster), "75_chr22_"), ]
mchrY = merging75[startsWith(as.character(merging75$cluster), "75_chrY_"), ]

mchr1$plotGroup = gsub( "^.*?75_chr1_","", mchr1$plotGroup)
mchr2$plotGroup = gsub( "^.*?75_chr2_","", mchr2$plotGroup)
mchr3$plotGroup = gsub( "^.*?75_chr3_","", mchr3$plotGroup)
mchr4$plotGroup = gsub( "^.*?75_chr4_","", mchr4$plotGroup)
mchr5$plotGroup = gsub( "^.*?75_chr5_","", mchr5$plotGroup)
mchr6$plotGroup = gsub( "^.*?75_chr6_","", mchr6$plotGroup)
mchr7$plotGroup = gsub( "^.*?75_chr7_","", mchr7$plotGroup)
mchr8$plotGroup = gsub( "^.*?75_chr8_","", mchr8$plotGroup)
mchr9$plotGroup = gsub( "^.*?75_chr9_","", mchr9$plotGroup)
mchr10$plotGroup = gsub( "^.*?75_chr10_","", mchr10$plotGroup)
mchr11$plotGroup = gsub( "^.*?75_chr11_","", mchr11$plotGroup)
mchr12$plotGroup = gsub( "^.*?75_chr12_","", mchr12$plotGroup)
mchr13$plotGroup = gsub( "^.*?75_chr13_","", mchr13$plotGroup)
mchr14$plotGroup = gsub( "^.*?75_chr14_","", mchr14$plotGroup)
mchr15$plotGroup = gsub( "^.*?75_chr15_","", mchr15$plotGroup)
mchr16$plotGroup = gsub( "^.*?75_chr16_","", mchr16$plotGroup)
mchr17$plotGroup = gsub( "^.*?75_chr17_","", mchr17$plotGroup)
mchr18$plotGroup = gsub( "^.*?75_chr18_","", mchr18$plotGroup)
mchr19$plotGroup = gsub( "^.*?75_chr19_","", mchr19$plotGroup)
mchr20$plotGroup = gsub( "^.*?75_chr20_","", mchr20$plotGroup)
mchr21$plotGroup = gsub( "^.*?75_chr21_","", mchr21$plotGroup)
mchr22$plotGroup = gsub( "^.*?75_chr22_","", mchr22$plotGroup)
mchrYplotGroup = gsub( "^.*?75_chrY_","", mchrY$plotGroup)

pchrY = ggplot(data= mchrY, aes(x = as.factor(plotGroup))) +
  geom_bar(aes(y = pairsOfPeaks), stat='identity', fill='lightblue', color='lightblue4') +
  geom_bar(aes(y = signifPairs), stat='identity', fill='pink', color='red') +
    ggtitle("Chromosome Y") + theme_classic() + theme(plot.title = element_text(face = "bold.italic", hjust = 0.5)) +
    xlab("Cluster ID") + ylab("Significant Pairs of Peaks")

pchr22 = ggplot(data= mchr22, aes(x = as.factor(plotGroup))) +
  geom_bar(aes(y = pairsOfPeaks), stat='identity', fill='lightblue', color='lightblue4') +
  geom_bar(aes(y = signifPairs), stat='identity', fill='pink', color='red') +
    ggtitle("Chromosome 22") + theme_classic() + theme(plot.title = element_text(face = "bold.italic", hjust = 0.5)) +
    xlab("Cluster ID") + ylab("Significant Pairs of Peaks")

pchr21 = ggplot(data= mchr21, aes(x = as.factor(plotGroup))) +
  geom_bar(aes(y = pairsOfPeaks), stat='identity', fill='lightblue', color='lightblue4') +
  geom_bar(aes(y = signifPairs), stat='identity', fill='pink', color='red') +
    ggtitle("Chromosome 21") + theme_classic() + theme(plot.title = element_text(face = "bold.italic", hjust = 0.5)) +
    xlab("Cluster ID") + ylab("Significant Pairs of Peaks")

pchr20 = ggplot(data= mchr20, aes(x = as.factor(plotGroup))) +
  geom_bar(aes(y = pairsOfPeaks), stat='identity', fill='lightblue', color='lightblue4') +
  geom_bar(aes(y = signifPairs), stat='identity', fill='pink', color='red') +
    ggtitle("Chromosome 20") + theme_classic() + theme(plot.title = element_text(face = "bold.italic", hjust = 0.5)) +
    xlab("Cluster ID") + ylab("Significant Pairs of Peaks")

pchr19 = ggplot(data= mchr19, aes(x = as.factor(plotGroup))) +
  geom_bar(aes(y = pairsOfPeaks), stat='identity', fill='lightblue', color='lightblue4') +
  geom_bar(aes(y = signifPairs), stat='identity', fill='pink', color='red') +
    ggtitle("Chromosome 19") + theme_classic() + theme(plot.title = element_text(face = "bold.italic", hjust = 0.5)) +
    xlab("Cluster ID") + ylab("Significant Pairs of Peaks")

pchr18 = ggplot(data= mchr18, aes(x = as.factor(plotGroup))) +
  geom_bar(aes(y = pairsOfPeaks), stat='identity', fill='lightblue', color='lightblue4') +
  geom_bar(aes(y = signifPairs), stat='identity', fill='pink', color='red') +
    ggtitle("Chromosome 18") + theme_classic() + theme(plot.title = element_text(face = "bold.italic", hjust = 0.5)) +
    xlab("Cluster ID") + ylab("Significant Pairs of Peaks")

pchr17 = ggplot(data= mchr17, aes(x = as.factor(plotGroup))) +
  geom_bar(aes(y = pairsOfPeaks), stat='identity', fill='lightblue', color='lightblue4') +
  geom_bar(aes(y = signifPairs), stat='identity', fill='pink', color='red') +
    ggtitle("Chromosome 17") + theme_classic() + theme(plot.title = element_text(face = "bold.italic", hjust = 0.5)) +
    xlab("Cluster ID") + ylab("Significant Pairs of Peaks")

pchr16 = ggplot(data= mchr16, aes(x = as.factor(plotGroup))) +
  geom_bar(aes(y = pairsOfPeaks), stat='identity', fill='lightblue', color='lightblue4') +
  geom_bar(aes(y = signifPairs), stat='identity', fill='pink', color='red') +
    ggtitle("Chromosome 16") + theme_classic() + theme(plot.title = element_text(face = "bold.italic", hjust = 0.5)) +
    xlab("Cluster ID") + ylab("Significant Pairs of Peaks")

pchr15 = ggplot(data= mchr15, aes(x = as.factor(plotGroup))) +
  geom_bar(aes(y = pairsOfPeaks), stat='identity', fill='lightblue', color='lightblue4') +
  geom_bar(aes(y = signifPairs), stat='identity', fill='pink', color='red') +
    ggtitle("Chromosome 15") + theme_classic() + theme(plot.title = element_text(face = "bold.italic", hjust = 0.5)) +
    xlab("Cluster ID") + ylab("Significant Pairs of Peaks")

pchr14 = ggplot(data= mchr14, aes(x = as.factor(plotGroup))) +
  geom_bar(aes(y = pairsOfPeaks), stat='identity', fill='lightblue', color='lightblue4') +
  geom_bar(aes(y = signifPairs), stat='identity', fill='pink', color='red') +
    ggtitle("Chromosome 14") + theme_classic() + theme(plot.title = element_text(face = "bold.italic", hjust = 0.5)) +
    xlab("Cluster ID") + ylab("Significant Pairs of Peaks")

pchr13 = ggplot(data= mchr13, aes(x = as.factor(plotGroup))) +
  geom_bar(aes(y = pairsOfPeaks), stat='identity', fill='lightblue', color='lightblue4') +
  geom_bar(aes(y = signifPairs), stat='identity', fill='pink', color='red') +
    ggtitle("Chromosome 13") + theme_classic() + theme(plot.title = element_text(face = "bold.italic", hjust = 0.5)) +
    xlab("Cluster ID") + ylab("Significant Pairs of Peaks")

pchr12 = ggplot(data= mchr12, aes(x = as.factor(plotGroup))) +
  geom_bar(aes(y = pairsOfPeaks), stat='identity', fill='lightblue', color='lightblue4') +
  geom_bar(aes(y = signifPairs), stat='identity', fill='pink', color='red') +
    ggtitle("Chromosome 12") + theme_classic() + theme(plot.title = element_text(face = "bold.italic", hjust = 0.5)) +
    xlab("Cluster ID") + ylab("Significant Pairs of Peaks")

pchr11 = ggplot(data= mchr11, aes(x = as.factor(plotGroup))) +
  geom_bar(aes(y = pairsOfPeaks), stat='identity', fill='lightblue', color='lightblue4') +
  geom_bar(aes(y = signifPairs), stat='identity', fill='pink', color='red') +
    ggtitle("Chromosome 11") + theme_classic() + theme(plot.title = element_text(face = "bold.italic", hjust = 0.5)) +
    xlab("Cluster ID") + ylab("Significant Pairs of Peaks")

pchr10 = ggplot(data= mchr10, aes(x = as.factor(plotGroup))) +
  geom_bar(aes(y = pairsOfPeaks), stat='identity', fill='lightblue', color='lightblue4') +
  geom_bar(aes(y = signifPairs), stat='identity', fill='pink', color='red') +
    ggtitle("Chromosome 10") + theme_classic() + theme(plot.title = element_text(face = "bold.italic", hjust = 0.5)) +
    xlab("Cluster ID") + ylab("Significant Pairs of Peaks")

pchr9 = ggplot(data= mchr9, aes(x = as.factor(plotGroup))) +
  geom_bar(aes(y = pairsOfPeaks), stat='identity', fill='lightblue', color='lightblue4') +
  geom_bar(aes(y = signifPairs), stat='identity', fill='pink', color='red') +
    ggtitle("Chromosome 9") + theme_classic() + theme(plot.title = element_text(face = "bold.italic", hjust = 0.5)) +
    xlab("Cluster ID") + ylab("Significant Pairs of Peaks")

pchr8 = ggplot(data= mchr8, aes(x = as.factor(plotGroup))) +
  geom_bar(aes(y = pairsOfPeaks), stat='identity', fill='lightblue', color='lightblue4') +
  geom_bar(aes(y = signifPairs), stat='identity', fill='pink', color='red') +
    ggtitle("Chromosome 8") + theme_classic() + theme(plot.title = element_text(face = "bold.italic", hjust = 0.5)) +
    xlab("Cluster ID") + ylab("Significant Pairs of Peaks")

pchr7 = ggplot(data= mchr7, aes(x = as.factor(plotGroup))) +
  geom_bar(aes(y = pairsOfPeaks), stat='identity', fill='lightblue', color='lightblue4') +
  geom_bar(aes(y = signifPairs), stat='identity', fill='pink', color='red') +
    ggtitle("Chromosome 7") + theme_classic() + theme(plot.title = element_text(face = "bold.italic", hjust = 0.5)) +
    xlab("Cluster ID") + ylab("Significant Pairs of Peaks")

pchr6 = ggplot(data= mchr6, aes(x = as.factor(plotGroup))) +
  geom_bar(aes(y = pairsOfPeaks), stat='identity', fill='lightblue', color='lightblue4') +
  geom_bar(aes(y = signifPairs), stat='identity', fill='pink', color='red') +
    ggtitle("Chromosome 6") + theme_classic() + theme(plot.title = element_text(face = "bold.italic", hjust = 0.5)) +
    xlab("Cluster ID") + ylab("Significant Pairs of Peaks")

pchr5 = ggplot(data= mchr5, aes(x = as.factor(plotGroup))) +
  geom_bar(aes(y = pairsOfPeaks), stat='identity', fill='lightblue', color='lightblue4') +
  geom_bar(aes(y = signifPairs), stat='identity', fill='pink', color='red') +
    ggtitle("Chromosome 5") + theme_classic() + theme(plot.title = element_text(face = "bold.italic", hjust = 0.5)) +
    xlab("Cluster ID") + ylab("Significant Pairs of Peaks")

pchr4 = ggplot(data= mchr4, aes(x = as.factor(plotGroup))) +
  geom_bar(aes(y = pairsOfPeaks), stat='identity', fill='lightblue', color='lightblue4') +
  geom_bar(aes(y = signifPairs), stat='identity', fill='pink', color='red') +
    ggtitle("Chromosome 4") + theme_classic() + theme(plot.title = element_text(face = "bold.italic", hjust = 0.5)) +
    xlab("Cluster ID") + ylab("Significant Pairs of Peaks")

pchr3 = ggplot(data= mchr3, aes(x = as.factor(plotGroup))) +
  geom_bar(aes(y = pairsOfPeaks), stat='identity', fill='lightblue', color='lightblue4') +
  geom_bar(aes(y = signifPairs), stat='identity', fill='pink', color='red') +
    ggtitle("Chromosome 3") + theme_classic() + theme(plot.title = element_text(face = "bold.italic", hjust = 0.5)) +
    xlab("Cluster ID") + ylab("Significant Pairs of Peaks")

pchr2 = ggplot(data= mchr2, aes(x = as.factor(plotGroup))) +
  geom_bar(aes(y = pairsOfPeaks), stat='identity', fill='lightblue', color='lightblue4') +
  geom_bar(aes(y = signifPairs), stat='identity', fill='pink', color='red') +
    ggtitle("Chromosome 2") + theme_classic() + theme(plot.title = element_text(face = "bold.italic", hjust = 0.5)) +
    xlab("Cluster ID") + ylab("Significant Pairs of Peaks")

pchr1 = ggplot(data= mchr1, aes(x = as.factor(plotGroup))) +
  geom_bar(aes(y = pairsOfPeaks), stat='identity', fill='lightblue', color='lightblue4') +
  geom_bar(aes(y = signifPairs), stat='identity', fill='pink', color='red') +
    ggtitle("Chromosome 1") + theme_classic() + theme(plot.title = element_text(face = "bold.italic", hjust = 0.5)) +
    xlab("Cluster ID") + ylab("Significant Pairs of Peaks")

pdf(file = "/sc/arion/projects/epigenAD/Bukola/merg75hist.pdf")

  grid.arrange(pchr1, pchr2, pchr3, pchr4, pchr5, pchr6, nrow = 3, ncol = 2)

  grid.arrange(pchr7, pchr8, pchr9, pchr10, pchr11, pchr12, nrow = 3, ncol = 2)

  grid.arrange(pchr13, pchr14, pchr15, pchr16, pchr17, pchr18, nrow = 3, ncol = 2)

  grid.arrange(pchr19, pchr20, pchr21, pchr22, pchrY, nrow = 3, ncol = 2)

dev.off()


pdf(file = "/sc/arion/projects/epigenAD/Bukola/forKey.pdf")
  pchr1 + theme(legend.position="top")

dev.off()


scp -r ajanab02@minerva.hpc.mssm.edu:/sc/arion/projects/epigenAD/Bukola/forKey.pdf /Users/bukola/Documents/finalplots

##### plot 4: CRD examples "75_chr4_103" and "75_chr21_35"

query = range(finalObject_75[finalObject_75$forComb == "75_chr21_35",])

library(EnsDb.Hsapiens.v86) # is this right???
ensdb = EnsDb.Hsapiens.v86


pdf(file = "/sc/arion/projects/epigenAD/Bukola/checkingTL.pdf")
plotDecorate( ensdb, treeList, treeListClusters_collapse, peakLocs2, query, showGenes = TRUE)
dev.off()

scp -r ajanab02@minerva.hpc.mssm.edu:/sc/arion/projects/epigenAD/Bukola/checkingTL.pdf /Users/bukola/Documents/finalplots

########################################################### Plot 5: Heatmap ####
input_mat = quantResid2 #our peaks matrix
cluster_peaks = finalObject_75[finalObject_75$forComb == "75_chr21_35",]$ID
clusterMat = input_mat[rownames(input_mat) %in% cluster_peaks,]

peakIDS = rownames(clusterMat) # peakID of the specific cluster
peakLocs2$ID = names(peakLocs2)
newCRD_subset = peakLocs2[peakLocs2$ID %in% peakIDS]
clusterMat_subset = clusterMat[rownames(clusterMat) %in% peakIDS,]

cov_mat <- cov(t(clusterMat_subset)) # this is for the number of peaks in a given cluster
diag(cov_mat) <- diag(cov_mat) + 1e-4
​
distance_df = matrix(0,ncol=length(peakIDS),nrow=length(peakIDS))
rownames(distance_df) = peakIDS
colnames(distance_df) = peakIDS
distance_df_melt = melt(distance_df)
dfDist_distance = distance(newCRD_subset[match(distance_df_melt$Var1,
  newCRD_subset$ID)], newCRD_subset[match(distance_df_melt$Var2,newCRD_subset$ID)])
distance_df_melt$distance = dfDist_distance
distance_df_melt = distance_df_melt[,c(1,2,4)]
distance_df = xtabs(distance~.,distance_df_melt)
# distance_df = xtabs(distance~., data = distance_df_melt)​

dist_matrix = distance_df

    # parameters to play around with:
distance_parameter = 0.5
s = 0.2
xmin <- 10e3
out <- (1-(xmin/dist_matrix)^s) * distance_parameter
    # parameters done

out[!is.finite(out)] <- 0
out[out < 0] <- 0
out # out is the rho matrix

GL <- glasso::glasso(cov_mat, out)
rownames(GL$wi) = rownames(cov_mat)
colnames(GL$wi) = rownames(cov_mat)
GL_data = melt(GL$wi)

distance_df_melt_temp = distance_df_melt
distance_df_melt_temp$ID = paste0(distance_df_melt_temp$Var1,"_",distance_df_melt_temp$Var2)
cor_mat = cor(t(clusterMat_subset))
cor_mat_melt = melt(cor_mat)
cor_mat_melt$ID = paste0(cor_mat_melt$Var1,"_",cor_mat_melt$Var2)
cor_mat_melt_merge = merge(cor_mat_melt,distance_df_melt_temp,by="ID")

GL_data$ID = paste0(GL_data$Var1, "_", GL_data$Var2)
GL_results = merge(GL_data, cor_mat_melt_merge, by = "ID")
GL_results$signif = 0
GL_results$signif[GL_results$value.x != 0] = 1


highlight_df = GL_results[GL_results$signif == 0,]

pdf(file = "/sc/arion/projects/epigenAD/Bukola/bcorr1.pdf")
ggplot(GL_results, aes(Var1, Var2, fill= value.y)) +
    geom_tile()+ scale_fill_gradient(low="white", high="green") +
    theme_classic() + theme(axis.text.x = element_text(angle=65, hjust=1)) +
    ggtitle("Before GLASSO: Chrom 4 Cluster 103") + xlab("Peaks in CRD") +
    ylab("Peaks in CRD")

    ggplot(GL_results, aes(Var1, Var2, fill= value.y)) +
        geom_tile()+ scale_fill_gradient(low="white", high="green") +
        theme_classic() + geom_tile(data = highlight_df, fill = "red") +
        theme(axis.text.x = element_text(angle=65, hjust=1)) +
        ggtitle("After GLASSO: Chrom 4 Cluster 103") + xlab("Peaks in CRD") +
        ylab("Peaks in CRD")

dev.off()

scp -r ajanab02@minerva.hpc.mssm.edu:/sc/arion/projects/epigenAD/Bukola/bcorr1.pdf /Users/bukola/Documents/finalplots

pdf(file = "/sc/arion/projects/epigenAD/Bukola/bcorr2.pdf", width = 15, height = 13)
ggplot(GL_results, aes(Var1, Var2, fill= value.y)) +
    geom_tile()+ scale_fill_gradient(low="white", high="green") +
    theme_classic() + ggtitle("Before GLASSO: Chrom 21 Cluster 35") + xlab("Peaks in CRD") +
    ylab("Peaks in CRD") + theme(axis.text.x = element_text(size = 15, angle=65, hjust=1, face = "bold"),
    axis.title.x = element_text(size = 15, face = "bold.italic"), axis.title.y = element_text(size= 15, face = "bold.italic"))

ggplot(GL_results, aes(Var1, Var2, fill= value.y)) +
        geom_tile()+ scale_fill_gradient(low="white", high="green") +
        theme_classic() + geom_tile(data = highlight_df, fill = "red") +
        ggtitle("After GLASSO: Chrom 21 Cluster 35") + xlab("Peaks in CRD") +
        ylab("Peaks in CRD") + theme(axis.text.x = element_text(size = 15, angle=65, hjust=1, face = "bold"),
        axis.title.x = element_text(size = 15, face = "bold.italic"), axis.title.y = element_text(size= 15, face = "bold.italic"))

dev.off()

scp -r ajanab02@minerva.hpc.mssm.edu:/sc/arion/projects/epigenAD/Bukola/bcorr2.pdf /Users/bukola/Documents/finalplots


#################################################################### PLOT 6: GL


mchr1$sigPairs = as.list(mchr1$sigPairs)
mchr2$sigPairs = as.list(mchr2$sigPairs)
mchr3$sigPairs = as.list(mchr3$sigPairs)
mchr4$sigPairs = as.list(mchr4$sigPairs)
mchr5$sigPairs = as.list(mchr5$sigPairs)
mchr6$sigPairs = as.list(mchr6$sigPairs)
mchr7$sigPairs = as.list(mchr7$sigPairs)
mchr8$sigPairs = as.list(mchr8$sigPairs)
mchr9$sigPairs = as.list(mchr9$sigPairs)
mchr10$sigPairs = as.list(mchr10$sigPairs)
mchr11$sigPairs = as.list(mchr11$sigPairs)
mchr12$sigPairs = as.list(mchr12$sigPairs)
mchr13$sigPairs = as.list(mchr13$sigPairs)
mchr14$sigPairs = as.list(mchr14$sigPairs)
mchr15$sigPairs = as.list(mchr15$sigPairs)
mchr16$sigPairs = as.list(mchr16$sigPairs)
mchr17$sigPairs = as.list(mchr17$sigPairs)
mchr18$sigPairs = as.list(mchr18$sigPairs)
mchr19$sigPairs = as.list(mchr19$sigPairs)
mchr20$sigPairs = as.list(mchr20$sigPairs)
mchr21$sigPairs = as.list(mchr21$sigPairs)
mchr22$sigPairs = as.list(mchr22$sigPairs)
mchrY$sigPairs = as.list(mchrY$sigPairs)

glY = ggplot(data= mchrY, aes(x = as.factor(plotGroup))) +
  geom_bar(aes(y = pairsOfPeaks), stat='identity', fill='lightblue', color='lightblue4') +
  geom_bar(aes(y = signifPairs), stat='identity', fill= 'pink', color='red') +
  geom_bar(aes(y = sigPairs), stat='identity', fill='mediumorchid1', color='mediumorchid4') +
    ggtitle("Chromosome Y") + theme_classic() + theme(plot.title = element_text(face = "bold.italic", hjust = 0.5)) +
    xlab("Cluster ID") + ylab("Significant Pairs of Peaks")

gl22 = ggplot(data= mchr22, aes(x = as.factor(plotGroup))) +
  geom_bar(aes(y = pairsOfPeaks), stat='identity', fill='lightblue', color='lightblue4') +
  geom_bar(aes(y = signifPairs), stat='identity', fill= 'pink', color='red') +
  geom_bar(aes(y = sigPairs), stat='identity', fill='mediumorchid1', color='mediumorchid4') +
    ggtitle("Chromosome 22") + theme_classic() + theme(plot.title = element_text(face = "bold.italic", hjust = 0.5)) +
    xlab("Cluster ID") + ylab("Significant Pairs of Peaks")

gl21 = ggplot(data= mchr21, aes(x = as.factor(plotGroup))) +
  geom_bar(aes(y = pairsOfPeaks), stat='identity', fill='lightblue', color='lightblue4') +
  geom_bar(aes(y = signifPairs), stat='identity', fill= 'pink', color='red') +
  geom_bar(aes(y = sigPairs), stat='identity', fill='mediumorchid1', color='mediumorchid4') +
    ggtitle("Chromosome 21") + theme_classic() + theme(plot.title = element_text(face = "bold.italic", hjust = 0.5)) +
    xlab("Cluster ID") + ylab("Significant Pairs of Peaks")

gl20 = ggplot(data= mchr20, aes(x = as.factor(plotGroup))) +
  geom_bar(aes(y = pairsOfPeaks), stat='identity', fill='lightblue', color='lightblue4') +
  geom_bar(aes(y = signifPairs), stat='identity', fill= 'pink', color='red') +
  geom_bar(aes(y = sigPairs), stat='identity', fill='mediumorchid1', color='mediumorchid4') +
    ggtitle("Chromosome 20") + theme_classic() + theme(plot.title = element_text(face = "bold.italic", hjust = 0.5)) +
    xlab("Cluster ID") + ylab("Significant Pairs of Peaks")

gl19 = ggplot(data= mchr19, aes(x = as.factor(plotGroup))) +
  geom_bar(aes(y = pairsOfPeaks), stat='identity', fill='lightblue', color='lightblue4') +
  geom_bar(aes(y = signifPairs), stat='identity', fill= 'pink', color='red') +
  geom_bar(aes(y = sigPairs), stat='identity', fill='mediumorchid1', color='mediumorchid4') +
    ggtitle("Chromosome 19") + theme_classic() + theme(plot.title = element_text(face = "bold.italic", hjust = 0.5)) +
    xlab("Cluster ID") + ylab("Significant Pairs of Peaks")

gl18 = ggplot(data= mchr18, aes(x = as.factor(plotGroup))) +
  geom_bar(aes(y = pairsOfPeaks), stat='identity', fill='lightblue', color='lightblue4') +
  geom_bar(aes(y = signifPairs), stat='identity', fill= 'pink', color='red') +
  geom_bar(aes(y = sigPairs), stat='identity', fill='mediumorchid1', color='mediumorchid4') +
    ggtitle("Chromosome 18") + theme_classic() + theme(plot.title = element_text(face = "bold.italic", hjust = 0.5)) +
    xlab("Cluster ID") + ylab("Significant Pairs of Peaks")

gl17 = ggplot(data= mchr17, aes(x = as.factor(plotGroup))) +
  geom_bar(aes(y = pairsOfPeaks), stat='identity', fill='lightblue', color='lightblue4') +
  geom_bar(aes(y = signifPairs), stat='identity', fill= 'pink', color='red') +
  geom_bar(aes(y = sigPairs), stat='identity', fill='mediumorchid1', color='mediumorchid4') +
    ggtitle("Chromosome 17") + theme_classic() + theme(plot.title = element_text(face = "bold.italic", hjust = 0.5)) +
    xlab("Cluster ID") + ylab("Significant Pairs of Peaks")

gl16 = ggplot(data= mchr16, aes(x = as.factor(plotGroup))) +
  geom_bar(aes(y = pairsOfPeaks), stat='identity', fill='lightblue', color='lightblue4') +
  geom_bar(aes(y = signifPairs), stat='identity', fill= 'pink', color='red') +
  geom_bar(aes(y = sigPairs), stat='identity', fill='mediumorchid1', color='mediumorchid4') +
    ggtitle("Chromosome 16") + theme_classic() + theme(plot.title = element_text(face = "bold.italic", hjust = 0.5)) +
    xlab("Cluster ID") + ylab("Significant Pairs of Peaks")

gl15 = ggplot(data= mchr15, aes(x = as.factor(plotGroup))) +
  geom_bar(aes(y = pairsOfPeaks), stat='identity', fill='lightblue', color='lightblue4') +
  geom_bar(aes(y = signifPairs), stat='identity', fill= 'pink', color='red') +
  geom_bar(aes(y = sigPairs), stat='identity', fill='mediumorchid1', color='mediumorchid4') +
    ggtitle("Chromosome 15") + theme_classic() + theme(plot.title = element_text(face = "bold.italic", hjust = 0.5)) +
    xlab("Cluster ID") + ylab("Significant Pairs of Peaks")

gl14 = ggplot(data= mchr14, aes(x = as.factor(plotGroup))) +
  geom_bar(aes(y = pairsOfPeaks), stat='identity', fill='lightblue', color='lightblue4') +
  geom_bar(aes(y = signifPairs), stat='identity', fill= 'pink', color='red') +
  geom_bar(aes(y = sigPairs), stat='identity', fill='mediumorchid1', color='mediumorchid4') +
    ggtitle("Chromosome 14") + theme_classic() + theme(plot.title = element_text(face = "bold.italic", hjust = 0.5)) +
    xlab("Cluster ID") + ylab("Significant Pairs of Peaks")

gl13 = ggplot(data= mchr13, aes(x = as.factor(plotGroup))) +
  geom_bar(aes(y = pairsOfPeaks), stat='identity', fill='lightblue', color='lightblue4') +
  geom_bar(aes(y = signifPairs), stat='identity', fill= 'pink', color='red') +
  geom_bar(aes(y = sigPairs), stat='identity', fill='mediumorchid1', color='mediumorchid4') +
    ggtitle("Chromosome 13") + theme_classic() + theme(plot.title = element_text(face = "bold.italic", hjust = 0.5)) +
    xlab("Cluster ID") + ylab("Significant Pairs of Peaks")

gl12 = ggplot(data= mchr12, aes(x = as.factor(plotGroup))) +
  geom_bar(aes(y = pairsOfPeaks), stat='identity', fill='lightblue', color='lightblue4') +
  geom_bar(aes(y = signifPairs), stat='identity', fill= 'pink', color='red') +
  geom_bar(aes(y = sigPairs), stat='identity', fill='mediumorchid1', color='mediumorchid4') +
    ggtitle("Chromosome 12") + theme_classic() + theme(plot.title = element_text(face = "bold.italic", hjust = 0.5)) +
    xlab("Cluster ID") + ylab("Significant Pairs of Peaks")

gl11 = ggplot(data= mchr11, aes(x = as.factor(plotGroup))) +
  geom_bar(aes(y = pairsOfPeaks), stat='identity', fill='lightblue', color='lightblue4') +
  geom_bar(aes(y = signifPairs), stat='identity', fill= 'pink', color='red') +
  geom_bar(aes(y = sigPairs), stat='identity', fill='mediumorchid1', color='mediumorchid4') +
    ggtitle("Chromosome 11") + theme_classic() + theme(plot.title = element_text(face = "bold.italic", hjust = 0.5)) +
    xlab("Cluster ID") + ylab("Significant Pairs of Peaks")

gl10 = ggplot(data= mchr10, aes(x = as.factor(plotGroup))) +
  geom_bar(aes(y = pairsOfPeaks), stat='identity', fill='lightblue', color='lightblue4') +
  geom_bar(aes(y = signifPairs), stat='identity', fill= 'pink', color='red') +
  geom_bar(aes(y = sigPairs), stat='identity', fill='mediumorchid1', color='mediumorchid4') +
    ggtitle("Chromosome 10") + theme_classic() + theme(plot.title = element_text(face = "bold.italic", hjust = 0.5)) +
    xlab("Cluster ID") + ylab("Significant Pairs of Peaks")

gl9 = ggplot(data= mchr9, aes(x = as.factor(plotGroup))) +
  geom_bar(aes(y = pairsOfPeaks), stat='identity', fill='lightblue', color='lightblue4') +
  geom_bar(aes(y = signifPairs), stat='identity', fill= 'pink', color='red') +
  geom_bar(aes(y = sigPairs), stat='identity', fill='mediumorchid1', color='mediumorchid4') +
    ggtitle("Chromosome 9") + theme_classic() + theme(plot.title = element_text(face = "bold.italic", hjust = 0.5)) +
    xlab("Cluster ID") + ylab("Significant Pairs of Peaks")

gl8 = ggplot(data= mchr8, aes(x = as.factor(plotGroup))) +
  geom_bar(aes(y = pairsOfPeaks), stat='identity', fill='lightblue', color='lightblue4') +
  geom_bar(aes(y = signifPairs), stat='identity', fill= 'pink', color='red') +
  geom_bar(aes(y = sigPairs), stat='identity', fill='mediumorchid1', color='mediumorchid4') +
    ggtitle("Chromosome 8") + theme_classic() + theme(plot.title = element_text(face = "bold.italic", hjust = 0.5)) +
    xlab("Cluster ID") + ylab("Significant Pairs of Peaks")

gl7 = ggplot(data= mchr7, aes(x = as.factor(plotGroup))) +
  geom_bar(aes(y = pairsOfPeaks), stat='identity', fill='lightblue', color='lightblue4') +
  geom_bar(aes(y = signifPairs), stat='identity', fill= 'pink', color='red') +
  geom_bar(aes(y = sigPairs), stat='identity', fill='mediumorchid1', color='mediumorchid4') +
    ggtitle("Chromosome 7") + theme_classic() + theme(plot.title = element_text(face = "bold.italic", hjust = 0.5)) +
    xlab("Cluster ID") + ylab("Significant Pairs of Peaks")

gl6 = ggplot(data= mchr6, aes(x = as.factor(plotGroup))) +
  geom_bar(aes(y = pairsOfPeaks), stat='identity', fill='lightblue', color='lightblue4') +
  geom_bar(aes(y = signifPairs), stat='identity', fill= 'pink', color='red') +
  geom_bar(aes(y = sigPairs), stat='identity', fill='mediumorchid1', color='mediumorchid4') +
    ggtitle("Chromosome 6") + theme_classic() + theme(plot.title = element_text(face = "bold.italic", hjust = 0.5)) +
    xlab("Cluster ID") + ylab("Significant Pairs of Peaks")

gl5 = ggplot(data= mchr5, aes(x = as.factor(plotGroup))) +
  geom_bar(aes(y = pairsOfPeaks), stat='identity', fill='lightblue', color='lightblue4') +
  geom_bar(aes(y = signifPairs), stat='identity', fill= 'pink', color='red') +
  geom_bar(aes(y = sigPairs), stat='identity', fill='mediumorchid1', color='mediumorchid4') +
    ggtitle("Chromosome 5") + theme_classic() + theme(plot.title = element_text(face = "bold.italic", hjust = 0.5)) +
    xlab("Cluster ID") + ylab("Significant Pairs of Peaks")

gl4 = ggplot(data= mchr4, aes(x = as.factor(plotGroup))) +
  geom_bar(aes(y = pairsOfPeaks), stat='identity', fill='lightblue', color='lightblue4') +
  geom_bar(aes(y = signifPairs), stat='identity', fill= 'pink', color='red') +
  geom_bar(aes(y = sigPairs), stat='identity', fill='mediumorchid1', color='mediumorchid4') +
    ggtitle("Chromosome 4") + theme_classic() + theme(plot.title = element_text(face = "bold.italic", hjust = 0.5)) +
    xlab("Cluster ID") + ylab("Significant Pairs of Peaks")

gl3 = ggplot(data= mchr3, aes(x = as.factor(plotGroup))) +
  geom_bar(aes(y = pairsOfPeaks), stat='identity', fill='lightblue', color='lightblue4') +
  geom_bar(aes(y = signifPairs), stat='identity', fill= 'pink', color='red') +
  geom_bar(aes(y = sigPairs), stat='identity', fill='mediumorchid1', color='mediumorchid4') +
    ggtitle("Chromosome 3") + theme_classic() + theme(plot.title = element_text(face = "bold.italic", hjust = 0.5)) +
    xlab("Cluster ID") + ylab("Significant Pairs of Peaks")

gl2 = ggplot(data= mchr2, aes(x = as.factor(plotGroup))) +
  geom_bar(aes(y = pairsOfPeaks), stat='identity', fill='lightblue', color='lightblue4') +
  geom_bar(aes(y = signifPairs), stat='identity', fill= 'pink', color='red') +
  geom_bar(aes(y = sigPairs), stat='identity', fill='mediumorchid1', color='mediumorchid4') +
    ggtitle("Chromosome 2") + theme_classic() + theme(plot.title = element_text(face = "bold.italic", hjust = 0.5)) +
    xlab("Cluster ID") + ylab("Significant Pairs of Peaks")

gl1 = ggplot(data= mchr1, aes(x = as.factor(plotGroup))) +
  geom_bar(aes(y = pairsOfPeaks), stat='identity', fill='lightblue', color='lightblue4') +
  geom_bar(aes(y = signifPairs), stat='identity', fill= 'pink', color='red') +
  geom_bar(aes(y = sigPairs), stat='identity', fill='mediumorchid1', color='mediumorchid4') +
    ggtitle("Chromosome 1") + theme_classic() + theme(plot.title = element_text(face = "bold.italic", hjust = 0.5)) +
    xlab("Cluster ID") + ylab("Significant Pairs of Peaks")


pdf(file = "/sc/arion/projects/epigenAD/Bukola/checkGL.pdf")

  grid.arrange(gl1, gl2, gl3, gl4, gl5, gl6, nrow = 3, ncol = 2)

  grid.arrange(gl7, gl8, gl9, gl10, gl11, gl12, nrow = 3, ncol = 2)

  grid.arrange(gl13, gl14, gl15, gl16, gl17, gl18, nrow = 3, ncol = 2)

  grid.arrange(gl19, gl20, gl21, gl22, glY, nrow = 3, ncol = 2)

dev.off()

scp -r ajanab02@minerva.hpc.mssm.edu:/sc/arion/projects/epigenAD/Bukola/checkGL.pdf /Users/bukola/Documents/finalplots




################################################################################
GL_df_75p = GL_df
GL_df_75p$sigPairs = NA

for(i in 1:uniqueN(GL_df_75p$cluster)){
  clus = unique(GL_df_75p$cluster)[i]
  subGL = GL_df_75p[GL_df_75p$cluster == unique(GL_df_75p$cluster)[i], ]
  countAll = sum(subGL$signif)
  GL_df_75p[GL_df_75p$cluster == unique(GL_df_75p$cluster)[i], "sigPairs"] = countAll
}

forPlotGL75 = GL_df_75p[, c(12, 13)]
forPlotGL75 = unique(forPlotGL75)
forPlotGL75$forComb = forPlotGL75$cluster
merging75 = merge(forPlotGL75, clustersIn75, by = "forComb")

merging75 = as.data.frame(merging75)
load75 = merging75[c(2,3,5,6)]
load75 = data.frame(lapply(load75, as.character), stringsAsFactors=FALSE)

write.csv(load75, file= "/sc/arion/projects/epigenAD/Bukola/merging75.csv")
scp -r ajanab02@minerva.hpc.mssm.edu:/sc/arion/projects/epigenAD/Bukola/merging75.csv /Users/bukola/Documents/finalplots
### junk code
#  geom_smooth(data = dfFinalObject_10, method = "nls", se = FALSE,
#  method.args = list(formula = y~a*exp(b/x), start=list(a=1, b=0.1)))
  stat_smooth(method = 'nls', formula = y ~ a*exp(b *x), aes(colour = "black"), se = FALSE, start = list(a=1,b=1))



LLM1 = lm(I(log(width)) ~ mean_abs_corr, data= dfFinalObject_10)

  ggplot(GL_results, aes(distance,value.y, color = factor(signif)))+geom_point()+xlab("Distance (bp)")+ylab("Correlation")+theme_classic() +
    ggtitle("Comparing Cor & Signif over Distance ") + labs(color = "Significance")
