## Original Glasso Code for simplify499.R
## Bukola Ajanaku
## Saved March 10, 2022

### Defining variables for GLASSO: *********************************************
input_mat = quantResid2 #our peaks matrix
cluster_peaks = finalObject[finalObject$forComb == "80_chr1_41",]$ID
clusterMat = input_mat[rownames(input_mat) %in% cluster_peaks,]

## GLASSO CODE: ################################################################
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

## GRAPHING ####################################################################
# 1: Not really important
ggplot(GL_data, aes(Var1, Var2, fill= value)) +
    geom_tile()+ scale_fill_gradient(low="white", high="green") +
    theme_classic() + geom_tile(data = highlight_df, fill = "red")

# 2: Medium Importance
plot(sort(distance_df_melt$distance),type= "b")

# 3: IMPORTANT
ggplot(GL_results, aes(distance,value.y, color = factor(signif)))+geom_point()+xlab("Distance (bp)")+ylab("Correlation")+theme_classic() +
  ggtitle("Comparing Cor & Signif over Distance ") + labs(color = "Significance")

# 4: IMPORTANT
highlight_df = GL_results[GL_results$signif == 0,]

ggplot(GL_results, aes(Var1, Var2, fill= value.y)) +
    geom_tile()+ scale_fill_gradient(low="white", high="green") +
    theme_classic() + geom_tile(data = highlight_df, fill = "red")

########### for changing s:

pdf(file = "/sc/arion/projects/epigenAD/Bukola/corrMats/80_chr1_41/s0.2.pdf")
    ggplot(GL_results, aes(Var1, Var2, fill= value.y)) +
        geom_tile()+ scale_fill_gradient(low="white", high="green") +
        theme_classic()

    highlight_df = GL_results[GL_results$signif == 0,]
    ggplot(GL_results, aes(Var1, Var2, fill= value.y)) +
        geom_tile()+ scale_fill_gradient(low="white", high="green") +
        theme_classic() + geom_tile(data = highlight_df, fill = "red")

dev.off()

# scp -r ajanab02@minerva.hpc.mssm.edu:/sc/arion/projects/epigenAD/Bukola/corrMats/80_chr1_41 /Users/bukola/Documents/plotsforKiran
