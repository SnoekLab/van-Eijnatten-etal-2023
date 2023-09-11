library(ggplot2)
library(randomForest)
library(RcppArmadillo)
library(cowplot)
library(dplyr)

# Import snp and fpkm data. 
load(file="C:/Users/Gebruiker/Documents/Elegans_internship/Data/obj_FourP_2022.out") 
load(file="C:/Users/Gebruiker/Documents/Elegans_internship/Data/obj_fpkm.out")

# Filter for mean samples fpkm > 1/32 and more than 20 samples have > 0 fpkm.
fpkm.selc <- apply(log2(fpkm$fpkm),1,mean)> -5 & apply(fpkm$fpkm>0,1,sum)>20 

# ge.mat contains all filtered samples with genes as observations and samples as features.
ge.mat <- fpkm$fpkm[fpkm.selc,]

# Some genes are transcribed more than others, but we are interested in relative differences. Take log2 ratio of mean.
rge.mat <- log2((ge.mat+1)/apply(ge.mat+1,1,mean))

# Do pca on rge.mat, and get the loadings of the samples on the first two principles components.
pco <- prcomp((rge.mat))
pco1 <- -pco$rotation[,1]
pco2 <- pco$rotation[,2]

# Select all rows for which the maximal log2 ratio is -1 < or > 1 (corresponds to a 2 fold change relative to the mean)
use.rge <- rge.mat#rge.mat[apply(abs(rge.mat),1,max)>1,]

# Select only non-parent strains.
use.lines <- c(1:198,207)
pco1.lines <- pco1[use.lines]

# Store names of all genes which are not filtered out.
gene_names <- fpkm$genes[fpkm.selc]

pval_int_spline <- readRDS("~/Elegans_internship/Data/spline_pvals/pval_int_spline.RDS")
pval_int_perm_spline <- readRDS("~/Elegans_internship/Data/spline_pvals/pval_int_spline_perm.RDS")

# Get frame with most significant pvalues per gene for each chromosome for notime model not permuted.
marker_int_spline1 <- c()
marker_int_spline2 <- c()
marker_int_spline3 <- c()
marker_int_spline4 <- c()
marker_int_spline5 <- c()
marker_int_splinex <- c()


for(i in 1:ncol(pval_int_spline)){
  marker_int_spline1 <- c(marker_int_spline1,min(pval_int_spline[FourP$snp_info$ChrID=="I",i]))
}

for(i in 1:ncol(pval_int_spline)){
  marker_int_spline2 <- c(marker_int_spline2,min(pval_int_spline[FourP$snp_info$ChrID=="II",i]))
}

for(i in 1:ncol(pval_int_spline)){
  marker_int_spline3 <- c(marker_int_spline3,min(pval_int_spline[FourP$snp_info$ChrID=="III",i]))
}

for(i in 1:ncol(pval_int_spline)){
  marker_int_spline4 <- c(marker_int_spline4,min(pval_int_spline[FourP$snp_info$ChrID=="IV",i]))
}

for(i in 1:ncol(pval_int_spline)){
  marker_int_spline5 <- c(marker_int_spline5,min(pval_int_spline[FourP$snp_info$ChrID=="V",i]))
}

for(i in 1:ncol(pval_int_spline)){
  marker_int_splinex <- c(marker_int_splinex,min(pval_int_spline[FourP$snp_info$ChrID=="X",i]))
}

marker_int_spline1 <- data.frame(marker = marker_int_spline1,chr = rep("I",length(marker_int_spline1)))
marker_int_spline2 <- data.frame(marker = marker_int_spline2,chr = rep("II",length(marker_int_spline2)))
marker_int_spline3 <- data.frame(marker = marker_int_spline3,chr = rep("III",length(marker_int_spline3)))
marker_int_spline4 <- data.frame(marker = marker_int_spline4,chr = rep("IV",length(marker_int_spline4)))
marker_int_spline5 <- data.frame(marker = marker_int_spline5,chr = rep("V",length(marker_int_spline5)))
marker_int_splinex <- data.frame(marker = marker_int_splinex,chr = rep("X",length(marker_int_splinex)))

min_p_per_chr <- rbind(marker_int_spline1,marker_int_spline2,marker_int_spline3,marker_int_spline4,marker_int_spline5,marker_int_splinex)

# Do same for permuted dataset.

marker_int_spline_perm1 <- c()
marker_int_spline_perm2 <- c()
marker_int_spline_perm3 <- c()
marker_int_spline_perm4 <- c()
marker_int_spline_perm5 <- c()
marker_int_spline_permx <- c()


for(i in 1:ncol(pval_int_perm_spline)){
  marker_int_spline_perm1 <- c(marker_int_spline_perm1,min(pval_int_perm_spline[FourP$snp_info$ChrID=="I",i]))
}

for(i in 1:ncol(pval_int_perm_spline)){
  marker_int_spline_perm2 <- c(marker_int_spline_perm2,min(pval_int_perm_spline[FourP$snp_info$ChrID=="II",i]))
}

for(i in 1:ncol(pval_int_perm_spline)){
  marker_int_spline_perm3 <- c(marker_int_spline_perm3,min(pval_int_perm_spline[FourP$snp_info$ChrID=="III",i]))
}

for(i in 1:ncol(pval_int_perm_spline)){
  marker_int_spline_perm4 <- c(marker_int_spline_perm4,min(pval_int_perm_spline[FourP$snp_info$ChrID=="IV",i]))
}

for(i in 1:ncol(pval_int_perm_spline)){
  marker_int_spline_perm5 <- c(marker_int_spline_perm5,min(pval_int_perm_spline[FourP$snp_info$ChrID=="V",i]))
}

for(i in 1:ncol(pval_int_perm_spline)){
  marker_int_spline_permx <- c(marker_int_spline_permx,min(pval_int_perm_spline[FourP$snp_info$ChrID=="X",i]))
}

marker_int_spline_perm1 <- data.frame(marker = marker_int_spline_perm1,chr = rep("I",length(marker_int_spline_perm1)))
marker_int_spline_perm2 <- data.frame(marker = marker_int_spline_perm2,chr = rep("II",length(marker_int_spline_perm2)))
marker_int_spline_perm3 <- data.frame(marker = marker_int_spline_perm3,chr = rep("III",length(marker_int_spline_perm3)))
marker_int_spline_perm4 <- data.frame(marker = marker_int_spline_perm4,chr = rep("IV",length(marker_int_spline_perm4)))
marker_int_spline_perm5 <- data.frame(marker = marker_int_spline_perm5,chr = rep("V",length(marker_int_spline_perm5)))
marker_int_spline_permx <- data.frame(marker = marker_int_spline_permx,chr = rep("X",length(marker_int_spline_permx)))

min_p_per_chr_perm <- rbind(marker_int_spline_perm1,marker_int_spline_perm2,marker_int_spline_perm3,marker_int_spline_perm4,marker_int_spline_perm5,marker_int_spline_permx)

cutoff_int_spline <- 4.85
sum(-log10(min_p_per_chr_perm[,1])>cutoff_int_spline)/sum(-log10(min_p_per_chr[,1])>cutoff_int_spline)

cutoff_int_spline01 <-4.47
sum(-log10(min_p_per_chr_perm[,1])>cutoff_int_spline01)/sum(-log10(min_p_per_chr[,1])>cutoff_int_spline01)

cutoff_test <- c()

for(i in 0:1000){
  cutoff_test[i+1] <- sum(-log10(min_p_per_chr_perm[,1])>i/100)/sum(-log10(min_p_per_chr[,1])>i/100)
}

plot(seq(0,10,by = 0.01),cutoff_test,xlab = "Cutoff",ylab="FDR")

sum(-log10(pval_int_perm_spline)>3)
sum(-log10(pval_int_spline)>3)










