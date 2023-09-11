### This script calculates the cutoffs for the interaction model, interaction variable, that ensure an FDR of 0.05 
### and 0.1 using a dataset in which the gene expression values are permuted between the samples (with a different 
### permutation for each gene).

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

pval_int_time_intlm <- readRDS("~/Elegans_internship/Data/new_ttest_pvals/pval_int_time_intlm.RDS")
ol_pval_int_time_intlm_perm <- readRDS("~/Elegans_internship/Data/new_ttest_pvals/pval_int_time_intlm_oldperm.RDS")
pval_int_time_intlm_perm <-readRDS("~/Elegans_internship/Data/new_perm_one_per_gene/pval_int_time_intlm_one_perm.RDS")


# Get frame with most significant pvalues per gene for each chromosome for the interaction model, interaction
# variable, not permuted.
min_p_int_int_chr1 <- c()
min_p_int_int_chr2 <- c()
min_p_int_int_chr3 <- c()
min_p_int_int_chr4 <- c()
min_p_int_int_chr5 <- c()
min_p_int_int_chrx <- c()


for(i in 1:ncol(pval_int_time_intlm)){
  min_p_int_int_chr1 <- c(min_p_int_int_chr1,min(pval_int_time_intlm[FourP$snp_info$ChrID=="I",i]))
}

for(i in 1:ncol(pval_int_time_intlm)){
  min_p_int_int_chr2 <- c(min_p_int_int_chr2,min(pval_int_time_intlm[FourP$snp_info$ChrID=="II",i]))
}

for(i in 1:ncol(pval_int_time_intlm)){
  min_p_int_int_chr3 <- c(min_p_int_int_chr3,min(pval_int_time_intlm[FourP$snp_info$ChrID=="III",i]))
}

for(i in 1:ncol(pval_int_time_intlm)){
  min_p_int_int_chr4 <- c(min_p_int_int_chr4,min(pval_int_time_intlm[FourP$snp_info$ChrID=="IV",i]))
}

for(i in 1:ncol(pval_int_time_intlm)){
  min_p_int_int_chr5 <- c(min_p_int_int_chr5,min(pval_int_time_intlm[FourP$snp_info$ChrID=="V",i]))
}

for(i in 1:ncol(pval_int_time_intlm)){
  min_p_int_int_chrx <- c(min_p_int_int_chrx,min(pval_int_time_intlm[FourP$snp_info$ChrID=="X",i]))
}

min_p_int_int_chr1 <- data.frame(marker = min_p_int_int_chr1,chr = rep("I",length(min_p_int_int_chr1)))
min_p_int_int_chr2 <- data.frame(marker = min_p_int_int_chr2,chr = rep("II",length(min_p_int_int_chr2)))
min_p_int_int_chr3 <- data.frame(marker = min_p_int_int_chr3,chr = rep("III",length(min_p_int_int_chr3)))
min_p_int_int_chr4 <- data.frame(marker = min_p_int_int_chr4,chr = rep("IV",length(min_p_int_int_chr4)))
min_p_int_int_chr5 <- data.frame(marker = min_p_int_int_chr5,chr = rep("V",length(min_p_int_int_chr5)))
min_p_int_int_chrx <- data.frame(marker = min_p_int_int_chrx,chr = rep("X",length(min_p_int_int_chrx)))

min_p_per_chr <- rbind(min_p_int_int_chr1,min_p_int_int_chr2,min_p_int_int_chr3,min_p_int_int_chr4,min_p_int_int_chr5,min_p_int_int_chrx)

# Do same for permuted dataset.

min_p_int_int_chr_perm1 <- c()
min_p_int_int_chr_perm2 <- c()
min_p_int_int_chr_perm3 <- c()
min_p_int_int_chr_perm4 <- c()
min_p_int_int_chr_perm5 <- c()
min_p_int_int_chr_permx <- c()


for(i in 1:ncol(pval_int_time_intlm_perm)){
  min_p_int_int_chr_perm1 <- c(min_p_int_int_chr_perm1,min(pval_int_time_intlm_perm[FourP$snp_info$ChrID=="I",i]))
}

for(i in 1:ncol(pval_int_time_intlm_perm)){
  min_p_int_int_chr_perm2 <- c(min_p_int_int_chr_perm2,min(pval_int_time_intlm_perm[FourP$snp_info$ChrID=="II",i]))
}

for(i in 1:ncol(pval_int_time_intlm_perm)){
  min_p_int_int_chr_perm3 <- c(min_p_int_int_chr_perm3,min(pval_int_time_intlm_perm[FourP$snp_info$ChrID=="III",i]))
}

for(i in 1:ncol(pval_int_time_intlm_perm)){
  min_p_int_int_chr_perm4 <- c(min_p_int_int_chr_perm4,min(pval_int_time_intlm_perm[FourP$snp_info$ChrID=="IV",i]))
}

for(i in 1:ncol(pval_int_time_intlm_perm)){
  min_p_int_int_chr_perm5 <- c(min_p_int_int_chr_perm5,min(pval_int_time_intlm_perm[FourP$snp_info$ChrID=="V",i]))
}

for(i in 1:ncol(pval_int_time_intlm_perm)){
  min_p_int_int_chr_permx <- c(min_p_int_int_chr_permx,min(pval_int_time_intlm_perm[FourP$snp_info$ChrID=="X",i]))
}

min_p_int_int_chr_perm1 <- data.frame(marker = min_p_int_int_chr_perm1,chr = rep("I",length(min_p_int_int_chr_perm1)))
min_p_int_int_chr_perm2 <- data.frame(marker = min_p_int_int_chr_perm2,chr = rep("II",length(min_p_int_int_chr_perm2)))
min_p_int_int_chr_perm3 <- data.frame(marker = min_p_int_int_chr_perm3,chr = rep("III",length(min_p_int_int_chr_perm3)))
min_p_int_int_chr_perm4 <- data.frame(marker = min_p_int_int_chr_perm4,chr = rep("IV",length(min_p_int_int_chr_perm4)))
min_p_int_int_chr_perm5 <- data.frame(marker = min_p_int_int_chr_perm5,chr = rep("V",length(min_p_int_int_chr_perm5)))
min_p_int_int_chr_permx <- data.frame(marker = min_p_int_int_chr_permx,chr = rep("X",length(min_p_int_int_chr_permx)))

min_p_per_chr_perm <- rbind(min_p_int_int_chr_perm1,min_p_int_int_chr_perm2,min_p_int_int_chr_perm3,min_p_int_int_chr_perm4,min_p_int_int_chr_perm5,min_p_int_int_chr_permx)

#Calculate cutoffs that ensures 0.05 and 0.1 FDR
cutoff_int_int <- 5.60
sum(-log10(min_p_per_chr_perm[,1])>cutoff_int_int)/sum(-log10(min_p_per_chr[,1])>cutoff_int_int)
## Now that we've used correct samples its impossible to threshold at 0.05.

cutoff_int_int01 <- 4.52
sum(-log10(min_p_per_chr_perm[,1])>cutoff_int_int01)/sum(-log10(min_p_per_chr[,1])>cutoff_int_int01)

cutoff_int_int02 <- 3.68
sum(-log10(min_p_per_chr_perm[,1])>cutoff_int_int02)/sum(-log10(min_p_per_chr[,1])>cutoff_int_int02)

# Intelligently choose cutoff.
test <- 6
sum(-log10(min_p_per_chr_perm[,1])>test)/sum(-log10(min_p_per_chr[,1])>test)

cutoff_test <- c()

for(i in 0:1000){
  cutoff_test[i+1] <- sum(-log10(min_p_per_chr_perm[,1])>i/100)/sum(-log10(min_p_per_chr[,1])>i/100)
}

plot(seq(0,10,by = 0.01),cutoff_test,xlab = "Cutoff",ylab="FDR",ylim = c(0,0.2))

par(mfrow = c(1,2))
plot(apply(-log10(pval_int_time_intlm_perm),2,max),ylab = "P-value new perm")
plot(apply(-log10(ol_pval_int_time_intlm_perm),2,max),ylab = "P-value old  perm")











