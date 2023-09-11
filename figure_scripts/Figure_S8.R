.libPaths("~/R/x86_64-pc-linux-gnu-library/3.6/")

library(ggplot2)
library(cowplot)
library(dplyr)
library(viridis)
library(tidyverse)
library(lemon)

# Import snp and fpkm data. 
load(file="~/Documents/elegans/data/obj_FourP_2022.out") 
load(file="~/Documents/elegans/data/obj_fpkm.out")

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

# Load in p-values of the genotype variables for models without development, with additive development and 
# with an interaction with development. Also load in p-value of the interaction term in the interaction model.
pval_gen_notime_lm <- readRDS("~/Documents/elegans/pvals_correct/pval_gen_notime_lm.RDS")
pval_gen_time_lm <- readRDS("~/Documents/elegans/map_with_mean_dev_expr/pval_additive_marker_mean_dev_age.RDS")
pval_gen_time_intlm <- readRDS("~/Documents/elegans/map_with_mean_dev_expr/pval_gen_time_int_mean_dev.RDS")
pval_int_time_intlm <- readRDS("~/Documents/elegans/map_with_mean_dev_expr/pval_int_time_int_mean_dev.RDS")

# Find all eQTLs for models without development. eQTLs are defined as the marker with the highest -log10(p-value),
# provided this highest value is higher than the cutoff we calculated with the permutation. We call a maximum of 
# one eQTL per chromosome. With the below code we make a dataframe with the position on the genome of the
# markers, a factor whether the eQTL is cis or trans, and the chromosome ID.

cutoff_notime <- 4.06

min_p_no_time_pos_chr1 <- c()
min_p_no_time_pos_chr2 <- c()
min_p_no_time_pos_chr3 <- c()
min_p_no_time_pos_chr4 <- c()
min_p_no_time_pos_chr5 <- c()
min_p_no_time_pos_chrx <- c()

ct_gen_notime_chr1 <- c()
ct_gen_notime_chr2 <- c()
ct_gen_notime_chr3 <- c()
ct_gen_notime_chr4 <- c()
ct_gen_notime_chr5 <- c()
ct_gen_notime_chrx <- c()

# Save all genes which are affected by the large hotspot on chromosome X.
hotspot <- c()
#saveRDS(hotspot,"hotspot_index.RDS")

for(i in 1:ncol(pval_gen_notime_lm)){
  if(max(-log10(pval_gen_notime_lm[FourP$snp_info$ChrID=="I",i]))>cutoff_notime){
    min_p_no_time_pos_chr1 <- c(min_p_no_time_pos_chr1,FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="I"][which.min(pval_gen_notime_lm[FourP$snp_info$ChrID=="I",i])])
    loc <- gsub(".*:","",fpkm$locus[fpkm.selc][i])
    middle <- (as.numeric(sub("(^[^-]+)-.*", "\\1", loc)) + as.numeric(sub('.+-(.+)', '\\1', loc))) / 2 / 1000000
    gene_chr <- strsplit(as.character(fpkm$locus[fpkm.selc][i]), "[:]")[[1]][1]
    ifelse(abs(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="I"][which.min(pval_gen_notime_lm[FourP$snp_info$ChrID=="I",i])] - middle) < 2 & gene_chr =="I",ct_gen_notime_chr1 <- c(ct_gen_notime_chr1,"cis"),ct_gen_notime_chr1 <- c(ct_gen_notime_chr1,"trans"))
  }
}

for(i in 1:ncol(pval_gen_notime_lm)){
  if(max(-log10(pval_gen_notime_lm[FourP$snp_info$ChrID=="II",i]))>cutoff_notime){
    min_p_no_time_pos_chr2 <- c(min_p_no_time_pos_chr2,FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="II"][which.min(pval_gen_notime_lm[FourP$snp_info$ChrID=="II",i])])
    loc <- gsub(".*:","",fpkm$locus[fpkm.selc][i])
    middle <- (as.numeric(sub("(^[^-]+)-.*", "\\1", loc)) + as.numeric(sub('.+-(.+)', '\\1', loc))) / 2 / 1000000
    gene_chr <- strsplit(as.character(fpkm$locus[fpkm.selc][i]), "[:]")[[1]][1]
    ifelse(abs(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="II"][which.min(pval_gen_notime_lm[FourP$snp_info$ChrID=="II",i])] - middle) < 2 & gene_chr =="II",ct_gen_notime_chr2 <- c(ct_gen_notime_chr2,"cis"),ct_gen_notime_chr2 <- c(ct_gen_notime_chr2,"trans"))
  }
}

for(i in 1:ncol(pval_gen_notime_lm)){
  if(max(-log10(pval_gen_notime_lm[FourP$snp_info$ChrID=="III",i]))>cutoff_notime){
    min_p_no_time_pos_chr3 <- c(min_p_no_time_pos_chr3,FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="III"][which.min(pval_gen_notime_lm[FourP$snp_info$ChrID=="III",i])])
    loc <- gsub(".*:","",fpkm$locus[fpkm.selc][i])
    middle <- (as.numeric(sub("(^[^-]+)-.*", "\\1", loc)) + as.numeric(sub('.+-(.+)', '\\1', loc))) / 2 / 1000000
    gene_chr <- strsplit(as.character(fpkm$locus[fpkm.selc][i]), "[:]")[[1]][1]
    ifelse(abs(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="III"][which.min(pval_gen_notime_lm[FourP$snp_info$ChrID=="III",i])] - middle) < 2 & gene_chr =="III",ct_gen_notime_chr3 <- c(ct_gen_notime_chr3,"cis"),ct_gen_notime_chr3 <- c(ct_gen_notime_chr3,"trans"))  }
}

for(i in 1:ncol(pval_gen_notime_lm)){
  if(max(-log10(pval_gen_notime_lm[FourP$snp_info$ChrID=="IV",i]))>cutoff_notime){
    min_p_no_time_pos_chr4 <- c(min_p_no_time_pos_chr4,FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="IV"][which.min(pval_gen_notime_lm[FourP$snp_info$ChrID=="IV",i])])
    loc <- gsub(".*:","",fpkm$locus[fpkm.selc][i])
    middle <- (as.numeric(sub("(^[^-]+)-.*", "\\1", loc)) + as.numeric(sub('.+-(.+)', '\\1', loc))) / 2 / 1000000
    gene_chr <- strsplit(as.character(fpkm$locus[fpkm.selc][i]), "[:]")[[1]][1]
    ifelse(abs(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="IV"][which.min(pval_gen_notime_lm[FourP$snp_info$ChrID=="IV",i])] - middle) < 2 & gene_chr =="IV",ct_gen_notime_chr4 <- c(ct_gen_notime_chr4,"cis"),ct_gen_notime_chr4 <- c(ct_gen_notime_chr4,"trans"))  }
}

for(i in 1:ncol(pval_gen_notime_lm)){
  if(max(-log10(pval_gen_notime_lm[FourP$snp_info$ChrID=="V",i]))>cutoff_notime){
    min_p_no_time_pos_chr5 <- c(min_p_no_time_pos_chr5,FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="V"][which.min(pval_gen_notime_lm[FourP$snp_info$ChrID=="V",i])])
    loc <- gsub(".*:","",fpkm$locus[fpkm.selc][i])
    middle <- (as.numeric(sub("(^[^-]+)-.*", "\\1", loc)) + as.numeric(sub('.+-(.+)', '\\1', loc))) / 2 / 1000000
    gene_chr <- strsplit(as.character(fpkm$locus[fpkm.selc][i]), "[:]")[[1]][1]
    ifelse(abs(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="V"][which.min(pval_gen_notime_lm[FourP$snp_info$ChrID=="V",i])] - middle) < 2 & gene_chr =="V",ct_gen_notime_chr5 <- c(ct_gen_notime_chr5,"cis"),ct_gen_notime_chr5 <- c(ct_gen_notime_chr5,"trans"))
  }
}

for(i in 1:ncol(pval_gen_notime_lm)){
  if(max(-log10(pval_gen_notime_lm[FourP$snp_info$ChrID=="X",i]))>cutoff_notime){
    min_p_no_time_pos_chrx <- c(min_p_no_time_pos_chrx,FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="X"][which.min(pval_gen_notime_lm[FourP$snp_info$ChrID=="X",i])])
    loc <- gsub(".*:","",fpkm$locus[fpkm.selc][i])
    middle <- (as.numeric(sub("(^[^-]+)-.*", "\\1", loc)) + as.numeric(sub('.+-(.+)', '\\1', loc))) / 2 / 1000000
    gene_chr <- strsplit(as.character(fpkm$locus[fpkm.selc][i]), "[:]")[[1]][1]
    ifelse(abs(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="X"][which.min(pval_gen_notime_lm[FourP$snp_info$ChrID=="X",i])] - middle) < 2 & gene_chr =="X",ct_gen_notime_chrx <- c(ct_gen_notime_chrx,"cis"),ct_gen_notime_chrx <- c(ct_gen_notime_chrx,"trans"))
    if(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="X"][which.min(pval_gen_notime_lm[FourP$snp_info$ChrID=="X",i])]==3.404162){
      hotspot <- c(hotspot,i)
    }
  }
}

min_p_no_time_pos_chr1 <- data.frame(marker = min_p_no_time_pos_chr1,ct = ct_gen_notime_chr1, chr = rep("I",length(min_p_no_time_pos_chr1)))
min_p_no_time_pos_chr2 <- data.frame(marker = min_p_no_time_pos_chr2,ct = ct_gen_notime_chr2,chr = rep("II",length(min_p_no_time_pos_chr2)))
min_p_no_time_pos_chr3 <- data.frame(marker = min_p_no_time_pos_chr3,ct = ct_gen_notime_chr3,chr = rep("III",length(min_p_no_time_pos_chr3)))
min_p_no_time_pos_chr4 <- data.frame(marker = min_p_no_time_pos_chr4,ct = ct_gen_notime_chr4,chr = rep("IV",length(min_p_no_time_pos_chr4)))
min_p_no_time_pos_chr5 <- data.frame(marker = min_p_no_time_pos_chr5,ct = ct_gen_notime_chr5,chr = rep("V",length(min_p_no_time_pos_chr5)))
min_p_no_time_pos_chrx <- data.frame(marker = min_p_no_time_pos_chrx,ct = ct_gen_notime_chrx,chr = rep("X",length(min_p_no_time_pos_chrx)))

min_p_max_per_chr <- rbind(min_p_no_time_pos_chr1,min_p_no_time_pos_chr2,min_p_no_time_pos_chr3,min_p_no_time_pos_chr4,min_p_no_time_pos_chr5,min_p_no_time_pos_chrx)



# Now do the same for the genotype variable in the additive model.

cutoff_gen_time <- 4.46

min_p_time_pos_chr1 <- c()
min_p_time_pos_chr2 <- c()
min_p_time_pos_chr3 <- c()
min_p_time_pos_chr4 <- c()
min_p_time_pos_chr5 <- c()
min_p_time_pos_chrx <- c()

ct_gen_time_chr1 <- c()
ct_gen_time_chr2 <- c()
ct_gen_time_chr3 <- c()
ct_gen_time_chr4 <- c()
ct_gen_time_chr5 <- c()
ct_gen_time_chrx <- c()

for(i in 1:ncol(pval_gen_time_lm)){
  if(max(-log10(pval_gen_time_lm[FourP$snp_info$ChrID=="I",i]))>cutoff_gen_time){
    min_p_time_pos_chr1 <- c(min_p_time_pos_chr1,FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="I"][which.min(pval_gen_time_lm[FourP$snp_info$ChrID=="I",i])])
    loc <- gsub(".*:","",fpkm$locus[fpkm.selc][i])
    middle <- (as.numeric(sub("(^[^-]+)-.*", "\\1", loc)) + as.numeric(sub('.+-(.+)', '\\1', loc))) / 2 / 1000000
    gene_chr <- strsplit(as.character(fpkm$locus[fpkm.selc][i]), "[:]")[[1]][1]
    ifelse(abs(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="I"][which.min(pval_gen_time_lm[FourP$snp_info$ChrID=="I",i])] - middle) < 2 & gene_chr =="I",ct_gen_time_chr1 <- c(ct_gen_time_chr1,"cis"),ct_gen_time_chr1 <- c(ct_gen_time_chr1,"trans"))
  }
}

for(i in 1:ncol(pval_gen_time_lm)){
  if(max(-log10(pval_gen_time_lm[FourP$snp_info$ChrID=="II",i]))>cutoff_gen_time){
    min_p_time_pos_chr2 <- c(min_p_time_pos_chr2,FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="II"][which.min(pval_gen_time_lm[FourP$snp_info$ChrID=="II",i])])
    loc <- gsub(".*:","",fpkm$locus[fpkm.selc][i])
    middle <- (as.numeric(sub("(^[^-]+)-.*", "\\1", loc)) + as.numeric(sub('.+-(.+)', '\\1', loc))) / 2 / 1000000
    gene_chr <- strsplit(as.character(fpkm$locus[fpkm.selc][i]), "[:]")[[1]][1]
    ifelse(abs(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="II"][which.min(pval_gen_time_lm[FourP$snp_info$ChrID=="II",i])] - middle) < 2 & gene_chr =="II",ct_gen_time_chr2 <- c(ct_gen_time_chr2,"cis"),ct_gen_time_chr2 <- c(ct_gen_time_chr2,"trans"))
  }
}

for(i in 1:ncol(pval_gen_time_lm)){
  if(max(-log10(pval_gen_time_lm[FourP$snp_info$ChrID=="III",i]))>cutoff_gen_time){
    min_p_time_pos_chr3 <- c(min_p_time_pos_chr3,FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="III"][which.min(pval_gen_time_lm[FourP$snp_info$ChrID=="III",i])])
    loc <- gsub(".*:","",fpkm$locus[fpkm.selc][i])
    middle <- (as.numeric(sub("(^[^-]+)-.*", "\\1", loc)) + as.numeric(sub('.+-(.+)', '\\1', loc))) / 2 / 1000000
    gene_chr <- strsplit(as.character(fpkm$locus[fpkm.selc][i]), "[:]")[[1]][1]
    ifelse(abs(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="III"][which.min(pval_gen_time_lm[FourP$snp_info$ChrID=="III",i])] - middle) < 2 & gene_chr =="III",ct_gen_time_chr3 <- c(ct_gen_time_chr3,"cis"),ct_gen_time_chr3 <- c(ct_gen_time_chr3,"trans"))
  }
}

for(i in 1:ncol(pval_gen_time_lm)){
  if(max(-log10(pval_gen_time_lm[FourP$snp_info$ChrID=="IV",i]))>cutoff_gen_time){
    min_p_time_pos_chr4 <- c(min_p_time_pos_chr4,FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="IV"][which.min(pval_gen_time_lm[FourP$snp_info$ChrID=="IV",i])])
    loc <- gsub(".*:","",fpkm$locus[fpkm.selc][i])
    middle <- (as.numeric(sub("(^[^-]+)-.*", "\\1", loc)) + as.numeric(sub('.+-(.+)', '\\1', loc))) / 2 / 1000000
    gene_chr <- strsplit(as.character(fpkm$locus[fpkm.selc][i]), "[:]")[[1]][1]
    ifelse(abs(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="IV"][which.min(pval_gen_time_lm[FourP$snp_info$ChrID=="IV",i])] - middle) < 2 & gene_chr =="IV",ct_gen_time_chr4 <- c(ct_gen_time_chr4,"cis"),ct_gen_time_chr4 <- c(ct_gen_time_chr4,"trans"))
  }
}

for(i in 1:ncol(pval_gen_time_lm)){
  if(max(-log10(pval_gen_time_lm[FourP$snp_info$ChrID=="V",i]))>cutoff_gen_time){
    min_p_time_pos_chr5 <- c(min_p_time_pos_chr5,FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="V"][which.min(pval_gen_time_lm[FourP$snp_info$ChrID=="V",i])])
    loc <- gsub(".*:","",fpkm$locus[fpkm.selc][i])
    middle <- (as.numeric(sub("(^[^-]+)-.*", "\\1", loc)) + as.numeric(sub('.+-(.+)', '\\1', loc))) / 2 / 1000000
    gene_chr <- strsplit(as.character(fpkm$locus[fpkm.selc][i]), "[:]")[[1]][1]
    ifelse(abs(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="V"][which.min(pval_gen_time_lm[FourP$snp_info$ChrID=="V",i])] - middle) < 2 & gene_chr =="V",ct_gen_time_chr5 <- c(ct_gen_time_chr5,"cis"),ct_gen_time_chr5 <- c(ct_gen_time_chr5,"trans"))
  }
}

for(i in 1:ncol(pval_gen_time_lm)){
  if(max(-log10(pval_gen_time_lm[FourP$snp_info$ChrID=="X",i]))>cutoff_gen_time){
    min_p_time_pos_chrx <- c(min_p_time_pos_chrx,FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="X"][which.min(pval_gen_time_lm[FourP$snp_info$ChrID=="X",i])])
    loc <- gsub(".*:","",fpkm$locus[fpkm.selc][i])
    middle <- (as.numeric(sub("(^[^-]+)-.*", "\\1", loc)) + as.numeric(sub('.+-(.+)', '\\1', loc))) / 2 / 1000000
    gene_chr <- strsplit(as.character(fpkm$locus[fpkm.selc][i]), "[:]")[[1]][1]
    ifelse(abs(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="X"][which.min(pval_gen_time_lm[FourP$snp_info$ChrID=="X",i])] - middle) < 2 & gene_chr =="X",ct_gen_time_chrx <- c(ct_gen_time_chrx,"cis"),ct_gen_time_chrx <- c(ct_gen_time_chrx,"trans"))
  }
}

min_p_time_pos_chr1 <- data.frame(marker = min_p_time_pos_chr1,ct = ct_gen_time_chr1,chr = rep("I",length(min_p_time_pos_chr1)))
min_p_time_pos_chr2 <- data.frame(marker = min_p_time_pos_chr2,ct = ct_gen_time_chr2,chr = rep("II",length(min_p_time_pos_chr2)))
min_p_time_pos_chr3 <- data.frame(marker = min_p_time_pos_chr3,ct = ct_gen_time_chr3,chr = rep("III",length(min_p_time_pos_chr3)))
min_p_time_pos_chr4 <- data.frame(marker = min_p_time_pos_chr4,ct = ct_gen_time_chr4,chr = rep("IV",length(min_p_time_pos_chr4)))
min_p_time_pos_chr5 <- data.frame(marker = min_p_time_pos_chr5,ct = ct_gen_time_chr5,chr = rep("V",length(min_p_time_pos_chr5)))
min_p_time_pos_chrx <- data.frame(marker = min_p_time_pos_chrx,ct = ct_gen_time_chrx,chr = rep("X",length(min_p_time_pos_chrx)))

min_p_max_per_chr_time <- rbind(min_p_time_pos_chr1,min_p_time_pos_chr2,min_p_time_pos_chr3,min_p_time_pos_chr4,min_p_time_pos_chr5,min_p_time_pos_chrx)


# Now do the same for the genotype variable in the interaction model.

cutoff_gen_int <- 4.52

min_p_int_pos_chr1 <- c()
min_p_int_pos_chr2 <- c()
min_p_int_pos_chr3 <- c()
min_p_int_pos_chr4 <- c()
min_p_int_pos_chr5 <- c()
min_p_int_pos_chrx <- c()

ct_gen_int_chr1 <- c()
ct_gen_int_chr2 <- c()
ct_gen_int_chr3 <- c()
ct_gen_int_chr4 <- c()
ct_gen_int_chr5 <- c()
ct_gen_int_chrx <- c()

for(i in 1:ncol(pval_gen_time_intlm)){
  if(max(-log10(pval_gen_time_intlm[FourP$snp_info$ChrID=="I",i]))>cutoff_gen_int){
    min_p_int_pos_chr1 <- c(min_p_int_pos_chr1,FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="I"][which.min(pval_gen_time_intlm[FourP$snp_info$ChrID=="I",i])])
    loc <- gsub(".*:","",fpkm$locus[fpkm.selc][i])
    middle <- (as.numeric(sub("(^[^-]+)-.*", "\\1", loc)) + as.numeric(sub('.+-(.+)', '\\1', loc))) / 2 / 1000000
    gene_chr <- strsplit(as.character(fpkm$locus[fpkm.selc][i]), "[:]")[[1]][1]
    ifelse(abs(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="I"][which.min(pval_gen_time_intlm[FourP$snp_info$ChrID=="I",i])] - middle) < 2 & gene_chr =="I",ct_gen_int_chr1 <- c(ct_gen_int_chr1,"cis"),ct_gen_int_chr1 <- c(ct_gen_int_chr1,"trans"))
  }
}

for(i in 1:ncol(pval_gen_time_intlm)){
  if(max(-log10(pval_gen_time_intlm[FourP$snp_info$ChrID=="II",i]))>cutoff_gen_int){
    min_p_int_pos_chr2 <- c(min_p_int_pos_chr2,FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="II"][which.min(pval_gen_time_intlm[FourP$snp_info$ChrID=="II",i])])
    loc <- gsub(".*:","",fpkm$locus[fpkm.selc][i])
    middle <- (as.numeric(sub("(^[^-]+)-.*", "\\1", loc)) + as.numeric(sub('.+-(.+)', '\\1', loc))) / 2 / 1000000
    gene_chr <- strsplit(as.character(fpkm$locus[fpkm.selc][i]), "[:]")[[1]][1]
    ifelse(abs(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="II"][which.min(pval_gen_time_intlm[FourP$snp_info$ChrID=="II",i])] - middle) < 2 & gene_chr =="II",ct_gen_int_chr2 <- c(ct_gen_int_chr2,"cis"),ct_gen_int_chr2 <- c(ct_gen_int_chr2,"trans"))
  }
}

for(i in 1:ncol(pval_gen_time_intlm)){
  if(max(-log10(pval_gen_time_intlm[FourP$snp_info$ChrID=="III",i]))>cutoff_gen_int){
    min_p_int_pos_chr3 <- c(min_p_int_pos_chr3,FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="III"][which.min(pval_gen_time_intlm[FourP$snp_info$ChrID=="III",i])])
    loc <- gsub(".*:","",fpkm$locus[fpkm.selc][i])
    middle <- (as.numeric(sub("(^[^-]+)-.*", "\\1", loc)) + as.numeric(sub('.+-(.+)', '\\1', loc))) / 2 / 1000000
    gene_chr <- strsplit(as.character(fpkm$locus[fpkm.selc][i]), "[:]")[[1]][1]
    ifelse(abs(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="III"][which.min(pval_gen_time_intlm[FourP$snp_info$ChrID=="III",i])] - middle) < 2 & gene_chr =="III",ct_gen_int_chr3 <- c(ct_gen_int_chr3,"cis"),ct_gen_int_chr3 <- c(ct_gen_int_chr3,"trans"))
  }
}

for(i in 1:ncol(pval_gen_time_intlm)){
  if(max(-log10(pval_gen_time_intlm[FourP$snp_info$ChrID=="IV",i]))>cutoff_gen_int){
    min_p_int_pos_chr4 <- c(min_p_int_pos_chr4,FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="IV"][which.min(pval_gen_time_intlm[FourP$snp_info$ChrID=="IV",i])])
    loc <- gsub(".*:","",fpkm$locus[fpkm.selc][i])
    middle <- (as.numeric(sub("(^[^-]+)-.*", "\\1", loc)) + as.numeric(sub('.+-(.+)', '\\1', loc))) / 2 / 1000000
    gene_chr <- strsplit(as.character(fpkm$locus[fpkm.selc][i]), "[:]")[[1]][1]
    ifelse(abs(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="IV"][which.min(pval_gen_time_intlm[FourP$snp_info$ChrID=="IV",i])] - middle) < 2 & gene_chr =="IV",ct_gen_int_chr4 <- c(ct_gen_int_chr4,"cis"),ct_gen_int_chr4 <- c(ct_gen_int_chr4,"trans"))
  }
}

for(i in 1:ncol(pval_gen_time_intlm)){
  if(max(-log10(pval_gen_time_intlm[FourP$snp_info$ChrID=="V",i]))>cutoff_gen_int){
    min_p_int_pos_chr5 <- c(min_p_int_pos_chr5,FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="V"][which.min(pval_gen_time_intlm[FourP$snp_info$ChrID=="V",i])])
    loc <- gsub(".*:","",fpkm$locus[fpkm.selc][i])
    middle <- (as.numeric(sub("(^[^-]+)-.*", "\\1", loc)) + as.numeric(sub('.+-(.+)', '\\1', loc))) / 2 / 1000000
    gene_chr <- strsplit(as.character(fpkm$locus[fpkm.selc][i]), "[:]")[[1]][1]
    ifelse(abs(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="V"][which.min(pval_gen_time_intlm[FourP$snp_info$ChrID=="V",i])] - middle) < 2 & gene_chr =="V",ct_gen_int_chr5 <- c(ct_gen_int_chr5,"cis"),ct_gen_int_chr5 <- c(ct_gen_int_chr5,"trans"))
  }
}

for(i in 1:ncol(pval_gen_time_intlm)){
  if(max(-log10(pval_gen_time_intlm[FourP$snp_info$ChrID=="X",i]))>cutoff_gen_int){
    min_p_int_pos_chrx <- c(min_p_int_pos_chrx,FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="X"][which.min(pval_gen_time_intlm[FourP$snp_info$ChrID=="X",i])])
    loc <- gsub(".*:","",fpkm$locus[fpkm.selc][i])
    middle <- (as.numeric(sub("(^[^-]+)-.*", "\\1", loc)) + as.numeric(sub('.+-(.+)', '\\1', loc))) / 2 / 1000000
    gene_chr <- strsplit(as.character(fpkm$locus[fpkm.selc][i]), "[:]")[[1]][1]
    ifelse(abs(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="X"][which.min(pval_gen_time_intlm[FourP$snp_info$ChrID=="X",i])] - middle) < 2 & gene_chr =="X",ct_gen_int_chrx <- c(ct_gen_int_chrx,"cis"),ct_gen_int_chrx <- c(ct_gen_int_chrx,"trans"))
  }
}

min_p_int_pos_chr1 <- data.frame(marker = min_p_int_pos_chr1,ct = ct_gen_int_chr1,chr = rep("I",length(min_p_int_pos_chr1)))
min_p_int_pos_chr2 <- data.frame(marker = min_p_int_pos_chr2,ct = ct_gen_int_chr2,chr = rep("II",length(min_p_int_pos_chr2)))
min_p_int_pos_chr3 <- data.frame(marker = min_p_int_pos_chr3,ct = ct_gen_int_chr3,chr = rep("III",length(min_p_int_pos_chr3)))
min_p_int_pos_chr4 <- data.frame(marker = min_p_int_pos_chr4,ct = ct_gen_int_chr4,chr = rep("IV",length(min_p_int_pos_chr4)))
min_p_int_pos_chr5 <- data.frame(marker = min_p_int_pos_chr5,ct = ct_gen_int_chr5,chr = rep("V",length(min_p_int_pos_chr5)))
min_p_int_pos_chrx <- data.frame(marker = min_p_int_pos_chrx,ct = ct_gen_int_chrx,chr = rep("X",length(min_p_int_pos_chrx)))

min_p_max_per_chr_int <- rbind(min_p_int_pos_chr1,min_p_int_pos_chr2,min_p_int_pos_chr3,min_p_int_pos_chr4,min_p_int_pos_chr5,min_p_int_pos_chrx)

# Now do the same for the model without development, with the 0.1 FDR cutoff (earlier we used 0.05).

cutoff_notime01 <- 3.52

min_p_notime_01_pos_chr1 <- c()
min_p_notime_01_pos_chr2 <- c()
min_p_notime_01_pos_chr3 <- c()
min_p_notime_01_pos_chr4 <- c()
min_p_notime_01_pos_chr5 <- c()
min_p_notime_01_pos_chrx <- c()

ct_gen_notime01_chr1 <- c()
ct_gen_notime01_chr2 <- c()
ct_gen_notime01_chr3 <- c()
ct_gen_notime01_chr4 <- c()
ct_gen_notime01_chr5 <- c()
ct_gen_notime01_chrx <- c()

hotspot01 <- c()
#saveRDS(hotspot,"hotspot_index.RDS")

for(i in 1:ncol(pval_gen_notime_lm)){
  if(max(-log10(pval_gen_notime_lm[FourP$snp_info$ChrID=="I",i]))>cutoff_notime01){
    min_p_notime_01_pos_chr1 <- c(min_p_notime_01_pos_chr1,FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="I"][which.min(pval_gen_notime_lm[FourP$snp_info$ChrID=="I",i])])
    loc <- gsub(".*:","",fpkm$locus[fpkm.selc][i])
    middle <- (as.numeric(sub("(^[^-]+)-.*", "\\1", loc)) + as.numeric(sub('.+-(.+)', '\\1', loc))) / 2 / 1000000
    gene_chr <- strsplit(as.character(fpkm$locus[fpkm.selc][i]), "[:]")[[1]][1]
    ifelse(abs(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="I"][which.min(pval_gen_notime_lm[FourP$snp_info$ChrID=="I",i])] - middle) < 2 & gene_chr =="I",ct_gen_notime01_chr1 <- c(ct_gen_notime01_chr1,"cis"),ct_gen_notime01_chr1 <- c(ct_gen_notime01_chr1,"trans"))
  }
}

for(i in 1:ncol(pval_gen_notime_lm)){
  if(max(-log10(pval_gen_notime_lm[FourP$snp_info$ChrID=="II",i]))>cutoff_notime01){
    min_p_notime_01_pos_chr2 <- c(min_p_notime_01_pos_chr2,FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="II"][which.min(pval_gen_notime_lm[FourP$snp_info$ChrID=="II",i])])
    loc <- gsub(".*:","",fpkm$locus[fpkm.selc][i])
    middle <- (as.numeric(sub("(^[^-]+)-.*", "\\1", loc)) + as.numeric(sub('.+-(.+)', '\\1', loc))) / 2 / 1000000
    gene_chr <- strsplit(as.character(fpkm$locus[fpkm.selc][i]), "[:]")[[1]][1]
    ifelse(abs(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="II"][which.min(pval_gen_notime_lm[FourP$snp_info$ChrID=="II",i])] - middle) < 2 & gene_chr =="II",ct_gen_notime01_chr2 <- c(ct_gen_notime01_chr2,"cis"),ct_gen_notime01_chr2 <- c(ct_gen_notime01_chr2,"trans"))
  }
}

for(i in 1:ncol(pval_gen_notime_lm)){
  if(max(-log10(pval_gen_notime_lm[FourP$snp_info$ChrID=="III",i]))>cutoff_notime01){
    min_p_notime_01_pos_chr3 <- c(min_p_notime_01_pos_chr3,FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="III"][which.min(pval_gen_notime_lm[FourP$snp_info$ChrID=="III",i])])
    loc <- gsub(".*:","",fpkm$locus[fpkm.selc][i])
    middle <- (as.numeric(sub("(^[^-]+)-.*", "\\1", loc)) + as.numeric(sub('.+-(.+)', '\\1', loc))) / 2 / 1000000
    gene_chr <- strsplit(as.character(fpkm$locus[fpkm.selc][i]), "[:]")[[1]][1]
    ifelse(abs(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="III"][which.min(pval_gen_notime_lm[FourP$snp_info$ChrID=="III",i])] - middle) < 2 & gene_chr =="III",ct_gen_notime01_chr3 <- c(ct_gen_notime01_chr3,"cis"),ct_gen_notime01_chr3 <- c(ct_gen_notime01_chr3,"trans"))
  }
}

for(i in 1:ncol(pval_gen_notime_lm)){
  if(max(-log10(pval_gen_notime_lm[FourP$snp_info$ChrID=="IV",i]))>cutoff_notime01){
    min_p_notime_01_pos_chr4 <- c(min_p_notime_01_pos_chr4,FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="IV"][which.min(pval_gen_notime_lm[FourP$snp_info$ChrID=="IV",i])])
    loc <- gsub(".*:","",fpkm$locus[fpkm.selc][i])
    middle <- (as.numeric(sub("(^[^-]+)-.*", "\\1", loc)) + as.numeric(sub('.+-(.+)', '\\1', loc))) / 2 / 1000000
    gene_chr <- strsplit(as.character(fpkm$locus[fpkm.selc][i]), "[:]")[[1]][1]
    ifelse(abs(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="IV"][which.min(pval_gen_notime_lm[FourP$snp_info$ChrID=="IV",i])] - middle) < 2 & gene_chr =="IV",ct_gen_notime01_chr4 <- c(ct_gen_notime01_chr4,"cis"),ct_gen_notime01_chr4 <- c(ct_gen_notime01_chr4,"trans"))
  }
}

for(i in 1:ncol(pval_gen_notime_lm)){
  if(max(-log10(pval_gen_notime_lm[FourP$snp_info$ChrID=="V",i]))>cutoff_notime01){
    min_p_notime_01_pos_chr5 <- c(min_p_notime_01_pos_chr5,FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="V"][which.min(pval_gen_notime_lm[FourP$snp_info$ChrID=="V",i])])
    loc <- gsub(".*:","",fpkm$locus[fpkm.selc][i])
    middle <- (as.numeric(sub("(^[^-]+)-.*", "\\1", loc)) + as.numeric(sub('.+-(.+)', '\\1', loc))) / 2 / 1000000
    gene_chr <- strsplit(as.character(fpkm$locus[fpkm.selc][i]), "[:]")[[1]][1]
    ifelse(abs(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="V"][which.min(pval_gen_notime_lm[FourP$snp_info$ChrID=="V",i])] - middle) < 2 & gene_chr =="V",ct_gen_notime01_chr5 <- c(ct_gen_notime01_chr5,"cis"),ct_gen_notime01_chr5 <- c(ct_gen_notime01_chr5,"trans"))
  }
}

for(i in 1:ncol(pval_gen_notime_lm)){
  if(max(-log10(pval_gen_notime_lm[FourP$snp_info$ChrID=="X",i]))>cutoff_notime01){
    min_p_notime_01_pos_chrx <- c(min_p_notime_01_pos_chrx,FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="X"][which.min(pval_gen_notime_lm[FourP$snp_info$ChrID=="X",i])])
    loc <- gsub(".*:","",fpkm$locus[fpkm.selc][i])
    middle <- (as.numeric(sub("(^[^-]+)-.*", "\\1", loc)) + as.numeric(sub('.+-(.+)', '\\1', loc))) / 2 / 1000000
    gene_chr <- strsplit(as.character(fpkm$locus[fpkm.selc][i]), "[:]")[[1]][1]
    ifelse(abs(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="X"][which.min(pval_gen_notime_lm[FourP$snp_info$ChrID=="X",i])] - middle) < 2 & gene_chr =="X",ct_gen_notime01_chrx <- c(ct_gen_notime01_chrx,"cis"),ct_gen_notime01_chrx <- c(ct_gen_notime01_chrx,"trans"))
    if(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="X"][which.min(pval_gen_notime_lm[FourP$snp_info$ChrID=="X",i])]==3.404162){
      hotspot01 <- c(hotspot01,i)
    }
  }
}

min_p_notime_01_pos_chr1 <- data.frame(marker = min_p_notime_01_pos_chr1,ct = ct_gen_notime01_chr1,chr = rep("I",length(min_p_notime_01_pos_chr1)))
min_p_notime_01_pos_chr2 <- data.frame(marker = min_p_notime_01_pos_chr2,ct = ct_gen_notime01_chr2,chr = rep("II",length(min_p_notime_01_pos_chr2)))
min_p_notime_01_pos_chr3 <- data.frame(marker = min_p_notime_01_pos_chr3,ct = ct_gen_notime01_chr3,chr = rep("III",length(min_p_notime_01_pos_chr3)))
min_p_notime_01_pos_chr4 <- data.frame(marker = min_p_notime_01_pos_chr4,ct = ct_gen_notime01_chr4,chr = rep("IV",length(min_p_notime_01_pos_chr4)))
min_p_notime_01_pos_chr5 <- data.frame(marker = min_p_notime_01_pos_chr5,ct = ct_gen_notime01_chr5,chr = rep("V",length(min_p_notime_01_pos_chr5)))
min_p_notime_01_pos_chrx <- data.frame(marker = min_p_notime_01_pos_chrx,ct = ct_gen_notime01_chrx,chr = rep("X",length(min_p_notime_01_pos_chrx)))

min_p_max_per_chr01 <- rbind(min_p_notime_01_pos_chr1,min_p_notime_01_pos_chr2,min_p_notime_01_pos_chr3,min_p_notime_01_pos_chr4,min_p_notime_01_pos_chr5,min_p_notime_01_pos_chrx)

# Do the same for the genotype variable in the additive model.

cutoff_gen_time01 <- 3.95

min_p_time_01_pos_chr1 <- c()
min_p_time_01_pos_chr2 <- c()
min_p_time_01_pos_chr3 <- c()
min_p_time_01_pos_chr4 <- c()
min_p_time_01_pos_chr5 <- c()
min_p_time_01_pos_chrx <- c()

ct_gen_time01_chr1 <- c()
ct_gen_time01_chr2 <- c()
ct_gen_time01_chr3 <- c()
ct_gen_time01_chr4 <- c()
ct_gen_time01_chr5 <- c()
ct_gen_time01_chrx <- c()

for(i in 1:ncol(pval_gen_time_lm)){
  if(max(-log10(pval_gen_time_lm[FourP$snp_info$ChrID=="I",i]))>cutoff_gen_time01){
    min_p_time_01_pos_chr1 <- c(min_p_time_01_pos_chr1,FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="I"][which.min(pval_gen_time_lm[FourP$snp_info$ChrID=="I",i])])
    loc <- gsub(".*:","",fpkm$locus[fpkm.selc][i])
    middle <- (as.numeric(sub("(^[^-]+)-.*", "\\1", loc)) + as.numeric(sub('.+-(.+)', '\\1', loc))) / 2 / 1000000
    gene_chr <- strsplit(as.character(fpkm$locus[fpkm.selc][i]), "[:]")[[1]][1]
    ifelse(abs(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="I"][which.min(pval_gen_time_lm[FourP$snp_info$ChrID=="I",i])] - middle) < 2 & gene_chr =="I",ct_gen_time01_chr1 <- c(ct_gen_time01_chr1,"cis"),ct_gen_time01_chr1 <- c(ct_gen_time01_chr1,"trans"))
  }
}

for(i in 1:ncol(pval_gen_time_lm)){
  if(max(-log10(pval_gen_time_lm[FourP$snp_info$ChrID=="II",i]))>cutoff_gen_time01){
    min_p_time_01_pos_chr2 <- c(min_p_time_01_pos_chr2,FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="II"][which.min(pval_gen_time_lm[FourP$snp_info$ChrID=="II",i])])
    loc <- gsub(".*:","",fpkm$locus[fpkm.selc][i])
    middle <- (as.numeric(sub("(^[^-]+)-.*", "\\1", loc)) + as.numeric(sub('.+-(.+)', '\\1', loc))) / 2 / 1000000
    gene_chr <- strsplit(as.character(fpkm$locus[fpkm.selc][i]), "[:]")[[1]][1]
    ifelse(abs(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="II"][which.min(pval_gen_time_lm[FourP$snp_info$ChrID=="II",i])] - middle) < 2 & gene_chr =="II",ct_gen_time01_chr2 <- c(ct_gen_time01_chr2,"cis"),ct_gen_time01_chr2 <- c(ct_gen_time01_chr2,"trans"))
  }
}

for(i in 1:ncol(pval_gen_time_lm)){
  if(max(-log10(pval_gen_time_lm[FourP$snp_info$ChrID=="III",i]))>cutoff_gen_time01){
    min_p_time_01_pos_chr3 <- c(min_p_time_01_pos_chr3,FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="III"][which.min(pval_gen_time_lm[FourP$snp_info$ChrID=="III",i])])
    loc <- gsub(".*:","",fpkm$locus[fpkm.selc][i])
    middle <- (as.numeric(sub("(^[^-]+)-.*", "\\1", loc)) + as.numeric(sub('.+-(.+)', '\\1', loc))) / 2 / 1000000
    gene_chr <- strsplit(as.character(fpkm$locus[fpkm.selc][i]), "[:]")[[1]][1]
    ifelse(abs(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="III"][which.min(pval_gen_time_lm[FourP$snp_info$ChrID=="III",i])] - middle) < 2 & gene_chr =="III",ct_gen_time01_chr3 <- c(ct_gen_time01_chr3,"cis"),ct_gen_time01_chr3 <- c(ct_gen_time01_chr3,"trans"))
  }
}

for(i in 1:ncol(pval_gen_time_lm)){
  if(max(-log10(pval_gen_time_lm[FourP$snp_info$ChrID=="IV",i]))>cutoff_gen_time01){
    min_p_time_01_pos_chr4 <- c(min_p_time_01_pos_chr4,FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="IV"][which.min(pval_gen_time_lm[FourP$snp_info$ChrID=="IV",i])])
    loc <- gsub(".*:","",fpkm$locus[fpkm.selc][i])
    middle <- (as.numeric(sub("(^[^-]+)-.*", "\\1", loc)) + as.numeric(sub('.+-(.+)', '\\1', loc))) / 2 / 1000000
    gene_chr <- strsplit(as.character(fpkm$locus[fpkm.selc][i]), "[:]")[[1]][1]
    ifelse(abs(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="IV"][which.min(pval_gen_time_lm[FourP$snp_info$ChrID=="IV",i])] - middle) < 2 & gene_chr =="IV",ct_gen_time01_chr4 <- c(ct_gen_time01_chr4,"cis"),ct_gen_time01_chr4 <- c(ct_gen_time01_chr4,"trans"))
  }
}

for(i in 1:ncol(pval_gen_time_lm)){
  if(max(-log10(pval_gen_time_lm[FourP$snp_info$ChrID=="V",i]))>cutoff_gen_time01){
    min_p_time_01_pos_chr5 <- c(min_p_time_01_pos_chr5,FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="V"][which.min(pval_gen_time_lm[FourP$snp_info$ChrID=="V",i])])
    loc <- gsub(".*:","",fpkm$locus[fpkm.selc][i])
    middle <- (as.numeric(sub("(^[^-]+)-.*", "\\1", loc)) + as.numeric(sub('.+-(.+)', '\\1', loc))) / 2 / 1000000
    gene_chr <- strsplit(as.character(fpkm$locus[fpkm.selc][i]), "[:]")[[1]][1]
    ifelse(abs(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="V"][which.min(pval_gen_time_lm[FourP$snp_info$ChrID=="V",i])] - middle) < 2 & gene_chr =="V",ct_gen_time01_chr5 <- c(ct_gen_time01_chr5,"cis"),ct_gen_time01_chr5 <- c(ct_gen_time01_chr5,"trans"))
  }
}

for(i in 1:ncol(pval_gen_time_lm)){
  if(max(-log10(pval_gen_time_lm[FourP$snp_info$ChrID=="X",i]))>cutoff_gen_time01){
    min_p_time_01_pos_chrx <- c(min_p_time_01_pos_chrx,FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="X"][which.min(pval_gen_time_lm[FourP$snp_info$ChrID=="X",i])])
    loc <- gsub(".*:","",fpkm$locus[fpkm.selc][i])
    middle <- (as.numeric(sub("(^[^-]+)-.*", "\\1", loc)) + as.numeric(sub('.+-(.+)', '\\1', loc))) / 2 / 1000000
    gene_chr <- strsplit(as.character(fpkm$locus[fpkm.selc][i]), "[:]")[[1]][1]
    ifelse(abs(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="X"][which.min(pval_gen_time_lm[FourP$snp_info$ChrID=="X",i])] - middle) < 2 & gene_chr =="X",ct_gen_time01_chrx <- c(ct_gen_time01_chrx,"cis"),ct_gen_time01_chrx <- c(ct_gen_time01_chrx,"trans"))
  }
}

min_p_time_01_pos_chr1 <- data.frame(marker = min_p_time_01_pos_chr1,ct = ct_gen_time01_chr1,chr = rep("I",length(min_p_time_01_pos_chr1)))
min_p_time_01_pos_chr2 <- data.frame(marker = min_p_time_01_pos_chr2,ct = ct_gen_time01_chr2,chr = rep("II",length(min_p_time_01_pos_chr2)))
min_p_time_01_pos_chr3 <- data.frame(marker = min_p_time_01_pos_chr3,ct = ct_gen_time01_chr3,chr = rep("III",length(min_p_time_01_pos_chr3)))
min_p_time_01_pos_chr4 <- data.frame(marker = min_p_time_01_pos_chr4,ct = ct_gen_time01_chr4,chr = rep("IV",length(min_p_time_01_pos_chr4)))
min_p_time_01_pos_chr5 <- data.frame(marker = min_p_time_01_pos_chr5,ct = ct_gen_time01_chr5,chr = rep("V",length(min_p_time_01_pos_chr5)))
min_p_time_01_pos_chrx <- data.frame(marker = min_p_time_01_pos_chrx,ct = ct_gen_time01_chrx,chr = rep("X",length(min_p_time_01_pos_chrx)))

min_p_max_per_chr_time01 <- rbind(min_p_time_01_pos_chr1,min_p_time_01_pos_chr2,min_p_time_01_pos_chr3,min_p_time_01_pos_chr4,min_p_time_01_pos_chr5,min_p_time_01_pos_chrx)

# Do the same for the genotype variable in the interaction model.

cutoff_gen_int01 <- 3.91

min_p_int01_pos_chr1 <- c()
min_p_int01_pos_chr2 <- c()
min_p_int01_pos_chr3 <- c()
min_p_int01_pos_chr4 <- c()
min_p_int01_pos_chr5 <- c()
min_p_int01_pos_chrx <- c()

ct_gen_int01_chr1 <- c()
ct_gen_int01_chr2 <- c()
ct_gen_int01_chr3 <- c()
ct_gen_int01_chr4 <- c()
ct_gen_int01_chr5 <- c()
ct_gen_int01_chrx <- c()

for(i in 1:ncol(pval_gen_time_intlm)){
  if(max(-log10(pval_gen_time_intlm[FourP$snp_info$ChrID=="I",i]))>cutoff_gen_int01){
    min_p_int01_pos_chr1 <- c(min_p_int01_pos_chr1,FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="I"][which.min(pval_gen_time_intlm[FourP$snp_info$ChrID=="I",i])])
    loc <- gsub(".*:","",fpkm$locus[fpkm.selc][i])
    middle <- (as.numeric(sub("(^[^-]+)-.*", "\\1", loc)) + as.numeric(sub('.+-(.+)', '\\1', loc))) / 2 / 1000000
    gene_chr <- strsplit(as.character(fpkm$locus[fpkm.selc][i]), "[:]")[[1]][1]
    ifelse(abs(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="I"][which.min(pval_gen_time_intlm[FourP$snp_info$ChrID=="I",i])] - middle) < 2 & gene_chr =="I",ct_gen_int01_chr1 <- c(ct_gen_int01_chr1,"cis"),ct_gen_int01_chr1 <- c(ct_gen_int01_chr1,"trans"))
  }
}

for(i in 1:ncol(pval_gen_time_intlm)){
  if(max(-log10(pval_gen_time_intlm[FourP$snp_info$ChrID=="II",i]))>cutoff_gen_int01){
    min_p_int01_pos_chr2 <- c(min_p_int01_pos_chr2,FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="II"][which.min(pval_gen_time_intlm[FourP$snp_info$ChrID=="II",i])])
    loc <- gsub(".*:","",fpkm$locus[fpkm.selc][i])
    middle <- (as.numeric(sub("(^[^-]+)-.*", "\\1", loc)) + as.numeric(sub('.+-(.+)', '\\1', loc))) / 2 / 1000000
    gene_chr <- strsplit(as.character(fpkm$locus[fpkm.selc][i]), "[:]")[[1]][1]
    ifelse(abs(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="II"][which.min(pval_gen_time_intlm[FourP$snp_info$ChrID=="II",i])] - middle) < 2 & gene_chr =="II",ct_gen_int01_chr2 <- c(ct_gen_int01_chr2,"cis"),ct_gen_int01_chr2 <- c(ct_gen_int01_chr2,"trans"))
  }
}

for(i in 1:ncol(pval_gen_time_intlm)){
  if(max(-log10(pval_gen_time_intlm[FourP$snp_info$ChrID=="III",i]))>cutoff_gen_int01){
    min_p_int01_pos_chr3 <- c(min_p_int01_pos_chr3,FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="III"][which.min(pval_gen_time_intlm[FourP$snp_info$ChrID=="III",i])])
    loc <- gsub(".*:","",fpkm$locus[fpkm.selc][i])
    middle <- (as.numeric(sub("(^[^-]+)-.*", "\\1", loc)) + as.numeric(sub('.+-(.+)', '\\1', loc))) / 2 / 1000000
    gene_chr <- strsplit(as.character(fpkm$locus[fpkm.selc][i]), "[:]")[[1]][1]
    ifelse(abs(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="III"][which.min(pval_gen_time_intlm[FourP$snp_info$ChrID=="III",i])] - middle) < 2 & gene_chr =="III",ct_gen_int01_chr3 <- c(ct_gen_int01_chr3,"cis"),ct_gen_int01_chr3 <- c(ct_gen_int01_chr3,"trans"))
  }
}

for(i in 1:ncol(pval_gen_time_intlm)){
  if(max(-log10(pval_gen_time_intlm[FourP$snp_info$ChrID=="IV",i]))>cutoff_gen_int01){
    min_p_int01_pos_chr4 <- c(min_p_int01_pos_chr4,FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="IV"][which.min(pval_gen_time_intlm[FourP$snp_info$ChrID=="IV",i])])
    loc <- gsub(".*:","",fpkm$locus[fpkm.selc][i])
    middle <- (as.numeric(sub("(^[^-]+)-.*", "\\1", loc)) + as.numeric(sub('.+-(.+)', '\\1', loc))) / 2 / 1000000
    gene_chr <- strsplit(as.character(fpkm$locus[fpkm.selc][i]), "[:]")[[1]][1]
    ifelse(abs(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="IV"][which.min(pval_gen_time_intlm[FourP$snp_info$ChrID=="IV",i])] - middle) < 2 & gene_chr =="IV",ct_gen_int01_chr4 <- c(ct_gen_int01_chr4,"cis"),ct_gen_int01_chr4 <- c(ct_gen_int01_chr4,"trans"))
  }
}

for(i in 1:ncol(pval_gen_time_intlm)){
  if(max(-log10(pval_gen_time_intlm[FourP$snp_info$ChrID=="V",i]))>cutoff_gen_int01){
    min_p_int01_pos_chr5 <- c(min_p_int01_pos_chr5,FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="V"][which.min(pval_gen_time_intlm[FourP$snp_info$ChrID=="V",i])])
    loc <- gsub(".*:","",fpkm$locus[fpkm.selc][i])
    middle <- (as.numeric(sub("(^[^-]+)-.*", "\\1", loc)) + as.numeric(sub('.+-(.+)', '\\1', loc))) / 2 / 1000000
    gene_chr <- strsplit(as.character(fpkm$locus[fpkm.selc][i]), "[:]")[[1]][1]
    ifelse(abs(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="V"][which.min(pval_gen_time_intlm[FourP$snp_info$ChrID=="V",i])] - middle) < 2 & gene_chr =="V",ct_gen_int01_chr5 <- c(ct_gen_int01_chr5,"cis"),ct_gen_int01_chr5 <- c(ct_gen_int01_chr5,"trans"))
  }
}

for(i in 1:ncol(pval_gen_time_intlm)){
  if(max(-log10(pval_gen_time_intlm[FourP$snp_info$ChrID=="X",i]))>cutoff_gen_int01){
    min_p_int01_pos_chrx <- c(min_p_int01_pos_chrx,FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="X"][which.min(pval_gen_time_intlm[FourP$snp_info$ChrID=="X",i])])
    loc <- gsub(".*:","",fpkm$locus[fpkm.selc][i])
    middle <- (as.numeric(sub("(^[^-]+)-.*", "\\1", loc)) + as.numeric(sub('.+-(.+)', '\\1', loc))) / 2 / 1000000
    gene_chr <- strsplit(as.character(fpkm$locus[fpkm.selc][i]), "[:]")[[1]][1]
    ifelse(abs(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="X"][which.min(pval_gen_time_intlm[FourP$snp_info$ChrID=="X",i])] - middle) < 2 & gene_chr =="X",ct_gen_int01_chrx <- c(ct_gen_int01_chrx,"cis"),ct_gen_int01_chrx <- c(ct_gen_int01_chrx,"trans"))
  }
}

min_p_int01_pos_chr1 <- data.frame(marker = min_p_int01_pos_chr1,ct = ct_gen_int01_chr1,chr = rep("I",length(min_p_int01_pos_chr1)))
min_p_int01_pos_chr2 <- data.frame(marker = min_p_int01_pos_chr2,ct = ct_gen_int01_chr2,chr = rep("II",length(min_p_int01_pos_chr2)))
min_p_int01_pos_chr3 <- data.frame(marker = min_p_int01_pos_chr3,ct = ct_gen_int01_chr3,chr = rep("III",length(min_p_int01_pos_chr3)))
min_p_int01_pos_chr4 <- data.frame(marker = min_p_int01_pos_chr4,ct = ct_gen_int01_chr4,chr = rep("IV",length(min_p_int01_pos_chr4)))
min_p_int01_pos_chr5 <- data.frame(marker = min_p_int01_pos_chr5,ct = ct_gen_int01_chr5,chr = rep("V",length(min_p_int01_pos_chr5)))
min_p_int01_pos_chrx <- data.frame(marker = min_p_int01_pos_chrx,ct = ct_gen_int01_chrx,chr = rep("X",length(min_p_int01_pos_chrx)))

min_p_max_per_chr_int01 <- rbind(min_p_int01_pos_chr1,min_p_int01_pos_chr2,min_p_int01_pos_chr3,min_p_int01_pos_chr4,min_p_int01_pos_chr5,min_p_int01_pos_chrx)

# Add a factor called model which specifies to which models the frames correspond, and one called sig which specifies the FDR.
min_p_max_per_chr <- min_p_max_per_chr %>% mutate(model = as.factor(rep("No time",nrow(min_p_max_per_chr))),sig = as.factor(rep(0.5,nrow(min_p_max_per_chr))))
min_p_max_per_chr01 <- min_p_max_per_chr01 %>% mutate(model = as.factor(rep("No time",nrow(min_p_max_per_chr01))),sig = as.factor(rep(1,nrow(min_p_max_per_chr01))))

min_p_max_per_chr_time <- min_p_max_per_chr_time %>% mutate(model = as.factor(rep("Time",nrow(min_p_max_per_chr_time))),sig = as.factor(rep(0.5,nrow(min_p_max_per_chr_time))))
min_p_max_per_chr_time01 <- min_p_max_per_chr_time01 %>% mutate(model = as.factor(rep("Time",nrow(min_p_max_per_chr_time01))),sig = as.factor(rep(1,nrow(min_p_max_per_chr_time01))))

min_p_max_per_chr_int <- min_p_max_per_chr_int %>% mutate(model = as.factor(rep("Interaction",nrow(min_p_max_per_chr_int))),sig = as.factor(rep(0.5,nrow(min_p_max_per_chr_int))))
min_p_max_per_chr_int01 <- min_p_max_per_chr_int01 %>% mutate(model = as.factor(rep("Interaction",nrow(min_p_max_per_chr_int01))),sig = as.factor(rep(1,nrow(min_p_max_per_chr_int01))))

# Bind the 0.05 FDR eqTLs and the 0.1 FDR eQTLs in their own respective frames.
hist_frame <- rbind(min_p_max_per_chr,min_p_max_per_chr_time,min_p_max_per_chr_int)
hist_frame01 <- rbind(min_p_max_per_chr01,min_p_max_per_chr_time01,min_p_max_per_chr_int01)


# Now collect in a similar manner the eQTLs with a significant interaction in the interaction model (FDR = 0.05).

cutoff_int_int <- 5.49#6

min_p_int_int_chr1 <- c()
min_p_int_int_chr2 <- c()
min_p_int_int_chr3 <- c()
min_p_int_int_chr4 <- c()
min_p_int_int_chr5 <- c()
min_p_int_int_chrx <- c()

ct_int_int_chr1 <- c()
ct_int_int_chr2 <- c()
ct_int_int_chr3 <- c()
ct_int_int_chr4 <- c()
ct_int_int_chr5 <- c()
ct_int_int_chrx <- c()



for(i in 1:ncol(pval_int_time_intlm)){
  if(max(-log10(pval_int_time_intlm[FourP$snp_info$ChrID=="I",i]))>cutoff_int_int){
    min_p_int_int_chr1 <- c(min_p_int_int_chr1,FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="I"][which.min(pval_int_time_intlm[FourP$snp_info$ChrID=="I",i])])
    loc <- gsub(".*:","",fpkm$locus[fpkm.selc][i])
    middle <- (as.numeric(sub("(^[^-]+)-.*", "\\1", loc)) + as.numeric(sub('.+-(.+)', '\\1', loc))) / 2 / 1000000
    gene_chr <- strsplit(as.character(fpkm$locus[fpkm.selc][i]), "[:]")[[1]][1]
    ifelse(abs(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="I"][which.min(pval_int_time_intlm[FourP$snp_info$ChrID=="I",i])] - middle) < 2 & gene_chr =="I",ct_int_int_chr1 <- c(ct_int_int_chr1,"cis"),ct_int_int_chr1 <- c(ct_int_int_chr1,"trans"))
  }
}

for(i in 1:ncol(pval_int_time_intlm)){
  if(max(-log10(pval_int_time_intlm[FourP$snp_info$ChrID=="II",i]))>cutoff_int_int){
    min_p_int_int_chr2 <- c(min_p_int_int_chr2,FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="II"][which.min(pval_int_time_intlm[FourP$snp_info$ChrID=="II",i])])
    loc <- gsub(".*:","",fpkm$locus[fpkm.selc][i])
    middle <- (as.numeric(sub("(^[^-]+)-.*", "\\1", loc)) + as.numeric(sub('.+-(.+)', '\\1', loc))) / 2 / 1000000
    gene_chr <- strsplit(as.character(fpkm$locus[fpkm.selc][i]), "[:]")[[1]][1]
    ifelse(abs(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="II"][which.min(pval_int_time_intlm[FourP$snp_info$ChrID=="II",i])] - middle) < 2 & gene_chr =="II",ct_int_int_chr2 <- c(ct_int_int_chr2,"cis"),ct_int_int_chr2 <- c(ct_int_int_chr2,"trans"))
  }
}

for(i in 1:ncol(pval_int_time_intlm)){
  if(max(-log10(pval_int_time_intlm[FourP$snp_info$ChrID=="III",i]))>cutoff_int_int){
    min_p_int_int_chr3 <- c(min_p_int_int_chr3,FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="III"][which.min(pval_int_time_intlm[FourP$snp_info$ChrID=="III",i])])
    loc <- gsub(".*:","",fpkm$locus[fpkm.selc][i])
    middle <- (as.numeric(sub("(^[^-]+)-.*", "\\1", loc)) + as.numeric(sub('.+-(.+)', '\\1', loc))) / 2 / 1000000
    gene_chr <- strsplit(as.character(fpkm$locus[fpkm.selc][i]), "[:]")[[1]][1]
    ifelse(abs(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="III"][which.min(pval_int_time_intlm[FourP$snp_info$ChrID=="III",i])] - middle) < 2 & gene_chr =="III",ct_int_int_chr3 <- c(ct_int_int_chr3,"cis"),ct_int_int_chr3 <- c(ct_int_int_chr3,"trans"))
  }
}

for(i in 1:ncol(pval_int_time_intlm)){
  if(max(-log10(pval_int_time_intlm[FourP$snp_info$ChrID=="IV",i]))>cutoff_int_int){
    min_p_int_int_chr4 <- c(min_p_int_int_chr4,FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="IV"][which.min(pval_int_time_intlm[FourP$snp_info$ChrID=="IV",i])])
    loc <- gsub(".*:","",fpkm$locus[fpkm.selc][i])
    middle <- (as.numeric(sub("(^[^-]+)-.*", "\\1", loc)) + as.numeric(sub('.+-(.+)', '\\1', loc))) / 2 / 1000000
    gene_chr <- strsplit(as.character(fpkm$locus[fpkm.selc][i]), "[:]")[[1]][1]
    ifelse(abs(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="IV"][which.min(pval_int_time_intlm[FourP$snp_info$ChrID=="IV",i])] - middle) < 2 & gene_chr =="IV",ct_int_int_chr4 <- c(ct_int_int_chr4,"cis"),ct_int_int_chr4 <- c(ct_int_int_chr4,"trans"))
  }
}

for(i in 1:ncol(pval_int_time_intlm)){
  if(max(-log10(pval_int_time_intlm[FourP$snp_info$ChrID=="V",i]))>cutoff_int_int){
    min_p_int_int_chr5 <- c(min_p_int_int_chr5,FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="V"][which.min(pval_int_time_intlm[FourP$snp_info$ChrID=="V",i])])
    loc <- gsub(".*:","",fpkm$locus[fpkm.selc][i])
    middle <- (as.numeric(sub("(^[^-]+)-.*", "\\1", loc)) + as.numeric(sub('.+-(.+)', '\\1', loc))) / 2 / 1000000
    gene_chr <- strsplit(as.character(fpkm$locus[fpkm.selc][i]), "[:]")[[1]][1]
    ifelse(abs(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="V"][which.min(pval_int_time_intlm[FourP$snp_info$ChrID=="V",i])] - middle) < 2 & gene_chr =="V",ct_int_int_chr5 <- c(ct_int_int_chr5,"cis"),ct_int_int_chr5 <- c(ct_int_int_chr5,"trans"))
  }
}

for(i in 1:ncol(pval_int_time_intlm)){
  if(max(-log10(pval_int_time_intlm[FourP$snp_info$ChrID=="X",i]))>cutoff_int_int){
    min_p_int_int_chrx <- c(min_p_int_int_chrx,FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="X"][which.min(pval_int_time_intlm[FourP$snp_info$ChrID=="X",i])])
    loc <- gsub(".*:","",fpkm$locus[fpkm.selc][i])
    middle <- (as.numeric(sub("(^[^-]+)-.*", "\\1", loc)) + as.numeric(sub('.+-(.+)', '\\1', loc))) / 2 / 1000000
    gene_chr <- strsplit(as.character(fpkm$locus[fpkm.selc][i]), "[:]")[[1]][1]
    ifelse(abs(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="X"][which.min(pval_int_time_intlm[FourP$snp_info$ChrID=="X",i])] - middle) < 2 & gene_chr =="X",ct_int_int_chrx <- c(ct_int_int_chrx,"cis"),ct_int_int_chrx <- c(ct_int_int_chrx,"trans"))
  }
}

min_p_int_int_chr1 <- data.frame(marker = min_p_int_int_chr1,ct = ct_int_int_chr1, chr = rep("I",length(min_p_int_int_chr1)))
min_p_int_int_chr2 <- data.frame(marker = min_p_int_int_chr2,ct = ct_int_int_chr2,chr = rep("II",length(min_p_int_int_chr2)))
min_p_int_int_chr3 <- data.frame(marker = min_p_int_int_chr3,ct = ct_int_int_chr3,chr = rep("III",length(min_p_int_int_chr3)))
min_p_int_int_chr4 <- data.frame(marker = min_p_int_int_chr4,ct = ct_int_int_chr4,chr = rep("IV",length(min_p_int_int_chr4)))
min_p_int_int_chr5 <- data.frame(marker = min_p_int_int_chr5,ct = ct_int_int_chr5,chr = rep("V",length(min_p_int_int_chr5)))
min_p_int_int_chrx <- data.frame(marker = min_p_int_int_chrx,ct = ct_int_int_chrx,chr = rep("X",length(min_p_int_int_chrx)))

min_p_max_int_chr_int <- rbind(min_p_int_int_chr1,min_p_int_int_chr2,min_p_int_int_chr3,min_p_int_int_chr4,min_p_int_int_chr5,min_p_int_int_chrx)

# Do the same for FDR = 0.1.

cutoff_int_int01 <- 4.72

min_p_int_int01_chr1 <- c()
min_p_int_int01_chr2 <- c()
min_p_int_int01_chr3 <- c()
min_p_int_int01_chr4 <- c()
min_p_int_int01_chr5 <- c()
min_p_int_int01_chrx <- c()

ct_int_int01_chr1 <- c()
ct_int_int01_chr2 <- c()
ct_int_int01_chr3 <- c()
ct_int_int01_chr4 <- c()
ct_int_int01_chr5 <- c()
ct_int_int01_chrx <- c()

for(i in 1:ncol(pval_int_time_intlm)){
  if(max(-log10(pval_int_time_intlm[FourP$snp_info$ChrID=="I",i]))>cutoff_int_int01){
    min_p_int_int01_chr1 <- c(min_p_int_int01_chr1,FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="I"][which.min(pval_int_time_intlm[FourP$snp_info$ChrID=="I",i])])
    loc <- gsub(".*:","",fpkm$locus[fpkm.selc][i])
    middle <- (as.numeric(sub("(^[^-]+)-.*", "\\1", loc)) + as.numeric(sub('.+-(.+)', '\\1', loc))) / 2 / 1000000
    gene_chr <- strsplit(as.character(fpkm$locus[fpkm.selc][i]), "[:]")[[1]][1]
    ifelse(abs(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="I"][which.min(pval_int_time_intlm[FourP$snp_info$ChrID=="I",i])] - middle) < 2 & gene_chr =="I",ct_int_int01_chr1 <- c(ct_int_int01_chr1,"cis"),ct_int_int01_chr1 <- c(ct_int_int01_chr1,"trans"))
  }
}

for(i in 1:ncol(pval_int_time_intlm)){
  if(max(-log10(pval_int_time_intlm[FourP$snp_info$ChrID=="II",i]))>cutoff_int_int01){
    min_p_int_int01_chr2 <- c(min_p_int_int01_chr2,FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="II"][which.min(pval_int_time_intlm[FourP$snp_info$ChrID=="II",i])])
    loc <- gsub(".*:","",fpkm$locus[fpkm.selc][i])
    middle <- (as.numeric(sub("(^[^-]+)-.*", "\\1", loc)) + as.numeric(sub('.+-(.+)', '\\1', loc))) / 2 / 1000000
    gene_chr <- strsplit(as.character(fpkm$locus[fpkm.selc][i]), "[:]")[[1]][1]
    ifelse(abs(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="II"][which.min(pval_int_time_intlm[FourP$snp_info$ChrID=="II",i])] - middle) < 2 & gene_chr =="II",ct_int_int01_chr2 <- c(ct_int_int01_chr2,"cis"),ct_int_int01_chr2 <- c(ct_int_int01_chr2,"trans"))
  }
}

for(i in 1:ncol(pval_int_time_intlm)){
  if(max(-log10(pval_int_time_intlm[FourP$snp_info$ChrID=="III",i]))>cutoff_int_int01){
    min_p_int_int01_chr3 <- c(min_p_int_int01_chr3,FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="III"][which.min(pval_int_time_intlm[FourP$snp_info$ChrID=="III",i])])
    loc <- gsub(".*:","",fpkm$locus[fpkm.selc][i])
    middle <- (as.numeric(sub("(^[^-]+)-.*", "\\1", loc)) + as.numeric(sub('.+-(.+)', '\\1', loc))) / 2 / 1000000
    gene_chr <- strsplit(as.character(fpkm$locus[fpkm.selc][i]), "[:]")[[1]][1]
    ifelse(abs(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="III"][which.min(pval_int_time_intlm[FourP$snp_info$ChrID=="III",i])] - middle) < 2 & gene_chr =="III",ct_int_int01_chr3 <- c(ct_int_int01_chr3,"cis"),ct_int_int01_chr3 <- c(ct_int_int01_chr3,"trans"))
  }
}

for(i in 1:ncol(pval_int_time_intlm)){
  if(max(-log10(pval_int_time_intlm[FourP$snp_info$ChrID=="IV",i]))>cutoff_int_int01){
    min_p_int_int01_chr4 <- c(min_p_int_int01_chr4,FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="IV"][which.min(pval_int_time_intlm[FourP$snp_info$ChrID=="IV",i])])
    loc <- gsub(".*:","",fpkm$locus[fpkm.selc][i])
    middle <- (as.numeric(sub("(^[^-]+)-.*", "\\1", loc)) + as.numeric(sub('.+-(.+)', '\\1', loc))) / 2 / 1000000
    gene_chr <- strsplit(as.character(fpkm$locus[fpkm.selc][i]), "[:]")[[1]][1]
    ifelse(abs(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="IV"][which.min(pval_int_time_intlm[FourP$snp_info$ChrID=="IV",i])] - middle) < 2 & gene_chr =="IV",ct_int_int01_chr4 <- c(ct_int_int01_chr4,"cis"),ct_int_int01_chr4 <- c(ct_int_int01_chr4,"trans"))
  }
}

for(i in 1:ncol(pval_int_time_intlm)){
  if(max(-log10(pval_int_time_intlm[FourP$snp_info$ChrID=="V",i]))>cutoff_int_int01){
    min_p_int_int01_chr5 <- c(min_p_int_int01_chr5,FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="V"][which.min(pval_int_time_intlm[FourP$snp_info$ChrID=="V",i])])
    loc <- gsub(".*:","",fpkm$locus[fpkm.selc][i])
    middle <- (as.numeric(sub("(^[^-]+)-.*", "\\1", loc)) + as.numeric(sub('.+-(.+)', '\\1', loc))) / 2 / 1000000
    gene_chr <- strsplit(as.character(fpkm$locus[fpkm.selc][i]), "[:]")[[1]][1]
    ifelse(abs(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="V"][which.min(pval_int_time_intlm[FourP$snp_info$ChrID=="V",i])] - middle) < 2 & gene_chr =="V",ct_int_int01_chr5 <- c(ct_int_int01_chr5,"cis"),ct_int_int01_chr5 <- c(ct_int_int01_chr5,"trans"))
  }
}

for(i in 1:ncol(pval_int_time_intlm)){
  if(max(-log10(pval_int_time_intlm[FourP$snp_info$ChrID=="X",i]))>cutoff_int_int01){
    min_p_int_int01_chrx <- c(min_p_int_int01_chrx,FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="X"][which.min(pval_int_time_intlm[FourP$snp_info$ChrID=="X",i])])
    loc <- gsub(".*:","",fpkm$locus[fpkm.selc][i])
    middle <- (as.numeric(sub("(^[^-]+)-.*", "\\1", loc)) + as.numeric(sub('.+-(.+)', '\\1', loc))) / 2 / 1000000
    gene_chr <- strsplit(as.character(fpkm$locus[fpkm.selc][i]), "[:]")[[1]][1]
    ifelse(abs(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="X"][which.min(pval_int_time_intlm[FourP$snp_info$ChrID=="X",i])] - middle) < 2 & gene_chr =="X",ct_int_int01_chrx <- c(ct_int_int01_chrx,"cis"),ct_int_int01_chrx <- c(ct_int_int01_chrx,"trans"))
  }
}

min_p_int_int01_chr1 <- data.frame(marker = min_p_int_int01_chr1,ct = ct_int_int01_chr1, chr = rep("I",length(min_p_int_int01_chr1)))
min_p_int_int01_chr2 <- data.frame(marker = min_p_int_int01_chr2,ct = ct_int_int01_chr2,chr = rep("II",length(min_p_int_int01_chr2)))
min_p_int_int01_chr3 <- data.frame(marker = min_p_int_int01_chr3,ct = ct_int_int01_chr3,chr = rep("III",length(min_p_int_int01_chr3)))
min_p_int_int01_chr4 <- data.frame(marker = min_p_int_int01_chr4,ct = ct_int_int01_chr4,chr = rep("IV",length(min_p_int_int01_chr4)))
min_p_int_int01_chr5 <- data.frame(marker = min_p_int_int01_chr5,ct = ct_int_int01_chr5,chr = rep("V",length(min_p_int_int01_chr5)))
min_p_int_int01_chrx <- data.frame(marker = min_p_int_int01_chrx,ct = ct_int_int01_chrx,chr = rep("X",length(min_p_int_int01_chrx)))

min_p_max_int_chr_int01 <- rbind(min_p_int_int01_chr1,min_p_int_int01_chr2,min_p_int_int01_chr3,min_p_int_int01_chr4,min_p_int_int01_chr5,min_p_int_int01_chrx)

# Add the model and sig factors to these frames as well.
min_p_max_int_chr_int <- min_p_max_int_chr_int %>% mutate(model = as.factor(rep("Interaction term",nrow(min_p_max_int_chr_int))),sig = as.factor(rep(0.5,nrow(min_p_max_int_chr_int))))
min_p_max_int_chr_int01 <- min_p_max_int_chr_int01 %>% mutate(model = as.factor(rep("Interaction term",nrow(min_p_max_int_chr_int01))),sig = as.factor(rep(1,nrow(min_p_max_int_chr_int01))))

# Bind the genotype and interaction variables into frames.
hist_frame_int <- rbind(min_p_max_per_chr_int,min_p_max_int_chr_int)
hist_frame_int01 <- rbind(min_p_max_per_chr_int01,min_p_max_int_chr_int01)

# Recode factor model for better interpretibility.
hist_frame$model <- fct_recode(hist_frame$model,"IMM" = "Interaction", "SMM" = "No time", "AAM" = "Time")

hist_frame01$model <- fct_recode(hist_frame01$model,"IMM" = "Interaction", "SMM" = "No time", "AAM" = "Time")

hist_frame_int$model <- fct_recode(hist_frame_int$model,"IMM" = "Interaction", "IMI" = "Interaction term")

hist_frame_int01$model <- fct_recode(hist_frame_int01$model,"IMM" = "Interaction", "IMI" = "Interaction term")

# Bind interaction term eQTLs to hist_frame, resulting in eqtl_frame and eqtl_01_frame.
eqtl_frame <- rbind(hist_frame,hist_frame_int[hist_frame_int$model == "IMI",])
eqtl_01_frame <- rbind(hist_frame01,hist_frame_int01[hist_frame_int01$model ==  "IMI",])

levels(eqtl_frame$model) <- c("SMM","AAM","IMM", "IMI")
levels(eqtl_01_frame$model) <- c("SMM","AAM","IMM", "IMI")

# Prepare a separate frame for the cis/trans plot.
ct_frame <- rbind(hist_frame,hist_frame_int[hist_frame_int$model == "IMI",])

# Pad with 0s
zero_eqtl <- data.frame(model = c(rep("SMM",6),rep("AAM",6),rep("IMM",6),rep("IMI",6)),
                        sig = rep(0.5,24),
                        marker = rep(0,24),
                        ct = rep("cis",24),
                        chr = rep(c("I",'II',"III","IV","V","X"),4))

eqtl_frame <- rbind(eqtl_frame,zero_eqtl)
ct_frame <- rbind(ct_frame,zero_eqtl)

legend <- get_legend(ggplot(eqtl_frame, aes(x = marker,fill = model)) +
                       geom_histogram(data = eqtl_01_frame, aes(x = marker),fill = "grey", bins = 200) +
                       geom_histogram(bins = 200) +
                       ylab("eQTL count") +
                       xlab("Position on chromosome in Mbp") +
                       facet_grid(factor(model,levels = c("SMM","AAM","IMM", "IMI")) ~chr,scale = "free_x") +
                       coord_cartesian(ylim = c(0,150)) +
                       scale_fill_manual(values = c("blue3","#FF00FF","#1AE4B6FF","red3")) +
                       theme_minimal() +
                       theme(panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank(),strip.text.y = element_blank(),legend.position = "top",text = element_text(size = 35)) +
                       guides(fill=guide_legend(title="Model")))


ssm <- ggplot(eqtl_frame[eqtl_frame$model == "SMM",], aes(x = marker,fill = model)) +
  geom_histogram(data = eqtl_01_frame[eqtl_01_frame$model == "SMM",], aes(x = marker),fill = "grey", bins = 200) +
  geom_histogram(bins = 200) +
  ylab(" ") +
  xlab("Position on chromosome in Mbp") +
  facet_grid(~chr,scale = "free_x",space = "free_x") +
  coord_cartesian(ylim = c(0,600)) +
  scale_fill_manual(values = c("blue3")) +
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank(),strip.text.y = element_blank(),legend.position = "none",text = element_text(size = 35),axis.line = element_line(size = 0.2),axis.ticks= element_line(size = 0.2),axis.title.x = element_blank(),panel.spacing.x = unit(6, "mm")) +
  panel_border(color = "black",size = 0.2)+
  scale_x_continuous(expand =c(0,0),breaks = c(5,10,15,20))#+


ad <- ggplot(eqtl_frame[eqtl_frame$model == "AAM",], aes(x = marker,fill = model)) +
  geom_histogram(data = eqtl_01_frame[eqtl_01_frame$model == "AMM",], aes(x = marker),fill = "grey", bins = 200) +
  geom_histogram(bins = 200) +
  ylab("eQTL count") +
  xlab("Position on chromosome in Mbp") +
  facet_grid(~chr,scale = "free_x",space = "free_x") +
  coord_cartesian(ylim = c(0,600)) +
  scale_fill_manual(values = c("#FF00FF")) +
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank(),strip.text.y = element_blank(),legend.position = "none",text = element_text(size = 35),axis.line = element_line(size = 0.2),axis.ticks= element_line(size = 0.2),axis.title.x = element_blank(),panel.spacing.x = unit(6, "mm")) +
  panel_border(color = "black",size = 0.2)+
  scale_x_continuous(expand =c(0,0),breaks = c(5,10,15,20))

intm <- ggplot(eqtl_frame[eqtl_frame$model %in% c("IMM", "IMI"),], aes(x = marker,fill = model)) +
  geom_histogram(data = eqtl_01_frame[eqtl_01_frame$model %in% c("IMM", "IMI"),], aes(x = marker),fill = "grey", bins = 200) +
  geom_histogram(bins = 200) +
  ylab(" ") +
  xlab("Position on chromosome in Mbp") +
  facet_grid(~chr,scale = "free_x",space = "free_x") +
  coord_cartesian(ylim = c(0,600)) +
  scale_fill_manual(values = c("#1AE4B6FF","red3")) +
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank(),strip.text.y = element_blank(),legend.position = "none",text = element_text(size = 35),axis.line = element_line(size = 0.2),axis.ticks= element_line(size = 0.2),panel.spacing.x = unit(6, "mm"))+
  panel_border(color = "black",size = 0.2)+
  scale_x_continuous(expand =c(0,0),,breaks = c(5,10,15,20))

plt <- plot_grid(legend,ssm,ad,intm,nrow = 4,rel_heights = c(0.1,1,1,1),labels = c("","A","B","C"),align = "h",axis = "tblr",label_size = 23)

png(filename = "plots_new_perm/one_perm_hist_mean_dev_age.png",width = 1500, height = 1200)
plt
dev.off()

pdf("pdfs_figures/Figure_S8.pdf",width = 20, height = 17)
plt
dev.off()


