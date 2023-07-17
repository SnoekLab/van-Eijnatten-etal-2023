### This script produces a histogram of the distribution of eQTLs detected with linear models over the genome.

library(ggplot2)
library(cowplot)
library(viridis)
library(profvis)
library(ggvenn)

# Import marker and fpkm data. 
load(file="~/Documents/elegans/data/obj_FourP_2022.out") 
load(file="~/Documents/elegans/data/obj_fpkm.out")

# Filter for genes with mean fpkm > 1/32 and more than 20 samples have > 0 fpkm.
fpkm.selc <- apply(log2(fpkm$fpkm),1,mean)> -5 & apply(fpkm$fpkm>0,1,sum)>20 

# use.rge contains all filtered samples with genes as observations and samples as features.
use.rge <- fpkm$fpkm[fpkm.selc,]

# Take log2 ratio of mean.
ruse.rge <- log2((use.rge+1)/apply(use.rge+1,1,mean))

# Do pca on ruse.rge, and get the loadings of the samples on the first two principles components.
pco <- prcomp((ruse.rge))
pco1 <- -pco$rotation[,1]
pco2 <- pco$rotation[,2]

use.rge <- ruse.rge

# Select only non-parent strains.
use.lines <- c(1:198,207)
pco1.lines <- pco1[use.lines]

# Store names of all genes which are not filtered out.
gene_names <- as.character(fpkm$genes[fpkm.selc])

### Get eQTLs 


get_eqtl <- function(gene,chr,pvals,cutoff,model,sig){
  index <- min(which(gene_name == ""))
  gene_name[index] <<- gene_names[gene]
  pos[index] <<- FourP$snp_info$SNP_position[FourP$snp_info$ChrID==chr][which.max(pvals)]
  loc <- gsub(".*:","",fpkm$locus[fpkm.selc][gene])
  middle <- (as.numeric(sub("(^[^-]+)-.*", "\\1", loc)) + as.numeric(sub('.+-(.+)', '\\1', loc))) / 2 / 1000000
  gene_pos[index] <<- middle
  gene_chr <- strsplit(as.character(fpkm$locus[fpkm.selc][gene]), "[:]")[[1]][1]
  gene.chr[index] <<- gene_chr
  ct[index] <<- ifelse(abs(FourP$snp_info$SNP_position[FourP$snp_info$ChrID==chr][which.max(pvals)] - middle) < 2 & gene_chr ==chr,ct <- c("cis"),ct <- c("trans"))
  marker_chr[index] <<- chr
  gene_index[index] <<- gene
}

map_function <- function(gene,cutoff,pvals,model,sig){
  pvals_chr1 <- pvals[FourP$snp_info$ChrID=="I"]
  pvals_chr2 <- pvals[FourP$snp_info$ChrID=="II"]
  pvals_chr3 <- pvals[FourP$snp_info$ChrID=="III"]
  pvals_chr4 <- pvals[FourP$snp_info$ChrID=="IV"]
  pvals_chr5 <- pvals[FourP$snp_info$ChrID=="V"]
  pvals_chrx <- pvals[FourP$snp_info$ChrID=="X"]
  if(max(pvals_chr1)>cutoff){
    get_eqtl(gene,"I",pvals_chr1,cutoff,model,sig)
  }
  if(max(pvals_chr2)>cutoff){
    get_eqtl(gene,"II",pvals_chr2,cutoff,model,sig)
  }
  if(max(pvals_chr3)>cutoff){
    get_eqtl(gene,"III",pvals_chr3,cutoff,model,sig)
  }
  if(max(pvals_chr4)>cutoff){
    get_eqtl(gene,"IV",pvals_chr4,cutoff,model,sig)
  }
  if(max(pvals_chr5)>cutoff){
    get_eqtl(gene,"V",pvals_chr5,cutoff,model,sig)
  }
  if(max(pvals_chrx)>cutoff){
    get_eqtl(gene,"X",pvals_chrx,cutoff,model,sig)
  }
}

cutoff_notime <- 4.06

### eQTLs SMM
pval_gen_notime_lm <- readRDS("~/Documents/elegans/pvals_correct/pval_gen_notime_lm.RDS")

gene_name <- character(length = 18605)
pos <- numeric(length = 18605)
gene_pos <- numeric(length = 18605)
gene.chr <- character(length = 18605)
ct <- character(length = 18605)
marker_chr <- character(length = 18605)
gene_index <- numeric(length = 18605)

sapply(1:ncol(pval_gen_notime_lm),function(x) map_function(x,cutoff_notime,-log10(pval_gen_notime_lm[,x]),"Single marker model",0.05))

eqtl.smm <- data.frame(gene_name = gene_name, marker_pos = pos, gene_pos = gene_pos,ct = ct, marker_chr = marker_chr, gene_chr = gene.chr,gene_index = gene_index)

### eQTLs AAM

cutoff_dev <- 4.44

pval_gen_time_lm <- readRDS("~/Documents/elegans/pvals_correct/pval_gen_time_lm.RDS")

gene_name <- character(length = 8666)
pos <- numeric(length = 8666)
gene_pos <- numeric(length = 8666)
gene.chr <- character(length = 8666)
ct <- character(length = 8666)
marker_chr <- character(length = 8666)
gene_index <- numeric(length = 8666)

sapply(1:ncol(pval_gen_time_lm),function(x) map_function(x,cutoff_dev,-log10(pval_gen_time_lm[,x]),"Additive model",0.05))

eqtl.am <- data.frame(gene_name = gene_name, marker_pos = pos, gene_pos = gene_pos,ct = ct, marker_chr = marker_chr, gene_chr = gene.chr,gene_index = gene_index)

### eQTLs IMM

cutoff_gen_int <- 4.60

pval_gen_time_intlm <- readRDS("~/Documents/elegans/pvals_correct/pval_gen_time_intlm.RDS")

gene_name <- character(length = 6886)
pos <- numeric(length = 6886)
gene_pos <- numeric(length = 6886)
gene.chr <- character(length = 6886)
ct <- character(length = 6886)
marker_chr <- character(length = 6886)
gene_index <- numeric(length = 6886)

sapply(1:ncol(pval_gen_time_intlm),function(x) map_function(x,cutoff_gen_int,-log10(pval_gen_time_intlm[,x]),"Interaction model, marker",0.05))

eqtl.im.gen <- data.frame(gene_name = gene_name, marker_pos = pos, gene_pos = gene_pos,ct = ct, marker_chr = marker_chr, gene_chr = gene.chr,gene_index = gene_index)

### eQTLs IMI

cutoff_int_int <- 5.60

pval_int_time_intlm <- readRDS("~/Documents/elegans/pvals_correct/pval_int_time_intlm.RDS")

gene_name <- character(length = 1161)
pos <- numeric(length = 1161)
gene_pos <- numeric(length = 1161)
gene.chr <- character(length = 1161)
ct <- character(length = 1161)
marker_chr <- character(length = 1161)
gene_index <- numeric(length = 1161)

sapply(1:ncol(pval_int_time_intlm),function(x) map_function(x,cutoff_int_int,-log10(pval_int_time_intlm[,x]),"Interaction model, interaction",0.05))

eqtl.im.int <- data.frame(gene_name = gene_name, marker_pos = pos, gene_pos = gene_pos,ct = ct, marker_chr = marker_chr, gene_chr = gene.chr,gene_index = gene_index)

### Put all eQTLs in list (gene/marker combo's)

eqtl.list <- list()

eqtl.list[[1]] <- paste0(eqtl.smm$gene_index,eqtl.smm$marker_chr)
eqtl.list[[2]] <- paste0(eqtl.am$gene_index,eqtl.am$marker_chr)
eqtl.list[[3]] <- paste0(eqtl.im.gen$gene_index,eqtl.im.gen$marker_chr)
eqtl.list[[4]] <- paste0(eqtl.im.int$gene_index,eqtl.im.int$marker_chr)

names(eqtl.list) <- c("SMM", "AM", "IMM","IMI")

### Use list to make venn diagram
venn.plot <- ggvenn(eqtl.list,c("SMM", "AM", "IMM","IMI")
                    ,text_size = 8,
                    set_name_size = 10)

png(filename = "plots_new_perm/venn_general.png",width = 1000,height = 800)
venn.plot
dev.off()

pdf("pdfs_figures/Figure_S6.pdf",width = 13,height = 11)
venn.plot
dev.off()


