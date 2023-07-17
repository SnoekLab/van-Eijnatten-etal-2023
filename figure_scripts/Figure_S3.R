
library(ggplot2)
library(cowplot)
library(tidyverse)
library(viridis)

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

# Find developmental marker genes with significant expression in our dataset.
dev_genes <- read.table(file="~/Documents/elegans/data/Dev_clust1_genes.txt",header = T)
dev_genes <- as.vector(as.character(dev_genes$genes))
dev_genes <- dev_genes[dev_genes %in% gene_names]

# Define matrix with gene expression for these marker genes.
dev_frame <- t(use.rge[gene_names %in% dev_genes,use.lines])

# Mean center and scale in the sample direction.
dev_frame <- apply(dev_frame,2,scale)

# Define vector with the mean gene expression of the developmental marker genes for each sample.
means_gen <- rowMeans(dev_frame)

# Turn into dataframe.
dev_frame <- as.data.frame(dev_frame)
colnames(dev_frame) <- dev_genes[dev_genes%in%gene_names]

# use raptor ages for ordering mpRILs. Using Cel_larval as reference with prior of 48 with sd 1.

raptor_ages <- readRDS("RAPToR_on_mpRILs/raptor_ages_noprior_cel_larval.rds")

dev_frame <- dev_frame[order(raptor_ages),]

dev_frame <- dev_frame %>% mutate(mpRIL = as.factor(names(pco1.lines)[order(raptor_ages)]))

dev_frame$mpRIL <- factor(dev_frame$mpRIL, levels = names(pco1.lines)[order(raptor_ages)])

dev_frame_long <- dev_frame %>% pivot_longer(cols = -c(54), names_to = "Gene",values_to = "Expression")

raptor_hm <- dev_frame_long %>% ggplot(aes(x = mpRIL,y = Gene,fill = Expression)) +
  geom_tile() + 
  scale_x_discrete(breaks = dev_frame$mpRIL[c(T,F,F,F,F)]) + 
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),text = element_text(size = 30)) +
  theme(axis.text.y = element_text(size=10, hjust=1))+
  scale_fill_viridis(option = "A",name = "Z-score")+
  xlab("mpRILs sorted by raptor age estimates") +
  ylab("Developmental marker gene")

png("plots_w_cutoffs/hm_raptor.png",height = 650, width = 1000)
raptor_hm
dev.off()

pdf("pdfs_figures/Figure_S3.pdf",width = 15, height = 9)
raptor_hm
dev.off()

