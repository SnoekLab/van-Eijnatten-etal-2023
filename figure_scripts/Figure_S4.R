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

# Select developmental marker genes from cluster 1, Snoek et al 2014 which are significantly expressed in our dataset.
dev_genes <- read.table(file="~/Documents/elegans/data/Dev_clust1_genes.txt",header = T)
dev_genes <- as.vector(as.character(dev_genes$genes))
dev_genes <- dev_genes[dev_genes %in% gene_names]

dev_frame <- t(use.rge[gene_names %in% dev_genes,use.lines])

# Calculate mean correlation of expression dev_genes with pc1.
mean(apply(dev_frame,2,function (x) cor(x,pco1.lines)))

# Mean center and scale data in the sample direction
dev_frame <- apply(dev_frame,2,scale)

dev_frame <- as.data.frame(dev_frame)
colnames(dev_frame) <- dev_genes[dev_genes%in%gene_names]

# Correlation after scaling in column direction.
mean(apply(dev_frame,2,function (x) cor(x,pco1.lines)))

# Order samples by their loading on the first principle component.
dev_frame <- dev_frame[order(pco1.lines),]

# Include sample names as factor in frame.
dev_frame <- dev_frame %>% mutate(mpRIL = as.factor(names(sort(pco1.lines))))

# Sort this factor as well.
dev_frame$mpRIL <- factor(dev_frame$mpRIL, levels = names(sort(pco1.lines)))

# Convert frame to long form.
dev_frame_long <- dev_frame %>% pivot_longer(cols = -c(54), names_to = "Gene",values_to = "Gexpr")

# Make figure
heatmap_supplement <- dev_frame_long %>% ggplot(aes(x = mpRIL,y = Gene,fill = Gexpr)) +
  geom_tile() + 
  scale_x_discrete(breaks = dev_frame$mpRIL[c(T,F,F,F,F)]) + 
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),text = element_text(size = 30)) +
  theme(axis.text.y = element_text(size=10, hjust=1))+
  scale_fill_viridis(option = "A")+
  xlab("mpRILs sorted by loading on PC1") +
  ylab("Developmental marker gene")

png("plots_new_perm/heatmap_supplement.png",width =  1000, height = 650)
heatmap_supplement
dev.off()

pdf("pdfs_figures/Figure_S4.pdf",width = 15, height =9)
heatmap_supplement
dev.off()



