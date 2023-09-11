
library(ggplot2)
library(cowplot)
library(dplyr)
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

# Define with the samples in the same order as in the colomns of the datset (and the pca).
gtype <- c(rep("mpRIL",198),rep(c("JU1511","JU1926","JU1931","JU1941"),each=2),"mpRIL")

# Define frame with loadings and samples for plotting.
gtype_frame <- data.frame(pc1 = pco1, pc2 = pco2, genotype = gtype)

# Make it so that genotype variable is a factor.
gtype_frame$genotype <- factor(gtype_frame$genotype, levels = c("mpRIL", "JU1511", "JU1926", "JU1931", "JU1941"))

# Produce plot.
gtype.pca <- ggplot(data = gtype_frame,aes(x = pc1, y = pc2,col = genotype)) +
  geom_point(size =3) +
  scale_color_viridis(discrete = T,direction = -1) +
  theme_minimal() +
  theme(text = element_text(size = 20))+
  xlab("PC1 48%")+
  ylab("PC2 5%")+
  guides(col=guide_legend(title="Genotype"))

pdf("pdfs_figures/Figure_S5.pdf",width = 8, height= 6)
gtype.pca
dev.off()



