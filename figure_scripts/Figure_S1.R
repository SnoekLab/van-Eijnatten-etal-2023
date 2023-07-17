library(ggplot2)
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

rownames(pco$x) <- gene_names

### Read in developmental marker genes, cluster 1 Snoek et al 2014
dev_genes <- read.table(file="~/Documents/elegans/data/Dev_clust1_genes.txt",header = T)
dev_genes <- as.vector(as.character(dev_genes$genes))

# Find developmental marker genes with significant expression in our dataset.
dev_genes <- dev_genes[dev_genes %in% gene_names]

### Make frames with scores of genes on pcs
dev_scores <- data.frame(pco$x[rownames(pco$x)%in%dev_genes,])

### Make factor to distinguish developmental marker genes from other genes.
dev_fac <- ifelse(gene_names %in% dev_genes,"Dev gene","Other gene")

### Plot
dev_scores_plot <- ggplot(data.frame(pco$x),aes(x = -PC1, y = PC2, col = dev_fac))+
  geom_point(size = 3,alpha=0.8)+
  theme_minimal()+
  xlab("PC1")+
  ylab("PC2")+
  theme(text = element_text(size = 25),panel.border = element_rect(fill = NA,linewidth = 0.2))+
  scale_color_discrete(type = c("black","lightgreen"),name = "Gene type")

png(filename = "plots_new_perm/dev_scores_plot.png",width = 800,height = 500)
dev_scores_plot
dev.off()

pdf("pdfs_figures/Figure_S1.pdf",width = 9,height = 4.5)
dev_scores_plot
dev.off()






