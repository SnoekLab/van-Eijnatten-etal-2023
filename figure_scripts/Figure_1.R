library(ggplot2)
library(RcppArmadillo)
library(cowplot)
library(tidyverse)
library(viridis)


## Figure 1: Developmental age estimation using PCA. 
    ## A) PCA 
    ## B) Gene expression of the 53 developmental indicator genes from Snoek et al., 2014 
    ##    with expression in this dataset  


# Import marker and fpkm data. 
load(file="~/Documents/elegans/data/obj_FourP_2022.out") # Is this one used in this script? 
load(file="~/Documents/elegans/data/obj_fpkm.out") # should we add this as supplement?

# Filter for mean samples fpkm > 1/32 and more than 20 samples have > 0 fpkm.
fpkm.selc <- apply(log2(fpkm$fpkm),1,mean)> -5 & apply(fpkm$fpkm>0,1,sum)>20 

# use.rge contains all filtered genes with genes as observations and samples as features.
use.rge <- fpkm$fpkm[fpkm.selc,]

#Take log2 ratio of mean.
ruse.rge <- log2((use.rge+1)/apply(use.rge+1,1,mean))

### Panel A) --------------------------------------------------------------------
# Do pca on ruse.rge, and get the loadings of the samples on the first two principal components.
pco <- prcomp((ruse.rge))
pco1 <- -pco$rotation[,1]
pco2 <- pco$rotation[,2]

use.rge <- ruse.rge

# Select only non-parent strains.
use.lines <- c(1:198,207)
pco1.lines <- pco1[use.lines]

# Store names of all genes which are not filtered out.
gene_names <- as.character(fpkm$genes[fpkm.selc])

# Define a dataframe with the expression of vit-2 and loadings on pc1 and pc2 for each of the non-parental samples.
gexpr <- c(use.rge[fpkm$genes[fpkm.selc] == "vit-2",use.lines])
gexpr <- scale(gexpr)
gtype <- c(rep("mpRIL",199))
gene <- rep("vit-2",199)

cor(gexpr,pco1.lines)

topl <- data.frame(pco1[use.lines],pco2[use.lines],gtype,gexpr,gene)

colnames(topl)[1:2] <- c("pco1","pco2")

# Plot expression vit-2 over PC1.
vit_2_plot <- ggplot(topl,aes(pco1,pco2,col=gexpr))+
  geom_point(size = 5)+
  scale_color_gradientn(colors=plasma(256)) +
  xlab("PC1 (48%)") +
  ylab("PC2 (5%)")+
  theme_cowplot()+
  guides(col=guide_colorbar(title="Z-score vit-2 levels",frame.colour = "black",frame.linewidth = 0.2))+
  theme(text = element_text(size = 40),
        legend.text = element_text(size = 40),
        axis.text.y = element_text(angle = 90, size = 40,vjust = 0.1),
        axis.text.x = element_text(size = 40,hjust = 0.1),
        legend.position = "top",
        panel.border = element_rect(size=0.2,color="black"),
        axis.line = element_blank(),
        legend.justification = "center",
        legend.key.width = unit(2, 'cm'))

vit_2_plot

# Select developmental marker genes from cluster 1, Snoek et al 2014, which are significantly expressed in our dataset.
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

# Sort levels of this factor as well.
dev_frame$mpRIL <- factor(dev_frame$mpRIL, levels = names(sort(pco1.lines)))

# Convert frame to long form.
dev_frame_long <- dev_frame %>% pivot_longer(cols = -c(54), names_to = "Gene",values_to = "Gexpr")

pc1load <- topl[as.character(dev_frame_long$mpRIL),"pco1"]
for.line <- data.frame(dev_frame_long,pc1load) 

### Panel B) --------------------------------------------------------------------
### Make line figure

line.fig <- ggplot(for.line,aes(x = pc1load,y = Gexpr, col= Gene)) +
  geom_point(col="gray",size=1,alpha=0.2)+
  geom_smooth(se=F,size=0.2) +
  xlab("mpRILs loading on PC1") +
  ylab("Z-score gene expression")+
  theme_cowplot()+
  theme(axis.line = element_blank(),
        panel.border = element_rect(size=0.2,color="black",fill=NA),
        panel.grid.major = element_line(size = 0.2,color = "grey"),
        text = element_text(size = 40),
        axis.text.x = element_text(size = 40),
        axis.text.y = element_text(size = 40),
        legend.position = "none")
line.fig

# Paste subfigures together.
pc_plot <- plot_grid(vit_2_plot,line.fig,align = "v",ncol = 2,labels = c("A","B"),label_size = 35)

# Save
png("plots_new_perm/vit2_and_lineplot.png",width = 1700,height=750)
pc_plot
dev.off()

pdf("pdfs_figures/Figure_1.pdf",width = 30, height = 14)
pc_plot
dev.off()









