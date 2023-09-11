library(readxl)
library(tidyr)
library(dplyr)
library(ggplot2)
library(viridis)
library(cowplot)

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

# Read in data of Hendriks et al (publicly available)
gene_classes <- read_excel("~/Documents/elegans/data/oscillating_genes_Hendriks_2014.xlsx")

# Read in wormbase identifiers and gene names of our gene expression matrix. 
wormbase_selc <- read.table("~/Documents/elegans/data/accession_fpkmselc.txt", header = TRUE, sep = "\t",stringsAsFactors = FALSE)
colnames(wormbase_selc) <- c("Gene_names", "Wormbase_accession")
wormbase_selc <- wormbase_selc[wormbase_selc$Wormbase_accession!="Multiple entries : 2",]
colnames(gene_classes) <- gene_classes[3,]
gene_classes <- gene_classes[4:nrow(gene_classes),]

### Select all genes annotated as rising by Hendricks et al
rising_frame <- gene_classes[gene_classes$class == "rising",]
rising_frame <- rising_frame[,-c(5,6)]

# Vector index_rising gathers all genes from our filtered gene expression matrix that are identified as 
# rising by Hendriks et al.
index_rising <- c()

for(i in 1:length(gene_names)){
  split <- strsplit(gene_names[i],",")
  in_ref <- wormbase_selc[wormbase_selc$Gene_names %in% split[[1]],2]%in%rising_frame$`Gene WB ID`
  if(sum(in_ref) > 0){
    index_rising <- c(index_rising,i) 
  }
}

### Make figure S2

sorted_frame_rise <- t(use.rge[gene_names %in% gene_names[index_rising],use.lines])

# Mean center and scale in the sample direction.
sorted_frame_rise <- apply(sorted_frame_rise,2,scale)

# Mean center and scale data in the sample direction
# sorted_frame_rise <- apply(sorted_frame_rise,2,function(x) scale(x))
sorted_frame_rise <- as.data.frame(sorted_frame_rise)
colnames(sorted_frame_rise) <- gene_names[gene_names %in% gene_names[index_rising]]

# Correlation after scaling in column direction.
mean(apply(sorted_frame_rise,2,function (x) cor(x,pco1.lines)))

# Order samples by their loading on the first principle component.

slope <- function(gexpr){
  return(lm(gexpr~pco1.lines)$coefficients[2])
}

slopes <- apply(sorted_frame_rise,2,slope)

sorted_frame_rise <- sorted_frame_rise[order(pco1.lines),]

sorted_frame_rise <- sorted_frame_rise[,order(slopes,decreasing = T)]

# Include sample names as factor in frame.
sorted_frame_rise <- sorted_frame_rise[, !duplicated(colnames(sorted_frame_rise))] %>% mutate(mpRIL = as.factor(names(sort(pco1.lines))))

# Sort this factor as well.
sorted_frame_rise$mpRIL <- factor(sorted_frame_rise$mpRIL, levels = names(sort(pco1.lines)))

# Convert frame to long form.
sorted_frame_rise_long <- sorted_frame_rise %>% pivot_longer(cols = -c(ncol(sorted_frame_rise)), names_to = "Gene",values_to = "Gexpr")

# Sort this factor as well.
sorted_frame_rise_long$Gene <- factor(sorted_frame_rise_long$Gene, levels = names(sort(slopes[!duplicated(names(slopes))],decreasing = T)))

# Produce heatmap.
heatmap_sorted_rise <- sorted_frame_rise_long %>% ggplot(aes(x = mpRIL,y = Gene,fill = Gexpr)) +
  geom_tile() + 
  scale_x_discrete(breaks = sorted_frame_rise_long$mpRIL[c(T,F,F,F,F)]) + 
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),text = element_text(size = 25)) +
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank())+

  scale_fill_viridis(option = "A",name = "Z-score")+
  xlab("mpRILs sorted by loading on PC1") +
  ylab("Rising genes Hendriks et al")

png(filename = "plots_new_perm/rising_hendriks_sorted.png",width  = 1000, height = 625)
heatmap_sorted_rise
dev.off()

pdf("pdfs_figures/Figure_S2.pdf",width  = 10, height = 6.25)
heatmap_sorted_rise
dev.off()





