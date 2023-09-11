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

km <- readRDS("km_hendriks_fluc.RDS")

### Make heatmaps of km clusters

for(i in 1:10){
  clus_frame <- t(use.rge[gene_names %in% gene_names[rownames(use.rge)%in%names(km$cluster)[km$cluster==i]],use.lines])
  
  # Mean center and scale data in the sample direction
  clus_frame <- apply(clus_frame,2,function(x) scale(x))
  clus_frame <- as.data.frame(clus_frame)
  colnames(clus_frame) <- gene_names[gene_names %in% gene_names[rownames(use.rge)%in%names(km$cluster)[km$cluster==i]]]
  
  # Correlation after scaling in column direction.
  mean(apply(clus_frame,2,function (x) cor(x,pco1.lines)))
  
  # Order samples by their loading on the first principle component.
  clus_frame <- clus_frame[order(pco1.lines),]
  
  # Include sample names as factor in frame.
  clus_frame <- clus_frame[, !duplicated(colnames(clus_frame))] %>% mutate(mpRIL = as.factor(names(sort(pco1.lines))))
  
  # Sort this factor as well.
  clus_frame$mpRIL <- factor(clus_frame$mpRIL, levels = names(sort(pco1.lines)))
  
  # Convert frame to long form.
  clus_frame_long <- clus_frame %>% pivot_longer(cols = -c(ncol(clus_frame)), names_to = "Gene",values_to = "Gexpr")
  
  
  # Produce heatmap.
  assign(paste0("hm_cluster_",i),clus_frame_long %>% ggplot(aes(x = mpRIL,y = Gene,fill = Gexpr)) +
           geom_tile() + 
           scale_x_discrete(breaks = clus_frame_long$mpRIL[c(T,F,F,F,F)]) + 
           theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),text = element_text(size = 45)) +
           theme(axis.text.y = element_blank(),axis.ticks.y = element_blank())+
           scale_fill_viridis(option = "A")+
           xlab("mpRILs sorted by loading on pc1") +
           ylab("Genes")+
           theme(legend.position = "none",plot.title = element_text(hjust = 0.5)) +
           ggtitle(paste("cluster",i)))
}


heatmap.fig <- plot_grid(hm_cluster_1, hm_cluster_2, hm_cluster_3, hm_cluster_4,hm_cluster_5, hm_cluster_6,hm_cluster_7, hm_cluster_8,hm_cluster_9, hm_cluster_10,nrow = 5)

#### Make lineplots

for(i in 1:10){
  line_frame <- t(use.rge[gene_names %in% gene_names[rownames(use.rge)%in%names(km$cluster)[km$cluster==i]],use.lines])
  
  # Mean center and scale data in the sample direction
  # line_frame <- apply(line_frame,2,scale)
  line_frame <- as.data.frame(line_frame)
  colnames(line_frame) <- gene_names[gene_names %in% gene_names[rownames(use.rge)%in%names(km$cluster)[km$cluster==i]]]
  
  # Correlation after scaling in column direction.
  mean(apply(line_frame,2,function (x) cor(x,pco1.lines)))
  
  # Order samples by their loading on the first principle component.
  line_frame <- line_frame[order(pco1.lines),]
  
  # Include sample names as factor in frame.
  line_frame <- line_frame[, !duplicated(colnames(line_frame))] %>% mutate(pc1load = pco1.lines[order(pco1.lines)])
  
  # Convert frame to long form.
  line_frame_long <- line_frame %>% pivot_longer(cols = -c("pc1load"), names_to = "Gene",values_to = "Gexpr")
  
  assign(paste0("line_cluster_",i),ggplot(line_frame_long,aes(x = pc1load,y = Gexpr)) +
           geom_point(col="gray",size=1,alpha=0.2)+
           geom_smooth(se=F,size=0.2) +
           theme_minimal()+
           ylab("Normalized expression")+ 
           xlab("Developmental age (PC1)")+ 
           theme(text = element_text(size = 45),plot.title = element_text(hjust = 0.5))+
           ggtitle(paste("cluster",i)))
  
}


lineplot.fig <- plot_grid(line_cluster_1, line_cluster_2, line_cluster_3, line_cluster_4,line_cluster_5, line_cluster_6,line_cluster_7, line_cluster_8,line_cluster_9, line_cluster_10,nrow = 5)


pdf("pdfs_figures/Figure_S17.pdf",width = 50, height = 42)
plot_grid(heatmap.fig,lineplot.fig,ncol = 2)
dev.off()


