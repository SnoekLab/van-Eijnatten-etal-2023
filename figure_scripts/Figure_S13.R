
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

# Read in p-values IMI model
pval_int_time_intlm <- readRDS("~/Documents/elegans/pvals_correct/pval_int_time_intlm.RDS")

# Initialise vector with names of genes we want to plot.
genes <- c("spsb-2","rpl-34","wrm-1")

# Get position and chromosome of most significant marker for these genes.
markers <- sapply(genes,function(x) which.min(pval_int_time_intlm[FourP$snp_info$ChrID == "V",gene_names == x]))
marker_pos <- as.character(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="V"][markers]*10^6)
marker_chr <- c("V","V","V")

# Define function to plot each gene.

plot_gene <- function(dev,row,col,chr){
  plot_frame <- data.frame(dev,use.rge[row,use.lines],t(unname(round(FourP$snp_mat[FourP$snp_info$ChrID=="V",][as.numeric(col),use.lines]*2,0)/2)), chr = chr)
  names(plot_frame) <- c("dev", "gexpr","gt")
  name <- gsub('\\.', '-', gene_names[row])
  plt <- plot_grid(ggplot(plot_frame,aes(dev,gexpr,col=as.factor(gt)))+
                     geom_point(size =3)+
                     geom_smooth(data = plot_frame[plot_frame$gt%in%c(0,1),],method="lm",size=1.5)+
                     scale_color_manual(values = c("0"="#377eb8","0.5"="#4daf4a","1"="#e41a1c"),name = "Genotype")+
                     #scale_color_viridis(option="H",discrete = T)+
                     scale_fill_manual(values = c("0" = "#377eb8", "1" = "#e41a1c"),guide = "none")+
                     theme_cowplot() +
                     ggtitle(paste0(paste("Gene", name),"\n",paste("Marker",FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="V"][col],"Mb","chr",chr))) +
                     xlab("Estimated age (PC1)") + 
                     ylab("Center log ratio gene expression")+
                     scale_size(name = "Gene expression")+
                     theme(text = element_text(size = 45),axis.text = element_text(size =45),legend.position = "none"),ggplot(plot_frame[plot_frame$gt != 0.5,],aes(y = gexpr,fill = as.factor(gt)))+
                     geom_boxplot()+
                     theme_void()+
                     scale_fill_manual(values = c("#377eb8","#e41a1c")) +
                     theme(legend.position = "none") +
                     ylab("Normalized expression") +
                     theme(legend.title = element_blank()),rel_widths = c(6,1),align = "h",axis="rtbl")+
    panel_border(size = 0.2,color = "black")
  #ggsave(paste0("spline_plots/sig_marker/",name,"chr",chr,".png"),plot=plt,device = "png")
  return(plt)
}

plot_frame <- data.frame(dev = pco1.lines,use.rge[which(gene_names==genes[1]),use.lines],t(unname(round(FourP$snp_mat[as.numeric(markers[1]),use.lines]*2,0)/2)), chr = marker_chr[1])
names(plot_frame) <- c("dev", "gexpr","gt")
legend <- get_legend(ggplot(plot_frame,aes(dev,gexpr,col=as.factor(gt)))+
                       geom_point(size =3)+
                       geom_smooth(data = plot_frame[plot_frame$gt%in%c(0,1),],method="lm",size=1.5)+
                       scale_color_manual(values = c("0"="#377eb8","0.5"="#4daf4a","1"="#e41a1c"),name = "Genotype")+
                       scale_fill_manual(values = c("0" = "blue", "1" = "red"),guide = "none")+
                       theme_cowplot() +
                       xlab("Developmental age") + 
                       ylab("Log2 ratio gene expression")+
                       scale_size(name = "Gene expression")+
                       theme(text = element_text(size = 50),legend.position = "top",legend.direction = "horizontal",legend.justification="center" ,legend.box.just = "bottom"))

# Call functions and paste plots together using plot_grid()
gene_1 <- plot_gene(pco1.lines,which(gene_names==genes[1]),markers[1],marker_chr[1])
gene_1
gene_2 <- plot_gene(pco1.lines,which(gene_names==genes[2]),markers[2],marker_chr[2])
gene_2
gene_3 <- plot_gene(pco1.lines,which(gene_names==genes[3]),markers[3],marker_chr[3])
gene_3

### Make  dataframe with genotypes and PC1 loadings at interactions hotspot locus
to.pl.20.440891 <- data.frame(t(unname(round(FourP$snp_mat[FourP$snp_info$SNP_position==20.440891,use.lines]*2,0)/2)),pco1.lines)
colnames(to.pl.20.440891) <- c("Genotype","Developmental_age")


### Function to get p-values for PC1~genotype linear model
map_pc <- function(gen,pc1){
  ok <- summary(lm(pc1~gen))
  pval <- as.numeric(ok$coefficients[2,4])
  return(-log10(pval))
}

map_pc(as.numeric(FourP$snp_mat[which(FourP$snp_info$SNP_position==20.440891),use.lines]),pco1.lines)

fig_20.440891 <- ggplot(to.pl.20.440891[to.pl.20.440891$Genotype != 0.5,],aes(y = Developmental_age,fill = as.factor(Genotype)))+
  geom_boxplot()+
  scale_fill_manual(values = c("0" = "#377eb8", "1" = "#e41a1c"),guide = "none") +
  theme_minimal() +
  ylab("Developmental age") +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 50)) +
  ggtitle("20.440891 Mb, chr V\n -log10(Pval) = 3.2")+
  panel_border(size = 0.2,color = "black")

png(filename = "plots_new_perm/int_hotspot_examples.png",height = 2000,width = 2100)
plot_grid(legend,plot_grid(gene_1,gene_2,gene_3,fig_20.440891,nrow = 2,align = "h",labels = c("A","B","C","D"),label_size =40,axis="rtbl"),nrow = 2,rel_heights = c(1,10))+
  theme(plot.margin = unit(c(20,20,20,20), "points"))
dev.off()

pdf("pdfs_figures/Figure_S13.pdf",height = 30,width = 31.5)
plot_grid(legend,plot_grid(gene_1,gene_2,gene_3,fig_20.440891,nrow = 2,align = "h",labels = c("A","B","C","D"),label_size =40,axis="rtbl"),nrow = 2,rel_heights = c(1,10))+
  theme(plot.margin = unit(c(20,20,20,20), "points"))
dev.off()

