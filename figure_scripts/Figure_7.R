
library(ggplot2)
library(cowplot)
library(dplyr)
library(viridis)
library(ggpattern)

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

# Read in p-values SMM and AAM.
pval_gen_notime_lm <- readRDS("~/Documents/elegans/pvals_correct/pval_gen_notime_lm.RDS")
pval_gen_time_lm <- readRDS("~/Documents/elegans/pvals_correct/pval_gen_time_lm.RDS")

# Initialise vector with names of genes we want to plot.
genes <- c("sams-5","nspd-3")

which.min(pval_gen_notime_lm[FourP$snp_info$ChrID=="V",gene_names=="sams-5"])
which.min(pval_gen_time_lm[FourP$snp_info$ChrID=="V",gene_names=="sams-5"])

which.min(pval_gen_notime_lm[FourP$snp_info$ChrID=="X",gene_names=="nspd-3"])
which.min(pval_gen_time_lm[FourP$snp_info$ChrID=="X",gene_names=="nspd-3"])

### Store most significant markers on chromosome X for these genes according to AAM and SMM.
markers_notime <- c(which.min(pval_gen_notime_lm[FourP$snp_info$ChrID=="V",gene_names=="sams-5"]),which.min(pval_gen_notime_lm[FourP$snp_info$ChrID=="X",gene_names=="nspd-3"]))
markers_time <- c(which.min(pval_gen_time_lm[FourP$snp_info$ChrID=="V",gene_names=="sams-5"]),which.min(pval_gen_time_lm[FourP$snp_info$ChrID=="X",gene_names=="nspd-3"]))

### Get position of these markers.
marker_notime_pos <- c(as.character(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="V"][markers_notime[1]]*10^6),
                       as.character(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="X"][markers_notime[2]]*10^6))
marker_time_pos <- c(as.character(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="V"][markers_time[1]]*10^6),
                     as.character(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="X"][markers_time[2]]*10^6))

marker_chr <- c("V","X")


# Specify function for plotting the best marker found with the single marker model.
plot_gene <- function(dev,row,col,chr){
  plot_frame <- data.frame(dev,use.rge[row,use.lines],t(unname(round(FourP$snp_mat[FourP$snp_info$ChrID==chr,][as.numeric(col),use.lines]*2,0)/2)), chr = chr)
  names(plot_frame) <- c("dev", "gexpr","gt")
  name <- gsub('\\.', '-', gene_names[row])
  plt <- plot_grid(ggplot(plot_frame,aes(dev,gexpr,col=as.factor(gt)))+
                     geom_point(size =3,alpha = 0.4)+
                     geom_smooth(data = plot_frame[plot_frame$gt%in%c(0,1),],method="lm",size=1.5, linetype = "dashed")+
                     scale_color_manual(values = c("0"="#377eb8","0.5"="#4daf4a","1"="#e41a1c"),name = "Genotype")+
                     #scale_color_viridis(option="H",discrete = T)+
                     scale_fill_manual(values = c("0" = "#377eb8", "1" = "#e41a1c"),guide = "none")+
                     theme_cowplot() +
                     ggtitle(paste0(paste("Gene",name),"\n",paste("Marker position",FourP$snp_info$SNP_position[FourP$snp_info$ChrID==chr][col],"Mb","chr",chr))) +
                     xlab("Estimated age (PC1)") + 
                     ylab("Gene expression")+
                     scale_size(name = "Gene expression")+
                     theme(text = element_text(size = 29),legend.position = "none",axis.text.x = element_text(size = 29),axis.text.y = element_text(size = 29)),ggplot(plot_frame[plot_frame$gt != 0.5,],aes(y = gexpr,fill = as.factor(gt)))+
                     geom_boxplot()+
                     theme_void()+
                     scale_fill_manual(values = c("#377eb8","#e41a1c")) +
                     theme(legend.position = "none") +
                     ylab("Normalized expression") +
                     theme(legend.title = element_blank()),rel_widths = c(6,1),align = "h",axis="rtbl")
  #ggsave(paste0("spline_plots/sig_marker/",name,"chr",chr,".png"),plot=plt,device = "png")
  return(plt)
}

# Specify function for plotting the best marker found with the additive model.
plot_gene_time <- function(dev,row,col,chr){
  plot_frame <- data.frame(dev,use.rge[row,use.lines],t(unname(round(FourP$snp_mat[FourP$snp_info$ChrID==chr,][as.numeric(col),use.lines]*2,0)/2)), chr = chr)
  names(plot_frame) <- c("dev", "gexpr","gt")
  name <- gsub('\\.', '-', gene_names[row])
  plt <- plot_grid(ggplot(plot_frame,aes(dev,gexpr,col=as.factor(gt)))+
                     geom_point(size =3)+
                     geom_smooth(data = plot_frame[plot_frame$gt%in%c(0,1),],method="lm",size=1.5)+
                     scale_color_manual(values = c("0"="#377eb8","0.5"="#4daf4a","1"="#e41a1c"),name = "Genotype")+
                     #scale_color_viridis(option="H",discrete = T)+
                     scale_fill_manual(values = c("0" = "#377eb8", "1" = "#e41a1c"),guide = "none")+
                     theme_cowplot() +
                     ggtitle(paste0(paste("Gene",name),"\n",paste("Marker position",FourP$snp_info$SNP_position[FourP$snp_info$ChrID==chr][col],"Mb","chr",chr))) +
                     xlab("Estimated age (PC1)") + 
                     ylab("Gene expression")+
                     scale_size(name = "Gene expression")+
                     theme(text = element_text(size = 29),legend.position = "none",axis.text.x = element_text(size = 29),axis.text.y = element_text(size = 29)),ggplot(plot_frame[plot_frame$gt != 0.5,],aes(y = gexpr,fill = as.factor(gt)))+
                     geom_boxplot_pattern(pattern_density = 0.5,pattern_color = "white",pattern_fill= "white")+
                     theme_void()+
                     scale_fill_manual(values = c("#377eb8","#e41a1c")) +
                     theme(legend.position = "none") +
                     ylab("Normalized expression") +
                     theme(legend.title = element_blank()),rel_widths = c(6,1),align = "h",axis="rtbl")
  #ggsave(paste0("spline_plots/sig_marker/",name,"chr",chr,".png"),plot=plt,device = "png")
  return(plt)
}

# Extract legend we can paste to figure with plot_grid()
plot_frame <- data.frame(dev = pco1.lines,use.rge[which(gene_names==genes[1]),use.lines],t(unname(round(FourP$snp_mat[FourP$snp_info$ChrID==marker_chr[1],][as.numeric(markers_notime[1]),use.lines]*2,0)/2)), chr = marker_chr[1])
names(plot_frame) <- c("dev", "gexpr","gt")
legend <- get_legend(ggplot(plot_frame,aes(dev,gexpr,col=as.factor(gt)))+
                       geom_point(size =3)+
                       geom_smooth(data = plot_frame[plot_frame$gt%in%c(0,1),],method="lm",size=1.5)+
                       scale_color_manual(values = c("0"="#377eb8","0.5"="#4daf4a","1"="#e41a1c"),name = "Genotype")+
                       scale_fill_manual(values = c("0" = "blue", "1" = "red"),guide = "none")+
                       theme_cowplot() +
                       xlab("Developmental age") + 
                       ylab("Gene expression")+
                       scale_size(name = "Gene expression")+
                       theme(text = element_text(size = 32),legend.position = "top",
                             plot.title = element_text(size = 29),legend.direction = "horizontal",legend.justification="center" ,legend.box.just = "bottom"))

topleft_notime <- plot_gene(pco1.lines,which(gene_names == genes[1]),markers_notime[1],marker_chr[1])
topleft_time <- plot_gene_time(pco1.lines,which(gene_names == genes[1]),markers_time[1],marker_chr[1])

botright_notime <- plot_gene(pco1.lines,which(gene_names == genes[2]),markers_notime[2],marker_chr[2])
botright_time <- plot_gene_time(pco1.lines,which(gene_names == genes[2]),markers_time[2],marker_chr[2])

title_left <- ggdraw() + draw_label("AAM-only eQTL", fontface='bold', size = 33)
title_right <- ggdraw() + draw_label("SMM-only eQTL", fontface='bold', size = 33)

title_left_col <- ggdraw()+draw_label("Best marker SMM", fontface='bold', size = 33)
title_right_col <- ggdraw()+draw_label("Best marker AAM", fontface='bold', size = 33)

# gene_name_left <- ggdraw() + draw_label(genes[1], fontface='bold', size = 30)
# gene_name_right <- ggdraw() + draw_label(genes[2], fontface='bold', size = 30)

# legend in middle
png(filename = "plots_new_perm/marker_change.png",width = 1500, height = 1600)
plot_grid(plot_grid(title_left_col,title_right_col,nrow =1),plot_grid(plot_grid(title_left,plot_grid(topleft_notime+panel_border(colour = "black"),topleft_time+panel_border(colour = "black"),nrow =1,labels = c("A","B"),label_size =30),nrow =2,rel_heights = c(0.1,1)),legend,plot_grid(title_right,plot_grid(botright_notime+panel_border(colour = "black"),botright_time+panel_border(colour = "black"),labels = c("C","D"),label_size = 30, nrow = 1),nrow =2, rel_heights = c(0.1,1)),nrow =3,rel_heights = c(1,0.1,1)),nrow = 2, rel_heights = c(0.1,1))
dev.off()

pdf("pdfs_figures/Figure_7.pdf",width = 19, height = 20)
plot_grid(plot_grid(title_left_col,title_right_col,nrow =1),plot_grid(plot_grid(title_left,plot_grid(topleft_notime+panel_border(colour = "black"),topleft_time+panel_border(colour = "black"),nrow =1,labels = c("A","B"),label_size =30),nrow =2,rel_heights = c(0.1,1)),legend,plot_grid(title_right,plot_grid(botright_notime+panel_border(colour = "black"),botright_time+panel_border(colour = "black"),labels = c("C","D"),label_size = 30, nrow = 1),nrow =2, rel_heights = c(0.1,1)),nrow =3,rel_heights = c(1,0.1,1)),nrow = 2, rel_heights = c(0.1,1))
dev.off()







