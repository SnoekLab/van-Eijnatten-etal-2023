
library(ggplot2)
library(cowplot)
library(dplyr)
library(viridis)
library(lemon)

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

# Read in p-values.
pval_int_spline <- readRDS("~/Documents/elegans/redo_test_and_perm/pval_int_spline_2model.RDS")

# Initialise vector with names of genes we want to plot.
genes <- c("bli-3","ZK1225.4","21ur-15524,ZK1053.4,ZK1053.6","Y20F4.5")

# Get position and chromosome of most significant marker for these genes.
gene_index <- which(gene_names%in%genes)
markers <- sapply(gene_index,function(x) which.min(pval_int_spline[,x]))
marker_pos <- as.character(FourP$snp_info$SNP_position[markers]*10^6)
marker_chr <- as.character(FourP$snp_info$ChrID[markers])

# Define function for plotting.
expr_time_spline <- function(dev,row,col,chr){
  plot_frame <- data.frame(dev,use.rge[row,use.lines],t(unname(round(FourP$snp_mat[as.numeric(col),use.lines]*2,0)/2)), chr = chr)
  names(plot_frame) <- c("dev", "gexpr","gt")
  name <- gsub('\\.', '-', gene_names[row])
  spline_predict <- data.frame(prediction = predict(lm(use.rge[row,use.lines]~as.numeric(FourP$snp_mat[as.numeric(col),use.lines])*ns(dev,knots = quantile(dev)[2:4])),interval = "confidence"))
  names(spline_predict) <- c("prediction","lower","upper")
  plot_frame <- cbind(plot_frame,spline_predict)
  plt <- plot_grid(ggplot(plot_frame,aes(dev,gexpr,col=as.factor(gt)))+
                     geom_point(size =3)+
                     geom_line(data = plot_frame[plot_frame$gt%in%c(0,1),],aes(x = dev, y = prediction, color = as.factor(gt)),size = 1.5)+
                     geom_ribbon(data = plot_frame[plot_frame$gt%in%c(0,1),],aes(x = dev,ymin=lower,ymax = upper,fill = as.factor(gt)), alpha = 0.05)+
                     scale_color_manual(values = c("0"="#377eb8","0.5"="#4daf4a","1"="#e41a1c"),name = "Genotype")+
                     #scale_color_viridis(option="H",discrete = T)+
                     scale_fill_manual(values = c("0" = "#377eb8", "1" = "#e41a1c"),guide = "none")+
                     theme_cowplot() +
                     ggtitle(paste0(paste("Gene", name),"\n",paste("Marker",FourP$snp_info$SNP_position[col],"Mb","chr",chr))) +
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

# Call function and paste plots together using plot_grid()
gene_1 <- expr_time_spline(pco1.lines,gene_index[1],markers[1],marker_chr[1])
gene_1
gene_2 <- expr_time_spline(pco1.lines,gene_index[2],markers[2],marker_chr[2])
gene_2
gene_3 <- expr_time_spline(pco1.lines,gene_index[3],markers[3],marker_chr[3])
gene_3
gene_4 <- expr_time_spline(pco1.lines,gene_index[4],markers[4],marker_chr[4])
gene_4

legend_frame  <- data.frame(pco1.lines,use.rge[gene_index[1],use.lines],t(unname(round(FourP$snp_mat[as.numeric(markers[1]),use.lines]*2,0)/2)), chr = marker_chr[1])
names(legend_frame) <- c("dev", "gexpr","gt")
spline_predict <- data.frame(prediction = predict(lm(use.rge[gene_index[1],use.lines]~as.numeric(FourP$snp_mat[as.numeric(markers[1]),use.lines])*ns(pco1.lines,knots = quantile(pco1.lines)[2:4])),interval = "confidence"))
names(spline_predict) <- c("prediction","lower","upper")
legend_frame <- cbind(legend_frame,spline_predict)

# Extract legend we want to paste to the figures with get_legend()
legend <- get_legend(ggplot(legend_frame,aes(dev,gexpr,col=as.factor(gt)))+
                       geom_point(size =3)+
                       geom_line(data = legend_frame[legend_frame$gt%in%c(0,1),],aes(x = dev, y = prediction, color = as.factor(gt)),size = 1.5)+
                       geom_ribbon(data = legend_frame[legend_frame$gt%in%c(0,1),],aes(x = dev,ymin=lower,ymax = upper,fill = as.factor(gt)), alpha = 0.05)+
                       scale_color_manual(values = c("0"="#377eb8","0.5"="#4daf4a","1"="#e41a1c"),name = "Genotype")+
                       scale_fill_manual(values = c("0" = "blue", "1" = "red"),guide = "none")+
                       theme_cowplot() +
                       ggtitle(paste0(paste("Gene", genes[1]),",",paste(" marker",FourP$snp_info$SNP_position[markers[1]],"chr",marker_chr[1]))) +
                       xlab("Developmental age") + 
                       ylab("Log2 ratio gene expression")+
                       scale_size(name = "Gene expression")+
                       theme(text = element_text(size = 50),legend.position = "top",
                             plot.title = element_text(size = 20),legend.direction = "horizontal",legend.justification="center" ,legend.box.just = "bottom"))

# Produce figure.

png(filename = "plots_new_perm/categ_int_spline.png",height = 2000,width =2000)
plot_grid(legend,plot_grid(gene_2,gene_4,gene_1,gene_3,nrow = 2,align = "h",labels = c("A","B","C","D"),
                    axis="rtbl",label_size = 40),nrow=2,rel_heights = c(1,5),align = "v",axis="tb",hjust = 0.5)+
  theme(plot.margin = unit(c(20,20,20,20), "points"))
dev.off()

pdf("pdfs_figures/Figure_A1_3.pdf",width = 27, height = 27)
plot_grid(legend,plot_grid(gene_2,gene_4,gene_1,gene_3,nrow = 2,align = "h",labels = c("A","B","C","D"),
                           axis="rtbl",label_size = 40),nrow=2,rel_heights = c(1,5),align = "v",axis="tb",hjust = 0.5)+
  theme(plot.margin = unit(c(20,20,20,20), "points"))
dev.off()









