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

# Read in p-values spline, interaction term.
pval_int_spline <- readRDS("~/Documents/elegans/redo_test_and_perm/pval_int_spline_2model.RDS")

# Initialise vector with names of genes we want to plot.
genes <- c("epn-1","F47B7.2","T23F2.9,bus-8","dpy-23","lrp-1")

# Get position and chromosome of most significant marker for this gene.
gene_index <- which(gene_names%in%genes)
markers <- sapply(gene_index,function(x) which.min(pval_int_spline[,x]))
marker_pos <- as.character(FourP$snp_info$SNP_position[markers]*10^6)
marker_chr <- as.character(FourP$snp_info$ChrID[markers])

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
                     ggtitle(paste0(paste("Gene", name),",",paste(" marker",FourP$snp_info$SNP_position[col],"Mb","chr",chr))) +
                     xlab("Estimated age (PC1)") + 
                     ylab("Center log ratio gene expression")+
                     scale_size(name = "Gene expression")+
                     theme(text = element_text(size = 27),axis.text = element_text(size =45),legend.position = "none",plot.title = element_text(size = 20)),ggplot(plot_frame[plot_frame$gt != 0.5,],aes(y = gexpr,fill = as.factor(gt)))+
                     geom_boxplot()+
                     theme_void()+
                     scale_fill_manual(values = c("#377eb8","#e41a1c")) +
                     theme(legend.position = "none") +
                     ylab("Normalized expression") +
                     theme(legend.title = element_blank()),rel_widths = c(6,1),align = "h",axis="rtbl")
  #ggsave(paste0("spline_plots/sig_marker/",name,"chr",chr,".png"),plot=plt,device = "png")
  return(plt)
}

# Call functions and paste plots together using plot_grid()
gene_1 <- expr_time_spline(pco1.lines,gene_index[1],markers[1],marker_chr[1])
gene_1
gene_2 <- expr_time_spline(pco1.lines,gene_index[2],markers[2],marker_chr[2])
gene_2
gene_3 <- expr_time_spline(pco1.lines,gene_index[3],markers[3],marker_chr[3])
gene_3
gene_4 <- expr_time_spline(pco1.lines,gene_index[4],markers[4],marker_chr[4])
gene_4
gene_5 <- expr_time_spline(pco1.lines,gene_index[5],markers[5],marker_chr[5])
gene_5

### Get legend
legend_frame  <- data.frame(pco1.lines,use.rge[gene_index[1],use.lines],t(unname(round(FourP$snp_mat[as.numeric(markers[1]),use.lines]*2,0)/2)), chr = marker_chr[1])
names(legend_frame) <- c("dev", "gexpr","gt")
spline_predict <- data.frame(prediction = predict(lm(use.rge[gene_index[1],use.lines]~as.numeric(FourP$snp_mat[as.numeric(markers[1]),use.lines])*ns(pco1.lines,knots = quantile(pco1.lines)[2:4])),interval = "confidence"))
names(spline_predict) <- c("prediction","lower","upper")
legend_frame <- cbind(legend_frame,spline_predict)

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
                       theme(text = element_text(size = 27),legend.position = "right",
                             plot.title = element_text(size = 20),legend.direction = "vertical",legend.justification="center"))

### Get prediction line and confidence intervals according to natural spline model in frames for each gene/marker combination 

gene_1_frame <- data.frame(dev = pco1.lines,gexpr = use.rge[gene_index[1],use.lines], gt = c(t(unname(round(FourP$snp_mat[as.numeric(markers[1]),use.lines]*2,0)/2))), chr = marker_chr[1],gene = genes[1])
gene_1_predict <- data.frame(prediction = predict(lm(use.rge[gene_index[1],use.lines]~as.numeric(FourP$snp_mat[as.numeric(markers[1]),use.lines])*ns(pco1.lines,knots = quantile(pco1.lines)[2:4])),interval = "confidence"))
names(gene_1_predict) <- c("prediction","lower","upper")
gene_1_frame <- cbind(gene_1_frame,gene_1_predict)                        
                    
gene_2_frame <- data.frame(dev = pco1.lines,gexpr = use.rge[gene_index[2],use.lines], gt = c(t(unname(round(FourP$snp_mat[as.numeric(markers[2]),use.lines]*2,0)/2))), chr = marker_chr[2],gene = genes[2])
gene_2_predict <- data.frame(prediction = predict(lm(use.rge[gene_index[2],use.lines]~as.numeric(FourP$snp_mat[as.numeric(markers[2]),use.lines])*ns(pco1.lines,knots = quantile(pco1.lines)[2:4])),interval = "confidence"))
names(gene_2_predict) <- c("prediction","lower","upper")
gene_2_frame <- cbind(gene_2_frame,gene_2_predict)        

gene_3_frame <- data.frame(dev = pco1.lines,gexpr = use.rge[gene_index[3],use.lines], gt = c(t(unname(round(FourP$snp_mat[as.numeric(markers[3]),use.lines]*2,0)/2))), chr = marker_chr[3],gene = genes[3])
gene_3_predict <- data.frame(prediction = predict(lm(use.rge[gene_index[3],use.lines]~as.numeric(FourP$snp_mat[as.numeric(markers[3]),use.lines])*ns(pco1.lines,knots = quantile(pco1.lines)[2:4])),interval = "confidence"))
names(gene_3_predict) <- c("prediction","lower","upper")
gene_3_frame <- cbind(gene_3_frame,gene_3_predict)   

gene_4_frame <- data.frame(dev = pco1.lines,gexpr = use.rge[gene_index[4],use.lines], gt = c(t(unname(round(FourP$snp_mat[as.numeric(markers[4]),use.lines]*2,0)/2))), chr = marker_chr[4],gene = genes[4])
gene_4_predict <- data.frame(prediction = predict(lm(use.rge[gene_index[4],use.lines]~as.numeric(FourP$snp_mat[as.numeric(markers[4]),use.lines])*ns(pco1.lines,knots = quantile(pco1.lines)[2:4])),interval = "confidence"))
names(gene_4_predict) <- c("prediction","lower","upper")
gene_4_frame <- cbind(gene_4_frame,gene_4_predict)   

gene_5_frame <- data.frame(dev = pco1.lines,gexpr = use.rge[gene_index[5],use.lines], gt = c(t(unname(round(FourP$snp_mat[as.numeric(markers[5]),use.lines]*2,0)/2))), chr = marker_chr[5],gene = genes[5])
gene_5_predict <- data.frame(prediction = predict(lm(use.rge[gene_index[5],use.lines]~as.numeric(FourP$snp_mat[as.numeric(markers[5]),use.lines])*ns(pco1.lines,knots = quantile(pco1.lines)[2:4])),interval = "confidence"))
names(gene_5_predict) <- c("prediction","lower","upper")
gene_5_frame <- cbind(gene_5_frame,gene_5_predict)   

### Bind in single frame
Condensed_frame <- rbind(gene_1_frame,gene_2_frame,gene_3_frame,gene_4_frame,gene_5_frame)

### Filter out heterozygous
Condensed_frame <- Condensed_frame[Condensed_frame$gt%in%c(0,1),] %>% 
  mutate(uniq_id=paste(gt,gene,sep=""))

Condensed_frame$gene <- factor(Condensed_frame$gene, levels = genes)

mapping <- Condensed_frame %>% distinct(gt,uniq_id)
cols <- c("#377eb8","#e41a1c")
names(cols) = unique(mapping$gt)
plot_cols = cols[as.factor(mapping$gt)]
names(plot_cols) = mapping$uniq_id
Condensed_frame$gt <-as.factor(Condensed_frame$gt)


# Plot with different colors per gene. Extract to separate legends and add using cowplot.

plot_col_blue <- c("#accbff","#92bbff","#78aaff","#649eff","#4188ff")
names(plot_col_blue) <- genes

# Get legend blue genes
legend_blue <- get_legend(ggplot(Condensed_frame)+
  geom_line(data = Condensed_frame[Condensed_frame$gt%in%c(0),],aes(x = dev, y = prediction, color = as.factor(gene)),size = 1.5)+
  geom_ribbon(data = Condensed_frame[Condensed_frame$gt%in%c(0),],aes(x = dev,ymin=lower,ymax = upper,fill = as.factor(uniq_id)), alpha = 0.05)+
  scale_color_manual(values = plot_col_blue,name = "Allele 0")+
  scale_fill_manual(values = plot_cols,guide = "none")+
  theme_cowplot() +
  xlab("Estimated age (PC1)") + 
  ylab("Log2 ratio gene expression")+
  scale_size(name = "Gene expression")+
  theme(text = element_text(size = 27)))

plot_col_red <- c("#ff0000","#d70000","#c60000","#b70000","#9b0000")
names(plot_col_red) <- genes

# Get legend red genes.
legend_red <- get_legend(ggplot(Condensed_frame)+
  geom_line(data = Condensed_frame[Condensed_frame$gt%in%c(1),],aes(x = dev, y = prediction, color = as.factor(gene)),size = 1.5)+
  geom_ribbon(data = Condensed_frame[Condensed_frame$gt%in%c(1),],aes(x = dev,ymin=lower,ymax = upper,fill = as.factor(uniq_id)), alpha = 0.05)+
  scale_color_manual(values = plot_col_red ,name = "Allele 1")+
  scale_fill_manual(values = plot_cols,guide = "none")+
  theme_cowplot() +
  xlab("Estimated age (PC1)") + 
  ylab("Log2 ratio gene expression")+
  scale_size(name = "Gene expression")+
  theme(text = element_text(size = 27)))

plot_cols2 <- c("#accbff","#ff0000","#92bbff","#d70000","#78aaff","#c60000","#649eff","#b70000","#4188ff","#9b0000")
names(plot_cols2) <- mapping$uniq_id

cool.plot <- plot_grid(ggplot(Condensed_frame)+
  geom_line(data = Condensed_frame[Condensed_frame$gt%in%c(0,1),],aes(x = dev, y = prediction, color = as.factor(uniq_id)),size = 1.5)+
  geom_ribbon(data = Condensed_frame[Condensed_frame$gt%in%c(0,1),],aes(x = dev,ymin=lower,ymax = upper,fill = as.factor(uniq_id)), alpha = 0.05)+
  scale_color_manual(values = plot_cols2,guide = "none")+
  scale_fill_manual(values = plot_cols2,guide = "none")+
  theme_cowplot() +
  xlab("Estimated age (PC1)") + 
  ylab("Center log ratio gene expression")+
  scale_size(name = "Gene expression")+
  theme(text = element_text(size = 27),panel.border = element_rect(linewidth = 0.2,fill = NA),
        axis.text = element_text(size =27)),
  plot_grid(legend_blue,legend_red,nrow = 2), nrow=1,rel_widths = c(8,3))

png(filename = "plots_new_perm/similar_int.png",width = 1100,height = 600)
cool.plot
dev.off()

pdf("pdfs_figures/Figure_A2_2.pdf",width = 14, height = 8)
cool.plot
dev.off()











