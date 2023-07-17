
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

# Get gene expression and genotype for the gene/marker we are interested in.
gno <- grep("C04F12.6,C36F7.5",fpkm$genes[fpkm.selc])
mpgexp <- use.rge[gno,use.lines]
mrk <- as.numeric(round(FourP$snp_mat["9710411",use.lines]*2,0)/2)

# Run a linear additive age model.
model <- summary(lm(mpgexp~pco1.lines+mrk))

# Remove developmental variation using the slope over development as calculated in model.
plot_frame <- data.frame(dev = pco1.lines, gexpr_nodev = mpgexp-pco1.lines*model$coefficients[2,1], gexpr = mpgexp,gt = mrk)

# Produce plot without removing gene expression.
normal_plot <- plot_grid(ggplot(plot_frame,aes(dev,gexpr,col=as.factor(gt)))+
            geom_point(size =3)+
            geom_smooth(data = plot_frame[plot_frame$gt%in%c(0,1),],method="lm",size=1.5)+
            scale_color_manual(values = c("0"="#377eb8","0.5"="#4daf4a","1"="#e41a1c"),name = "Genotype")+
            scale_fill_manual(values = c("0" = "#377eb8", "1" = "#e41a1c"),guide = "none")+
            theme_cowplot() +
            ggtitle(paste0(paste("Gene", "C04F12.6,C36F7.5"),"\n",paste("Marker","9710411","Mb","chr","I"))) +
            xlab("Estimated age (PC1)") + 
            ylab("Center log ratio gene expression")+
            scale_size(name = "Gene expression")+
            theme(text = element_text(size = 35),axis.text = element_text(size =45),legend.position = "none",legend.title = element_text(hjust = 0.5)),ggplot(plot_frame[plot_frame$gt != 0.5,],aes(y = gexpr,fill = as.factor(gt)))+
            geom_boxplot()+
            theme_void()+
            scale_fill_manual(values = c("#377eb8","#e41a1c")) +
            theme(legend.position = "none") +
            ylab("Normalized expression") +
            theme(legend.title = element_blank()),rel_widths = c(6,1),align = "h",axis="rtbl")

# After removing gene expression.
remove_dev <- plot_grid(ggplot(plot_frame,aes(dev,gexpr_nodev,col=as.factor(gt)))+
                          geom_point(size =3)+
                          geom_smooth(data = plot_frame[plot_frame$gt%in%c(0,1),],method="lm",size=1.5)+
                          scale_color_manual(values = c("0"="#377eb8","0.5"="#4daf4a","1"="#e41a1c"),name = "Genotype")+
                          scale_fill_manual(values = c("0" = "#377eb8", "1" = "#e41a1c"),guide = "none")+
                          theme_cowplot() +
                          ggtitle(paste0(paste("Gene", "C04F12.6,C36F7.5"),"\n",paste("Marker","9710411","Mb","chr","I"))) +
                          xlab("Estimated age (PC1)") + 
                          ylab("Center log ratio gene expression")+
                          scale_size(name = "Gene expression")+
                          theme(text = element_text(size = 35),axis.text = element_text(size =45),legend.position = "none",legend.title = element_text(hjust = 0.5)),ggplot(plot_frame[plot_frame$gt != 0.5,],aes(y = gexpr_nodev,fill = as.factor(gt)))+
                          geom_boxplot()+
                          theme_void()+
                          scale_fill_manual(values = c("#377eb8","#e41a1c")) +
                          theme(legend.position = "none") +
                          ylab("Normalized expression") +
                          theme(legend.title = element_blank()),rel_widths = c(6,1),align = "h",axis="rtbl")

# Extract legend.
legend <- get_legend(ggplot(plot_frame,aes(dev,gexpr,col=as.factor(gt)))+
                       geom_point(size =3)+
                       geom_smooth(data = plot_frame[plot_frame$gt%in%c(0,1),],method="lm",size=1.5)+
                       scale_color_manual(values = c("0"="#377eb8","0.5"="#4daf4a","1"="#e41a1c"),name = "Genotype")+
                       scale_fill_manual(values = c("0" = "blue", "1" = "red"),guide = "none")+
                       theme_cowplot() +
                       xlab("Developmental age") + 
                       ylab("Log2 ratio gene expression")+
                       scale_size(name = "Gene expression")+
                       theme(text = element_text(size = 40),legend.position = "top",
                             legend.direction = "horizontal",legend.justification="center" ,legend.box.just = "bottom"))

# Produce plot.

png(filename = "plots_new_perm/remove_dev.png",height = 750, width = 1500)
plot_grid(legend,plot_grid(normal_plot,remove_dev,ncol = 2,align = "v",labels = c("A","B"),
                    axis="rtbl",label_size = 30),nrow=2,rel_heights = c(1,5),hjust = 0.5)
dev.off()

pdf("pdfs_figures/Figure_S10.pdf",width = 22,height = 11)
plot_grid(legend,plot_grid(normal_plot,remove_dev,ncol = 2,align = "v",labels = c("A","B"),
                           axis="rtbl",label_size = 30),nrow=2,rel_heights = c(1,5),hjust = 0.5)
dev.off()


