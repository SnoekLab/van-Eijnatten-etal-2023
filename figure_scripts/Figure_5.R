library(ggplot2)
library(cowplot)
library(tidyverse)
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

# Read in p-values SMM and AAM.
pval_gen_notime_lm <- readRDS("~/Documents/elegans/pvals_correct/pval_gen_notime_lm.RDS")
pval_gen_time_lm <- readRDS("~/Documents/elegans/pvals_correct/pval_gen_time_lm.RDS")

# Collect lowest p-value per chromosome into frame for single marker model.

cutoff_notime <- 4.06 

min_p_no_time_pos_chr1 <- c()
min_p_no_time_pos_chr2 <- c()
min_p_no_time_pos_chr3 <- c()
min_p_no_time_pos_chr4 <- c()
min_p_no_time_pos_chr5 <- c()
min_p_no_time_pos_chrx <- c()

for(i in 1:ncol(pval_gen_notime_lm)){
  min_p_no_time_pos_chr1 <- c(min_p_no_time_pos_chr1,min(pval_gen_notime_lm[FourP$snp_info$ChrID=="I",i]))
}

for(i in 1:ncol(pval_gen_notime_lm)){
  min_p_no_time_pos_chr2 <- c(min_p_no_time_pos_chr2,min(pval_gen_notime_lm[FourP$snp_info$ChrID=="II",i])) 
}

for(i in 1:ncol(pval_gen_notime_lm)){
  min_p_no_time_pos_chr3 <- c(min_p_no_time_pos_chr3,min(pval_gen_notime_lm[FourP$snp_info$ChrID=="III",i]))
}

for(i in 1:ncol(pval_gen_notime_lm)){
  min_p_no_time_pos_chr4 <- c(min_p_no_time_pos_chr4,min(pval_gen_notime_lm[FourP$snp_info$ChrID=="IV",i]))
}

for(i in 1:ncol(pval_gen_notime_lm)){
  min_p_no_time_pos_chr5 <- c(min_p_no_time_pos_chr5,min(pval_gen_notime_lm[FourP$snp_info$ChrID=="V",i]))
}

for(i in 1:ncol(pval_gen_notime_lm)){
  min_p_no_time_pos_chrx <- c(min_p_no_time_pos_chrx,min(pval_gen_notime_lm[FourP$snp_info$ChrID=="X",i]))
}

min_p_no_time_pos_chr1 <- data.frame(pval = min_p_no_time_pos_chr1)
min_p_no_time_pos_chr2 <- data.frame(pval = min_p_no_time_pos_chr2)
min_p_no_time_pos_chr3 <- data.frame(pval = min_p_no_time_pos_chr3)
min_p_no_time_pos_chr4 <- data.frame(pval = min_p_no_time_pos_chr4)
min_p_no_time_pos_chr5 <- data.frame(pval = min_p_no_time_pos_chr5)
min_p_no_time_pos_chrx <- data.frame(pval = min_p_no_time_pos_chrx)

min_p_max_per_chr <- rbind(min_p_no_time_pos_chr1,min_p_no_time_pos_chr2,min_p_no_time_pos_chr3,min_p_no_time_pos_chr4,min_p_no_time_pos_chr5,min_p_no_time_pos_chrx)

# Collect lowest p-value for additive development model into frame.

cutoff_gen_time <- 4.44 

min_p_time_pos_chr1 <- c()
min_p_time_pos_chr2 <- c()
min_p_time_pos_chr3 <- c()
min_p_time_pos_chr4 <- c()
min_p_time_pos_chr5 <- c()
min_p_time_pos_chrx <- c()

for(i in 1:ncol(pval_gen_time_lm)){
  min_p_time_pos_chr1 <- c(min_p_time_pos_chr1,min(pval_gen_time_lm[FourP$snp_info$ChrID=="I",i]))
}

for(i in 1:ncol(pval_gen_time_lm)){
  min_p_time_pos_chr2 <- c(min_p_time_pos_chr2,min(pval_gen_time_lm[FourP$snp_info$ChrID=="II",i]))
}

for(i in 1:ncol(pval_gen_time_lm)){
  min_p_time_pos_chr3 <- c(min_p_time_pos_chr3,min(pval_gen_time_lm[FourP$snp_info$ChrID=="III",i]))
}

for(i in 1:ncol(pval_gen_time_lm)){
  min_p_time_pos_chr4 <- c(min_p_time_pos_chr4,min(pval_gen_time_lm[FourP$snp_info$ChrID=="IV",i]))
}

for(i in 1:ncol(pval_gen_time_lm)){
  min_p_time_pos_chr5 <- c(min_p_time_pos_chr5,min(pval_gen_time_lm[FourP$snp_info$ChrID=="V",i]))
}

for(i in 1:ncol(pval_gen_time_lm)){
  min_p_time_pos_chrx <- c(min_p_time_pos_chrx,min(pval_gen_time_lm[FourP$snp_info$ChrID=="X",i]))
}

min_p_time_pos_chr1 <- data.frame(marker = min_p_time_pos_chr1)
min_p_time_pos_chr2 <- data.frame(marker = min_p_time_pos_chr2)
min_p_time_pos_chr3 <- data.frame(marker = min_p_time_pos_chr3)
min_p_time_pos_chr4 <- data.frame(marker = min_p_time_pos_chr4)
min_p_time_pos_chr5 <- data.frame(marker = min_p_time_pos_chr5)
min_p_time_pos_chrx <- data.frame(marker = min_p_time_pos_chrx)

min_p_max_per_chr_time <- rbind(min_p_time_pos_chr1,min_p_time_pos_chr2,min_p_time_pos_chr3,min_p_time_pos_chr4,min_p_time_pos_chr5,min_p_time_pos_chrx)

# Define frame with pvals for both models.
scatter_frame <- data.frame(notime = -log10(min_p_max_per_chr$pval), time = -log10(min_p_max_per_chr_time$marker))

# Define additional variable which markers whether the p-value is significant for the both models, only single 
# marker, only additive development or for neither model.
scatter_frame <- scatter_frame %>% mutate(quadrant = case_when(notime >= cutoff_notime & time >= cutoff_gen_time   ~ "Both significant",
                                                               notime < cutoff_notime & time >= cutoff_gen_time  ~ "AAM-only eQTLs",
                                                               notime < cutoff_notime & time < cutoff_gen_time ~ "Neither significant",
                                                               TRUE ~ "SMM-only eQTLs"))

scatter_frame$quadrant <- factor(scatter_frame$quadrant, levels = c("Neither significant", "SMM-only eQTLs", "AAM-only eQTLs", "Both significant"))

# Check how many are significant only when the developmental age variable is included in the model.
dim(scatter_frame[scatter_frame$quadrant == "AAM-only eQTLs",])
dim(scatter_frame[scatter_frame$quadrant == "SMM-only eQTLs",])
dim(scatter_frame[scatter_frame$quadrant == "Neither significant",])
dim(scatter_frame[scatter_frame$quadrant == "Both significant",])

### Number of developmental effect eQTLs
dim(scatter_frame[scatter_frame$time-scatter_frame$notime > 10,])

# produce the scatter without axis limits.
scatter_full <- ggplot(data = scatter_frame,aes(x=notime,y=time)) + 
  geom_point(alpha = 0.6)+
  geom_abline(linewidth = 1) + 
  geom_abline(intercept = 10,slope = 1,linewidth = 1) +
  geom_segment(aes(x=0,y=10,xend =10,yend =10), col = "red",linewidth = 3 ) + 
  geom_segment(aes(x =10,y =0,xend =10,yend =10), col = "red",linewidth = 3) +
  theme_minimal() + 
  xlab("-log10(P) SMM") + 
  ylab("-log10(P) AAM") + 
  theme(text = element_text(size = 40),
        panel.border = element_rect(size=0.2,color="black",fill=NA),
        axis.text.x = element_text(size = 40),
        axis.text.y = element_text(size = 40),
        axis.ticks.x = element_line(size = 0.2),
        axis.ticks.y = element_line(size=0.2),
        panel.grid.major = element_line(size = 0.2,color = "grey")) +
  scale_x_continuous(expand =c(0,0), limits = c(0,152),breaks = c(50,100,150))+
  scale_y_continuous(expand =c(0,0), limits = c(0,152),breaks = c(50,100,150))

#Get legend

legend <- get_legend(ggplot(data = scatter_frame,aes(x=notime,y=time,col = quadrant)) +
                       geom_point(alpha=0.3) +
                       geom_abline() + 
                       theme_minimal() +
                       theme(text = element_text(size = 37),legend.position = "top",legend.justification = "center",legend.direction = "horizontal") +
                       xlab("-log10(P) SMM") + 
                       ylab("-log10(P) AAM") +
                       xlim(0,10)+ylim(0,10) +
                       geom_segment(aes(x=0,y=10,xend =10,yend =10), col = "red") + 
                       geom_segment(aes(x = 10,y = 0,xend = 10,yend = 10), col = "red") +
                       geom_hline(yintercept = cutoff_gen_time,col = "black") +
                       scale_color_viridis(discrete = T,option = "H",direction = -1) +
                       guides(colour = guide_legend(override.aes = list(size=10),nrow=2,byrow=TRUE,title="Quadrant")))


# Produce zoomed in scatter.
scatter_quadrant <- ggplot(data = scatter_frame,aes(x=notime,y=time,col = quadrant)) +
  geom_point(alpha=0.3) +
  theme_minimal() +
  theme(text = element_text(size = 40)) +
  xlab("-log10(P) SMM") + 
  ylab("-log10(P) AAM") +
  geom_vline(xintercept = cutoff_notime ,col = "black") + 
  geom_hline(yintercept = cutoff_gen_time,col = "black") +
  scale_color_viridis(discrete = T,option = "H",direction = -1) +
  theme(legend.position = "none",
        panel.border = element_rect(size=0.2,color="black",fill=NA),
        axis.text.x = element_text(size = 40),
        axis.text.y = element_text(size = 40),
        axis.ticks.x = element_line(size = 0.2),
        axis.ticks.y = element_line(size=0.2)) +
  scale_x_continuous(expand =c(0,0),limits = c(0,10),breaks = c(2.5,5,7.5,10))+
  scale_y_continuous(expand =c(0,0),limits = c(0,10),breaks = c(2.5,5,7.5,10))

scatter.plot <- plot_grid(scatter_full,scatter_quadrant,align = "v",ncol = 2,labels = c("A","B"),label_size = 35)+ 
  theme(plot.margin = unit(c(70,50,50,0), "points"))

png("plots_new_perm/scatter_plot.png",width = 1700,height=850)
scatter.plot
dev.off()

pdf("pdfs_figures/Figure_5.pdf",width = 23, height = 12)
scatter.plot
dev.off()











