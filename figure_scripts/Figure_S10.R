library(ggplot2)
library(cowplot)
library(dplyr)
library(viridis)
library(lemon)
library(tidyr)

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

# Get eQTLs for SMM and AAM

cutoff_notime <- 4.06

min_p_no_time_pos_chr1 <- c()
min_p_no_time_pos_chr2 <- c()
min_p_no_time_pos_chr3 <- c()
min_p_no_time_pos_chr4 <- c()
min_p_no_time_pos_chr5 <- c()
min_p_no_time_pos_chrx <- c()

gene_1_notime <-c()
gene_2_notime <-c()
gene_3_notime <-c()
gene_4_notime <-c()
gene_5_notime <-c()
gene_x_notime <-c()

pval_1_notime <-c()
pval_2_notime <-c()
pval_3_notime <-c()
pval_4_notime <-c()
pval_5_notime <-c()
pval_x_notime <-c()

position_notime_1 <-c()
position_notime_2 <-c()
position_notime_3 <-c()
position_notime_4 <-c()
position_notime_5 <-c()
position_notime_x <-c()

for(i in 1:ncol(pval_gen_notime_lm)){
  min_p_no_time_pos_chr1 <- c(min_p_no_time_pos_chr1,which.min(pval_gen_notime_lm[FourP$snp_info$ChrID=="I",i]))
  gene_1_notime <- c(gene_1_notime,i)
  position_notime_1 <- c(position_notime_1,FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="I"][which.min(pval_gen_notime_lm[FourP$snp_info$ChrID=="I",i])])
  pval_1_notime <- c(pval_1_notime,min(pval_gen_notime_lm[FourP$snp_info$ChrID=="I",i]))
}

for(i in 1:ncol(pval_gen_notime_lm)){
  min_p_no_time_pos_chr2 <- c(min_p_no_time_pos_chr2,which.min(pval_gen_notime_lm[FourP$snp_info$ChrID=="II",i])) 
  gene_2_notime <- c(gene_2_notime,i)
  position_notime_2 <- c(position_notime_2,FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="II"][which.min(pval_gen_notime_lm[FourP$snp_info$ChrID=="II",i])])
  pval_2_notime <- c(pval_2_notime,min(pval_gen_notime_lm[FourP$snp_info$ChrID=="II",i]))
}

for(i in 1:ncol(pval_gen_notime_lm)){
  min_p_no_time_pos_chr3 <- c(min_p_no_time_pos_chr3,which.min(pval_gen_notime_lm[FourP$snp_info$ChrID=="III",i]))
  position_notime_3 <- c(position_notime_3,FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="III"][which.min(pval_gen_notime_lm[FourP$snp_info$ChrID=="III",i])])
  pval_3_notime <- c(pval_3_notime,min(pval_gen_notime_lm[FourP$snp_info$ChrID=="III",i]))
  gene_3_notime <- c(gene_3_notime,i)
}

for(i in 1:ncol(pval_gen_notime_lm)){
  min_p_no_time_pos_chr4 <- c(min_p_no_time_pos_chr4,which.min(pval_gen_notime_lm[FourP$snp_info$ChrID=="IV",i]))
  position_notime_4 <- c(position_notime_4,FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="IV"][which.min(pval_gen_notime_lm[FourP$snp_info$ChrID=="IV",i])])
  gene_4_notime <- c(gene_4_notime,i)
  pval_4_notime <- c(pval_4_notime,min(pval_gen_notime_lm[FourP$snp_info$ChrID=="IV",i]))
}

for(i in 1:ncol(pval_gen_notime_lm)){
  min_p_no_time_pos_chr5 <- c(min_p_no_time_pos_chr5,which.min(pval_gen_notime_lm[FourP$snp_info$ChrID=="V",i]))
  position_notime_5 <- c(position_notime_5,FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="V"][which.min(pval_gen_notime_lm[FourP$snp_info$ChrID=="V",i])])
  gene_5_notime <- c(gene_5_notime,i)
  pval_5_notime <- c(pval_5_notime,min(pval_gen_notime_lm[FourP$snp_info$ChrID=="V",i]))
}

for(i in 1:ncol(pval_gen_notime_lm)){
  min_p_no_time_pos_chrx <- c(min_p_no_time_pos_chrx,which.min(pval_gen_notime_lm[FourP$snp_info$ChrID=="X",i]))
  position_notime_x <- c(position_notime_x,FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="X"][which.min(pval_gen_notime_lm[FourP$snp_info$ChrID=="X",i])])
  gene_x_notime <- c(gene_x_notime,i)
  pval_x_notime <- c(pval_x_notime,min(pval_gen_notime_lm[FourP$snp_info$ChrID=="X",i]))
}

min_p_no_time_pos_chr1 <- data.frame(marker = min_p_no_time_pos_chr1, gene = gene_1_notime, chr = rep("I",ncol(pval_gen_notime_lm)),pval = pval_1_notime,position = position_notime_1)
min_p_no_time_pos_chr2 <- data.frame(marker = min_p_no_time_pos_chr2, gene = gene_2_notime, chr = rep("II",ncol(pval_gen_notime_lm)),pval = pval_2_notime,position = position_notime_2)
min_p_no_time_pos_chr3 <- data.frame(marker = min_p_no_time_pos_chr3, gene = gene_3_notime, chr = rep("III",ncol(pval_gen_notime_lm)),pval = pval_3_notime,position = position_notime_3)
min_p_no_time_pos_chr4 <- data.frame(marker = min_p_no_time_pos_chr4, gene = gene_4_notime, chr = rep("IV",ncol(pval_gen_notime_lm)),pval = pval_4_notime,position = position_notime_4)
min_p_no_time_pos_chr5 <- data.frame(marker = min_p_no_time_pos_chr5, gene = gene_5_notime, chr = rep("V",ncol(pval_gen_notime_lm)),pval = pval_5_notime,position = position_notime_5)
min_p_no_time_pos_chrx <- data.frame(marker = min_p_no_time_pos_chrx, gene = gene_x_notime, chr = rep("X",ncol(pval_gen_notime_lm)),pval = pval_x_notime,position = position_notime_x)

marker_notime <- rbind(min_p_no_time_pos_chr1,min_p_no_time_pos_chr2,min_p_no_time_pos_chr3,min_p_no_time_pos_chr4,min_p_no_time_pos_chr5,min_p_no_time_pos_chrx)


# Collect lowest p-value for additive development model into frame.

cutoff_gen_time <- 4.44

min_p_time_pos_chr1 <- c()
min_p_time_pos_chr2 <- c()
min_p_time_pos_chr3 <- c()
min_p_time_pos_chr4 <- c()
min_p_time_pos_chr5 <- c()
min_p_time_pos_chrx <- c()

gene_1_time <-c()
gene_2_time <-c()
gene_3_time <-c()
gene_4_time <-c()
gene_5_time <-c()
gene_x_time <-c()

pval_1_time <-c()
pval_2_time <-c()
pval_3_time <-c()
pval_4_time <-c()
pval_5_time <-c()
pval_x_time <-c()

position_time_1 <-c()
position_time_2 <-c()
position_time_3 <-c()
position_time_4 <-c()
position_time_5 <-c()
position_time_x <-c()

for(i in 1:ncol(pval_gen_time_lm)){
  min_p_time_pos_chr1 <- c(min_p_time_pos_chr1,which.min(pval_gen_time_lm[FourP$snp_info$ChrID=="I",i]))
  gene_1_time <-c(gene_1_time,i)
  pval_1_time <- c(pval_1_time,min(pval_gen_time_lm[FourP$snp_info$ChrID=="I",i]))
  position_time_1 <- c(position_time_1,FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="I"][which.min(pval_gen_time_lm[FourP$snp_info$ChrID=="I",i])])
}

for(i in 1:ncol(pval_gen_time_lm)){
  min_p_time_pos_chr2 <- c(min_p_time_pos_chr2,which.min(pval_gen_time_lm[FourP$snp_info$ChrID=="II",i]))
  gene_2_time <-c(gene_2_time,i)
  pval_2_time <- c(pval_2_time,min(pval_gen_time_lm[FourP$snp_info$ChrID=="II",i]))
  position_time_2 <- c(position_time_2,FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="II"][which.min(pval_gen_time_lm[FourP$snp_info$ChrID=="II",i])])
}

for(i in 1:ncol(pval_gen_time_lm)){
  min_p_time_pos_chr3 <- c(min_p_time_pos_chr3,which.min(pval_gen_time_lm[FourP$snp_info$ChrID=="III",i]))
  gene_3_time <-c(gene_3_time,i)
  pval_3_time <- c(pval_3_time,min(pval_gen_time_lm[FourP$snp_info$ChrID=="III",i]))
  position_time_3 <- c(position_time_3,FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="III"][which.min(pval_gen_time_lm[FourP$snp_info$ChrID=="III",i])])
}

for(i in 1:ncol(pval_gen_time_lm)){
  min_p_time_pos_chr4 <- c(min_p_time_pos_chr4,which.min(pval_gen_time_lm[FourP$snp_info$ChrID=="IV",i]))
  gene_4_time <-c(gene_4_time,i)
  pval_4_time <- c(pval_4_time,min(pval_gen_time_lm[FourP$snp_info$ChrID=="IV",i]))
  position_time_4 <- c(position_time_4,FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="IV"][which.min(pval_gen_time_lm[FourP$snp_info$ChrID=="IV",i])])
}

for(i in 1:ncol(pval_gen_time_lm)){
  min_p_time_pos_chr5 <- c(min_p_time_pos_chr5,which.min(pval_gen_time_lm[FourP$snp_info$ChrID=="V",i]))
  gene_5_time <-c(gene_5_time,i)
  pval_5_time <- c(pval_5_time,min(pval_gen_time_lm[FourP$snp_info$ChrID=="V",i]))
  position_time_5 <- c(position_time_5,FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="V"][which.min(pval_gen_time_lm[FourP$snp_info$ChrID=="V",i])])
}

for(i in 1:ncol(pval_gen_time_lm)){
  min_p_time_pos_chrx <- c(min_p_time_pos_chrx,which.min(pval_gen_time_lm[FourP$snp_info$ChrID=="X",i]))
  gene_x_time <-c(gene_x_time,i)
  pval_x_time <- c(pval_x_time,min(pval_gen_time_lm[FourP$snp_info$ChrID=="X",i]))
  position_time_x <- c(position_time_x,FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="X"][which.min(pval_gen_time_lm[FourP$snp_info$ChrID=="X",i])])
}

min_p_time_pos_chr1 <- data.frame(marker = min_p_time_pos_chr1, gene = gene_1_time, chr = rep("I",ncol(pval_gen_time_lm)),pval = pval_1_time,position = position_time_1)
min_p_time_pos_chr2 <- data.frame(marker = min_p_time_pos_chr2, gene = gene_2_time, chr = rep("II",ncol(pval_gen_time_lm)),pval = pval_2_time,position = position_time_2)
min_p_time_pos_chr3 <- data.frame(marker = min_p_time_pos_chr3, gene = gene_3_time, chr = rep("III",ncol(pval_gen_time_lm)),pval = pval_3_time,position = position_time_3)
min_p_time_pos_chr4 <- data.frame(marker = min_p_time_pos_chr4, gene = gene_4_time, chr = rep("IV",ncol(pval_gen_time_lm)),pval = pval_4_time,position = position_time_4)
min_p_time_pos_chr5 <- data.frame(marker = min_p_time_pos_chr5, gene = gene_5_time, chr = rep("V",ncol(pval_gen_time_lm)),pval = pval_5_time,position = position_time_5)
min_p_time_pos_chrx <- data.frame(marker = min_p_time_pos_chrx, gene = gene_x_time, chr = rep("X",ncol(pval_gen_time_lm)),pval = pval_x_time,position = position_time_x)

marker_time <- rbind(min_p_time_pos_chr1,min_p_time_pos_chr2,min_p_time_pos_chr3,min_p_time_pos_chr4,min_p_time_pos_chr5,min_p_time_pos_chrx)


# Define frame with pvals for both models.
scatter_frame <- data.frame(pval_notime = -log10(marker_notime$pval), pval_time = -log10(marker_time$pval), position_notime = marker_notime$position,position_time = marker_time$position,chr = marker_notime$chr)

# Define additional variable which markers whether the p-value is significant for the both models, only single 
# marker, only additive development or for neither model.
scatter_frame <- scatter_frame %>% mutate(quadrant = case_when(pval_notime >= cutoff_notime & pval_time >= cutoff_gen_time   ~ "Both significant",
                                                               pval_notime < cutoff_notime & pval_time >= cutoff_gen_time  ~ "Only w/ time significant",
                                                               pval_notime < cutoff_notime & pval_time < cutoff_gen_time ~ "Neither significant",
                                                               TRUE ~ "Only w/o time significant"))


### Make extra eQTLs with -log10(P) = 0 to get the x-axis range right

zero_eqtl_1 <- data.frame(pval_notime = rep(0,18),
                        pval_time = rep(0,18),
                        position_notime = rep(0,18),
                        position_time = rep(2,18),
                        chr = rep(c("I",'II',"III","IV","V","X"),3),
                        quadrant = c(rep("Only w/ time significant",6),rep("Only w/o time significant",6),rep("Both significant",6)))
zero_eqtl_2 <- data.frame(pval_notime = rep(0,18),
                          pval_time = rep(0,18),
                          position_notime = rep(2,18),
                          position_time = rep(0,18),
                          chr = rep(c("I",'II',"III","IV","V","X"),3),
                          quadrant = c(rep("Only w/ time significant",6),rep("Only w/o time significant",6),rep("Both significant",6)))

chr_3_eqtl_1 <- data.frame(pval_notime = rep(0,6),
                         pval_time = rep(0,6),
                         position_notime = rep(0,6),
                         position_time = c(tail(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="I"])[5],tail(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="II"])[5],tail(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="III"])[5],tail(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="IV"])[5],tail(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="V"])[5],tail(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="X"])[5]),
                         chr = c("I","II","III","IV","V","X"),
                         quadrant = rep("Only w/o time significant",6))

chr_3_eqtl_2 <- data.frame(pval_notime = rep(0,6),
                         pval_time = rep(0,6),
                         position_notime = c(tail(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="I"])[5],tail(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="II"])[5],tail(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="III"])[5],tail(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="IV"])[5],tail(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="V"])[5],tail(FourP$snp_info$SNP_position[FourP$snp_info$ChrID=="X"])[5]),
                         position_time = rep(0,6),
                         chr = c("I","II","III","IV","V","X"),
                         quadrant = rep("Only w/o time significant",6))

# Bind to eQTL frame
scatter_frame <- rbind(scatter_frame,zero_eqtl_1,zero_eqtl_2,chr_3_eqtl_1,chr_3_eqtl_2)

legend_frame <- scatter_frame %>% pivot_longer(cols = c("position_notime","position_time"),names_to = "Model")

legend <- get_legend(ggplot(legend_frame[legend_frame$quadrant=="Only w/ time significant",],aes(x = value,fill = Model)) +
             geom_histogram(bins = 200) +
             ylab("eQTL count") +
             xlab("Position on chromosome in Mbp") +
             facet_grid(Model~chr,scale = "free_x") +
             theme_minimal() +
             theme(panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank(),strip.text.y = element_blank(),text = element_text(size = 35),axis.line = element_line(size = 0.2),axis.ticks= element_line(size = 0.2),axis.title.x = element_blank(),panel.spacing.x = unit(6, "mm"),legend.position = "top",legend.direction = "horizontal") +
             panel_border(color = "black",size = 0.2)+
             scale_x_continuous(expand =c(0,0))+
             scale_fill_manual(values =  c("blue","red"),labels = c("Single marker model","Additive model")))

topleft_hist <- plot_grid(legend,plot_grid(topleft_notime,topleft_time,nrow = 2),nrow =2, rel_heights = c(0.1,1))

### Make histogram of SMM-only eQTLs
botright_notime <- ggplot(scatter_frame[abs(scatter_frame$position_notime-scatter_frame$position_time)>1&scatter_frame$quadrant=="Only w/o time significant",], aes(x = position_notime)) +
  geom_histogram(bins = 200,fill = "blue") +
  ylab("eQTL count") +
  xlab("Position on chromosome in Mbp") +
  facet_grid(~chr,scale = "free_x",space = "free_x") +
  scale_fill_manual(values = c("#7A0403FF")) +
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank(),strip.text.y = element_blank(),legend.position = "none",text = element_text(size = 35),axis.line = element_line(size = 0.2),axis.ticks= element_line(size = 0.2),axis.title.x = element_blank(),panel.spacing.x = unit(12, "mm")) +
  panel_border(color = "black",size = 0.2)+
  scale_x_continuous(expand =c(0,0)) + 
  scale_y_continuous(breaks = c(0,300,600,900))

### Make histogram of best marker on same chromosome for same genes according to AAM.
botright_time <- ggplot(scatter_frame[abs(scatter_frame$position_notime-scatter_frame$position_time)>1&scatter_frame$quadrant=="Only w/o time significant",], aes(x = position_time)) +
  geom_histogram(bins = 200,fill = "red") +
  ylab("eQTL count") +
  xlab("Position on chromosome in Mbp") +
  facet_grid(~chr,scale = "free_x",space = "free_x") +
  scale_fill_manual(values = c("blue3")) +
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank(),strip.text.y = element_blank(),legend.position = "none",text = element_text(size = 35),axis.line = element_line(size = 0.2),axis.ticks= element_line(size = 0.2),panel.spacing.x = unit(12, "mm")) +
  panel_border(color = "black",size = 0.2) + 
  scale_x_continuous(expand =c(0,0))

### Paste figures together
botright_hist <- plot_grid(legend,plot_grid(botright_notime,botright_time,nrow = 2,align = "vh",axis ="lrtb"),nrow =2, rel_heights = c(0.1,1))

png(filename = "plots_new_perm/hist_botright.png",width = 1500,height = 1000)
botright_hist
dev.off()

pdf("pdfs_figures/Figure_S10.pdf",width = 15, height = 10)
botright_hist
dev.off()








