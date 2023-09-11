library(forcats)
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

# Read in p-values SMM and AAM.
pval_gen_notime_lm <- readRDS("~/Documents/elegans/pvals_correct/pval_gen_notime_lm.RDS")
pval_gen_time_lm <- readRDS("~/Documents/elegans/pvals_correct/pval_gen_time_lm.RDS")

### Make a dataframe with the lowest p-values per chromosome for each gene.

pval_notime_chr1 <- c()
pval_notime_chr2 <- c()
pval_notime_chr3 <- c()
pval_notime_chr4 <- c()
pval_notime_chr5 <- c()
pval_notime_chrx <- c()

gene_notime_chr1 <- c()
gene_notime_chr2 <- c()
gene_notime_chr3 <- c()
gene_notime_chr4 <- c()
gene_notime_chr5 <- c()
gene_notime_chrx <- c()

marker_notime_chr1 <- c()
marker_notime_chr2 <- c()
marker_notime_chr3 <- c()
marker_notime_chr4 <- c()
marker_notime_chr5 <- c()
marker_notime_chrx <- c()

for(i in 1:ncol(pval_gen_notime_lm)){
  pval_notime_chr1 <- c(pval_notime_chr1,min(pval_gen_notime_lm[FourP$snp_info$ChrID=="I",i]))
  gene_notime_chr1 <- c(gene_notime_chr1,i)
  marker_notime_chr1 <- c(marker_notime_chr1,which.min(pval_gen_notime_lm[FourP$snp_info$ChrID=="I",i]))
}

for(i in 1:ncol(pval_gen_notime_lm)){
  pval_notime_chr2 <- c(pval_notime_chr2,min(pval_gen_notime_lm[FourP$snp_info$ChrID=="II",i])) 
  gene_notime_chr2 <- c(gene_notime_chr2,i)
  marker_notime_chr2 <- c(marker_notime_chr2,which.min(pval_gen_notime_lm[FourP$snp_info$ChrID=="II",i]))
}

for(i in 1:ncol(pval_gen_notime_lm)){
  pval_notime_chr3 <- c(pval_notime_chr3,min(pval_gen_notime_lm[FourP$snp_info$ChrID=="III",i]))
  gene_notime_chr3 <- c(gene_notime_chr3,i)
  marker_notime_chr3 <- c(marker_notime_chr3,which.min(pval_gen_notime_lm[FourP$snp_info$ChrID=="III",i]))
}

for(i in 1:ncol(pval_gen_notime_lm)){
  pval_notime_chr4 <- c(pval_notime_chr4,min(pval_gen_notime_lm[FourP$snp_info$ChrID=="IV",i]))
  gene_notime_chr4 <- c(gene_notime_chr4,i)
  marker_notime_chr4 <- c(marker_notime_chr4,which.min(pval_gen_notime_lm[FourP$snp_info$ChrID=="IV",i]))
}

for(i in 1:ncol(pval_gen_notime_lm)){
  pval_notime_chr5 <- c(pval_notime_chr5,min(pval_gen_notime_lm[FourP$snp_info$ChrID=="V",i]))
  gene_notime_chr5 <- c(gene_notime_chr5,i)
  marker_notime_chr5 <- c(marker_notime_chr5,which.min(pval_gen_notime_lm[FourP$snp_info$ChrID=="V",i]))
}

for(i in 1:ncol(pval_gen_notime_lm)){
  pval_notime_chrx <- c(pval_notime_chrx,min(pval_gen_notime_lm[FourP$snp_info$ChrID=="X",i]))
  gene_notime_chrx <- c(gene_notime_chrx,i)
  marker_notime_chrx <- c(marker_notime_chrx,which.min(pval_gen_notime_lm[FourP$snp_info$ChrID=="X",i]))
}

pval_notime_chr1 <- data.frame(pval = pval_notime_chr1, gene = gene_notime_chr1, marker = marker_notime_chr1,chr = rep("I",ncol(pval_gen_notime_lm)))
pval_notime_chr2 <- data.frame(pval = pval_notime_chr2, gene = gene_notime_chr2, marker = marker_notime_chr2,chr = rep("II",ncol(pval_gen_notime_lm)))
pval_notime_chr3 <- data.frame(pval = pval_notime_chr3, gene = gene_notime_chr3, marker = marker_notime_chr3,chr = rep("III",ncol(pval_gen_notime_lm)))
pval_notime_chr4 <- data.frame(pval = pval_notime_chr4, gene = gene_notime_chr4, marker = marker_notime_chr4,chr = rep("IV",ncol(pval_gen_notime_lm)))
pval_notime_chr5 <- data.frame(pval = pval_notime_chr5, gene = gene_notime_chr5, marker = marker_notime_chr5,chr = rep("V",ncol(pval_gen_notime_lm)))
pval_notime_chrx <- data.frame(pval = pval_notime_chrx, gene = gene_notime_chrx, marker = marker_notime_chrx,chr = rep("X",ncol(pval_gen_notime_lm)))

min_p_max_per_chr <- rbind(pval_notime_chr1,pval_notime_chr2,pval_notime_chr3,pval_notime_chr4,pval_notime_chr5,pval_notime_chrx)

# Collect lowest p-value for additive development model into frame.

pval_time_chr1 <- c()
pval_time_chr2 <- c()
pval_time_chr3 <- c()
pval_time_chr4 <- c()
pval_time_chr5 <- c()
pval_time_chrx <- c()

gene_time_chr1 <- c()
gene_time_chr2 <- c()
gene_time_chr3 <- c()
gene_time_chr4 <- c()
gene_time_chr5 <- c()
gene_time_chrx <- c()

marker_time_chr1 <- c()
marker_time_chr2 <- c()
marker_time_chr3 <- c()
marker_time_chr4 <- c()
marker_time_chr5 <- c()
marker_time_chrx <- c()

for(i in 1:ncol(pval_gen_time_lm)){
  pval_time_chr1 <- c(pval_time_chr1,min(pval_gen_time_lm[FourP$snp_info$ChrID=="I",i]))
  gene_time_chr1 <- c(gene_time_chr1,i)
  marker_time_chr1 <- c(marker_time_chr1,which.min(pval_gen_time_lm[FourP$snp_info$ChrID=="I",i]))
}

for(i in 1:ncol(pval_gen_time_lm)){
  pval_time_chr2 <- c(pval_time_chr2,min(pval_gen_time_lm[FourP$snp_info$ChrID=="II",i]))
  gene_time_chr2 <- c(gene_time_chr2,i)
  marker_time_chr2 <- c(marker_time_chr2,which.min(pval_gen_time_lm[FourP$snp_info$ChrID=="II",i]))
}

for(i in 1:ncol(pval_gen_time_lm)){
  pval_time_chr3 <- c(pval_time_chr3,min(pval_gen_time_lm[FourP$snp_info$ChrID=="III",i]))
  gene_time_chr3 <- c(gene_time_chr3,i)
  marker_time_chr3 <- c(marker_time_chr3,which.min(pval_gen_time_lm[FourP$snp_info$ChrID=="III",i]))
}

for(i in 1:ncol(pval_gen_time_lm)){
  pval_time_chr4 <- c(pval_time_chr4,min(pval_gen_time_lm[FourP$snp_info$ChrID=="IV",i]))
  gene_time_chr4 <- c(gene_time_chr4,i)
  marker_time_chr4 <- c(marker_time_chr4,which.min(pval_gen_time_lm[FourP$snp_info$ChrID=="IV",i]))
}

for(i in 1:ncol(pval_gen_time_lm)){
  pval_time_chr5 <- c(pval_time_chr5,min(pval_gen_time_lm[FourP$snp_info$ChrID=="V",i]))
  gene_time_chr5 <- c(gene_time_chr5,i)
  marker_time_chr5 <- c(marker_time_chr5,which.min(pval_gen_time_lm[FourP$snp_info$ChrID=="V",i]))
}

for(i in 1:ncol(pval_gen_time_lm)){
  pval_time_chrx <- c(pval_time_chrx,min(pval_gen_time_lm[FourP$snp_info$ChrID=="X",i]))
  gene_time_chrx <- c(gene_time_chrx,i)
  marker_time_chrx <- c(marker_time_chrx,which.min(pval_gen_time_lm[FourP$snp_info$ChrID=="X",i]))
}

pval_time_chr1 <- data.frame(pval_time = pval_time_chr1, gene_time = gene_time_chr1, marker_time = marker_time_chr1,chr_time = rep("I",ncol(pval_gen_time_lm)))
pval_time_chr2 <- data.frame(pval_time = pval_time_chr2, gene_time = gene_time_chr2, marker_time = marker_time_chr2,chr_time = rep("II",ncol(pval_gen_time_lm)))
pval_time_chr3 <- data.frame(pval_time = pval_time_chr3, gene_time = gene_time_chr3, marker_time = marker_time_chr3,chr_time = rep("III",ncol(pval_gen_time_lm)))
pval_time_chr4 <- data.frame(pval_time = pval_time_chr4, gene_time = gene_time_chr4, marker_time = marker_time_chr4,chr_time = rep("IV",ncol(pval_gen_time_lm)))
pval_time_chr5 <- data.frame(pval_time = pval_time_chr5, gene_time = gene_time_chr5, marker_time = marker_time_chr5,chr_time = rep("V",ncol(pval_gen_time_lm)))
pval_time_chrx <- data.frame(pval_time = pval_time_chrx, gene_time = gene_time_chrx, marker_time = marker_time_chrx,chr_time = rep("X",ncol(pval_gen_time_lm)))

min_p_max_per_chr_time <- rbind(pval_time_chr1,pval_time_chr2,pval_time_chr3,pval_time_chr4,pval_time_chr5,pval_time_chrx)

# Define frame with pvals for both models.

min_p_max_per_chr$pval <- -log10(min_p_max_per_chr$pval)
min_p_max_per_chr_time$pval_time <- -log10(min_p_max_per_chr_time$pval_time)

eqtlframe <- cbind(min_p_max_per_chr,min_p_max_per_chr_time)

### Make frame with only strong developmental effect eQTLs

very_sig_dev_frame <- eqtlframe[eqtlframe$pval_time-eqtlframe$pval > 10,]

# Partial eta squared without correcting for dev age.

smh <- function(mpgexp,mrk){
  model <- aov(lm(mpgexp~mrk))
  msqmrk <- sum(model$effects[2]^2)
  msqres <- sum(model$effects[c(-1,-2)]^2) # /197
  return(msqmrk/(msqmrk+msqres))
}

# Partial eta squared additive age model

amh <- function(mpgexp,mrk){
  model <- aov(lm(mpgexp~pco1.lines+mrk))
  msqmrk <- sum(model$effects[3]^2)
  msqres <- sum(model$effects[c(-1,-2,-3)]^2) # /197
  return(msqmrk/(msqmrk+msqres))
}


# Calculate partial eta squared for strong developmental effect eQTLs.
nsh_smh <- sapply(1:nrow(very_sig_dev_frame),function(x) smh(mpgexp = fpkm$fpkm[fpkm.selc,use.lines][very_sig_dev_frame[x,"gene"],], mrk = t(unname(FourP$snp_mat[FourP$snp_info$ChrID==very_sig_dev_frame[x,"chr"],use.lines][very_sig_dev_frame[x,"marker"],]))))
amh_smh <- sapply(1:nrow(very_sig_dev_frame),function(x) amh(mpgexp = fpkm$fpkm[fpkm.selc,use.lines][very_sig_dev_frame[x,"gene"],], mrk = t(unname(FourP$snp_mat[FourP$snp_info$ChrID==very_sig_dev_frame[x,"chr_time"],use.lines][very_sig_dev_frame[x,"marker_time"],]))))

nsh_smh <- data.frame(nsh = nsh_smh, Model = rep("Single marker model",length(nsh_smh)))
amh_smh <- data.frame(nsh = amh_smh, Model = rep("Additive model",length(amh_smh)))

nsh_frame <- rbind(nsh_smh,amh_smh)

### Make frame from which we can select the SMM-only and AAM-only eQTLs

cutoff_notime <- 4.06
cutoff_gen_time <- 4.44

scatter_frame <- eqtlframe %>% mutate(quadrant = case_when(pval >= cutoff_notime & pval_time >= cutoff_gen_time   ~ "Both significant",
                                                               pval < cutoff_notime & pval_time >= cutoff_gen_time  ~ "Only w/ time significant",
                                                               pval < cutoff_notime & pval_time < cutoff_gen_time ~ "Neither significant",
                                                               TRUE ~ "Only w/o time significant"))
# Partial eta squared for AAM-only eQTLs

topleft <- scatter_frame[scatter_frame$quadrant == "Only w/ time significant",]

nsh_topleft <- sapply(1:nrow(topleft),function(x) smh(mpgexp = fpkm$fpkm[fpkm.selc,use.lines][topleft[x,"gene"],], mrk = t(unname(FourP$snp_mat[FourP$snp_info$ChrID==topleft[x,"chr"],use.lines][topleft[x,"marker"],]))))
amh_topleft <- sapply(1:nrow(topleft),function(x) amh(mpgexp = fpkm$fpkm[fpkm.selc,use.lines][topleft[x,"gene"],], mrk = t(unname(FourP$snp_mat[FourP$snp_info$ChrID==topleft[x,"chr_time"],use.lines][topleft[x,"marker_time"],]))))

nsh_topleft <- data.frame(nsh = nsh_topleft, Model = rep("Single marker model",length(nsh_topleft)))
amh_topleft <- data.frame(nsh = amh_topleft, Model = rep("Additive model",length(amh_topleft)))

topleft_frame <- rbind(nsh_topleft,amh_topleft)

# Partial eta squared for SMM-only eQTLs.

botright <- scatter_frame[scatter_frame$quadrant == "Only w/o time significant",]

nsh_botright <- sapply(1:nrow(botright),function(x) smh(mpgexp = fpkm$fpkm[fpkm.selc,use.lines][botright[x,"gene"],], mrk = t(unname(FourP$snp_mat[FourP$snp_info$ChrID==botright[x,"chr"],use.lines][botright[x,"marker"],]))))
amh_botright <- sapply(1:nrow(botright),function(x) amh(mpgexp = fpkm$fpkm[fpkm.selc,use.lines][botright[x,"gene"],], mrk = t(unname(FourP$snp_mat[FourP$snp_info$ChrID==botright[x,"chr_time"],use.lines][botright[x,"marker_time"],]))))

nsh_botright <- data.frame(nsh = nsh_botright, Model = rep("Single marker model",length(nsh_botright)))
amh_botright <- data.frame(nsh = amh_botright, Model = rep("Additive model",length(amh_botright)))

botright_frame <- rbind(nsh_botright,amh_botright)

### Put all into a single frame
paper_frame <- rbind(data.frame(nsh_frame,quadrant = rep("Ten orders of magnitude increase",nrow(nsh_frame))),data.frame(topleft_frame,quadrant = rep("Top left quadrant",nrow(topleft_frame))),data.frame(botright_frame,quadrant = rep("Bottom right quadrant",nrow(botright_frame))))
paper_frame <- paper_frame[paper_frame$Model != "Dev removed",]
paper_frame$quadrant <- factor(paper_frame$quadrant, levels = c("Ten orders of magnitude increase", "Bottom right quadrant","Top left quadrant"))
paper_frame$Model <- factor(paper_frame$Model,levels = c("Single marker model","Additive model"))

load("partial_eta.out")

paper_frame$Model <- fct_recode(paper_frame$Model,"SMM" = "Single marker model", "AAM" = "Additive model")
levels(paper_frame$Model)<- c("SMM","AAM")

### Make legend of figure
legend <- get_legend(paper_frame %>% ggplot(aes(x = nsh,fill = Model)) + 
                       geom_density(alpha = 0.5, position = "identity") + 
                       theme_minimal() + 
                       xlab("Narrow sense heritability") +
                       ylab("Density")+
                       facet_wrap(~quadrant,nrow = 2,ncol =2,scales = "free") +
                       scale_fill_viridis(discrete = T,option = "C")+
                       theme(text = element_text(size = 80),legend.direction = "horizontal",legend.position = "top"))

n.bot <- nrow(paper_frame[paper_frame$quadrant=="Bottom right quadrant",])/2
n.top <- nrow(paper_frame[paper_frame$quadrant=="Top left quadrant",])/2
n.sde <- nrow(paper_frame[paper_frame$quadrant=="Ten orders of magnitude increase",])/2

paper_fig_bot <- paper_frame[paper_frame$quadrant=="Bottom right quadrant",] %>% ggplot(aes(x = nsh,fill = Model)) + 
  geom_density(alpha = 0.5, position = "identity") + 
  theme_classic() + 
  xlab("Marker effect") +
  ylab("Density")+
  scale_fill_viridis(discrete = T,option = "C")+
  ggtitle(paste("SMM-only eQTLs, n =",n.bot))+
  theme(text = element_text(size = 60),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5,size = 55),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

paper_fig_top <- paper_frame[paper_frame$quadrant=="Top left quadrant",] %>% ggplot(aes(x = nsh,fill = Model)) + 
  geom_density(alpha = 0.5, position = "identity") + 
  theme_classic() + 
  scale_fill_viridis(discrete = T,option = "C")+
  ggtitle(paste("AAM-only eQTLs, n = ",n.top))+
  xlab("Marker effect") +
  ylim(c(0,31)) +
  theme(text = element_text(size = 60),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5,size = 55),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 


png(filename = "plots_new_perm/heritability_density.png",height = 1000,width = 2100)
plot_grid(legend,plot_grid(paper_fig_bot+panel_border(size = 0.2,color = "black"),paper_fig_top+panel_border(size = 0.2,color = "black"),nrow = 1,align = "h",labels = c("A","B","C"),label_size = 50,axis = "tblr"), nrow = 2, rel_heights = c(0.15,1))
dev.off()

pdf("pdfs_figures/Figure_8.pdf",width = 30, height = 17)
plot_grid(legend,plot_grid(paper_fig_bot+panel_border(size = 0.2,color = "black"),paper_fig_top+panel_border(size = 0.2,color = "black"),nrow = 1,align = "h",labels = c("A","B","C"),label_size = 50,axis = "tblr"), nrow = 2, rel_heights = c(0.15,1))
dev.off()






