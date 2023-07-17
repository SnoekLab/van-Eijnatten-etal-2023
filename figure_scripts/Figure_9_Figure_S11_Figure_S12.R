
library(heritability)
library(ggplot2)
library(cowplot)
library(viridis)
library(rrBLUP)
library(dplyr)
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

### Calculate kinship matrix using A.mat function from rrBLUP package
kin_mat <- A.mat(t(as.matrix(FourP$snp_mat[,use.lines])))

### Calculate h2 gene expression without covariates
h2 <- apply(use.rge[,use.lines],1,function(x) marker_h2(data.vector = x,
                                                            geno.vector = colnames(use.rge[,use.lines]),
                                                            K=kin_mat))

#saveRDS(h2,file = "h2_nodev.RDS")
#h2 <- readRDS("h2_nodev.RDS")

### Put h2 values into vector
h2.vec <- numeric(length=length(h2))

for(i in 1:length(h2)){
  h2.vec[i] <- h2[[i]]$h2
}

### Calculate h2 with PC1 as covariate.
capture.output(h2_pc1 <- apply(use.rge[,use.lines],1,function(x) marker_h2(data.vector = x,
                                                  geno.vector = colnames(use.rge[,use.lines]),
                                                  covariates = pco1.lines,
                                                  K=kin_mat)),file = nullfile())

#saveRDS(h2_pc1,file = "h2_dev.RDS")
#h2_pc1 <- readRDS("h2_dev.RDS")

### Put in vector
h2.vec.pc1 <- numeric(length=length(h2_pc1))

for(i in 1:length(h2_pc1)){
  h2.vec.pc1[i] <- h2_pc1[[i]]$h2
}


### Get eQTLs SMM

get_eqtl <- function(gene,chr,pvals,cutoff,model,sig){
  index <- min(which(gene_name == ""))
  gene_name[index] <<- gene_names[gene]
  pos[index] <<- FourP$snp_info$SNP_position[FourP$snp_info$ChrID==chr][which.max(pvals)]
  loc <- gsub(".*:","",fpkm$locus[fpkm.selc][gene])
  middle <- (as.numeric(sub("(^[^-]+)-.*", "\\1", loc)) + as.numeric(sub('.+-(.+)', '\\1', loc))) / 2 / 1000000
  gene_pos[index] <<- middle
  gene_chr <- strsplit(as.character(fpkm$locus[fpkm.selc][gene]), "[:]")[[1]][1]
  gene.chr[index] <<- gene_chr
  ct[index] <<- ifelse(abs(FourP$snp_info$SNP_position[FourP$snp_info$ChrID==chr][which.max(pvals)] - middle) < 2 & gene_chr ==chr,ct <- c("cis"),ct <- c("trans"))
  marker_chr[index] <<- chr
  gene_index[index] <<- gene
}

map_function <- function(gene,cutoff,pvals,model,sig){
  pvals_chr1 <- pvals[FourP$snp_info$ChrID=="I"]
  pvals_chr2 <- pvals[FourP$snp_info$ChrID=="II"]
  pvals_chr3 <- pvals[FourP$snp_info$ChrID=="III"]
  pvals_chr4 <- pvals[FourP$snp_info$ChrID=="IV"]
  pvals_chr5 <- pvals[FourP$snp_info$ChrID=="V"]
  pvals_chrx <- pvals[FourP$snp_info$ChrID=="X"]
  if(max(pvals_chr1)>cutoff){
    get_eqtl(gene,"I",pvals_chr1,cutoff,model,sig)
  }
  if(max(pvals_chr2)>cutoff){
    get_eqtl(gene,"II",pvals_chr2,cutoff,model,sig)
  }
  if(max(pvals_chr3)>cutoff){
    get_eqtl(gene,"III",pvals_chr3,cutoff,model,sig)
  }
  if(max(pvals_chr4)>cutoff){
    get_eqtl(gene,"IV",pvals_chr4,cutoff,model,sig)
  }
  if(max(pvals_chr5)>cutoff){
    get_eqtl(gene,"V",pvals_chr5,cutoff,model,sig)
  }
  if(max(pvals_chrx)>cutoff){
    get_eqtl(gene,"X",pvals_chrx,cutoff,model,sig)
  }
}

cutoff_notime <- 4.06

pval_gen_notime_lm <- readRDS("~/Documents/elegans/pvals_correct/pval_gen_notime_lm.RDS")

gene_name <- character(length = 18605)
pos <- numeric(length = 18605)
gene_pos <- numeric(length = 18605)
gene.chr <- character(length = 18605)
ct <- character(length = 18605)
marker_chr <- character(length = 18605)
gene_index <- numeric(length = 18605)

sapply(1:ncol(pval_gen_notime_lm),function(x) map_function(x,cutoff_notime,-log10(pval_gen_notime_lm[,x]),"Single marker model",0.05))

eqtl.smm <- data.frame(gene_name = gene_name, marker_pos = pos, gene_pos = gene_pos,ct = ct, marker_chr = marker_chr, gene_chr = gene.chr,gene_index = gene_index)

### Get eQTLs AAM

cutoff_dev <- 4.44

pval_gen_time_lm <- readRDS("~/Documents/elegans/pvals_correct/pval_gen_time_lm.RDS")

gene_name <- character(length = 8666)
pos <- numeric(length = 8666)
gene_pos <- numeric(length = 8666)
gene.chr <- character(length = 8666)
ct <- character(length = 8666)
marker_chr <- character(length = 8666)
gene_index <- numeric(length = 8666)

sapply(1:ncol(pval_gen_time_lm),function(x) map_function(x,cutoff_dev,-log10(pval_gen_time_lm[,x]),"Additive model",0.05))

eqtl.am <- data.frame(gene_name = gene_name, marker_pos = pos, gene_pos = gene_pos,ct = ct, marker_chr = marker_chr, gene_chr = gene.chr,gene_index = gene_index)

### Get SMM-only eQTLs
sm.eqtls <- anti_join(eqtl.smm[,c("gene_index","marker_chr")],eqtl.am[,c("gene_index","marker_chr")])

### Store in vector for each gene whether it has a SMM-only eQTL

genes.sm.eqtls <- logical(length=12029)

for(i in 1:12029){
  if(i %in% sm.eqtls$gene_index){
    genes.sm.eqtls[i] <- T
  }
  else{
    genes.sm.eqtls[i] <-F
  }
}

### Store in vector for each gene whether it has a AAM-only eQTL

dyn.eqtls <- anti_join(eqtl.am[,c("gene_index","marker_chr")],eqtl.smm[,c("gene_index","marker_chr")])

genes.dyn.eqtls <- logical(length=12029)

for(i in 1:12029){
  if(i %in% dyn.eqtls$gene_index){
    genes.dyn.eqtls[i] <- T
  }
  else{
    genes.dyn.eqtls[i] <-F
  }
}

### Make vector describing what types of eQTL each gene has

eqtls.vec <- character(length = 12029)

for(i in 1:12029){
  if(genes.sm.eqtls[i] == T & genes.dyn.eqtls[i] == F){
    eqtls.vec[i] <- "SMM-only eQTL but no AAM-only eQTL"
  }
  if(genes.sm.eqtls[i] == F & genes.dyn.eqtls[i] == T){
    eqtls.vec[i] <- "AAM-only eQTL but no SMM-only eQTL"
  }
  if(genes.sm.eqtls[i] == T & genes.dyn.eqtls[i] == T){
    eqtls.vec[i] <- "SMM-only eQTL and AAM-only eQTL"
  }
  if(genes.sm.eqtls[i] == F & genes.dyn.eqtls[i] == F & !(i%in%eqtl.smm$gene_index)){
    eqtls.vec[i] <- "No eQTL"
  }
  if(genes.sm.eqtls[i] == F & genes.dyn.eqtls[i] == F & i%in%eqtl.smm$gene_index){
    eqtls.vec[i] <- "Only same eQTLs between models"
  }
}

### Make frame to plot figure S11

h2.frame.2 <- data.frame(h2_nodev = h2.vec, h2_dev = h2.vec.pc1, eqtl = eqtls.vec)
h2.frame.2$eqtl <- factor(h2.frame.2$eqtl,levels = c("No eQTL","Only same eQTLs between models","SMM-only eQTL but no AAM-only eQTL","AAM-only eQTL but no SMM-only eQTL","SMM-only eQTL and AAM-only eQTL"))

h2.plot.2 <- ggplot(h2.frame.2,aes(x = h2_nodev,y = h2_dev,col = eqtl)) +
  geom_point(alpha = 0.7,size = 3)+
  geom_abline(linewidth = 2)+
  theme_minimal()+
  xlab(expression(paste(h["2"]," without cofactor")))+
  ylab(expression(paste(h["2"]," with PC1 as cofactor")))+
  theme(text = element_text(size = 25),legend.position = "bottom",legend.direction = "vertical",
        axis.ticks.x = element_line(size = 0.2),
        axis.ticks.y = element_line(size=0.2),
        panel.grid.major = element_line(size = 0.2,color = "grey"),
        panel.grid.minor = element_line(size = 0.2, color = "grey"))+
  scale_color_viridis(discrete = T,option = "H",direction = -1)+
  guides(col=guide_legend("Type eQTL"))

png(filename = "plots_new_perm/h2_nodev_vs_dev_col_eqtl.png",width =600,height = 700)
h2.plot.2
dev.off()

pdf("pdfs_figures/Figure_S11.pdf",width = 8, height = 9)
h2.plot.2
dev.off()

### Make frame to plot figure S12A

# Collect number of eQTLs per gene according to AAM

am.eqtls <- numeric(length=12029)

for(i in 1:12029){
  am.eqtls[i] <- sum(eqtl.am$gene_index==i)
}

h2.frame.am <- data.frame(h2_nodev = h2.vec, h2_dev = h2.vec.pc1, eqtl = am.eqtls)

am.eqtl <- ggplot(h2.frame.am,aes(x = h2_nodev,y = h2_dev,col = eqtl)) +
  geom_point(alpha = 0.7,size =3)+
  geom_abline(linewidth = 2)+
  theme_minimal()+
  xlab(expression(paste(h["2"]," without cofactor")))+
  ylab(expression(paste(h["2"]," with PC1 as cofactor")))+
  ggtitle("AM")+
  theme(text = element_text(size = 35),legend.position = "none",
        axis.ticks.x = element_line(size = 0.2),
        axis.ticks.y = element_line(size=0.2),
        panel.grid.major = element_line(size = 0.2,color = "grey"),
        panel.grid.minor = element_line(size = 0.2, color = "grey"),
        plot.title = element_text(hjust = 0.5))+
  scale_color_viridis(option = "H")+
  guides(col=guide_legend("Number eQTL (AM)"))

# Do same for SMM, figure S12B

smm.eqtls <- numeric(length=12029)

for(i in 1:12029){
  smm.eqtls[i] <- sum(eqtl.smm$gene_index==i)
}

h2.frame.smm <- data.frame(h2_nodev = h2.vec, h2_dev = h2.vec.pc1, eqtl = smm.eqtls)

smm.eqtl <- ggplot(h2.frame.smm,aes(x = h2_nodev,y = h2_dev,col = eqtl)) +
  geom_point(alpha = 0.7,size =3)+
  geom_abline(linewidth = 2)+
  theme_minimal()+
  xlab(expression(paste(h["2"]," without cofactor")))+
  ylab(expression(paste(h["2"]," with PC1 as cofactor")))+
  ggtitle("SMM")+
  theme(text = element_text(size = 35),legend.position = "none",
        axis.ticks.x = element_line(size = 0.2),
        axis.ticks.y = element_line(size=0.2),
        panel.grid.major = element_line(size = 0.2,color = "grey"),
        panel.grid.minor = element_line(size = 0.2, color = "grey"),
        plot.title = element_text(hjust = 0.5))+
  scale_color_viridis(option = "H")+
  guides(col=guide_legend("Number eQTL (SMM)"))

### Make legend

legend<- get_legend(ggplot(h2.frame.smm,aes(x = h2_nodev,y = h2_dev,col = eqtl)) +
                      geom_point(alpha = 0.7,size =3)+
                      geom_abline(linewidth = 2)+
                      theme_minimal()+
                      xlab("H2 without cofactor")+
                      ylab("H2 with PC1 as cofactor")+
                      theme(text = element_text(size = 35),legend.position = "top",legend.direction = "horizontal")+
                      scale_color_viridis(option = "H")+
                      guides(col=guide_legend("Number eQTL")))

h2.plot.smm.am <- plot_grid(legend,plot_grid(smm.eqtl,am.eqtl,ncol = 2,rel_widths = c(1,1),axis = "tblr",align = "h",labels = c("A","B"),label_size =30),nrow = 2, rel_heights = c(1,10))

png(filename = "plots_new_perm/h2_smm_am.png",width = 1600, height = 850)
h2.plot.smm.am
dev.off()

pdf("pdfs_figures/Figure_S12.pdf",width = 20,height = 11)
h2.plot.smm.am
dev.off()

### Now figure 9 from the main text.

eqtl.nr.frame <- data.frame(h2_nopc1=h2.vec,h2_pc1=h2.vec.pc1,smm_eqtls = factor(smm.eqtls), am_eqtls = factor(am.eqtls))

legend<- get_legend(ggplot(eqtl.nr.frame,aes(y = h2_nopc1,fill = smm_eqtls))+
         geom_boxplot()+
         facet_grid(~smm_eqtls, switch = "x")+
         scale_fill_viridis(discrete = T,option = "H",direction = -1) +
         theme_minimal()+
         ylab("NSH without PC1 as cofactor")+
         theme(text = element_text(size =25),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               strip.background = element_blank(),
               strip.text.x = element_blank(),
               axis.text.x=element_blank(),legend.position = "top",legend.direction = "horizontal")+
         guides(fill=guide_legend("",nrow = 1)))

### Label number of genes in red (separate label frame for each plot)
label.frame.smm <- data.frame(
  label = as.numeric(table(smm.eqtls)),
  smm_eqtls = factor(c(0,1,2,3,4,5,6)),
  y     = c(-0.05,-0.05,-0.05,-0.05,-0.05,-0.05,-0.05),
  x = c(0,0,0,0,0,0,0)
)


### Make boxplot of nsh without covariate, facetted by number of eQTLs per gene according to SMM. 
smm_nsh_boxplot <- ggplot(eqtl.nr.frame,aes(y = h2_nopc1,fill = smm_eqtls))+
  geom_boxplot()+
  geom_text(data = label.frame.smm, mapping = aes(x = x,y = y, label = label),fontface = "bold",col = "red") +
  facet_grid(~smm_eqtls, switch = "x")+
  scale_fill_viridis(discrete = T,option = "H",direction = -1) +
  theme_minimal()+
  ylab(expression(paste(h["2"]," without PC1")))+
  xlab("Number of eQTLs in SMM")+
  scale_y_continuous(expand = c(0.05,0.05),limits = c(-0.1,1),breaks = c(0,0.25,0.5,0.75,1))+
  theme(text = element_text(size =17),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.text.x=element_blank(),
        legend.position = "none")

label.frame.am <- data.frame(
  label = as.numeric(table(am.eqtls)),
  am_eqtls = factor(c(0,1,2,3,4,5,6)),
  y  = c(-0.05,-0.05,-0.05,-0.05,-0.05,-0.05,-0.05),
  x = c(0,0,0,0,0,0,0)
)

### Make boxplot of nsh without covariate, facetted by number of eQTLs per gene according to AAM 
am_nsh_boxplot <- ggplot(eqtl.nr.frame,aes(y = h2_nopc1,fill = am_eqtls))+
  geom_boxplot()+
  geom_text(data = label.frame.am, mapping = aes(x = x,y = y, label = label),fontface = "bold",col = "red") +
  facet_grid(~am_eqtls, switch = "x")+
  scale_fill_viridis(discrete = T,option = "H",direction = -1) +
  theme_minimal()+
  ylab(expression(paste(h["2"]," without PC1")))+
  xlab("Number of eQTLs in AAM")+
  scale_y_continuous(expand = c(0.05,0.05),limits = c(-0.1,1),breaks = c(0,0.25,0.5,0.75,1))+
  theme(text = element_text(size =17),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.text.x=element_blank(),legend.position = "none")


### Make boxplot of nsh with PC1 as covariate, facetted by number of eQTLs per gene according to SMM. 
smm_nsh_pc1_boxplot <- ggplot(eqtl.nr.frame,aes(y = h2_pc1,fill = smm_eqtls))+
  geom_boxplot()+
  geom_text(data = label.frame.smm, mapping = aes(x = x,y = y, label = label),fontface = "bold",col = "red") +
  facet_grid(~smm_eqtls, switch = "x")+
  scale_fill_viridis(discrete = T,option = "H",direction = -1) +
  theme_minimal()+
  ylab(expression(paste(h["2"]," with PC1")))+
  xlab("Number of eQTLs in SMM")+
  scale_y_continuous(expand = c(0.05,0.05),limits = c(-0.1,1),breaks = c(0,0.25,0.5,0.75,1))+
  theme(text = element_text(size =17),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.text.x=element_blank(),legend.position = "none")

### Make boxplot of nsh with PC1 as covariate, facetted by number of eQTLs per gene according to AAM. 
am_nsh_pc1_boxplot <- ggplot(eqtl.nr.frame,aes(y = h2_pc1,fill = am_eqtls))+
  geom_boxplot()+
  geom_text(data = label.frame.am, mapping = aes(x = x,y = y, label = label),fontface = "bold",col = "red") +
  facet_grid(~am_eqtls, switch = "x")+
  scale_fill_viridis(discrete = T,option = "H",direction = -1) +
  theme_minimal()+
  ylab(expression(paste(h["2"]," with PC1")))+
  xlab("Number of eQTLs in AAM")+
  scale_y_continuous(expand = c(0.05,0.05),limits = c(-0.1,1),breaks = c(0,0.25,0.5,0.75,1))+
  theme(text = element_text(size =17),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.text.x=element_blank(),legend.position = "none")

png(filename = "plots_new_perm/nsh_cowplot.png",width =700,height = 550)
plot_grid(legend,plot_grid(plot_grid(smm_nsh_boxplot)+panel_border(color = "black",size = 0.2),plot_grid(smm_nsh_pc1_boxplot)+panel_border(color = "black",size = 0.2),plot_grid(am_nsh_boxplot)+panel_border(color = "black",size = 0.2),plot_grid(am_nsh_pc1_boxplot)+panel_border(color = "black",size = 0.2),nrow =2,align = "v",axis = "tblr",labels = c("A","B","C","D"),label_size = 15),nrow = 2,rel_heights = c(1,10))+theme(plot.margin = unit(c(20,20,20,20), "points"))
dev.off()

pdf("pdfs_figures/Figure_9.pdf",width = 10, height = 8)
plot_grid(legend,plot_grid(plot_grid(smm_nsh_boxplot)+panel_border(color = "black",size = 0.2),plot_grid(smm_nsh_pc1_boxplot)+panel_border(color = "black",size = 0.2),plot_grid(am_nsh_boxplot)+panel_border(color = "black",size = 0.2),plot_grid(am_nsh_pc1_boxplot)+panel_border(color = "black",size = 0.2),nrow =2,align = "v",axis = "tblr",labels = c("A","B","C","D"),label_size = 15),nrow = 2,rel_heights = c(1,10))+theme(plot.margin = unit(c(20,20,20,20), "points"))
dev.off()

