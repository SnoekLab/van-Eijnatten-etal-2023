library(ggplot2)
library(cowplot)
library(viridis)

### Read in data from Boeck et al (public)
expr.mat <- read.table("dev_stage_investigation/Data_time_resolved_transcriptome/Supplemental_Table_S2/Suppl_Table_2_expression_of_wormbase_genes_readcounts_dcpm_per_sample.txt",header = T)

load(file="~/Documents/elegans/data/obj_fpkm.out")

use.lines <- c(1:198,207)

### We need to match genes between their dataset and ours, despite different identifiers. Below are some steps to get as many matching genes as possible. 

genes.us <- as.character(fpkm$genes[-which(fpkm$genes=="-")])

split.genes <- sapply(genes.us,function(x) strsplit(x,","))

genes.us.split <- c()

for(i in 1:length(genes.us)){
  genes.us.split <- c(genes.us.split,split.genes[i][[1]])
}

genes.us.split<- unique(genes.us.split)

fpkm.split.genes <- matrix(0,ncol = 199, nrow = length(genes.us.split))

for(i in 1:length(genes.us.split)){
  gene.index <- grep(paste0("(^|,)",genes.us.split[i],"(,|$)"),fpkm$genes)
  if(length(gene.index)>1){
    lengths <- sapply(gene.index,function(x) length(fpkm$genes[x]))
    shortest <- which.min(lengths)
    gene.index <- gene.index[shortest]
  }
  fpkm.split.genes[i,] <- fpkm$fpkm[gene.index,use.lines]
}

rownames(fpkm.split.genes) <- genes.us.split

#write.table(genes.us.split,file = "total_genes_nofilter.txt",row.names = F,sep = ",")

### Translation table between gene names and sequence names
us.seq.names <- read.table("dev_stage_investigation/sequence_names_total_genes_nofilter.txt",header = T,sep = "\t")


for(i in 1:nrow(us.seq.names)){
  if(us.seq.names[i,2]=="not found"){
    us.seq.names[i,2] <- us.seq.names[i,1]
  }
}

### Get duplicated entries
double.id <- us.seq.names[44903:nrow(us.seq.names),]

### Remove from vector
us.seq.names <- us.seq.names[1:44902,]

for(i in 1:length(genes.us.split)){
  gene.seq <- us.seq.names[us.seq.names$Your.Input==genes.us.split[i],2]
  rownames(fpkm.split.genes)[i] <- gene.seq
}

doubles <-matrix(ncol = 199, nrow = 0)

for(i in unique(double.id$Your.Input)){
  seq.gene <- double.id[double.id$Your.Input==i,2]
  seq.1 <- seq.gene[1]
  rownames(fpkm.split.genes)[which(rownames(fpkm.split.genes)==i)]<-seq.1
  if(length(seq.gene)>1){
    seq.2 <-seq.gene[2]
    doubles <- rbind(doubles,fpkm.split.genes[genes.us.split==i,])
    rownames(doubles)[nrow(doubles)]<-seq.2
  }
}

fpkm.split.genes <- rbind(fpkm.split.genes,doubles)

genes_them <- expr.mat$WormbaseName

fpkm.split.both <- fpkm.split.genes[rownames(fpkm.split.genes)%in%genes_them,]

fpkm.split.both <- fpkm.split.both[!duplicated(rownames(fpkm.split.both)),]

genes.them.both <- genes_them[genes_them %in%rownames(fpkm.split.both)]

## We transform their counts into CPM

L4_a <- expr.mat$L4_counts[genes_them %in%rownames(fpkm.split.both)]
L4_a <- L4_a/sum(L4_a)*1e6

L4_b <- expr.mat$L4b_counts[genes_them %in%rownames(fpkm.split.both)]
L4_b <- L4_b/sum(L4_b)*1e6

YA_a <- expr.mat$YA_counts[genes_them %in%rownames(fpkm.split.both)]
YA_a <- YA_a/sum(YA_a)*1e6

YA_b <- expr.mat$N2_Yad.1_counts[genes_them %in%rownames(fpkm.split.both)]
YA_b <- YA_b/sum(YA_b)*1e6

fpkm.split.both <- fpkm.split.both[match(genes.them.both,rownames(fpkm.split.both)),]

#save(fpkm.split.both,file = "fpkm.split.both.out")

### Calculate correlations between our transcriptomes
cor.l4a <- apply(fpkm.split.both,2,function(x) cor(x,L4_a))
cor.l4b <- apply(fpkm.split.both,2,function(x) cor(x,L4_b))
cor.YAa <- apply(fpkm.split.both,2,function(x) cor(x,YA_a))
cor.YAb <- apply(fpkm.split.both,2,function(x) cor(x,YA_b))

### Get mean for the replicates
cor.L4 <- (cor.l4a+cor.l4b)/2
cor.YA <- (cor.YAa+cor.YAb)/2

### Plot
cor.L4 <- data.frame(correlation = cor.L4,stage = "L4")
cor.YA <- data.frame(correlation = cor.YA,stage = "YA")

plot.frame <- rbind(cor.L4,cor.YA)

box.cor.stage <- ggplot(plot.frame,aes(x = stage,y = correlation,fill = stage))+
  geom_boxplot()+
  theme_minimal()+
  ylab("Correlation with mpRILs")+
  xlab("")+
  ylim(0,1)+
  scale_fill_manual(values = c("blue","red")) +
  theme(text=element_text(size = 35),panel.border = element_rect(linewidth = 0.2, fill = NA))+
  guides(fill = guide_legend("Stage"))

png(filename = "plots_new_perm/box_cor.png",width = 850, height = 550)
box.cor.stage
dev.off()

pdf("pdfs_figures/Figure_S16.pdf",width = 10, height = 7)
box.cor.stage
dev.off()



















