library(readxl)
library(ggplot2)

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

# Import matrix w/ phenotypes, publicly available Snoek et al 2019
phenotypes <- read_excel("~/Documents/elegans/data/snoeketal_2019.xlsx",sheet = 7,col_names = T)

phenotypes <- phenotypes[use.lines,]

correlations <- apply(phenotypes[6:26],2,function(x) cor(as.numeric(x),pco1.lines, use="complete.obs"))
pheno <- colnames(phenotypes)[6:26]

plot.frame <- data.frame(pheno = pheno,cor = correlations)

pheno.plot <- ggplot(plot.frame,aes(x=pheno,y= correlations))+
  geom_bar(stat = "identity")+
  theme_minimal()+
  theme(axis.text.x=element_text(angle=90,vjust = 0.5, hjust=1),text = element_text(size =35),
        panel.border = element_rect(linewidth = 0.2,fill = NA))+
  ylim(-1,1)+
  xlab("Phenotype")+
  ylab("Correlation with PC1")

png(filename = "plots_new_perm/pheno_cor.png",width = 1000,height = 700)
pheno.plot
dev.off()

pdf("pdfs_figures/Figure_S18.pdf", width = 12, height = 8.5)
pheno.plot
dev.off()

