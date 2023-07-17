
library(ggplot2)
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

# Function which tests a linear model of pc1 ~ genotype
map_pc <- function(gen,pc1){
  ok <- summary(lm(pc1~gen))
  pval <- as.numeric(ok$coefficients[2,4])
  return(pval)
}

# For each gene, calculate p-value of genotype in single marker model and store these in pc1_gen.
pc1_gen <- apply(FourP$snp_mat[,use.lines],1,function (x) map_pc(as.numeric(x),pco1.lines))

# Prepare frame with p-values and marker position on the genome.
plot_frame <- data.frame(sig_gen = -log10(pc1_gen), pos = FourP$snp_info$SNP_position, chr = FourP$snp_info$ChrID)
plot_frame <- rbind(plot_frame,data.frame(sig_gen = 0,pos =0,chr = "II"))

# Produce plot.
plt <- ggplot(data = plot_frame,aes(x = pos, y = sig_gen,col = chr)) +
  geom_point()+
  theme_minimal()+
  theme(text = element_text(size = 35),
        legend.position="none",
        panel.spacing.x = unit(6, "mm"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.2))+
  facet_grid(~chr,scale = "free_x",space = "free_x")+
  xlab("Marker position in Mbp")+
  ylab("-log10 p-value") + 
  scale_color_viridis(discrete = T,direction = -1,option = "H") +
  geom_segment(x = 12, xend = 5, y = 6.218787, yend=  6.218787,lineend = "butt",
               linejoin = "mitre", size = 2, arrow = arrow(length = unit(5, "pt")),
               colour = "red", data = data.frame(chr = "X"))

png(filename = "plots_new_perm/pco1~markers.png",width = 1000,height = 700)
plt
dev.off()

pdf("pdfs_figures/Figure_S8.pdf",width = 15, height = 10)
plt
dev.off()









