library(ggplot2)
library(cowplot)
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

### Make dataframe with genotypes of the different worms at 3.405162 Mbp on chromsome X and the loadings on PC1. 
to.pl.3.404162 <- data.frame(t(unname(round(FourP$snp_mat[FourP$snp_info$SNP_position==3.404162,use.lines]*2,0)/2)),pco1.lines)
colnames(to.pl.3.404162) <- c("Genotype","Developmental_age")

### Function performs linear model of PC1~genotype 
map_pc <- function(gen,pc1){
  ok <- summary(lm(pc1~gen))
  pval <- as.numeric(ok$coefficients[2,4])
  return(pval)
}

### Get p-value
map_pc(as.numeric(FourP$snp_mat[which(FourP$snp_info$SNP_position==3.404162),use.lines]),pco1.lines)

### Make boxplot
fig_3.404162 <- ggplot(to.pl.3.404162[to.pl.3.404162$Genotype != 0.5,],aes(y = Developmental_age,fill = as.factor(Genotype)))+
                    geom_boxplot()+
                    scale_fill_manual(values = c("0" = "#377eb8", "1" = "#e41a1c"),guide = "none") +
                    theme_minimal() +
                    ylab("Developmental age") +
                    theme(axis.ticks.x = element_blank(),
                          axis.text.x = element_blank(),
                          panel.grid.major.x = element_blank(),
                          panel.grid.minor.x = element_blank(),
                          plot.title = element_text(hjust = 0.5),
                          text = element_text(size = 27)) +
                    ggtitle("3.404162 Mb, chr X\nPval = 6.04E-07")

### Do same for locus at 12.818318 Mbp on chromosome II.
to.pl.12.818318 <- data.frame(t(unname(round(FourP$snp_mat[FourP$snp_info$SNP_position==12.818318,use.lines]*2,0)/2)),pco1.lines)
colnames(to.pl.12.818318) <- c("Genotype","Developmental_age")

map_pc(as.numeric(FourP$snp_mat[which(FourP$snp_info$SNP_position==12.818318),use.lines]),pco1.lines)

fig_12.818318 <- ggplot(to.pl.12.818318[to.pl.12.818318$Genotype != 0.5,],aes(y = Developmental_age,fill = as.factor(Genotype)))+
  geom_boxplot()+
  scale_fill_manual(values = c("0" = "#377eb8", "1" = "#e41a1c"),guide = "none")+
  theme_minimal() +
  ylab("") +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 27)) +
  ggtitle("12.818318 Mb, chr II\nPval = 5.87E-05")

### And finally for locus at 4.042568Mbp on chromosome III.
to.pl.4.042568 <- data.frame(t(unname(round(FourP$snp_mat[FourP$snp_info$SNP_position==4.042568,use.lines]*2,0)/2)),pco1.lines)
colnames(to.pl.4.042568) <- c("Genotype","Developmental_age")

map_pc(as.numeric(FourP$snp_mat[which(FourP$snp_info$SNP_position==4.042568),use.lines]),pco1.lines)

fig_4.042568 <- ggplot(to.pl.4.042568[to.pl.4.042568$Genotype != 0.5,],aes(y = Developmental_age,fill = as.factor(Genotype)))+
  geom_boxplot()+
  scale_fill_manual(values = c("0" = "#377eb8", "1" = "#e41a1c"),guide = "none")+
  theme_minimal() +
  ylab("") +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 27)) +
  ggtitle("4.042568 Mb, chr III\nPval = 4.58E-05")

legend <- get_legend(ggplot(to.pl.4.042568[to.pl.4.042568$Genotype != 0.5,],aes(y = Developmental_age,fill = as.factor(Genotype)))+
             geom_boxplot()+
             scale_fill_manual(values = c("0" = "#377eb8", "1" = "#e41a1c"))+
             theme_minimal() +
             ylab("Developmental age") +
             theme(axis.ticks.x = element_blank(),
                   axis.text.x = element_blank(),
                   panel.grid.major.x = element_blank(),
                   panel.grid.minor.x = element_blank(),
                   plot.title = element_text(hjust = 0.5),
                   text = element_text(size = 30)) +
             guides(fill=guide_legend(title="Genotype")) +
             theme(legend.position = "top",legend.direction = "horizontal",legend.justification="center") +
             ggtitle("4.042568 Mb, chr III"))

png("plots_new_perm/dev_dis_hotspots.png",height = 700,width = 1150)
plot_grid(legend,plot_grid(fig_3.404162+panel_border(color = "black",size = 0.2),fig_12.818318+panel_border(color = "black",size = 0.2),fig_4.042568+panel_border(color = "black",size = 0.2),ncol =3, align = "h",labels = c("A","B","C"),label_size = 20),rel_heights = c(1,5),nrow = 2)
dev.off()

pdf("pdfs_figures/Figure_3.pdf",width = 17, height = 10)
plot_grid(legend,plot_grid(fig_3.404162+panel_border(color = "black",size = 0.2),fig_12.818318+panel_border(color = "black",size = 0.2),fig_4.042568+panel_border(color = "black",size = 0.2),ncol =3, align = "h",labels = c("A","B","C"),label_size = 20),rel_heights = c(1,5),nrow = 2)
dev.off()




