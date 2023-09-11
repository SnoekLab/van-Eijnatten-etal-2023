.libPaths("~/R/x86_64-pc-linux-gnu-library/3.6/")
library(ggplot2)
library(RcppArmadillo)

# Import snp and fpkm data. 
load(file="~/Documents/elegans/data/obj_FourP_2022.out") 
load(file="~/Documents/elegans/data/obj_fpkm.out")

# Filter for mean samples fpkm > 1/32 and more than 20 samples have > 0 fpkm.
fpkm.selc <- apply(log2(fpkm$fpkm),1,mean)> -5 & apply(fpkm$fpkm>0,1,sum)>20 

# ge.mat contains all filtered samples with genes as observations and samples as features.
ge.mat <- fpkm$fpkm[fpkm.selc,]

# Some genes are transcribed more than others, but we are interested in relative differences. Take log2 ratio of mean.
rge.mat <- log2((ge.mat+1)/apply(ge.mat+1,1,mean))

# Store names of all genes which are not filtered out.
gene_names <- as.character(fpkm$genes[fpkm.selc])

# Do pca on rge.mat, and get the loadings of the samples on the first two principles components.
pco <- prcomp((rge.mat))
pco1 <- -pco$rotation[,1]
pco2 <- pco$rotation[,2]

# Subset data  to exclude four parent strains.
use.lines <- c(1:198,207)
pco1.lines <- pco1[use.lines]

# Select all rows for which the maximal log2 ratio is -1 < or > 1 (corresponds to a 2 fold change relative to the mean)
use.rge <- rge.mat #rge.mat[apply(abs(rge.mat),1,max)>1,]

#Want to do test or full run?
ng <- nrow(use.rge)

# Import snp and fpkm data. 
load(file="~/Documents/elegans/data/obj_FourP_2022.out") 
load(file="~/Documents/elegans/data/obj_fpkm.out")
dev_genes <- read.table(file="~/Documents/elegans/data/Dev_clust1_genes.txt",header = T)
dev_genes <- as.vector(as.character(dev_genes$genes))

# Find developmental marker genes with significant expression in our dataset.
dev_genes <- dev_genes[dev_genes %in% gene_names]

# Define matrix with gene expression for these marker genes.
dev_frame <- t(use.rge[gene_names %in% dev_genes,use.lines])

# Mean center and scale in the sample direction.
dev_frame <- apply(dev_frame,2,scale)

# Define vector with the mean gene expression of the developmental marker genes for each sample.
means_gen <- rowMeans(dev_frame)

# Turn into dataframe.
dev_frame <- as.data.frame(dev_frame)
colnames(dev_frame) <- dev_genes[dev_genes%in%gene_names]

# Same but also gets p-values of interaction term.
fast.map.co.int <- function(expla,rge){
  ok <- fastLm(cbind(1,expla),as.numeric(rge))
  # ok <- fastLm(gxpr~xm)
  pval <- pt(abs(ok$coefficients)/ok$stderr,ok$df.residual,lower.tail = F) * 2
  pval_gen <- pval[2]
  pval_dev <- pval[3]
  pval_int <- pval[4]
  out <- as.numeric(c(pval_gen,pval_dev,pval_int,as.numeric(ok$coefficients[2]),as.numeric(ok$coefficients)[3]),as.numeric(ok$coefficients)[4])
  return(out)
}

# Same as above, but now store interaction p-value in 3rd column p_val_gene_i.
for(i in 1:ng){
  gexpr <- use.rge[i,use.lines]
  assign(paste0("p_val_gene_",i), apply(FourP$snp_mat[,use.lines],1,function (x) fast.map.co.int(data.frame(as.numeric(x),means_gen,as.numeric(x)*as.numeric(means_gen)),gexpr)))
}

pval_gen_time_intlm <- matrix(nrow =8933,ncol = ng)

# Put genotype p-value in single matrix.
for(j in 1:ng){
  name <-paste0("p_val_gene_",j)
  pval_gen_time_intlm[,j] <- get(name)[1,]
}

saveRDS(pval_gen_time_intlm, file = "~/Documents/elegans/map_with_mean_dev_expr/pval_gen_time_int_mean_dev.RDS")

pval_dev_time_intlm <- matrix(nrow =8933,ncol = ng)

# Put time p-value in single matrix
for(j in 1:ng){
  name <-paste0("p_val_gene_",j)
  pval_dev_time_intlm[,j] <- get(name)[2,]
}

saveRDS(pval_dev_time_intlm, file = "~/Documents/elegans/map_with_mean_dev_expr/pval_dev_time_int_mean_dev.RDS")

pval_int_time_intlm <- matrix(nrow =8933,ncol = ng)

#put Interaction p-value in single matrix.
for(j in 1:ng){
  name <-paste0("p_val_gene_",j)
  pval_int_time_intlm[,j] <- get(name)[3,]
}

saveRDS(pval_int_time_intlm, file = "~/Documents/elegans/map_with_mean_dev_expr/pval_int_time_int_mean_dev.RDS")





