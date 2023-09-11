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

# Do pca on rge.mat, and get the loadings of the samples on the first two principles components.
pco <- prcomp((rge.mat))
pco1 <- -pco$rotation[,1]
pco2 <- pco$rotation[,2]

# Select all rows for which the maximal log2 ratio is -1 < or > 1 (corresponds to a 2 fold change relative to the mean)
use.rge <- rge.mat #rge.mat[apply(abs(rge.mat),1,max)>1,]

#Want to do test or full run?
ng <- nrow(use.rge)

# Function which does a linear model of gene expression vs genotype + pco1 (no interaction term)
fast.map.co <- function(expla,rge){
  ok <- fastLm(cbind(1,expla),as.numeric(rge))
  pval <- pt(abs(ok$coefficients)/ok$stderr,ok$df.residual,lower.tail = F) * 2
  pval_gen <- pval[2]
  pval_dev <- pval[3]
  out <- as.numeric(c(pval_gen,pval_dev,as.numeric(ok$coefficients[2]),as.numeric(ok$coefficients[3]))) #,R2,R2.adj)
  return(out)
}

# Subset data  to exclude four parent strains.
use.lines <- c(1:198,207)
pco1.lines <- pco1[use.lines]

# Model with time. Time p-value now stored in the second column of p_val_gene_i.
for(i in 1:ng){
  gexpr <- sample(use.rge[i,use.lines])
  assign(paste0("p_val_gene_",i), apply(FourP$snp_mat[,use.lines],1,function (x) fast.map.co(data.frame(as.numeric(x),pco1.lines),gexpr)))
}

pval_gen_time_lm <- matrix(nrow =8933,ncol = ng)

# Collect genotype p-value of model including time into a single matrix.
for(j in 1:ng){
  name <-paste0("p_val_gene_",j)
  pval_gen_time_lm[,j] <- get(name)[1,]
}

saveRDS(pval_gen_time_lm, file = "~/Documents/elegans/pvals_and_effect_size/pval_gen_time_lm_one_perm.RDS") 

pval_dev_time_lm <- matrix(nrow =8933,ncol = ng)

# Put p-value time of model with time into a single matrix.
for(j in 1:ng){
  name <-paste0("p_val_gene_",j)
  pval_dev_time_lm[,j] <- get(name)[2,]
}

saveRDS(pval_dev_time_lm, file = "~/Documents/elegans/pvals_and_effect_size/pval_dev_time_lm_one_perm.RDS")
