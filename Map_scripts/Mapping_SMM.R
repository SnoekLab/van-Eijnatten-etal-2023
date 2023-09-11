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

# Function which tests a linear model of gene expression explained by genotype, calculates p-value by hand.
fast.map <- function(geno,rge){
  ok <- fastLm(rge~geno)
  pval <- pt(abs(ok$coefficients)/ok$stderr,ok$df.residual,lower.tail = F) * 2
  pval_gen <- pval[2]
  out <- as.numeric(c(pval_gen,ok$coefficients[2])) #,R2,R2.adj)
  return(out)
}

# Subset data  to exclude four parent strains.
use.lines <- c(1:198,207)
pco1.lines <- pco1[use.lines]

# For each gene, calculate p-value of genotype in single variable model. Store as vector in p_val_gen_i where i is the row number of gene in use.rge.
for(i in 1:ng){
  gexpr <- use.rge[i,use.lines]
  assign(paste0("p_val_gene_",i), apply(FourP$snp_mat[,use.lines],1,function (x) fast.map(as.numeric(x),gexpr)))
}

# Now I have to collect all these p-values into a single matrix. Define matrix below.
pval_gen_notime_lm <- matrix(nrow =8933, ncol = ng)
effect_ssm_model <- matrix(nrow = 8933, ncol = ng)


# Collect p-values and in column i put all the p-values of gene with row number i in use.rge.
for(j in 1:ng){
  name <-paste0("p_val_gene_",j)
  pval_gen_notime_lm[,j] <- get(name)[1,]
}

saveRDS(pval_gen_notime_lm, file = "~/Documents/elegans/pvals_and_effect_size/pval_ssm_model.RDS") 

for(j in 1:ng){
  name <-paste0("p_val_gene_",j)
  effect_ssm_model[,j] <- get(name)[2,]
}

saveRDS(effect_ssm_model, file = "~/Documents/elegans/pvals_and_effect_size/effect_ssm_model.RDS") 
