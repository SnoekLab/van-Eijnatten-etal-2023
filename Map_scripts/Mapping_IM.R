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

# Subset data  to exclude four parent strains.
use.lines <- c(1:198,207)
pco1.lines <- pco1[use.lines]

# Same but also gets p-values of interaction term.
fast.map.co.int <- function(expla,rge){
  ok <- fastLm(cbind(1,expla),as.numeric(rge))
  pval <- pt(abs(ok$coefficients)/ok$stderr,ok$df.residual,lower.tail = F) * 2
  pval_gen <- pval[2]
  pval_dev <- pval[3]
  pval_int <- pval[4]
  out <- as.numeric(c(pval_gen,pval_dev,pval_int,as.numeric(ok$coefficients[2]),as.numeric(ok$coefficients[3]),as.numeric(ok$coefficients)[4]))
  return(out)
}

# Same as above, but now store interaction p-value in 3rd column p_val_gene_i.
for(i in 1:ng){
  gexpr <- use.rge[i,use.lines]
  assign(paste0("p_val_gene_",i), apply(FourP$snp_mat[,use.lines],1,function (x) fast.map.co.int(data.frame(as.numeric(x),pco1.lines,as.numeric(x)*as.numeric(pco1.lines)),gexpr)))
}

pval_gen_time_intlm <- matrix(nrow =8933,ncol = ng)

# Put genotype p-value in single matrix.
for(j in 1:ng){
  name <-paste0("p_val_gene_",j)
  pval_gen_time_intlm[,j] <- get(name)[1,]
}

saveRDS(pval_gen_time_intlm, file = "~/Documents/elegans/pvals_and_effect_size/pval_gen_time_intlm.RDS")

pval_dev_time_intlm <- matrix(nrow =8933,ncol = ng)

# Put time p-value in single matrix
for(j in 1:ng){
  name <-paste0("p_val_gene_",j)
  pval_dev_time_intlm[,j] <- get(name)[2,]
}

saveRDS(pval_dev_time_intlm, file = "~/Documents/elegans/pvals_and_effect_size/pval_dev_time_intlm.RDS")

pval_int_time_intlm <- matrix(nrow =8933,ncol = ng)

#put Interaction p-value in single matrix.
for(j in 1:ng){
  name <-paste0("p_val_gene_",j)
  pval_int_time_intlm[,j] <- get(name)[3,]
}

saveRDS(pval_int_time_intlm, file = "~/Documents/elegans/pvals_and_effect_size/pval_int_time_intlm.RDS")

### Difference in intersect between the two models

effect_imm_model <- matrix(nrow =8933,ncol = ng)

for(j in 1:ng){
  name <-paste0("p_val_gene_",j)
  effect_imm_model[,j] <- get(name)[4,]
}

saveRDS(effect_imm_model, file = "~/Documents/elegans/pvals_and_effect_size/effect_imm_model.RDS")

### Slope of allele 0

slope_allele_0_int_model <- matrix(nrow =8933,ncol = ng)

for(j in 1:ng){
  name <-paste0("p_val_gene_",j)
  slope_allele_0_int_model[,j] <- get(name)[5,]
}

saveRDS(slope_allele_0_int_model, file = "~/Documents/elegans/pvals_and_effect_size/slope_allele_0_int_model.RDS")

### Difference in slope between allele 0 and 1. 

effect_imi_model <- matrix(nrow =8933,ncol = ng)

for(j in 1:ng){
  name <-paste0("p_val_gene_",j)
  effect_imi_model[,j] <- get(name)[6,]
}

saveRDS(effect_imi_model, file = "~/Documents/elegans/pvals_and_effect_size/effect_imi_model.RDS")






