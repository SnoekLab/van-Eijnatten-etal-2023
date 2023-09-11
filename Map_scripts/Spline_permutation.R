.libPaths("~/R/x86_64-pc-linux-gnu-library/3.6/")

library(splines)

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


lm_spline <- function(gexpr,marker,dev){
  lm_gen_dev_int <- lm(gexpr~marker*ns(dev,knots = quantile(dev)[2:4]))
  lm_gen_dev <- lm(gexpr~marker + ns(dev,knots = quantile(dev)[2:4]))
  lm_gen_int <- lm(gexpr~marker + I(marker*ns(dev,knots = quantile(dev)[2:4])))
  lm_dev_int <- lm(gexpr~ns(dev,knots = quantile(dev)[2:4])+I(marker*ns(dev,knots = quantile(dev)[2:4])))
  pval_gen <- anova(lm_gen_dev_int,lm_dev_int)$`Pr(>F)`[2]
  pval_ns_dev <- anova(lm_gen_dev_int,lm_gen_int)$`Pr(>F)`[2]
  pval_int <- anova(lm_gen_dev_int,lm_gen_dev)$`Pr(>F)`[2]
  out <- as.numeric(c(pval_gen,pval_ns_dev,pval_int))
  return(out)
}

for(i in 1:ng){
  gexpr <- sample(use.rge[i,use.lines])
  assign(paste0("p_val_gene_",i), apply(FourP$snp_mat[,use.lines],1,function (x) lm_spline(gexpr,as.numeric(x),pco1.lines)))
}

pval_gen_spline <- matrix(nrow =8933,ncol = ng)

# Put genotype p-value in single matrix.
for(j in 1:ng){
  name <-paste0("p_val_gene_",j)
  pval_gen_spline[,j] <- get(name)[1,]
}

saveRDS(pval_gen_spline, file = "~/Documents/elegans/redo_perm_one_perm_per_gene/pval_gen_spline_2model_one_perm.RDS")

pval_dev_spline <- matrix(nrow =8933,ncol = ng)

# Put time p-value in single matrix
for(j in 1:ng){
  name <-paste0("p_val_gene_",j)
  pval_dev_spline[,j] <- get(name)[2,]
}

saveRDS(pval_dev_spline, file = "~/Documents/elegans/redo_perm_one_perm_per_gene/pval_dev_spline_2model_one_perm.RDS")

pval_int_spline <- matrix(nrow =8933,ncol = ng)

# put Interaction p-value in single matrix.
for(j in 1:ng){
  name <-paste0("p_val_gene_",j)
  pval_int_spline[,j] <- get(name)[3,]
}

saveRDS(pval_int_spline, file = "~/Documents/elegans/redo_perm_one_perm_per_gene/pval_int_spline_2model_one_perm.RDS")









