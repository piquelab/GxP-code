# this script fits mash on all the preprocessed SCAIP FastQTL eQTL mapping results and saves the model to be run on chunked data
# based on https://stephenslab.github.io/mashr/articles/eQTL_outline.html
# based on SCAIP-genetic/mashr_eQTL/mashr-fit-model.R
# 9/16/25 MS

library(ashr)
library(mashr)

library(data.table)
library(dplyr)
library(doParallel)
cores <- as.integer(Sys.getenv("SLURM_STEP_TASKS_PER_NODE"))
registerDoParallel(cores = cores)
library(RhpcBLASctl)
blas_set_num_threads(1)


### Only run this section once to generate mashr object input
# 1. read in all the data
slope_files <-list.files(path=paste0("input/"),pattern=".*_slope.txt",full.name=T)
datasets <- gsub("_slope.txt","",list.files(path=paste0("input/"),pattern=".*_slope.txt"))
SE_files <- list.files(path=paste0("input/"),pattern=".*_SE.txt", full.name=T)
p_files <- list.files(path=paste0("input/"),pattern=".*_pvalue.txt", full.name=T)
lfsr_files <- list.files(path=paste0("input/"),pattern=".*lfsr.txt",full.name=T)

# 2. create list of files
ldf_slope <- lapply(slope_files, read.table, sep="\t", header=T)
ldf_SE <- lapply(SE_files, read.table, sep="\t", header=T)
ldf_p <- lapply(p_files, read.table, sep="\t", header=T)
ldf_lfsr <- lapply(lfsr_files, read.table, sep="\t", header=T)

# 3. merge them into one data frame:
slopes <- Reduce(function(x,y) merge(x = x, y = y, all=TRUE), ldf_slope)
SEs <- Reduce(function(x,y) merge(x = x, y = y, all=TRUE), ldf_SE)
pvalues <- Reduce(function(x,y) merge(x = x, y = y, all=TRUE), ldf_p)
lfsrs <- Reduce(function(x,y) merge(x = x, y = y, all=TRUE), ldf_lfsr)
                
# 4. clean up
slopes <- slopes[!duplicated(slopes$uniqID),]
rownames(slopes) <- slopes[,1]
slopes <- slopes[,-1]

SEs <- SEs[!duplicated(SEs$uniqID),]
rownames(SEs) <- SEs[,1]
SEs <- SEs[,-1]

pvalues <- pvalues[!duplicated(pvalues$uniqID),]
rownames(pvalues) <- pvalues[,1]
pvalues <- pvalues[,-1]
                
lfsrs <- lfsrs[!duplicated(lfsrs$uniqID),]
rownames(lfsrs) <- lfsrs[,1]
lfsrs <- lfsrs[,-1]
                
# 5. save as input
save(slopes, SEs, pvalues, lfsrs,file=paste0("./mashr_input.Rd"))


### Begin mashr model fitting #####
### 1. Read in and clean up input file
load(paste0("./mashr_input.Rd"))
                
# Remove rows where ANY column has NA or Inf values in the SE matrix
SEs <- SEs[!rowSums(is.finite(as.matrix(SEs)))<ncol(SEs),]
                
# Remove the same rows from pvalue and slope matrixes:
slopes <- slopes[rownames(SEs),]
pvalues <- pvalues[rownames(SEs),]
                
# Check:
sum(rowSums(is.na(as.matrix(slopes))>0))
sum(rowSums(is.na(as.matrix(pvalues))>0))
stopifnot(identical(rownames(slopes),rownames(SEs)))
stopifnot(identical(rownames(slopes),rownames(pvalues)))

# Save final SNP-gene list
write.table(rownames(pvalues),paste0("rownames_mash.txt"),row.names=F,col.names=F,quote=F)

# Sample 200k random tests:
set.seed(12)
sset <- sample(c(1:nrow(slopes)),200000)
data.temp = mash_set_data(as.matrix(slopes[sset,]), as.matrix(SEs[sset,]))
                
# calculate Vhat:
Vhat = estimate_null_correlation_simple(data.temp)
rm(data.temp)
data = mash_set_data(as.matrix(slopes), as.matrix(SEs),V=Vhat)
data.random <- mash_set_data(as.matrix(slopes[sset,]), as.matrix(SEs[sset,]),V=Vhat)

### 2. Set up the covariance matrices
# 2.1. select strong signals
m.1by1 = mash_1by1(data)
strong = get_significant_results(m.1by1, 0.05)

# 2.2. Obtain initial DATA-DRIVEN covariance matrices
# this step estimates the covariances of the observed data
U.pca = cov_pca(data,5,subset=strong)
print(names(U.pca)) 
                
# 2.3. Apply Extreme Deconvolution
# this step estimates the covariances of the actual underlying effects
U.ed = cov_ed(data, U.pca, subset=strong)

# 2.4. Obtain CANONICAL matrix:
U.c = cov_canonical(data.random)

# save all objects created up to now:
save(data,U.c,U.ed,Vhat, file="mashr-fit-model_objects.Rd")

### 3. Fit the model on sampled data
Sys.time()
m   = mash(data.random, c(U.c, U.ed), outputlevel=1)
Sys.time()
                
save(data,data.random,m,Vhat,file=paste0("mash-model-fit.Rd"))
                
# save the colnames:
write.table(colnames(pvalues),"TensorQTL_all_conditions-colnames.txt",sep="\n",col.names=F, row.names=F, quote=F)

### END


## # select a good number of chunks:
## tests <- nrow(data$Bhat)
## for(i in 1:round(tests/2)){
## div <- tests/i
## if(floor(div)==div){
##     print(paste0(i,"  ",div))
##         }
## }
