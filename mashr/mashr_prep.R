# Generate the input for mashr (pvalues,SEs)
# based on SCAIP-genetic/mashr_eQTL/mashr_prep.R
# 9/16/25 MS

# RUN:  /wsu/el7/groups/piquelab/R/4.1.0/bin/Rscript mashr_prep.R

library(ashr)
library(data.table)

eQTL_dir <- "/rs/rs_grp_gxp/RNAseq_analysis/GxP_20250730/eQTL_mapping/tensor/tensorqtl_output_cis-nominal_SV15_100kb/"

# read in, grab effect estimate and its standard error
m <- fread(paste0(eQTL_dir, dataset, "_all_nominal_pairs.txt.gz"),sep="\t")
colnames(m) <- c("phenotype_id","variant_id", "start_distance", "end_distance", 
                  "af", "ma_samples", "ma_count", "pval_nominal", "slope", "slope_se", "condition")

# add unique SNP identifier:
m$uniqID <- paste0(m$phenotype_id,"_",m$variant_id)

# get list of unique conditions
conditions <- unique(m$condition)

# process each condition
for(dataset in conditions) {
    # subset data for this condition
    m_sub <- m[condition == dataset]
    
    # calculate lfsr for this condition:
    ashr <- ash(m_sub$slope, m_sub$slope_se)
    m_sub$lfsr <- ashr$result$lfsr
    
    # save the slope, SE and pvals in separate files:
    slope <- m_sub[,c("uniqID","slope")]
    SE <- m_sub[,c("uniqID","slope_se")]
    pvals <- m_sub[,c("uniqID","pval_nominal")]
    lfsr <- m_sub[,c("uniqID","lfsr")]
    
    # rename columns
    colnames(slope)[2] <- dataset
    colnames(SE)[2] <- dataset
    colnames(pvals)[2] <- dataset
    colnames(lfsr)[2] <- dataset
    
    # write files
    write.table(slope, paste0("input/", dataset, "_slope.txt"), sep="\t", col.names=T, row.names=F, quote=F)
    write.table(SE, paste0("input/", dataset, "_SE.txt"), sep="\t", col.names=T, row.names=F, quote=F)
    write.table(pvals, paste0("input/", dataset, "_pvalue.txt"), sep="\t", col.names=T, row.names=F, quote=F)
    write.table(lfsr, paste0("input/", dataset, "_lfsr.txt"), sep="\t", col.names=T, row.names=F, quote=F)
}

# gzip all files at the end:
system("for file in input/*txt ; do gzip $file; done")

### END
