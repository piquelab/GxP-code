# GxP RNAseq Data Analysis
Environment from Cindy Kalita: `/wsu/home/groups/piquelab/cindy/tensorqtl.environment.yaml`   
Set up environment with:    
```{r}
conda activate tensorqtl_p3.11_env
export PATH=/wsu/el7/groups/piquelab/R/4.3.2/bin:${PATH}   
export LD_LIBRARY_PATH=/wsu/el7/groups/piquelab/R/4.3.2/lib64/R/lib:$LD_LIBRARY_PATH
```

# Scripts
- Outlier_removal.R  
- ComBat.R  
- Filter_expression.R  
- PCA_vis.R    
- Filtered_SVA.R  
- Limma.R  
- Iterative_limma_SVs.R  
- DEG_QQ.R  
- DEG_summary.R  
- Prep_QTL.R  
- Prep_genotypes.R  
- 01_tensorQTL.sh  
- 02_tensorQTL.sbatch  
- tensorQTL_analysis.R   

# Key Files and Folder Structure

### Input Data
```
/rs/rs_grp_gxp/RNAseq_analysis/GxP_01092026/
├── GxP_SamplesRemoved_01182026.RData          # Initial filtered data
├── GxP_Filtered_02052026.RData                 # Filtered & normalized by timepoint
└── GxP-eQTL_DNA_genotypes_filtered_SNPs_noChr.vcf.gz  # Genotype data
```

### Differential Expression Analysis
```
/rs/rs_grp_gxp/RNAseq_analysis/GxP_01092026/
├── DEG/
│   ├── {date}_BPA_T{6,24}.txt                 # BPA vs H2O results
│   ├── {date}_CRL_T{6,24}.txt                 # EtOH vs H2O results
│   ├── {date}_MBP_T{6,24}.txt                 # MBP vs EtOH results
│   ├── {date}_AllDEGs_FDR10.txt               # All DEGs combined
│   └── {date}_DEG_Summary_FDR10.txt           # DEG counts
├── DEG_SVs/
│   ├── {date}_BPA_T{6,24}_SV{0,1-2,...1-25}.txt
│   ├── {date}_CRL_T{6,24}_SV{0,1-2,...1-25}.txt
│   ├── {date}_MBP_T{6,24}_SV{0,1-2,...1-25}.txt
│   └── {date}_DEG_SVs_Summary_FDR10.txt       # DEG counts across SV iterations
├── Filtered_SVs/
│   ├── SV_tp{6,24}.RData                      # SVs for DEG analysis
│   └── SV_tp{6,24}.txt
└── PCA_analysis/
    ├── {date}_PCA_PC1-4_BothTimepoints_{covariate}.png
    └── {date}_PCA_VarianceExplained.png
```

### QTL Mapping
```
/rs/rs_grp_gxp/RNAseq_analysis/GxP_01092026/
├── GxP-eQTL_{BPA,H2O,MBP,EtOH}_T{6,24}_qnorm.bed.gz  # Expression BED files
├── GxP-eQTL_{BPA,H2O,MBP,EtOH}_T{6,24}_qnorm.sorted.bed.gz  # Sorted & indexed
├── Norm_SVs/
│   ├── SV_{BPA,H2O,MBP,EtOH}_T{6,24}.txt      # SVs by condition with dbGaP IDs
│   └── SV_{BPA_H2O,MBP_EtOH}_T{6,24}.RData    # Original pair-wise SVs (deprecated)
└── tensor/
    ├── GxP-eQTL_DNA_genotypes_filtered_SNPs_noChr.pgen  # PLINK format genotypes
    ├── GxP-eQTL_DNA_genotypes_filtered_SNPs_noChr.pvar
    ├── GxP-eQTL_DNA_genotypes_filtered_SNPs_noChr.psam
    ├── GxP-eQTL_genotypes_pruned.prune.in     # LD-pruned variants
    ├── genotype_pcs_tensorqtl.txt             # Top 5 genotype PCs
    ├── covariates/
    │   ├── SV0/
    │   │   └── covariates_{condition}_T{tp}.txt  # 3 genotype PCs only
    │   ├── SV1-2/
    │   │   └── covariates_{condition}_T{tp}.txt  # 3 genotype PCs + 2 SVs
    │   └── SV1-{3...20}/                       # Up to 3 genotype PCs + 20 SVs
    ├── SV0/
    │   ├── tensorqtl_conditions.tsv            # Condition table for this iteration
    │   └── {condition}_T{tp}_cis.cis_qtl.txt.gz  # tensorQTL results
    ├── SV1-2/
    │   ├── tensorqtl_conditions.tsv
    │   └── {condition}_T{tp}_cis.cis_qtl.txt.gz
    ├── SV1-{3...20}/                           # Additional SV iterations
    └── analysis/
        ├── {date}_QQ_tensorQTL_{condition}_T{tp}.png  # QQ plots
        └── {date}_eGene_Summary_FDR10.txt      # eGene counts across iterations
```

### QC Plots
```
/rs/rs_grp_gxp/RNAseq_analysis/GxP_01092026/
├── PCA_analysis/
│   └── {date}_PCA_*.png
└── QQ_plots/
    ├── {date}_QQ_DEG_{contrast}.png           # DEG QQ plots (timepoints overlaid)
    └── {date}_QQ_DEG_SVs_{contrast}_T{tp}.png # DEG_SVs QQ plots (iterations overlaid)
```

### Key Variables
- **Timepoints**: 6, 24 (hours)
- **Conditions**: BPA, H2O, MBP, EtOH
- **Treatment pairs**: BPA_H2O, MBP_EtOH
- **Contrasts**: BPA (BPA vs H2O), CRL (EtOH vs H2O), MBP (MBP vs EtOH)
- **SV iterations**: 0, 1-2, 1-3, ... 1-20, 1-25 (DEG analysis)
- **SV iterations**: 0, 1-2, 1-3, ... 1-20 (QTL analysis)
- **FDR threshold**: 0.10

## Outlier removal 
`Rscript Outlier_removal.R`  
Input: `/rs/rs_grp_gxp/RNAseq_analysis/GxP_20250730/GxP_774Samples_20250803.RData`  
Output: `GxP_SamplesRemoved_01182026.RData`  
Logged sample: `samples_removed_01182026.csv`  

We used diagonal linear discriminant analysis (DLDA) on log-transformed, normalized read counts, using significantly differentially expressed genes as features, 
to detect potential swaps between treatment samples and their matched controls. For each treatment contrast, DLDA scores were computed for each sample based on its 
expression profile across the selected gene set, including the treatment sample, its primary matched control, and a secondary matched control. We then calculated two differences 
for each cell line:  
    
∆DLDA1=DLDAtreatment-DLDAprimary control and  
∆DLDA2=DLDAprimary control-DLDAsecondary control.  
    
Cell lines were excluded from downstream analyses if both of the following criteria were met:  
	(1) the sign of ∆DLDA1 is opposite to the consensus direction observed across cell lines and is an outlier, defined as exceeding 1.5 times the interquartile range  
	(2) the sign of ∆DLDA2 is the same as the consensus direction observed across cell lines for ∆DLDA1.  
	
**The following samples were indicated for removal:**   
### BPA6
| Batch | dbGaPID |
|-------|---------|
| GxP12 | LP-102  |
| GxP12 | LP-104  |
| GxP37 | LP-052  |

### MBP6
| Batch | dbGaPID |
|-------|---------|
| GxP11 | LP-056  |
| GxP22 | LP-011  |
| GxP22 | LP-017  |
| GxP24 | LP-038  |

### MBP24
| Batch | dbGaPID |
|-------|---------|
| GxP21 | LP-012  |

### Sample IDs pulled from covariate sheet:
| Sample ID | Treatment  | CellLine | dbGaPID | Timepoint |
|-----------|------------|----------|---------|-----------|
| GxP12_02  | BPA_100nM  | H121R    | LP-102  | 6         |
| GxP12_03  | BPA_100nM  | H123M    | LP-104  | 6         |
| GxP12_05  | H2O        | H121R    | LP-102  | 6         |
| GxP12_06  | H2O        | H123M    | LP-104  | 6         |
| GxP37_02  | BPA_100nM  | H022T    | LP-052  | 6         |
| GxP37_05  | H2O        | H022T    | LP-052  | 6         |
| GxP11_08  | MBP_500nM  | H029R    | LP-056  | 6         |
| GxP11_11  | EtOH       | H029R    | LP-056  | 6         |
| GxP22_07  | MBP_500nM  | H097D    | LP-011  | 6         |
| GxP22_09  | MBP_500nM  | H107M    | LP-017  | 6         |
| GxP22_10  | EtOH       | H097D    | LP-011  | 6         |
| GxP22_12  | EtOH       | H107M    | LP-017  | 6         |
| GxP24_07  | MBP_500nM  | H001S    | LP-038  | 6         |
| GxP24_10  | EtOH       | H001S    | LP-038  | 6         |
| GxP21_20  | MBP_500nM  | H098M    | LP-012  | 24        |
| GxP21_23  | EtOH       | H098M    | LP-012  | 24        |  

## ComBat 
`Rscript ComBat.R`  
Input: `GxP_SamplesRemoved_01182026.RData`   
Output: `GxP_ComBatCorrected_01182026.RData`   

Remove known systematic batch effects. **Uncorrected data used for all downstream analysis.**   

Uncorrected: `m2`    
Sequencer effects removed: `m2.seq`   
Batch effects remove: `m2.batch`   

## Filter and process expression data
`Rscript Filter_expression.R`   
Input: `GxP_SamplesRemoved_01182026.RData`   
Output: `GxP_Filtered_02052026.RData`   

Separate by time point, filter for expressed genes, TMM normalize, filter for autosomal coding genes, voom normalize   
`f = ~Treatment + dbGaP_ID + trimmed_dClean.dFastq`    

## PCA QC Check
`Rscript PCA_vis.R`    
Input:`GxP_Filtered_02052026.RData`    
Output: `/PCA_analysis`   

Visualize first 4 PCs colored by key covariates.   

## Calcualte SVs for differential gene expression
`Rscript Filtered_SVA.R`    
Input: `GxP_Filtered_02052026.RData`    
Output: `/Filtered_SVs`    

Estimate and calculate SVs for differential gene expression. Processed by timepoint.       
Model: ~Treatment + dbGaP_ID + trimmed_dClean.dFastq     
Null: ~dbGaP_ID + trimmed_dClean.dFastq     

Estimated SVs for T6: 50      
Estimate SVs for T24: 52     

## Differential gene expression 
`Rscript Limma.R T6`    
`Rscript Limma.R T24`    
Input: `GxP_Filtered_02052026.RData`,   
Output: `/DEG`    
   
Run limma using: `f = ~Treatment + dbGaP_ID + trimmed_dClean.dFastq`     

## Iterate on voom and limma with up to 25 SVs
`sbatch -q primary -N 1 -n 1 --mem=64G -t 12:00:00 --wrap "Rscript Iterative_limma_SVs.R"`    
Input:`GxP_Filtered_02052026.RData`     
Output: `/DEG_SVs`     

Run limma using: `f = ~Treatment + dbGaP_ID + trimmed_dClean.dFastq + ...SVs`     

## Generate overlaid QQ plots for all limma results 
`Rscript DEG_QQ.R`    
Input: `GxP_Filtered_02052026.RData`,      
Output: `/QQ_Plots`      

Generate QQ plots for DEG and DEG_SVs results     

## Call DEGs ato 10% FDR and write summary files 
`Rscript DEG_summary.R`     
Input: `/DEG`, `/DEG_SVs`   
Ouptput:`/DEG/02052026_AllDEGs_FDR10.txt` - All DEGs with gene info   
        `DEG/02052026_DEG_Summary_FDR10.txt` - Summary counts for DEG   
        `DEG_SVs/[date]_DEG_SVs_Summary_FDR10.txt` - Summary counts for all SV iterations     

Call DEGs and create summary tables and new files with stratified gene info     

## Prep expression data and calculate SVs for QTL mapping
`Rscript Prep_QTL.R`    
Input: `GxP_Filtered_02052026.RData`    
Output: `GxP-eQTL_{condition}_T{tp}_qnorm.bed.gz`, `Norm_SV/SV_{pair}_T{tp}.RData`    

Regress out trimmed_dClean.dFastq, separate by treatment pairs, quantile normalize, calculate SVs, separate by condition, save bed files and SVs with dbgap info.        

Estimated SVs for BPA_6: 24    
Estimated SVs for MBP_6: 25   
Estimated SVs for BPA_24: 25    
Estimated SVs for MBP_24: 23   

## Prep genotype data and calculate PCs for QTL mapping 
```{r}
module load plink/2.0 
Rscript Prep_genotypes.R
```
Input: `GxP-eQTL_DNA_genotypes_filtered_SNPs_noChr.vcf.gz`   
Output: `/tensor`, `genotype_pcs_tensorqtl.txt`   

Filter genotype data, use plink formatting, perform LD pruning and generate genotype PCs   

## Run tensorQTL 
`bash 01_tensorQTL.sh`     
Input: `GxP-eQTL_{condition}_T{tp}_qnorm.sorted.bed.gz`, `Norm_SVs/SV_{condition}_T{tp}.txt`, `genotype_pcs_tensorqtl.txt`, `tensor/GxP-eQTL_DNA_genotypes_filtered_SNPs_noChr.[pgen|pvar|psam]`    
Output: `/tensor/SV{n}/`, `/tensor/covariates/SV{n}/`, `tensorqtl_conditions.tsv`    

Sort and index BED files, create covariate files combining top 3 genotype PCs with 0-20 SVs iteratively, generate condition tables for each SV iteration, submit tensorQTL job arrays.    

## tensorQTL batch submissions
02_tensorQTL.sbatch (submitted via 01_tensorQTL.sh)   
Input: PLINK genotypes, sorted BED files, covariate files    
Output: `/tensor/SV{n}/{condition}_T{tp}_cis.cis_qtl.txt.gz`   

Runs cis-eQTL mapping with 100kb window and FDR 0.1 for each condition-timepoint-SV iteration.   

## tensorQTL QQ plots and eGene calling 
`Rscript tensorQTL_analysis.R`    
Input: `/tensor/SV{n}/{condition}_T{tp}_cis.cis_qtl.txt.gz`   
Output: `/tensor/analysis/`, `{date}_QQ_tensorQTL_{condition}_T{tp}.png`, `{date}_eGene_Summary_FDR10.txt`   

Generate QQ plots for each condition-timepoint with all SV iterations overlaid.    
Create summary table showing number of eGenes at 10% FDR for each condition-timepoint-SV combination.   
