### File structure:

/rs/rs_grp_gxp/RNAseq_analysis/phthalates/
```
│
├── scripts/...
├── READ.ME <- you are here
├── GxP_Filtered.RData
├── GxP-eQTL_MBP_T6_int.bed.gz
├── GxP-eQTL_MBP_T24_int.bed.gz
├── GxP-eQTL_EtOH_T6_int.bed.gz
├── GxP-eQTL_EtOH_T24_int.bed.gz
├── PCs/...
├── starter_files
│   ├── Outlier_removal.R
│   ├── samples_removed_01182026.csv
│   ├── GxP_SamplesRemoved_01182026.RData
│   ├── GxP-eQTL_DNA_genotypes_filtered_SNPs_noChr.vcf.gz
│   ├── GxP-eQTL_DNA_genotypes_filtered_SNPs_noChr.pgen/pvar/psam
│   ├── genotype_pcs_tensorqtl.txt
│   └── read.me
├── DEG
│   ├── 04202026_MBP_T6.txt
│   ├── 04202026_MBP_T24.txt
│   ├── 04202026_DEG_Summary_FDR10.txt
│   ├── 04202026_AllDEGs_FDR10.txt
│   ├── GO_enrichment/...
│   ├── permutations/...
│   └── figures/...
├── tensor_cis
│   ├── logs...
│   ├── analysis/ 
│       ├── 04202026_eGene_Summary_FDR10.txt
│       └── ...
│   ├── covariates/ 
│       ├── PC0/...
│       └── ...
│   ├── PC0/ 
│       ├── tensorqtl_conditions.tsv
│       ├── EtOH_T6_cis.cis_qtl.txt.gz
│       └── ...
│   └── ...
├── tensor
│   ├── covariates/ 
│       └── PC1-15/...
│   ├── PC1-15/ 
│       ├── tensorqtl_conditions.tsv
│       ├── MBP_T6_chr1_chr1.bed.gz
│       ├── MBP_T6_chr1.cis_qtl_pairs.1.parquet
│       └── ...
│   ├── txt_files/ 
│           └── PC1-15
│               ├── MBP_T6_cis_nominal.txt
│               ├── MBP_T24_cis_nominal.txt
│               ├── EtOH_T6_cis_nominal.txt
│               └── EtOH_T24_cis_nominal.txt
│   ├── combined_nominal/ 
│           └── PC1-15
│               └── tensorqtl_nominal_PC1-15_combined.txt.gz
│   └── logs/...
```

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


### DEGs: 
```
Rscript scripts/Filter_expression.R
Rscript scripts/Limma.R T6
Rscript scripts/Limma.R T24
Rscript scripts/DEG_summary.R
Rscript scripts/DEG_visualization.R
Rscript scripts/GO_ORA.R MBP 6
Rscript scripts/GO_ORA.R MBP 24
Rscript scripts/Limma_permutations.R 6 100
Rscript scripts/Limma_permutations.R 24 100
Rscript scripts/QQ.R 6 100
Rscript scripts/QQ.R 24 100
```
### tensorQTL:
```
module load htslib
module load plink/2.0
Rscript scripts/Prep_genotypes.R

conda activate tensorqtl_p3.11_env
export LD_LIBRARY_PATH=/wsu/el7/groups/piquelab/R/4.3.2/lib64/R/lib:$LD_LIBRARY_PATH
export PATH=/wsu/el7/groups/piquelab/R/4.3.2/bin:$PATH

# Run cis mode with 100kB window 
bash scripts/01_tensorQTL.sh
#### Decided on 15 genotype PCs

# Run cis_nominal mode with 1mB window and 15 PCs
bash scripts/01_tensorQTL.sh
python3 scripts/convert_parquet.py
python3 scripts/combine_txt.py

Rscript scripts/tensorQTL_analysis.R
```
### mash 
```
Rscript check_mash_requirements.R <- OK

bash run_mash_pipeline.sh

Rscript scripts/classify_eQTLs_and_genes.R
Rscript scripts/mashr_visualization.R
Rscript scripts/select_top_reQTLs.R
Rscript scripts/save_mash_rds.R

sbatch -q primary -N1 -n1 --mem=32G -t 3:00:00 \
  --job-name=eqtl_vis \
  --parsable \
  --wrap="Rscript reQTL_visualization.R /rs/rs_grp_gxp/RNAseq_analysis/phthalates/mash_results/04212026_MBP_T6_top_reQTLs_tensorqtl.txt"
  
  sbatch -q primary -N1 -n1 --mem=32G -t 3:00:00 \
  --job-name=eqtl_vis \
  --parsable \
  --wrap="Rscript reQTL_visualization.R /rs/rs_grp_gxp/RNAseq_analysis/phthalates/mash_results/04212026_MBP_T24_top_reQTLs_tensorqtl.txt"
```

### mash posterior sampling

```
sbatch -q highmem -N1 -n1 --mem=150G -t 24:00:00 \
  --constraint=avx2 \
  --job-name=mash_posterior_samples \
  --output=logs/mash_posterior_samples_%j.log \
  --wrap="module purge && module load R && Rscript /rs/rs_grp_gxp/RNAseq_analysis/phthalates/scripts/mashr_posterior_samples.R"
  
  Submitted batch job 36479305
  
Rscript scripts/save_mash_rds.R
  
  sbatch -q express -N1 -n1 --mem=60G -t 12:00:00 \
  --constraint=avx2 \
  --job-name=classification \
  --output=logs/mash_posterior_samples_%j.log \
  --wrap="module purge && module load R && Rscript /rs/rs_grp_gxp/RNAseq_analysis/phthalates/scripts/classify_eQTLs_and_genes.R"
  
  sbatch -q express -N1 -n1 --mem=60G -t 12:00:00 \
  --constraint=avx2 \
  --job-name=visualization \
  --output=logs/mash_posterior_samples_%j.log \
  --wrap="module purge && module load R && Rscript /rs/rs_grp_gxp/RNAseq_analysis/phthalates/scripts/mashr_visualization.R"

  sbatch -q express -N1 -n1 --mem=60G -t 12:00:00 \
  --constraint=avx2 \
  --job-name=selection \
  --output=logs/mash_posterior_samples_%j.log \
  --wrap="module purge && module load R && Rscript /rs/rs_grp_gxp/RNAseq_analysis/phthalates/scripts/select_top_reQTLs.R"
  
  sbatch -q express -N1 -n1 --mem=60G -t 12:00:00 \
  --constraint=avx2 \
  --job-name=mash_comp \
  --output=logs/mash_posterior_samples_%j.log \
  --wrap="module purge && module load R && Rscript /rs/rs_grp_gxp/RNAseq_analysis/phthalates/scripts/compare_mash.R"
  
  
  sbatch -q express -N1 -n1 --mem=60G -t 12:00:00 \
  --constraint=avx2 \
  --job-name=posterior_analysis \
  --output=logs/mash_posterior_samples_%j.log \
  --wrap="module purge && module load R && Rscript /rs/rs_grp_gxp/RNAseq_analysis/phthalates/scripts/mashr_posterior_analyze.R"

```
  
  
