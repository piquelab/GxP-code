#!/bin/bash
# Job name
#SBATCH --job-name aligning_220307
# Submit to the primary QoS
#SBATCH -q primary
# Request one node
#SBATCH -N 1
# Total number of cores, in this example it will 1 node with 1 core each. 
#SBATCH -n 8
# Request memory
#SBATCH --mem=100G
# Request a node with avx2 instruction set 
#SBATCH --constraint=avx2
# Mail when the job begins, ends, fails, requeues 
#SBATCH --mail-type=ALL
# Where to send email alerts
#SBATCH --mail-user=hm8004@wayne.edu
# Create an output file that will be output_<jobid>.out 
#SBATCH -o output_%j.out
# Create an error file that will be error_<jobid>.out
#SBATCH -e errors_%j.err
# Set maximum time limit 
#SBATCH --time=1-01:00:00



module load hisat2/2.0.4
module load samtools/1.4


###
set=220307_NS500258_0517_AHJF3JBGXK

#Update these if copied from another directory
filePath=/wsu/home/groups/piquelab/OurData/Nextseq/gxp/${set}


genomeindex=/wsu/home/groups/piquelab/fungei/atac/hisat.ref/grch38_snp/genome_snp

###Align Reads###

hisat2 -p 8 -x ${genomeindex} -1 ${filePath}/${var}_R1_001.fastq.gz \
                              -2 ${filePath}/${var}_R2_001.fastq.gz \
      2> ${var}_aligned.bam.e | samtools view -b1 - > ${var}_aligned.bam



###Sort Reads###
samtools sort -@ 4 -T tmp_${var}_aligned.bam -o ${var}_sorted.bam ${var}_aligned.bam
samtools index ${var}_sorted.bam
samtools view -c ${var}_sorted.bam > ${var}_sorted_count.txt


###Quality Filter###
samtools view -b1 -q10 ${var}_sorted.bam > ${var}_quality.bam
samtools index ${var}_quality.bam
samtools view -c ${var}_quality.bam > ${var}_quality_count.txt


###Deduplication###
samtools rmdup ${var}_quality.bam ${var}_clean.bam
samtools index ${var}_clean.bam
samtools view -c ${var}_clean.bam > ${var}_clean_count.txt

echo ${var} >> finished.txt
