# GxP-code
This repository contains all scripts files of the GxP project. This repository has been superseeded by the (GxP)[github.com:piquelab/GxP.git]

# General report
This document contains all the analysis performed so far for the gxp project.
I processed the raw fastq files through a several steps, including raw reads alignment, quality check, sorting and  annotation, using a pipeline that incorporates various tools such as FastQC, HISAT2, Samtools, HTseq and DEseq2.

List of batches analysed:
GxP1: 

•	220324_NS500258_0518_AHJG33BGXK (considered as old sequences)

•	221020_NS500258_0531_AHJF3GBGXK (considered as new sequences)

Then merged 

 
Script used:
For each batch, it was created a folder dedicated for the analysis and results.
Inside each folder(batch), three folders were created; “bams” for the alignment, “counts” for reads counting, and “deseq” for the differential gene expression analysis.
First step, in the “bams” folder the following scripts were added and launched.
run.sh
```
set={Folder to the fastq files for the specific batch}
ls /wsu/home/groups/piquelab/OurData/Nextseq/GxP/ ${set} | while read fastq_file; do echo ${fastq_file%_R*}; done | uniq > names.txt
cat names.txt | while read var; do sbatch --export=var=${var} align.sh; done
```
#align.sh is the bash script file that is being called by run.sh

align.sh
```
#!/bin/bash
# Job name
#SBATCH --job-name alignment
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
set={Folder to the fastq files for the specific batch}

#Update this "path" if copied from another directory
filePath= /wsu/home/groups/piquelab/OurData/Nextseq/GxP/${set}

genomeindex=/wsu/home/groups/piquelab/fungei/atac/hisat.ref/grch38_snp_tran/genome_snp_tran

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
```

summarize.sh
```
#!/bin/bash
# Job name
#SBATCH --job-name alignmentsummary
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

#Making counts.csv
for i in `ls *sorted_count.txt`; do echo "$i"| cut -d_ -f1,2,3 >m ; cat m| tr "\n" " "; cat $i; done > sorted.txt
for i in `ls *quality_count.txt`; do echo "$i"| cut -d_ -f1,2,3 >m ; cat m| tr "\n" " "; cat $i; done > quality.txt
for i in `ls *clean_count.txt`; do echo "$i"| cut -d_ -f1,2,3 >m ; cat m| tr "\n" " "; cat $i; done > clean.txt

join sorted.txt quality.txt > tmp.txt
join tmp.txt clean.txt > all_counts.txt

rm m
rm tmp.txt
rm sorted.txt
rm quality.txt
rm clean.txt
```
#To make list of total number of reads processed by HTcnts
for i in `cat names.txt`; do echo -n "$i " >> total_reads.txt; head -1 ${i}_aligned.bam.e|cut -d ' ' -f1 >> total_reads.txt;done


#names.txt This file contains a list of all fastq files found and used for the alignment step. It is generated by run.sh
 
Second step, in the folder “counts” the following scripts were added and launched

run.sh
```
set={Folder to the fastq files for the specific batch}
ls /wsu/home/groups/piquelab/OurData/Nextseq/GxP/ ${set} | while read fastq_file; do echo ${fastq_file%_R*}; done | uniq > names.txt
cat names.txt | while read var; do sbatch --export=var=${var} HTcnt.sh; done
mkdir counts2 #all file containing counts will be move to this folder
```
#HTcnt.sh is the bash script file that is being called by run.sh

HTcnt.sh
```
#!/bin/bash
# Job name
#SBATCH --job-name readscounting
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
#SBATCH --mail-user=hx@wayne.edu
# Create an output file that will be output_<jobid>.out 
#SBATCH -o output_%j.out
# Create an error file that will be error_<jobid>.out
#SBATCH -e errors_%j.err
# Set maximum time limit 
#SBATCH --time=1-01:00:00

module unload python
module load anaconda3.python
source activate htseq

gtffile=/nfs/rprdata/ALOFT/AL/bams/Homo_sapiens.GRCh37.75.gtf.gz

htseq-count ../bams/${var}_clean.bam  $gtffile --stranded=reverse -f bam  > counts2/$var.cnts

echo ${var} >> Finished.txt
```
 
Move to the folder “counts2” which was previously created with run.sh script. Then add and launch the following script

summarize.sh
```
#!/bin/bash
# Job name
#SBATCH --job-name readscountingsummary
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


a=`head -n 1 ../names.txt`
## cd counts
cut -f 1 "$a".cnts >counts.txt
for i in ` cat ../names.txt`
do
	cut -f 2 "$i".cnts >tmp1
	cp counts.txt tmp2
	paste tmp2 tmp1 >counts.txt
done
# add header
# add header
echo "Genes" >tmp1.txt
`cat ../names.txt >tmp2.txt`
cat tmp1.txt tmp2.txt|paste -s -d '\t' >tmp3.txt
cp counts.txt tmp1.txt
cat tmp3.txt tmp1.txt >counts.txt
rm tmp1.txt
rm tmp2.txt
rm tmp3.txt
rm tmp1
rm tmp2
```
 

