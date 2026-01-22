#!/bin/bash
# Job name
#SBATCH --job-name counting_220324
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
module unload python
module load anaconda3.python
source activate htseq


gtffile=/wsu/home/hm/hm80/hm8004/piquelab/pbmc_handls/BAMs/Reference/Homo_sapiens.GRCh38.103.chr.gtf.gz

htseq-count ../bams/${var}_clean.bam  $gtffile --stranded=reverse -f bam  > counts2/$var.cnts

echo ${var} >> Finished.txt
