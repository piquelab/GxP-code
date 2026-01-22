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


#To make list of total number of reads processed by HTcnts
for i in `cat finished.txt`; do echo -n "$i " >> total_reads.txt; head -1 ${i}_aligned.bam.e|cut -d ' ' -f1 >> total_reads.txt;done
