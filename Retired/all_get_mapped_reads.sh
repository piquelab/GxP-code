#!/bin/bash
# Job name
#SBATCH --job-name rawreads_230712
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

folder="../gxp7/bams"  # Replace with the path to your BAM files folder

samples=("$folder"/GxP*quality*.bam)

total_mapped_reads=0
sample_count=0

output_file="all_mapped_reads.csv"  # Replace with the desired output CSV file name

echo "Sample,Mapped Reads" > "$output_file"

for sample in "${samples[@]}"; do
    mapped_reads=$(samtools idxstats "$sample" | awk '{mapped+=($3+$4)} END {print mapped}')
    echo "$(basename "$sample"),$mapped_reads" >> "$output_file"

    total_mapped_reads=$((total_mapped_reads + mapped_reads))
    sample_count=$((sample_count + 1))
done


average_mapped_reads=$((total_mapped_reads / sample_count))

echo "Average Mapped Reads: $average_mapped_reads" >> "$output_file"
echo $average_mapped_reads
