#!/bin/bash
# Job name
#SBATCH --job-name countingSummary_220307
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
