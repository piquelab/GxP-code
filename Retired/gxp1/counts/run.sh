set=220307_NS500258_0517_AHJF3JBGXK

ls /wsu/home/groups/piquelab/OurData/Nextseq/GxP/${set} | while read fastq_file; do echo ${fastq_file%_R*}; done | uniq > names.txt
cat names.txt | while read var; do sbatch --export=var=${var} HTcnt.sh; done
