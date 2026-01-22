samples <- read.table("/wsu/home/hm/hm80/hm8004/piquelab/Dany/gxp1/counts/names.txt", sep="\t")
samples <- as.character(samples[,1])

n=60612
matrix <- matrix(nrow=n)
df <- data.frame(matrix)
for (i in 1:length(samples)){
  read <- read.table(paste0(samples[i], ".cnts"), sep="\t")
  mapped <- sum(read[1:(nrow(read)-5),2])/sum(read[,2])
  HTseq_counted <- sum(read[,2])
  table_write <- t(c(samples[i],mapped, HTseq_counted))
  write.table(table_write, file="./HUVEC_alignment_HTseq.txt", sep="\t", append=TRUE, quote=FALSE, col.names=FALSE, row.names=FALSE)
  write.table(nrow(read), file="./HUVEC_nrow_HTSeq.txt", sep="\t", append=TRUE, quote=FALSE, col.names=FALSE, row.names=FALSE)
  df[,i] <- read[1:n,2]
  colnames(df)[i] <- samples[i]
}

rownames(df) <- read[1:n,1]

write.table(df, file="./HUVEC_HTseq_gene_counts.txt", quote=FALSE)

HTseq <- read.table("./HUVEC_alignment_HTseq.txt", sep="\t")
colnames(HTseq) <- c("barcode", "fraction_mapped", "HTseq_total_count")
write.table(HTseq, file="./HUVEC_alignment_HTseq.txt", sep="\t",quote=FALSE, row.names=FALSE)


