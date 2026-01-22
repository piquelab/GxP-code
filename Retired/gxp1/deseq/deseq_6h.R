require(ggplot2)
library(DESeq2)
library(qvalue)
library(annotables)
library(dplyr)
require(BiocParallel)
library(knitr)
library(ggrepel)

timestamp()

cores <- 8

platePrefix <- "cardiovascular"
ParallelSapply <- function(...,mc.cores=cores){
  simplify2array(mclapply(...,mc.cores=mc.cores))
}

## To run DESeq2 in parallel, using the
register(MulticoreParam(cores))

## Gene counts: this is our data for analysis
#countFile <- "/wsu/home/hm/hm80/hm8004/piquelab/Dany/gxp/counts1_new/counts2_new/counts.txt"
#countFile <- "/wsu/home/hm/hm80/hm8004/piquelab/Dany/gxp/counts1_new/counts2_new/HUVEC_HTseq_gene_counts.txt"
countFile <- "/wsu/home/hm/hm80/hm8004/piquelab/Dany/gxp-new/counts-combined-new-and-old/counts2/HUVEC_HTseq_gene_counts.txt"

#data <- read.table(countFile,sep="\t",as.is=T,header=1,row.names=1)
data <- read.table(countFile,sep=" ",as.is=T,header=1,row.names=1)
colnames(data) <- substring(colnames(data), 1, 7)

#samples_6h_DEP_500nM <- c("GxP1_03","GxP1_23","GxP1_68")
#samples_24h_DEP_1nM <- c("GxP1_73","GxP1_28","GxP1_08")
#samples_24h_MBP_500nM <- c("GxP1_53","GxP1_71")

#removing the missing samples
#missing_samples <-c ("GxP1_22","GxP1_11","GxP1_33","GxP1_32","GxP1_44","GxP1_82")
missing_samples <-c ("GxP1_22","GxP1_11","GxP1_33","GxP1_32","GxP1_44","GxP1_82", samples_6h_DEP_500nM)
data <- data[,!(colnames(data) %in% missing_samples)]

cov.file <- "/wsu/home/groups/piquelab/Dany/gxp/GxP1_metainfo_mod.csv"
cv <- read.table(cov.file, stringsAsFactors=F, header=T, comment="", sep=",")
colnames(cv)[1] <- "Samples"

cv$treatment <- sub("_","",substring(cv$Treatment,1,4))
cv$concentration <- substring(cv$Treatment,5,20)


#removing spaces
cv$treatment <- trimws(cv$treatment, which = "both")
cv$concentration <- trimws(cv$concentration, which = "both")

#Adding the percentage of reads mapped on exons
percent_file <- "/wsu/home/hm/hm80/hm8004/piquelab/Dany/gxp/counts1_new/counts2_new/HUVEC_alignment_HTseq.txt"
percent_mapped <- read.table(percent_file,sep="\t",as.is=T,header=1,row.names=1)

cv <- cv[!(cv$Samples %in% missing_samples),]
rownames(percent_mapped) <-  (substring(rownames(percent_mapped),1,7))
percent_mapped <- percent_mapped[cv$Samples,]

dim(percent_mapped)
dim(cv)
identical(cv$Samples, rownames(percent_mapped))
cv <- cbind(cv,percent_mapped)

#79 samples and 10 features/variables

kable(head(cv))

#Removing samples with ...
#samples_to_remove <- cv[cv$Treatment %in% c( "DEP_500nM"),]$Samples
#cv <- cv[!(cv$Samples %in% samples_to_remove),]
#data <- data[,!(colnames(data) %in% samples_to_remove)]


##################################################################
## annotation Table:
grch38
anno <- grch38 %>% select(ensgene,symbol,chr) %>% unique() %>% as.data.frame()
x <- c(1:22)
anno <- anno[grep(paste(x, collapse="|" ), anno$chr),]
head(anno)
data  <- data  %>% filter(rownames(data) %in% anno$ensgene)
#rownames(data) <- anno$ensgene
##################################################################


##################################################################
#### Manual conversion to R factor objects:

#### conc
cv$Treatment<- factor(cv$Treatment)
TreatmentLevels <- levels(cv$Treatment)
cat("#",TreatmentLevels,"\n")

#### Cell-line
cv$CellLine <- factor(cv$CellLine)
CellLineLevels <- levels(cv$CellLine)
cat("#",CellLineLevels,"\n")

#### time
cv$time <- factor(cv$Timepoint)
timeLevels <- levels(cv$time)
cat("#",timeLevels,"\n")

#### control
ConcentrationLevels <- c("H2O","EtOH")
cat("#",ConcentrationLevels,"\n")

#### Join
cv$Join <- factor(paste(cv$Treatment,cv$CellLine, sep="."))
JoinLevels <- levels(cv$Join)
cat("#",JoinLevels,"\n")



#### separating treatments and controls
TreatmentOnlyLevels <- TreatmentLevels[!(TreatmentLevels %in% ConcentrationLevels)]
cat("#",TreatmentOnlyLevels,"\n")
ControlOnlyLevels <- TreatmentLevels[(TreatmentLevels %in% ConcentrationLevels)]
cat("#",ControlOnlyLevels,"\n")

allColSamples <- paste(cv$Join, sep=".")
cat("#",allColSamples,"\n")

cv<- (cv[which(cv$Samples %in% colnames(data)),]) #So cv and data match

cv <- cv[order(cv$Samples),]
data <- data[,order(colnames(data))]
colnames(data)==cv$Samples #Check that files are in order


#removing the "Undeter"
data <- (data[,colnames(data)!="Undeter"])
colnames(data)==cv$Samples #Check that files are in order


#####################################################################
### Differential expression analysis

#filtering only 6 hours
data_6h <- data[colnames(data) %in% cv[cv$Timepoint=="6",]$Samples]
cv_6h <- cv[cv$Timepoint=="6",]

#### With filter

##### Pre-filtering genes
#Combine processed data into a DESeqDataSet and Remove genes with very low coverage/expression

###Merging the Controls(H2O and EtOH)

cv_for_OnlyTreatments <- cv_6h[!(cv_6h$Treatment %in% c("H2O","EtOH")),]
samples_for_OnlyTreatments <- cv_6h$Samples[!(cv_6h$Treatment %in% c("H2O","EtOH"))]
data_for_OnlyTreatments <- data_6h[,colnames(data_6h) %in% samples_for_OnlyTreatments]
#OnlyTreatments <- (data.frame((data_for_OnlyTreatments)))
OnlyTreatments <- data_for_OnlyTreatments

data_merged_6h <- OnlyTreatments

##H2O
cv_for_H2O <- cv_6h[cv_6h$Treatment %in% "H2O",]

Control <- data.frame()
CellLinesOrder <- list()
cv_merged <- data.frame()
for (cellLine in CellLineLevels) {
  samples_cellLine <- cv_for_H2O[(cv_for_H2O$CellLine %in% cellLine),]$Samples
  data_for_Control <- data_6h[,colnames(data_6h) %in% samples_cellLine]
  if (length(Control)==0) {
    Control <- rowSums(data_for_Control)
  }
  else{
    Control <- cbind(Control,data.frame(rowSums(data_for_Control)))
  }  
  CellLinesOrder <- append(CellLinesOrder,paste("H20",cellLine,sep ="."))
  cv_merged_1st <- cv_for_H2O[(cv_for_H2O$CellLine %in% cellLine),][1,] #keeping on the first one
  cv_merged <- rbind(cv_merged,cv_merged_1st)
}
colnames(Control) <- CellLinesOrder


data_merged_6h <- cbind(data_merged_6h,Control)
cv_merged_6h <- rbind(cv_for_OnlyTreatments, cv_merged)


##EtOH
cv_for_EtOH <- cv_6h[cv_6h$Treatment %in% "EtOH",]

Control <- data.frame()
CellLinesOrder <- list()
cv_merged <- data.frame()
for (cellLine in CellLineLevels) {
  samples_cellLine <- cv_for_EtOH[(cv_for_EtOH$CellLine %in% cellLine),]$Samples
  data_for_Control <- data_6h[,colnames(data_6h) %in% samples_cellLine]
  if (length(Control)==0) {
    Control <- rowSums(data_for_Control)
  }
  else{
    Control <- cbind(Control,data.frame(rowSums(data_for_Control)))
  }  
  CellLinesOrder <- append(CellLinesOrder,paste("EtOH",cellLine,sep ="."))
  cv_merged_1st <- cv_for_EtOH[(cv_for_EtOH$CellLine %in% cellLine),][1,] #keeping on the first one
  cv_merged <- rbind(cv_merged,cv_merged_1st)
}
colnames(Control) <- CellLinesOrder


data_merged_6h <- cbind(data_merged_6h,Control)
cv_merged_6h <- rbind(cv_merged_6h, cv_merged)



dim(data_merged_6h)

################################################################




###################
# Filtering genes #
###################

#Before running the DESeq, we performed different filters of genes by keeping genes that were relatively expressed in at least 50%  of the samples for each treatment and control contrast.

number_allgenes <- 0  #will be containing number of kept genes for each contrast
list_allgenes <- list() #will be containing list of kept genes for each contrast

#Removing the DEP_500nM treatment and corrisponding samples
TreatmentOnlyLevels <- TreatmentLevels[!(TreatmentLevels %in% c(ConcentrationLevels, "DEP_500nM"))]

for (i in 1:length(TreatmentOnlyLevels)){
  t <- TreatmentOnlyLevels[i]
  if (length(grep("BPA", t))!=0){ c <- "H2O"} else {c <- "EtOH"} #t=treatment and c= control
  #Collecting All samples for each contrast(treatment & control) --> s
  s <- cv_6h[cv_6h$Treatment %in% c(t, c),]$Samples 
  #filtering samples --> d
  d <- data_6h[,colnames(data_6h) %in% s]
  #filtering genes and make a new dataset
  keep <- rowSums(d>0) >= 0.5*ncol(d)
  d <- d[keep,]
  
  number_allgenes <- number_allgenes + nrow(d)
  print(paste0(t," and ", c, " ---> ", nrow(d), " genes"))
  list_allgenes <- append(list_allgenes, rownames(d))
}
#removing duplications in the list of collected genes
list_allgenes <- unique(list_allgenes) 
#filtering the whole dataset with the list of collected genes
data_filtered <- data_6h[rownames(data_6h) %in% list_allgenes,]


dim(data_filtered)

length(list_allgenes)

dim(data_6h)

cv_6h$Samples <- colnames(data_6h)

rownames(cv_6h) <- cv_6h$Samples 


###################
# The Deseq model #
###################

ddsFull <- DESeqDataSetFromMatrix(
    countData = round(data_6h),
    colData = cv_6h,
    design = ~ 0  + CellLine + fraction_mapped + Treatment ) #can add +cv$Ratio
#identical(colnames(ddsFull), cv_merged_6h$Samples)
ddsFull


## Fit the model on the whole plate

system.time(ddsFull <- DESeq(ddsFull,parallel=TRUE)) #kept the default fitType



res_BPA_1000nM <- results(ddsFull, contrast= c("Treatment","BPA_1000nM","H2O"))
res_BPA_100nM <- results(ddsFull, contrast= c("Treatment","BPA_100nM","H2O"))
res_BPA_500nM <- results(ddsFull, contrast= c("Treatment","BPA_500nM","H2O"))

res_DEP_1000nM <- results(ddsFull, contrast= c("Treatment","DEP_1000nM","EtOH"))
res_DEP_100nM <- results(ddsFull, contrast= c("Treatment","DEP_100nM","EtOH"))
res_DEP_1nM <- results(ddsFull, contrast= c("Treatment","DEP_1nM","EtOH"))
#res_DEP_500nM <- results(ddsFull, contrast= c("Treatment","DEP_500nM","EtOH"))

res_MBP_1000nM <- results(ddsFull, contrast= c("Treatment","MBP_1000nM","EtOH"))
res_MBP_100nM <- results(ddsFull, contrast= c("Treatment","MBP_100nM","EtOH"))
res_MBP_500nM <- results(ddsFull, contrast= c("Treatment","MBP_500nM","EtOH"))



#output summaries
cat("#","BPA_1000nM")
summary(res_BPA_1000nM)
cat("#","BPA_100nM")
summary(res_BPA_100nM)
cat("#","BPA_500nM")
summary(res_BPA_500nM) 
cat("#","DEP_1000nM")
summary(res_DEP_1000nM)
cat("#","DEP_100nM")
summary(res_DEP_100nM)
cat("#","DEP_1nM")
summary(res_DEP_1nM)
cat("#","DEP_500nM")
#summary(res_DEP_500nM)
cat("#","MBP_1000nM")
summary(res_MBP_1000nM)
cat("#","MBP_100nM")
summary(res_MBP_100nM)
cat("#","MBP_500nM")
summary(res_MBP_500nM) 


res_BPA_1000nM <- as.data.frame(res_BPA_1000nM)
res_BPA_100nM <- as.data.frame(res_BPA_100nM)
res_BPA_500nM <- as.data.frame(res_BPA_500nM)
res_DEP_1000nM <- as.data.frame(res_DEP_1000nM)
res_DEP_100nM <- as.data.frame(res_DEP_100nM)
res_DEP_1nM <- as.data.frame(res_DEP_1nM) 
#res_DEP_500nM <- as.data.frame(res_DEP_500nM)
res_MBP_1000nM <- as.data.frame(res_MBP_1000nM) 
res_MBP_100nM <- as.data.frame(res_MBP_100nM) 
res_MBP_500nM <- as.data.frame(res_MBP_500nM) 

res_000nM_sigBPA_1 <- res_BPA_1000nM  %>% filter(abs(log2FoldChange) >0.25 & padj <0.1)
res_BPA_1000nM_sig <- res_BPA_1000nM %>% filter(abs(log2FoldChange) >0.25 & padj <0.1)
res_BPA_100nM_sig <- res_BPA_100nM %>% filter(abs(log2FoldChange) >0.25 & padj <0.1)
res_BPA_500nM_sig <- res_BPA_500nM %>% filter(abs(log2FoldChange) >0.25 & padj <0.1)
res_DEP_1000nM_sig <- res_DEP_1000nM %>% filter(abs(log2FoldChange) >0.25 & padj <0.1)
res_DEP_100nM_sig <- res_DEP_100nM %>% filter(abs(log2FoldChange) >0.25 & padj <0.1)
res_DEP_1nM_sig <- res_DEP_1nM %>% filter(abs(log2FoldChange) >0.25 & padj <0.1)
#res_DEP_500nM_sig <- res_DEP_500nM %>% filter(abs(log2FoldChange) >0.25 & padj <0.1)
res_MBP_1000nM_sig <- res_MBP_1000nM %>% filter(abs(log2FoldChange) >0.25 & padj <0.1)
res_MBP_100nM_sig <- res_MBP_100nM %>% filter(abs(log2FoldChange) >0.25 & padj <0.1)
res_MBP_500nM_sig <- res_MBP_500nM %>% filter(abs(log2FoldChange) >0.25 & padj <0.1)



BPA_1000nM <- res_BPA_1000nM
BPA_100nM <- res_BPA_100nM
BPA_500nM <- res_BPA_500nM
DEP_1000nM <- res_DEP_1000nM
DEP_100nM <- res_DEP_100nM
DEP_1nM <- res_MBP_1000nM
DEP_500nM <- res_DEP_500nM
MBP_1000nM <- res_MBP_1000nM
MBP_100nM <- res_MBP_100nM
MBP_500nM <- res_MBP_500nM


my.pvalues_BPA_1000nM <- BPA_1000nM$pvalue
exp.pvalues_BPA_1000nM<-(rank(my.pvalues_BPA_1000nM, ties.method="first")+.5)/(length(my.pvalues_BPA_1000nM)+1)
my.pvalues_BPA_100nM <- BPA_100nM$pvalue
exp.pvalues_BPA_100nM<-(rank(my.pvalues_BPA_100nM, ties.method="first")+.5)/(length(my.pvalues_BPA_100nM)+1)
my.pvalues_BPA_500nM <- BPA_500nM$pvalue
exp.pvalues_BPA_500nM<-(rank(my.pvalues_BPA_500nM, ties.method="first")+.5)/(length(my.pvalues_BPA_500nM)+1)

my.pvalues_DEP_1000nM <- DEP_1000nM$pvalue
exp.pvalues_DEP_1000nM <-(rank(my.pvalues_DEP_1000nM, ties.method="first")+.5)/(length(my.pvalues_DEP_1000nM)+1)
my.pvalues_DEP_100nM <- DEP_100nM$pvalue
exp.pvalues_DEP_100nM <-(rank(my.pvalues_DEP_100nM, ties.method="first")+.5)/(length(my.pvalues_DEP_100nM)+1)
my.pvalues_DEP_1nM <- DEP_1nM$pvalue
exp.pvalues_DEP_1nM <-(rank(my.pvalues_DEP_1nM, ties.method="first")+.5)/(length(my.pvalues_DEP_1nM)+1)
#my.pvalues_DEP_500nM <- DEP_500nM$pvalue
#exp.pvalues_DEP_500nM <-(rank(my.pvalues_DEP_500nM, ties.method="first")+.5)/(length(my.pvalues_DEP_500nM)+1)

my.pvalues_MBP_1000nM <- MBP_1000nM$pvalue
exp.pvalues_MBP_1000nM <-(rank(my.pvalues_MBP_1000nM, ties.method="first")+.5)/(length(my.pvalues_MBP_1000nM)+1)
my.pvalues_MBP_100nM <- MBP_100nM$pvalue
exp.pvalues_MBP_100nM <-(rank(my.pvalues_MBP_100nM, ties.method="first")+.5)/(length(my.pvalues_MBP_100nM)+1)
my.pvalues_MBP_500nM <- MBP_500nM$pvalue
exp.pvalues_MBP_500nM <-(rank(my.pvalues_MBP_500nM, ties.method="first")+.5)/(length(my.pvalues_MBP_500nM)+1)


BPA_1000nM$group <- "BPA_1000nM vs H2O"
BPA_100nM$group <- "BPA_100nM vs H2O"
BPA_500nM$group <- "BPA_500nM vs H2O"
DEP_1000nM$group <- "DEP_1000nM vs EtOH"
DEP_100nM$group <- "DEP_100nM vs EtOH"
DEP_1nM$group <- "DEP_1nM vs EtOH"
#DEP_500nM$group <- "DEP_500nM vs EtOH"
MBP_1000nM$group <- "MBP_1000nM vs EtOH"
MBP_100nM$group <- "MBP_100nM vs EtOH"
MBP_500nM$group <- "MBP_500nM vs EtOH"

BPA_1000nM$exp.pvalues_BPA_1000nM <- exp.pvalues_BPA_1000nM
BPA_100nM$exp.pvalues_BPA_100nM <- exp.pvalues_BPA_100nM
BPA_500nM$exp.pvalues_BPA_500nM <- exp.pvalues_BPA_500nM
DEP_1000nM$exp.pvalues_DEP_1000nM <- exp.pvalues_DEP_1000nM
DEP_100nM$exp.pvalues_DEP_100nM <- exp.pvalues_DEP_100nM
DEP_1nM$exp.pvalues_DEP_1nM <- exp.pvalues_DEP_1nM
#DEP_500nM$exp.pvalues_DEP_500nM <- exp.pvalues_DEP_500nM
MBP_1000nM$exp.pvalues_MBP_1000nM <- exp.pvalues_MBP_1000nM
MBP_100nM$exp.pvalues_MBP_100nM <- exp.pvalues_MBP_100nM
MBP_500nM$exp.pvalues_MBP_500nM <- exp.pvalues_MBP_500nM



names(BPA_1000nM)[8]<- "exp.pvalues"
names(BPA_100nM)[8]<- "exp.pvalues"
names(BPA_500nM)[8]<- "exp.pvalues"
names(DEP_1000nM)[8]<- "exp.pvalues"
names(DEP_100nM)[8]<- "exp.pvalues"
names(DEP_1nM)[8]<- "exp.pvalues"
#names(DEP_500nM)[8]<- "exp.pvalues"
names(MBP_1000nM)[8]<- "exp.pvalues"
names(MBP_100nM)[8]<- "exp.pvalues"
names(MBP_500nM)[8]<- "exp.pvalues"

##### QQplots and PCA

res <- rbind(BPA_1000nM,BPA_100nM,BPA_500nM,DEP_1000nM,DEP_100nM,DEP_1nM,MBP_1000nM,MBP_100nM,MBP_500nM)
##names(res)[7] <- "symbol"

names(res)[7] <- "Treatment"

pdf("QQplot_Treatment_deseq_6h.pdf")
p1<- ggplot(res, aes(-log10(exp.pvalues), -log10(pvalue), colour=Treatment))+
  geom_point()+
  theme_classic()+
  geom_abline(intercept=0, slope=1)+
  labs(title="QQplot_Treatment_deseq_6h")+
  theme(legend.title= element_blank(), plot.title = element_text(hjust = 0.5), axis.title.x = element_text(size = rel(1)), axis.title.y = element_text(size = rel(1)), legend.text = element_text(size = rel(1)), title = element_text(size = rel(1)))
#labs(title= "DEGs induced by heat deactivated microbiome or antibiotic treated cells")
print(p1)
dev.off()

vsd <- vst(ddsFull, blind=FALSE)
pdf("PCA_6h.pdf")
pcaData <- plotPCA(vsd, intgroup=c("Treatment", "CellLine"), returnData=T)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pca <- ggplot(pcaData, aes(PC1, PC2, color=Treatment, shape=CellLine)) +
  geom_point(size=3) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  labs(title="PCA_6h")+
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  coord_fixed()
print(pca) 
dev.off()

##### Heatmaps

library(pheatmap)
library(gplots)

rld <- rlog(ddsFull)
head(assay(rld[1:5,1:5]))

sampleDists <- dist( t( assay(rld) ) )
#sampleDists
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$Treatment, rld$CellLine, sep="-" )
df <- as.data.frame(colData(rld)[,c("Treatment","CellLine")])
#if (nrow(df) > 100) stop("Too many rows for heatmap, who can read?!")
fontsize = 10 - nrow(df) / 15
pdf("htmap_anno_6h.pdf")
f <- pheatmap(sampleDistMatrix, annotation =df, main= "Timepoint_6h (sampleDists)",fontsize_row=fontsize, fontsize_col=fontsize)
print(f)
dev.off()

library("genefilter")
topVarGenes <- head(order(-rowVars(assay(rld))),20)

mat <- assay(rld)[ topVarGenes, ]
rownames(mat) <- anno$symbol[anno$ensgene %in% rownames(mat)]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld)[,c("Treatment","CellLine")])

fontsize = 10 - nrow(df) / 15
pdf("gene_clustering_6h.pdf")
ph <- pheatmap(mat, annotation_col=df, main="gene clustering (6h)",fontsize_col=fontsize, fontsize_row=fontsize, annotation_colors= )
print(ph)
dev.off()








# Load libraries
library(DESeq2)
library(reshape)
library(ggplot2)
library(ggrepel)
library(DEGreport)
library(RColorBrewer)
library(pheatmap)


#BPA


sigOE_ordered <- res_BPA_1000nM_sig[order(res_BPA_1000nM_sig$padj), ]
rownames(sigOE_ordered) <- anno$symbol[anno$ensgene %in% rownames(sigOE_ordered)]


top20_sigOE_genes <- rownames(sigOE_ordered[1:20, ])



ddsFull_BPA <- ddsFull[,ddsFull$treatment %in% c("BPA", "H2O")]
normalized_counts <- counts(ddsFull_BPA, normalized=T)

rownames(normalized_counts) <- anno$symbol[anno$ensgene %in% rownames(normalized_counts)]
top20_sigOE_norm <- normalized_counts[top20_sigOE_genes, ]

## use melt to modify the format of the data frame
melted_top20_sigOE <- data.frame(melt(top20_sigOE_norm))

## check the column header in the "melted" data frame
View(melted_top20_sigOE)

## add column names that make sense
colnames(melted_top20_sigOE) <- c("gene", "Samples", "normalized_counts")
cv_merged_BPA <- cv_merged_6h[cv_merged_6h$treatment %in% c("BPA","H2O"),]
dt <- merge(melted_top20_sigOE, cv_merged_BPA, by = "Samples")




pdf("DEG_BPA_dosages_TRAM2.pdf")
p <-ggplot((dt[dt$gene %in% "TRAM2",]),aes(x = Treatment, y = normalized_counts, group = CellLine, color = CellLine)) +
  geom_point(aes(x = factor(Treatment, level=c('H2O', 'BPA_100nM', 'BPA_500nM', 'BPA_1000nM')),
                 y = normalized_counts, color = CellLine )) +
  geom_line(alpha =0.7) +
  xlab("Treatments") +
  ylab("Normalized Counts") +
  ggtitle("gene: TRAM2") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title=element_text(hjust=0.5))
print(p)
dev.off()

pdf("DEG_BPA_dosages_CTTN.pdf")
p <-ggplot((dt[dt$gene %in% "CTTN",]),aes(x = Treatment, y = normalized_counts, group = CellLine, color = CellLine)) +
  geom_point(aes(x = factor(Treatment, level=c('H2O', 'BPA_100nM', 'BPA_500nM', 'BPA_1000nM')),
                 y = normalized_counts, color = CellLine )) +
  geom_line(alpha =0.7) +
  xlab("Treatments") +
  ylab("Normalized Counts") +
  ggtitle("gene: CTTN") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title=element_text(hjust=0.5))
print(p)
dev.off()

pdf("DEG_BPA_dosages_DYRK4.pdf")
p <-ggplot((dt[dt$gene %in% "DYRK4",]),aes(x = Treatment, y = normalized_counts, group = CellLine, color = CellLine)) +
  geom_point(aes(x = factor(Treatment, level=c('H2O', 'BPA_100nM', 'BPA_500nM', 'BPA_1000nM')),
                 y = normalized_counts, color = CellLine )) +
  geom_line(alpha =0.7) +
  xlab("Treatments") +
  ylab("Normalized Counts") +
  ggtitle("gene: DYRK4") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title=element_text(hjust=0.5))
print(p)
dev.off()

pdf("DEG_BPA_dosages_RBMS2.pdf")
p <-ggplot((dt[dt$gene %in% "RBMS2",]),aes(x = Treatment, y = normalized_counts, group = CellLine, color = CellLine)) +
  geom_point(aes(x = factor(Treatment, level=c('H2O', 'BPA_100nM', 'BPA_500nM', 'BPA_1000nM')),
                 y = normalized_counts, color = CellLine )) +
  geom_line(alpha =0.7) +
  xlab("Treatments") +
  ylab("Normalized Counts") +
  ggtitle("gene: RBMS2") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title=element_text(hjust=0.5))
print(p)
dev.off()

#DEP

sigOE_ordered <- res_DEP_1000nM_sig[order(res_DEP_1000nM_sig$padj), ]
rownames(sigOE_ordered) <- anno$symbol[anno$ensgene %in% rownames(sigOE_ordered)]

top20_sigOE_genes <- rownames(sigOE_ordered[1:20, ])



ddsFull_DEP <- ddsFull[,ddsFull$treatment %in% c("DEP", "EtOH")]
normalized_counts <- counts(ddsFull_DEP, normalized=T)

rownames(normalized_counts) <- anno$symbol[anno$ensgene %in% rownames(normalized_counts)]
top20_sigOE_norm <- normalized_counts[top20_sigOE_genes, ]

## use melt to modify the format of the data frame
melted_top20_sigOE <- data.frame(melt(top20_sigOE_norm))

## check the column header in the "melted" data frame
View(melted_top20_sigOE)

## add column names that make sense
colnames(melted_top20_sigOE) <- c("gene", "Samples", "normalized_counts")
cv_merged_DEP <- cv_merged_6h[cv_merged_6h$treatment %in% c("DEP","EtOH"),]
dt <- merge(melted_top20_sigOE, cv_merged_DEP, by = "Samples")

pdf("DEG_DEP_dosages_MAP4K1.pdf")
p <-ggplot((dt[dt$gene %in% "MAP4K1",]),aes(x = Treatment, y = normalized_counts, group = CellLine, color = CellLine)) +
  geom_point(aes(x = factor(Treatment, level=c('EtOH', 'DEP_1nM', 'DEP_100nM', 'DEP_1000nM')),
                 y = normalized_counts, color = CellLine )) +
  geom_line(alpha =0.7) +
  xlab("Treatments") +
  ylab("Normalized Counts") +
  ggtitle("gene: MAP4K1") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title=element_text(hjust=0.5))
print(p)
dev.off()


pdf("DEG_DEP_dosages_APBA1.pdf")
p <-ggplot((dt[dt$gene %in% "APBA1",]),aes(x = Treatment, y = normalized_counts, group = CellLine, color = CellLine)) +
  geom_point(aes(x = factor(Treatment, level=c('EtOH', 'DEP_1nM', 'DEP_100nM', 'DEP_1000nM')),
                 y = normalized_counts, color = CellLine )) +
  geom_line(alpha =0.7) +
  xlab("Treatments") +
  ylab("Normalized Counts") +
  ggtitle("gene: APBA1") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title=element_text(hjust=0.5))
print(p)
dev.off()


pdf("DEG_DEP_dosages_OAS2.pdf")
p <-ggplot((dt[dt$gene %in% "OAS2",]),aes(x = Treatment, y = normalized_counts, group = CellLine, color = CellLine)) +
  geom_point(aes(x = factor(Treatment, level=c('EtOH', 'DEP_1nM', 'DEP_100nM', 'DEP_1000nM')),
                 y = normalized_counts, color = CellLine )) +
  geom_line(alpha =0.7) +
  xlab("Treatments") +
  ylab("Normalized Counts") +
  ggtitle("gene: OAS2") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title=element_text(hjust=0.5))
print(p)
dev.off()



pdf("DEG_DEP_dosages_MSX2.pdf")
p <-ggplot((dt[dt$gene %in% "MSX2",]),aes(x = Treatment, y = normalized_counts, group = CellLine, color = CellLine)) +
  geom_point(aes(x = factor(Treatment, level=c('EtOH', 'DEP_1nM', 'DEP_100nM', 'DEP_1000nM')),
                 y = normalized_counts, color = CellLine )) +
  geom_line(alpha =0.7) +
  xlab("Treatments") +
  ylab("Normalized Counts") +
  ggtitle("gene: MSX2") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title=element_text(hjust=0.5))
print(p)
dev.off()


#MBP

sigOE_ordered <- res_MBP_100nM_sig[order(res_MBP_100nM_sig$padj), ]
rownames(sigOE_ordered) <- anno$symbol[anno$ensgene %in% rownames(sigOE_ordered)]

top20_sigOE_genes <- rownames(sigOE_ordered[1:20, ])



ddsFull_MBP <- ddsFull[,ddsFull$treatment %in% c("MBP", "EtOH")]
normalized_counts <- counts(ddsFull_MBP, normalized=T)

rownames(normalized_counts) <- anno$symbol[anno$ensgene %in% rownames(normalized_counts)]
top20_sigOE_norm <- normalized_counts[top20_sigOE_genes, ]

## use melt to modify the format of the data frame
melted_top20_sigOE <- data.frame(melt(top20_sigOE_norm))

## check the column header in the "melted" data frame
View(melted_top20_sigOE)

## add column names that make sense
colnames(melted_top20_sigOE) <- c("gene", "Samples", "normalized_counts")
cv_merged_MBP <- cv_merged_6h[cv_merged_6h$treatment %in% c("MBP","EtOH"),]
dt <- merge(melted_top20_sigOE, cv_merged_MBP, by = "Samples")

pdf("DEG_MBP_dosages_RPS20.pdf")
p <-ggplot((dt[dt$gene %in% "RPS20",]),aes(x = Treatment, y = normalized_counts, group = CellLine, color = CellLine)) +
  geom_point(aes(x = factor(Treatment, level=c('EtOH', 'MBP_100nM', 'MBP_500nM', 'MBP_1000nM')),
                 y = normalized_counts, color = CellLine )) +
  geom_line(alpha =0.7) +
  xlab("Treatments") +
  ylab("Normalized Counts") +
  ggtitle("gene: RPS20") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title=element_text(hjust=0.5))
print(p)
dev.off()


pdf("DEG_MBP_dosages_TMSB10.pdf")
p <-ggplot((dt[dt$gene %in% "TMSB10",]),aes(x = Treatment, y = normalized_counts, group = CellLine, color = CellLine)) +
  geom_point(aes(x = factor(Treatment, level=c('EtOH', 'MBP_100nM', 'MBP_500nM', 'MBP_1000nM')),
                 y = normalized_counts, color = CellLine )) +
  geom_line(alpha =0.7) +
  xlab("Treatments") +
  ylab("Normalized Counts") +
  ggtitle("gene: TMSB10") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title=element_text(hjust=0.5))
print(p)
dev.off()



pdf("DEG_MBP_dosages_CDH10.pdf")
p <-ggplot((dt[dt$gene %in% "CDH10",]),aes(x = Treatment, y = normalized_counts, group = CellLine, color = CellLine)) +
  geom_point(aes(x = factor(Treatment, level=c('EtOH', 'MBP_100nM', 'MBP_500nM', 'MBP_1000nM')),
                 y = normalized_counts, color = CellLine )) +
  geom_line(alpha =0.7) +
  xlab("Treatments") +
  ylab("Normalized Counts") +
  ggtitle("gene: CDH10") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title=element_text(hjust=0.5))
print(p)
dev.off()







pdf("DEG1_BPA_dosages.pdf")
p <- ggplot(data = (dt[dt$gene %in% "DYRK4",]), 
  mapping = aes(x = Treatment, y = normalized_counts, fill = CellLine)) +
  xlab("Treatments with dosage") +
  ylab("Normalized Counts") +
  ggtitle("gene: DYRK4") +
  geom_bar(stat = "identity", width = 1, position = position_dodge(0.8)) +
  geom_text(aes(label = paste(sprintf("%s", Samples)),
                hjust= 0.5,vjust=-2,group=CellLine),
            size = 2, position = position_dodge(width=.9)) +
  #theme_bw() +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title=element_text(hjust=0.5))
print(p)
dev.off()





pdf("DEG1_BPA_Samples.pdf")
p <- ggplot(data = (dt[dt$gene %in% "DYRK4",]), 
            mapping = aes(x = Samples, y = normalized_counts, fill = CellLine)) +
  xlab("Samples") +
  ylab("Normalized Counts") +
  ggtitle("gene: DYRK4") +
  geom_bar(stat = "identity", width = 0.5, position = position_dodge(0.8)) +
  geom_text(aes(label = paste(sprintf("%s", Treatment)),
                hjust= 0.5,vjust=-2,group=CellLine),
            size = 2, position = position_dodge(width=.9)) +
  #theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title=element_text(hjust=0.5))
print(p)
dev.off()


















##### volcano plots

library(ggpubr)
theme_set(theme_pubr())

lfc <- 0.25

#res_BPA_1000nM
de <- as.data.frame(res_BPA_1000nM)
de$diffexpressed <- "NO"
de$diffexpressed[de$log2FoldChange > lfc & de$padj < 0.1] <- "UP"
de$diffexpressed[de$log2FoldChange < lfc & de$padj < 0.1] <- "DOWN"

de$gene_symbol <- anno[anno$ensgene %in% rownames(res_BPA_1000nM),]$symbol
if (length(de$gene_symbol[de$diffexpressed != "NO"] )==0 ) {de["delabel"] = NA_character_} else {de$delabel[de$diffexpressed != "NO"] <- de$gene_symbol[de$diffexpressed != "NO"]}

p1 <- ggplot(data=de, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) + 
  geom_point(size = 0.7) + theme_minimal() + #geom_text_repel(size = 2.3,box.padding = 0.3, max.overlaps = Inf) + 
  labs(title = "BPA_1000nM_6h vs H2O") + 
  theme(plot.title.position = "plot", plot.title = element_text(hjust = 0.1),legend.title = element_text(size=6), legend.text = element_text(size=6))

#res_BPA_100nM
de <- as.data.frame(res_BPA_100nM)
de$diffexpressed <- "NO"
de$diffexpressed[de$log2FoldChange > lfc & de$padj < 0.1] <- "UP"
de$diffexpressed[de$log2FoldChange < -lfc & de$padj < 0.1] <- "DOWN"

de$gene_symbol <- anno[anno$ensgene %in% rownames(res_BPA_100nM),]$symbol
if (length(de$gene_symbol[de$diffexpressed != "NO"] )==0 ) {de["delabel"] = NA_character_} else {de$delabel[de$diffexpressed != "NO"] <- de$gene_symbol[de$diffexpressed != "NO"]}

p2 <- ggplot(data=de, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) + 
  geom_point(size = 0.7) + theme_minimal() + #geom_text_repel(size = 2.3,box.padding = 0.3, max.overlaps = Inf) + 
  labs(title = "BPA_100nM vs H2O") + 
  theme(plot.title.position = "plot", plot.title = element_text(hjust = 0.1),legend.title = element_text(size=6), legend.text = element_text(size=6))

#res_BPA_500nM
de <- as.data.frame(res_BPA_500nM)
de$diffexpressed <- "NO"
de$diffexpressed[de$log2FoldChange > lfc & de$padj < 0.1] <- "UP"
de$diffexpressed[de$log2FoldChange < -lfc & de$padj < 0.1] <- "DOWN"

de$gene_symbol <- anno[anno$ensgene %in% rownames(res_BPA_500nM),]$symbol
if (length(de$gene_symbol[de$diffexpressed != "NO"] )==0 ) {de["delabel"] = NA_character_} else {de$delabel[de$diffexpressed != "NO"] <- de$gene_symbol[de$diffexpressed != "NO"]}

p3 <- ggplot(data=de, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) + 
  geom_point(size = 0.7) + theme_minimal() + #geom_text_repel(size = 2.3,box.padding = 0.3, max.overlaps = Inf) + 
  labs(title = "BPA_500nM_6h vs H2O") + 
  theme(plot.title.position = "plot", plot.title = element_text(hjust = 0.1),legend.title = element_text(size=6), legend.text = element_text(size=6))

png("vplot_BPA_6h.png")
l = ggarrange(p1,p2,p3)
print(l)
dev.off()

########################


#res_DEP_1000nM
de <- as.data.frame(res_DEP_1000nM)
de$diffexpressed <- "NO"
de$diffexpressed[de$log2FoldChange > lfc & de$padj < 0.1] <- "UP"
de$diffexpressed[de$log2FoldChange < -lfc & de$padj < 0.1] <- "DOWN"

de$gene_symbol <- anno[anno$ensgene %in% rownames(res_DEP_1000nM),]$symbol
if (length(de$gene_symbol[de$diffexpressed != "NO"] )==0 ) {de["delabel"] = NA_character_} else {de$delabel[de$diffexpressed != "NO"] <- de$gene_symbol[de$diffexpressed != "NO"]}

#pdf("vplot_DEP_1000nM_6h.pdf");
p1 <- ggplot(data=de, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) + 
  geom_point(size = 0.7) + theme_minimal() + geom_text_repel(size = 2.3,box.padding = 0.3, max.overlaps = Inf) + 
  labs(title = "DEP_1000nM_6h vs EtOH") + 
  theme(plot.title.position = "plot", plot.title = element_text(hjust = 0.1),legend.title = element_text(size=5), legend.text = element_text(size=4))

#res_DEP_100nM
de <- as.data.frame(res_DEP_100nM)
de$diffexpressed <- "NO"
de$diffexpressed[de$log2FoldChange > lfc & de$padj < 0.1] <- "UP"
de$diffexpressed[de$log2FoldChange < -lfc & de$padj < 0.1] <- "DOWN"

de$gene_symbol <- anno[anno$ensgene %in% rownames(res_DEP_100nM),]$symbol
if (length(de$gene_symbol[de$diffexpressed != "NO"] )==0 ) {de["delabel"] = NA_character_} else {de$delabel[de$diffexpressed != "NO"] <- de$gene_symbol[de$diffexpressed != "NO"]}

p2 <- ggplot(data=de, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) + 
  geom_point(size = 0.7) + theme_minimal() + geom_text_repel(size = 2.3,box.padding = 0.3, max.overlaps = Inf) + 
  labs(title = "DEP_100nM_6h vs EtOH") + 
  theme(plot.title.position = "plot", plot.title = element_text(hjust = 0.1),legend.title = element_text(size=5), legend.text = element_text(size=4))

#res_DEP_1nM
de <- as.data.frame(res_DEP_1nM)
de$diffexpressed <- "NO"
de$diffexpressed[de$log2FoldChange > lfc & de$padj < 0.1] <- "UP"
de$diffexpressed[de$log2FoldChange < -lfc & de$padj < 0.1] <- "DOWN"

de$gene_symbol <- anno[anno$ensgene %in% rownames(res_DEP_1nM),]$symbol
if (length(de$gene_symbol[de$diffexpressed != "NO"] )==0 ) {de["delabel"] = NA_character_} else {de$delabel[de$diffexpressed != "NO"] <- de$gene_symbol[de$diffexpressed != "NO"]}

p3 <- ggplot(data=de, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) + 
  geom_point(size = 0.7) + theme_minimal() + geom_text_repel(size = 2.3,box.padding = 0.3, max.overlaps = Inf) + 
  labs(title = "DEP_1nM_6h vs EtOH") + 
  theme(plot.title.position = "plot", plot.title = element_text(hjust = 0.1),legend.title = element_text(size=5), legend.text = element_text(size=4))

#res_DEP_500nM
de <- as.data.frame(res_DEP_500nM)
de$diffexpressed <- "NO"
de$diffexpressed[de$log2FoldChange > lfc & de$padj < 0.1] <- "UP"
de$diffexpressed[de$log2FoldChange < -lfc & de$padj < 0.1] <- "DOWN"

de$gene_symbol <- anno[anno$ensgene %in% rownames(res_DEP_500nM),]$symbol
if (length(de$gene_symbol[de$diffexpressed != "NO"] )==0 ) {de["delabel"] = NA_character_} else {de$delabel[de$diffexpressed != "NO"] <- de$gene_symbol[de$diffexpressed != "NO"]}

p4 <- ggplot(data=de, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) + 
  geom_point(size = 0.7) + theme_minimal() + geom_text_repel(size = 2.3,box.padding = 0.3, max.overlaps = Inf) + 
  labs(title = "DEP_500nM_6h vs EtOH") + 
  theme(plot.title.position = "plot", plot.title = element_text(hjust = 0.1),legend.title = element_text(size=5), legend.text = element_text(size=4))

pdf("vplot_DEP_6h.pdf")
l = ggarrange(p1,p2,p3,p4)
print(l)
dev.off()

############################


#res_MBP_1000nM
de <- as.data.frame(res_MBP_1000nM)
de$diffexpressed <- "NO"
de$diffexpressed[de$log2FoldChange > lfc & de$padj < 0.1] <- "UP"
de$diffexpressed[de$log2FoldChange < -lfc & de$padj < 0.1] <- "DOWN"

de$gene_symbol <- anno[anno$ensgene %in% rownames(res_MBP_1000nM),]$symbol
if (length(de$gene_symbol[de$diffexpressed != "NO"] )==0 ) {de["delabel"] = NA_character_} else {de$delabel[de$diffexpressed != "NO"] <- de$gene_symbol[de$diffexpressed != "NO"]}

p1 <- ggplot(data=de, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) + 
  geom_point(size = 0.7) + theme_minimal() + geom_text_repel(size = 2.3,box.padding = 0.3, max.overlaps = Inf) + 
  labs(title = "MBP_1000nM_6h vs EtOH") + 
  theme(plot.title.position = "plot", plot.title = element_text(hjust = 0.1),legend.title = element_text(size=5), legend.text = element_text(size=4))

#res_MBP_100nM
de <- as.data.frame(res_MBP_100nM)
de$diffexpressed <- "NO"
de$diffexpressed[de$log2FoldChange > lfc & de$padj < 0.1] <- "UP"
de$diffexpressed[de$log2FoldChange < -lfc & de$padj < 0.1] <- "DOWN"

de$gene_symbol <- anno[anno$ensgene %in% rownames(res_MBP_100nM),]$symbol
if (length(de$gene_symbol[de$diffexpressed != "NO"] )==0 ) {de["delabel"] = NA_character_} else {de$delabel[de$diffexpressed != "NO"] <- de$gene_symbol[de$diffexpressed != "NO"]}

p2 <- ggplot(data=de, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) + 
  geom_point(size = 0.7) + theme_minimal() + geom_text_repel(size = 2.3,box.padding = 0.3, max.overlaps = Inf) + 
  labs(title = "MBP_100nM_6h vs EtOH") + 
  theme(plot.title.position = "plot", plot.title = element_text(hjust = 0.1),legend.title = element_text(size=5), legend.text = element_text(size=4))

#res_MBP_500nM
de <- as.data.frame(res_MBP_500nM)
de$diffexpressed <- "NO"
de$diffexpressed[de$log2FoldChange > lfc & de$padj < 0.1] <- "UP"
de$diffexpressed[de$log2FoldChange < -lfc & de$padj < 0.1] <- "DOWN"

de$gene_symbol <- anno[anno$ensgene %in% rownames(res_MBP_500nM),]$symbol
if (length(de$gene_symbol[de$diffexpressed != "NO"] )==0 ) {de["delabel"] = NA_character_} else {de$delabel[de$diffexpressed != "NO"] <- de$gene_symbol[de$diffexpressed != "NO"]}

p3 <- ggplot(data=de, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) +
  geom_point(size = 0.7) + theme_minimal() + geom_text_repel(size = 2.3,box.padding = 0.3, max.overlaps = Inf) + 
  labs(title = "MBP_500nM_6h vs EtOH") + 
  theme(plot.title.position = "plot", plot.title = element_text(hjust = 0.1),legend.title = element_text(size=5), legend.text = element_text(size=4))

pdf("vplot_MBP_6h.pdf")
l = ggarrange(p1,p2,p3)
print(l)
dev.off()








library(ggplot2)
library(gplots)
library(dplyr)
library(ggpubr)
library(ggpubr)
theme_set(theme_pubr())

plots_zscore= "plots_zscore_6h"

BPA_1000nM$zscore <- BPA_1000nM$log2FoldChange/BPA_1000nM$lfcSE
BPA_100nM$zscore <- BPA_100nM$log2FoldChange/BPA_100nM$lfcSE
BPA_500nM$zscore <- BPA_500nM$log2FoldChange/BPA_500nM$lfcSE
DEP_1000nM$zscore <- DEP_1000nM$log2FoldChange/DEP_1000nM$lfcSE
DEP_100nM$zscore <- DEP_100nM$log2FoldChange/DEP_100nM$lfcSE
DEP_1nM$zscore <- DEP_1nM$log2FoldChange/DEP_1nM$lfcSE
#DEP_500nM$zscore <- DEP_500nM$log2FoldChange/DEP_500nM$lfcSE
MBP_1000nM$zscore <- MBP_1000nM$log2FoldChange/MBP_1000nM$lfcSE
MBP_100nM$zscore <- MBP_100nM$log2FoldChange/MBP_100nM$lfcSE
MBP_500nM$zscore <- MBP_500nM$log2FoldChange/MBP_500nM$lfcSE



library(ggplot2)
library(gplots)
library(dplyr)
library(ggpubr)
library(ggpubr)
theme_set(theme_pubr())

plots_zscore= "plots_zscore_6h_BPA"

BPA_1000nM$zscore <- BPA_1000nM$log2FoldChange/BPA_1000nM$lfcSE
BPA_100nM$zscore <- BPA_100nM$log2FoldChange/BPA_100nM$lfcSE
BPA_500nM$zscore <- BPA_500nM$log2FoldChange/BPA_500nM$lfcSE
DEP_1000nM$zscore <- DEP_1000nM$log2FoldChange/DEP_1000nM$lfcSE
DEP_100nM$zscore <- DEP_100nM$log2FoldChange/DEP_100nM$lfcSE
DEP_1nM$zscore <- DEP_1nM$log2FoldChange/DEP_1nM$lfcSE
#DEP_500nM$zscore <- DEP_500nM$log2FoldChange/DEP_500nM$lfcSE
MBP_1000nM$zscore <- MBP_1000nM$log2FoldChange/MBP_1000nM$lfcSE
MBP_100nM$zscore <- MBP_100nM$log2FoldChange/MBP_100nM$lfcSE
MBP_500nM$zscore <- MBP_500nM$log2FoldChange/MBP_500nM$lfcSE




for (i in TreatmentOnlyLevels) {
  if (i!="DEP_500nM") {
    dfi <- get(i)
    dfi$gene_ID <- rownames(dfi)
    assign(as.character(i), dfi)
    pref_treat <- substr(i,1,3)
    treats <- TreatmentLevels[grep(pref_treat,TreatmentLevels)]
    for (j in treats) {
      if (i!=j) {
        if (j!="DEP_500nM") {
          dfj <- get(j)
          dfj$gene_ID <- rownames(dfj)
          assign(as.character(j), dfj)
          fname=paste( plots_zscore,'/',"zscore_", i, "_vs_", j, "_6h.png", sep="")
          png(fname)
          res = inner_join(dfi,dfj,by="gene_ID")
          #res$diffexpressed <- ifelse(res$padj.x < 0.1,
          #                            ifelse(res$pvalue.y < 0.1, yes = "BothSignificant", no = "SignificantForX"),
          #                            ifelse(res$pvalue.y < 0.1, yes = "SignificantForY", no = "NonSignificant"))
          
          res$diffexpressed="NS_Both"
          res$diffexpressed[res$padj.x<0.1 &res$padj.y>0.1] = "Diff.X"
          res$diffexpressed[res$padj.x>0.1 &res$padj.y<0.1] = "Diff.Y"
          res$diffexpressed[res$padj.x<0.1 &res$padj.y<0.1] = "Both_S"
          
          
          p <- ggplot(res,aes(x=zscore.x, zscore.y,col=diffexpressed)) + geom_point(alpha = 0.7) +
            geom_count()+
            geom_abline(intercept=0,slope=1) + labs(title = paste(i, " vs ", j, " 6h")) +
            #geom_smooth(method='lm') + stat_cor(method = "spearman") +
            xlab(paste("zscore for", i)) + ylab(paste("zscore for ",j)) +
            
            theme(plot.title.position = "plot",  plot.title = element_text(hjust = 0.5, size = 9),
                  legend.title = element_text(size=10), legend.text = element_text(size=7),
                  axis.text = element_text(size=8),axis.title = element_text(size = 10))
          print(p)
          dev.off()
        }
        
      }
    }
    
    
  }
  
}








