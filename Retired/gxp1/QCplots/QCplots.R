library(ggplot2)


### Plot1: RawReads vs Unique Reads
###########
RawReads <-read.table("/wsu/home/hm/hm80/hm8004/piquelab/Dany/gxp/bams/total_reads.txt",header=F, stringsAsFactors=F)

AllReads<-read.table("/wsu/home/hm/hm80/hm8004/piquelab/Dany/gxp/bams/all_counts.txt",header=F, stringsAsFactors=F)

colnames(RawReads) <- c("Plate_ID", "RawReads")

colnames(AllReads) <- c("Plate_ID", "Sorted", "Quality", "Clean")

RawReads$Plate_ID <- substring(RawReads$Plate_ID,1,7)                           
AllReads$Plate_ID <- substring(AllReads$Plate_ID,1,7)

##Bringing treatment information so plots can be organized by treatment
cov.file <- "/wsu/home/hm/hm80/hm8004/piquelab/Dany/gxp/GxP1_metainfo_mod.csv"
cv <- read.table(cov.file, stringsAsFactors=F, header=T, comment="", sep = ",")

colnames(cv)[1] <- "Plate_ID"
cv$treatment <- sub("_","",substring(cv$Treatment,1,4))
cv$concantration <- substring(cv$Treatment,5,20)

#removing spaces
cv$treatment <- trimws(cv$treatment, which = "both")  
cv$concantration <- trimws(cv$concantration, which = "both")   


temp <- merge(RawReads,AllReads, "Plate_ID")
PBMCMerged <- merge(temp,cv, "Plate_ID")
PBMCMerged$RawReads <- PBMCMerged$RawReads * 2
PBMCMerged$NonDup <- PBMCMerged$Clean / PBMCMerged$RawReads #Proportion of reads non-duplicates

pdf("./QC_Unique_vs_RawReads_All_treatment.pdf")
p <- ggplot(PBMCMerged, aes(RawReads, NonDup, color= Plate_ID))
PBP_Unique_vs_RawReads_All_Plate <- p + 
	geom_point(aes(shape= treatment)) + 
	labs(title = "Unique Read Percentage vs Raw Reads", x = "log10(Raw Read Count)", y = "Unique Read Percentatge") +
	theme(plot.title = element_text(hjust = 1, size = 10, color="blue", face = "bold.italic"))
									
print(PBP_Unique_vs_RawReads_All_Plate)  
dev.off()


#TIMEPOINTS 6hrs

cv_6hrs <- cv[grep("6",cv$Timepoint),]

temp <- merge(RawReads,AllReads, "Plate_ID")
PBMCMerged <- merge(temp,cv_6hrs, "Plate_ID")
PBMCMerged$RawReads <- PBMCMerged$RawReads * 2
PBMCMerged$NonDup <- PBMCMerged$Clean / PBMCMerged$RawReads #Proportion of reads non-duplicates

pdf("./QC_Unique_vs_RawReads_All_treatment_6hrs.pdf")
p <- ggplot(PBMCMerged, aes(RawReads, NonDup, color= Plate_ID))
PBP_Unique_vs_RawReads_All_Plate <- p + 
	geom_point(aes(shape= treatment)) + 
	labs(title = "Unique Read Percentage vs Raw Reads(6hrs)", x = "log10(Raw Read Count)", y = "Unique Read Percentatge") +
	theme(plot.title = element_text(hjust = 1, size = 10, color="blue", face = "bold.italic"))
									
print(PBP_Unique_vs_RawReads_All_Plate)  
dev.off()


#TIMEPOINTS 24hrs

cv_24hrs <- cv[grep("24",cv$Timepoint),]

temp <- merge(RawReads,AllReads, "Plate_ID")
PBMCMerged <- merge(temp,cv_24hrs, "Plate_ID")
PBMCMerged$RawReads <- PBMCMerged$RawReads * 2
PBMCMerged$NonDup <- PBMCMerged$Clean / PBMCMerged$RawReads #Proportion of reads non-duplicates


pdf("./QC_Unique_vs_RawReads_All_treatment_24hrs.pdf")
p <- ggplot(PBMCMerged, aes(RawReads, NonDup, color= Plate_ID))
PBP_Unique_vs_RawReads_All_Plate <- p + 
	geom_point(aes(shape= treatment)) + 
	labs(title = "Unique Read Percentage vs Raw Reads(24hrs)", x = "log10(Raw Read Count)", y = "Unique Read Percentatge") +
	theme(plot.title = element_text(hjust = 1, size = 10, color="blue", face = "bold.italic"))
									
print(PBP_Unique_vs_RawReads_All_Plate)  
dev.off()

#All

temp <- merge(RawReads,AllReads, "Plate_ID")
PBMCMerged <- merge(temp,cv, "Plate_ID")
PBMCMerged$RawReads <- PBMCMerged$RawReads * 2
PBMCMerged$NonDup <- PBMCMerged$Clean / PBMCMerged$RawReads #Proportion of reads non-duplicates


pdf("./QC_Unique_vs_RawReads_All_treatment_All.pdf")
p <- ggplot(PBMCMerged, aes(RawReads, NonDup, color= Plate_ID))
PBP_Unique_vs_RawReads_All_Plate <- p + 
	geom_point(aes(shape= treatment)) + 
	labs(title = "Unique Read Percentage vs Raw Reads(24hrs)", x = "log10(Raw Read Count)", y = "Unique Read Percentatge") +
	theme(plot.title = element_text(hjust = 1, size = 10, color="blue", face = "bold.italic"))
									
print(PBP_Unique_vs_RawReads_All_Plate)  
dev.off()


##########################################

counts <- AllReads

#colnames(counts) <- c("Plate_ID", "Sorted", "Quality", "Clean")
counts$Sorted.adj <- counts$Sorted - counts$Quality
counts$Quality.adj <- counts$Quality - counts$Clean

total_reads <- RawReads
colnames(total_reads) <- c("Plate_ID", "Total.Reads")

library(reshape2)
countsMelted <- melt(counts, id.vars = "Plate_ID", measure.vars = c( "Sorted.adj", "Quality.adj", "Clean"))

#countsMelted <- melt(counts, id.vars = "Plate_ID", measure.vars = c( "Quality.adj", "Clean"))

##library(ggplot2)

levels(countsMelted$variable)

countsMelted$variable <- factor(countsMelted$variable, levels = c("Clean", "Quality.adj", "Sorted.adj"))

#countsMelted$variable <- factor(countsMelted$variable, levels = c("Clean", "Quality.adj"))

countsMelted<-merge(countsMelted, total_reads, by="Plate_ID")

counts_Melted<-rbind(countsMelted[countsMelted$variable == "Clean",], countsMelted[countsMelted$variable == "Quality.adj",])
countsMelted$Total.Reads<-countsMelted$Total.Reads *2 #Because Hisat reports fragments instead of reads

pdf("./QC_Melted.pdf", height = 20, width = 20)
p <- ggplot(data = countsMelted, aes(x = Plate_ID, y = value, fill = factor(variable))) + 
  geom_bar(stat="identity", position= position_stack (reverse=TRUE)) +
  geom_point(data=countsMelted, aes(y = Total.Reads, x=Plate_ID))+ 
  theme(text = element_text(size=14,face= "bold.italic"))+
  coord_flip()

p
dev.off()

#######################################

tempTable <- merge(RawReads, AllReads, "Plate_ID")
tempTable <- merge(tempTable, cv, "Plate_ID")
tempTable$RawReads <- tempTable$RawReads * 2
tempTable$QvTR <- tempTable$Quality / tempTable$RawReads #proportion of quality reads to total reads
pdf("./QC_Quality_versus_RawReads.pdf")
p <- ggplot(tempTable, aes(RawReads, QvTR, color=Plate_ID))
QC_Quality_v_RawReads <- p + geom_point(aes(shape= treatment)) +
	labs(title = "Quality Read Percentage vs Raw Reads", x = "log10(Raw Read Count)", y = "Quality Read Percentage") +
	theme(plot.title = element_text(hjust = 1, size = 10, color="blue", face = "bold.italic"))

print(QC_Quality_v_RawReads)
dev.off()


######################################


tempTableCleanvQuality <- merge(RawReads, AllReads, "Plate_ID")
tempTableCleanvQuality <- merge(tempTableCleanvQuality, cv, "Plate_ID")
tempTableCleanvQuality$CleanVQuality <- tempTableCleanvQuality$Clean / tempTableCleanvQuality$Quality
pdf("./QC_Clean_Versus_QualityReads.pdf")
p <- ggplot(tempTableCleanvQuality, aes(Quality, CleanVQuality, color=Plate_ID))
QC_Clean_v_Quality <- p + geom_point(aes(shape = treatment)) + 
	labs(title = "Clean Versus Quality Reads", x = "log10(Quality Read Count)", y = "Clean Versus Quality Read Percentage")+
	theme(plot.title = element_text(hjust = 1, size = 10, color="blue", face = "bold.italic"))


print(QC_Clean_v_Quality)
dev.off()

