#Specify a working directory

getwd()
setwd("/Users/haniehroodashty/bin")
list.files()

#load HTseq_count data (two replicates for in vitro)


htseq_invivo_rep1 <- read.delim("result_HTSeq_count_invivo1.txt",header = FALSE)
htseq_invivo_rep2 <- read.delim("result_HTSeq_count_invivo2.txt",header = FALSE)

htseq_invivo = htseq_invivo_rep1

names(htseq_invivo) = c("gene_Id","reads")
m <- dim(htseq_invivo)[1]

#Remove summery text from the 5 last rows of HtSeq_count data
htseq_invivo<- htseq_invivo[1:(m-5),]


#Consider only genes which have reads
high_count <- subset(htseq_invivo, reads > 0)

#Find how many genes have zero,1,2, ... reads(frequency of genes)
read_freq_table = table(htseq_invivo$reads)

read_freq = as.data.frame(read_freq_table)
names(read_freq)[1] = 'reads'

#Find how many genes have 1,2,3,.. (frequency of genes)
high_count_freq_table= table(high_count$reads)
high_count_freq= as.data.frame(high_count_freq_table)

names(high_count_freq)[1]= 'Reads'

high_count_freq$Reads<-as.numeric(as.character(high_count_freq$Reads))

#Save data in corresponding data

saveRDS(high_count_freq, 'high_count_freq_Invivo_rep1.rds')

library(ggpubr)
library(ggplot2)
library(scales)

#ScatterPlot of Inviro replicate one or two
#Use high_count-freq_InvitroRep1 or high_count-freq_InvitroRep2 or high_count_freq_InvitroData as input data
high_count_freq = high_count_freq_Invivo_rep1
sp <- ggplot(high_count_freq,x = "Reads",y = "Freq" , aes(x = Reads, y = Freq)) + geom_point()

printplot <- sp +scale_y_continuous( name = "Frequency of genes\n" ,trans = 'log2',limits = c(1,3000)) + scale_x_continuous(name = "\nNumber of reads",trans='log2',limits = c(0.1,100000)) + theme(
  axis.text = element_text(size = 18), axis.title  = element_text(size = 20 , face = "bold")
) 
print(printplot)


#Compaire invivo replicate one and two (correlation of two replicates)
#Use htseq_invivo_rep1 anad htseq_invivo_rep2 as input data (data are after removing summery text)
htseq_invivoData1 <- readRDS('high_count_freq_Invitro_rep1.rds')
htseq_invivoData2 <- readRDS('high_count_freq_Invitro_rep2.rds')

#Change the second column of data for merging two data
names(htseq_invivoData1)[2]= 'Invivo1'
names(htseq_invivoData2)[2]= 'Invivo2'

mergeInvivo1Invivo2 <- merge(x = htseq_invivoData1, y = htseq_invivoData2, by = "gene_Id", all.x = TRUE)
mergeNew<- mergeInvivo1Invivo2[,2:3]+1

sp <- ggscatter(mergeNew ,x = "Invivo1",y = "Invivo2" , add = "reg.line" , add.params = list(color = "blue", fill = "lightgray"), 
                 conf.int = TRUE ,cor.coef = TRUE, cor.coef.size = 8 ,cor.method = "pearson")
print(sp +scale_y_continuous( name = "Invivo2\n" ,trans = 'log2',limits = c(1,600000)) + scale_x_continuous(name = "\nInvivo1",trans='log2',limits = c(1,100000))+ theme(
  axis.text = element_text(size = 18), axis.title  = element_text(size = 20 , face = "bold")
) )

cor.test(mergeNew$Invivo1,mergeNew$Invivo2 , method = "pearson")

# Pool two in vitro(in vivo) replicates if replicates have a good correlation
# First we need normalization CPM for both replicates
# Second add the number of reads rep1 to rep2 

htseq_invivo_rep1$reads <- (( htseq_invivo_rep1$reads ) / sum (htseq_invivo_rep1$reads) ) * 10 ^ 6
htseq_invivo_rep1$reads <- round (htseq_invivo_rep1$reads ,3)

# Save data in corresponding data for both replicates

saveRDS(htseq_invivo_rep1 ,'htseq_invivo_NormRep1.rds')
htseq_invivo_NormRep1 <- readRDS('htseq_invivo_NormRep1.rds')

saveRDS(htseq_invivo_rep2, 'htseq_invivo_NormRep2.rds')
htseq_invivo_NormRep2 <- readRDS('htseq_invivo_NormRep2.rds')

# DO the same as befor for making scatterPlot of invivo data
htseq_invivo$reads <- htseq_invivo_NormRep1$reads + htseq_invivo_NormRep2$reads

saveRDS(htseq_invivo, 'htseq_invivoData.rds')
htseq_invivoData <- readRDS('htseq_invivoData.rds')


