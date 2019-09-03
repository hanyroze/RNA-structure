#Specify a working directory

getwd()
setwd("/Users/haniehroodashty/bin")

#index fasta file
indexFa("danRer11.fa")
fa <- FaFile("danRer11.fa")

library(data.table)
library(ORFik)
library(GenomicAlignments)

#Load the aligned reads of in vivo(in vitro)

invivo1 <- GRanges(readGAlignments("NAI-N3_Oblong_S21_L008.bam"))
invivo2 <- GRanges(readGAlignments("NAI-N3_old_Oblong_S24_L008.bam"))

typeSample = invivo2

seqlevelsStyle(typeSample) <- "NCBI"

# Load the annotation file
library(GenomicFeatures)
txdb <- makeTxDbFromGFF (file = "annotation_danRer11_edited.gtf" )
seqlevelsStyle(txdb) <- "NCBI"

#Filter transcripts
txNames <- filterTranscripts(txdb,100,100,100)

# Get leaders
leaders <- fiveUTRsByTranscript(txdb,use.names=TRUE)[txNames]

# Get the CDS grouped by transcript
cds <- cdsBy(txdb,"tx", use.names = TRUE)[txNames]

# Get the trailers
trailers <- threeUTRsByTranscript(txdb,use.names=TRUE)[txNames]

# Get the exons grouped by transcript
tx <- exonsBy(txdb, by = "tx", use.names = TRUE)[txNames]

#Make coverage plot of one transcript by ggplot
tx_cvg <- coverageByTranscript(typeSample,tx,ignore.strand=TRUE)

a = IRanges::IntegerList(tx_cvg)
l <- lengths(tx_cvg,use.names = FALSE)
l <- rep.int(seq.int(length(l)), l)
res <- data.table(count= unlist(a) , genes= l)

#Give one gene (first or second ,..)
gene <- 1
test_gene<- res[genes==gene,]
test_gene[, position := cumsum(genes)]

#Plot one gene (count across position)
library(ggplot2)
plot <- ggplot(test_gene, aes(x = position, y = count)) + geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,
                                   vjust = 0.5)) + labs(title = paste("gene:", names(tx)[gene])) + ylab("Averaged counts") 

print(plot)

# Get windowCoverage in scale=100 for all transcripts (all width of transcripts are divided to 100 (so under 100 width is ignored))

seqlevelsStyle(tx) <- "NCBI"
coverage <- scaledWindowPositions(tx,typeSample)
coverageTran <- windowCoveragePlot(coverage = coverage,scoring ="transcriptNormalized")
print(coverageTran)

# Get windowCoverage in scale=100 for regions per transcripts
coverageLeader <- scaledWindowPositions(leaders,typeSample)
coverageLeader[, `:=` (fraction="In vivo replicate one",feature="transcripts")]

coverageCds <- scaledWindowPositions(cds,typeSample)
coverageCds[, `:=` (fraction="In vivo replicate one",feature="transcripts")]

coverageTrailers <- scaledWindowPositions(trailers,typeSample)
coverageTrailers[, `:=` (fraction="In vivo replicate one",feature="transcripts")]

coverage <- rbindlist(list(coverageLeader,coverageCds,coverageTrailers))
coverageRegions <- windowCoveragePlot(coverage = coverage ,scoring ="transcriptNormalized,")
print(coverageRegions)

# Multiple plot in one image

A<- ggarrange(coverageTran , coverageRegions,
              ncol = 1, nrow = 2,labels = c ("A","B"),
              heights = c(1,1))

# Normalization CPM  for both replicates in vivo(in vitro), htseq_invitro_rep1 is HTSeq_count for replicate one
names(txNames)=c ('gene_Id')
txNames$gene_Id <- sub("\\.\\d+$", "", txNames$gene_Id)# have the same  gene_Id for Invivo and txNames Data
htseq_invivo_rep1 <- merge(htseq_invivo_rep1, txNames, by='gene_Id')

coverageTran$score <- (( coverageTran$score ) / sum (htseq_invivo_rep1$reads) ) * 10 ^ 6
coverageTran$score <- round(coverageTran$score,3)
saveRDS(coverageTran, 'coverageTran_invivo_NormRep1.rds')

coverageLeader$score <- (( coverageLeader$score ) / sum (htseq_invivo_rep1$reads) ) * 10 ^ 6
coverageLeader$score <- round(coverageLeader$score,3)
saveRDS(coverageLeader, 'coverageLeader_invivo_NormRep1.rds')

coverageCds$score <- (( coverageCds$score ) / sum (htseq_invivo_rep1$reads) ) * 10 ^ 6
coverageCds$score <- round(coverageCds$score,3)
saveRDS(coverageCds, 'coverageCds_invivo_NormRep1.rds')

coverageTrailers$score <- (( coverageTrailers$score ) / sum (htseq_invivo_rep1$reads) ) * 10 ^ 6
coverageTrailers$score <- round(coverageTrailers$score,3)
saveRDS(coverageTrailers, 'coverageTrailer_invivo_NormRep1.rds')

# Do the Normalization CPM for replicate two as before
# Add two replicates, make coverage plot for regions per transcripts
# Make coverage plot for invivo data as before
coverageLeader_invivo$score <- coverageLeader_invivo_NormRep1$score + coverageLeader_invivo_NormRep2$score
saveRDS(coverageLeader_invivo, 'coverageLeader_invivoData.rds')

coverageCds_invivo$score <- coverageCds_invivo_NormRep1$score + coverageLeader_invivo_NormRep2$score
saveRDS(coverageCds_invivo, 'coverageCds_invivoData.rds')

coverageTrailer_invivo$score <- coverageTrailers_invivo_NormRep1$score + coverageTrailers_invivo_NormRep2$score
saveRDS(coverageTrailer_invivo, 'coverageTrailer_invivoData.rds')


##Number of reads for each gene from plot_htseq_count.r file

htseq_invivo_rep1 <- read.delim("result_HTSeq_count_NAI-N3_Oblong_S21_L008_yes.txt",header = FALSE)
htseq_invivo_rep <- read.delim("result_HTSeq_count_invitro_MCE_Oblong_yes.txt",header = FALSE)

htseqTypesample = htseq_invitro_rep1

names(htseqTypesample) = c("gene_Id","reads")
n <- dim(htseqTypesample)[1]
htseqTypesample <- htseqTypesample[1:(n-5),]


# Single gene plot for comparing  with vienna
htseqTypesample <- htseqTypesample [order(htseqTypesample$reads,decreasing = TRUE),]
pickThisone <- txNames[1]
# Get transcriptName from geneName by function geneToTxnames from functionMaster.R
sortedTranscripts <- geneToTxnames(txdb,htseqTypesample$gene_Id)

# Pick first transcript from sorted list of transcripts
pickThisone <- sortedTranscripts[1]

library(ORFik)

plotGene <- coveragePerTiling(tx[pickThisone],typeSample,as.data.table = TRUE)
singleGeneplot <- windowCoveragePlot(plotGene, scoring = "sum", title = pickThisone )
print(singleGeneplot)

# Vienna print top gene by function outputFasta from functionMaster.R
outputFasta(fa,tx,gene = pickThisone)


