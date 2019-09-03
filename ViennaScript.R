# Get vienna structure
# 1 get fasta file of gene we are looking at
# 2 send that to viennaRNA
# 3 look at structure

#Specify a working directory
getwd()
setwd("/Users/haniehroodashty/bin")
list.files()

# Load the annotation file
library(GenomicFeatures)
txdb <- makeTxDbFromGFF (file = "annotation_danRer11_edited.gtf" )

#index fasta file
library(Rsamtools)
library(ORFik)
indexFa("danRer11.fa")
fa <- FaFile("danRer11.fa")

seqlevelsStyle(txdb)  <- seqlevelsStyle(fa)
txNames <- filterTranscripts(txdb, 100, 100, 100)

# Get leaders
leaders <- fiveUTRsByTranscript(txdb,use.names=TRUE)[txNames]

# Get the CDS grouped by transcript
cds <- cdsBy(txdb,"tx", use.names = TRUE)[txNames]

# Get the trailers
trailers <- threeUTRsByTranscript(txdb,use.names=TRUE)[txNames]


#tx filtered by txNames
tx <- exonsBy(txdb, use.names = TRUE)[txNames]
tx <- tx[widthPerGroup(tx) > 300]
a <- extractTranscriptSeqs(fa, transcripts = tx)
writeXStringSet(a, filepath = "tx.fasta")

# Now do perl scripts to get minimum free energy 
library(pander)
system(p("/Adams_vienna_explain.sh -f tx.fasta -o ", getwd()))
# Now we have csvs folder which contains am.csv file (MFE of each position of transcripts)
# Now do analysis
setwd("/Users/haniehroodashty/bin/vienna/csvs")

library(data.table)
hits<- fread("am.csv" ,nrow = length(tx), header = FALSE , fill=TRUE )


# Check that all was made 
if(nrow(hits) != length(tx)) stop("did not create all")


h <- hits[,-1]
m <- setDT(melt(t(h)))

best <- m[, .(which.min(value)), by = Var2]
bestValue<- m[, .(score = min(value, na.rm = TRUE)), by = Var2]

saveRDS(best, 'best_filter.rds')
best <- readRDS('best.rds')

positions <- IRanges(best$V1, width = 1)

# Now create window coverage plot
# grl is the minimum free energy position per gene
# pmapfromtranscript means tx coord -> genomic coord

library(GenomicFeatures)

grl <- ORFik:::pmapFromTranscriptF(positions, tx[hits_filter$V1], removeEmpty = TRUE)# check this is correct


#In silico scatterplot
meansOfAllGenes <- data.table ( gene_Id = hits_filter[,1], Silico = rowMeans(hits_filter[,c(-1)], na.rm = T))

names(meansOfAllGenes)[1] = 'gene_Id'

meansOfAllGenes$gene_Id <- seq.int(nrow(meansOfAllGenes))


# Make window coveragePlot of for all transcripts 

txCoverage <- scaledWindowPositions(tx, grl)
txCoverage[, `:=` (fraction = "structures", feature = "transcripts")]

outName <- paste("structure_sum", ".pdf", sep="")
windowCoveragePlot(txCoverage, scoring = "sum")
outName <- paste("structure_zscore", ".pdf", sep="")
windowCoveragePlot(txCoverage, output = outName ,scoring = "zscore")
outName <- paste("structure_transcriptNorm", ".pdf", sep="")
windowCoveragePlot(txCoverage ,scoring = "transcriptNormalized")

# Make coverage plot for regions per transcript

coverageLeader <- scaledWindowPositions(leaders,grl)
coverageLeader[, `:=` (fraction="In silico",feature="transcripts")]

coverageCds <- scaledWindowPositions(cds,grl)
coverageCds[, `:=` (fraction="In silico",feature="transcripts")]

coverageTrailers <- scaledWindowPositions(trailers,grl)
coverageTrailers[, `:=` (fraction="In silico",feature="transcripts")]

coverage <- rbindlist(list(coverageLeader,coverageCds,coverageTrailers))
coverageRegions <- windowCoveragePlot(coverage = coverage ,scoring ="transcriptNormalized,")

# Coverage plot of all transcripts for Minimum free energy 

cov <- m[,-1]

colnames(cov) <- c("genes", "count")
cov <- cov[!is.na(cov$count),]
perGroup <- cov[,.N, by = genes]$N
if ( length(tx) != length(perGroup)) stop('you messed up the filtering of tx, they are not the same!')
cov[, ones := rep.int(1L, length(genes))]
cov[, position := cumsum(ones), by = genes]

cov$ones <- NULL

scaleTo = 100
scoring = "meanPos"

cov[, scalingFactor := ((scaleTo)/perGroup[genes])]
cov[, position := ceiling(scalingFactor * position)]

cov[position > scaleTo]$position <- scaleTo
groupFPF <- quote(list(genes, position))
res <- cov[, .(score = mean(count, na.rm = TRUE)), by = eval(groupFPF)]
res$feature <- 'transcripts'
res$fraction <- 'In silico'

outName <- paste("MFE_distribution", ".pdf", sep="")
coverageTx_plot_insilico<-windowCoveragePlot(res, scoring = "transcriptNormalized", colors = c("violetred3 ","deeppink4"))+ theme_bw(base_size = 15) +ylab("TranscriptNormalized MFE") + theme(legend.position = "none")
print(coverageTx_plot_insilico)


# Make leader,cds and trailer metacoverage with all free energy positions
# Per gene, find which position cds start, trailers start in transcript coordinate for tx 
names(cds) <- sub(x = names(cds), pattern = '_', replacement = '.')
names(tx) <- sub(x = names(tx), pattern = '_', replacement = '.')
names(trailers) <- sub(x = names(trailers), pattern = '_', replacement = '.')

cdsTrans <- ORFik:::asTX(cds,tx)
cdsStarts <- ORFik:::startSites(cdsTrans)
cdsStartsLong <- cdsStarts[cov$genes]

trailerTrans <- ORFik:::asTX(trailers,tx)
trailerStarts <- ORFik:::startSites(trailerTrans)
trailerStartsLong <- trailerStarts[cov$genes]

# Split hits per gene, into leader, cds, trailer

leaderhits <- cov[position < cdsStartsLong, ]

cdshits <- cov[position >= cdsStartsLong & position < trailerStartsLong, ]
cdshits$position <- cdshits$position - cdsStarts[cdshits$genes] + 1

# The position of cov here is not sacled here
trailerhits <- cov[position >= trailerStartsLong, ]
trailerhits$position <- trailerhits$position - trailerStarts[trailerhits$genes] + 1

perGroupLead <- leaderhits[,.N,by = genes]
perGroupCds <- cdshits[,.N,by = genes]
perGroupTrail <- trailerhits[,.N,by = genes]

scaleTo = 100
scoring = "meanPos"

leaderhits <- merge(x = leaderhits, y = perGroupLead, by = "genes")
leaderhits$scalingFactor <- scaleTo/leaderhits$N
leaderhits[, position := ceiling(scalingFactor * position)]
leaderhits$feature <- 'Leaders'

cdshits <- merge(x = cdshits, y = perGroupCds, by = "genes")
cdshits$scalingFactor <- scaleTo/cdshits$N
cdshits[, position := ceiling(scalingFactor * position)]
cdshits$feature <- 'Cds' 


trailerhits <- merge(x = trailerhits, y = perGroupTrail, by = "genes")
trailerhits$scalingFactor <- scaleTo/trailerhits$N
trailerhits[, position := ceiling(scalingFactor * position)]
trailerhits$feature <- 'Trailers'

fractionhits <- rbindlist(list(leaderhits,cdshits,trailerhits))
fractionhits[position > scaleTo]$position <- scaleTo
fractionhits$N<- NULL
fractionhits$scalingFactor <- NULL


groupFPF <- quote(list(genes, position, feature))

resPerTx <- fractionhits[, .(score = mean(count, na.rm = TRUE)), by = eval(groupFPF)]
resPerTx$fraction <- 'In silico'

# Make coverage plot for leader,cds and trailers of minimum free energy 
outName <- paste("MFE_distributionLeaderCdsTrailer", ".pdf", sep="")
coverageRegionTx_plot_insilico<- windowCoveragePlot(resPerTx ,scoring = 'transcriptNormalized', colors = c("violetred3 ","deeppink4"))+ theme_bw(base_size = 15) +ylab("TranscriptNormalized MFE") + theme(legend.position = "none")
print(coverageRegionTx_plot_insilico)


# Density plot of insilico data for best(most negative MFE) and mean score
# For best
BestScore_Position <- merge(x = best, y = bestValue, by = "Var2")

LeaderBestPosition <- BestScore_Position[V1 < cdsStarts, ]
LeaderBestPosition[, `:=` (feature="Leader")]

CdsBestPosition <- BestScore_Position[V1 >= cdsStarts & V1 < trailerStarts, ]
CdsBestPosition[, `:=` (feature="Cds")]

TrailerBestPosition<- BestScore_Position[V1 >= trailerStarts, ]
TrailerBestPosition[, `:=` (feature="Trailer")]

BestScoreFeature <- rbindlist(list(LeaderBestPosition,CdsBestPosition,TrailerBestPosition))
ggplot(BestScoreFeature,aes(x=score, fill = feature)) + geom_density(alpha=0.5)+ theme(
  axis.text = element_text(size = 18), axis.title  = element_text(size = 20 , face = "bold"),legend.title = element_text(size=20),legend.text= element_text(size=15)) 
saveRDS(BestScoreFeature, 'BestScoreFeature.rds')

# For mean score

LeaderMeanPerGene <- leaderhits[, .(MeanScore = mean(count, na.rm = TRUE)), by = genes] 
LeaderMeanPerGene[, `:=` (feature="Leader")]

CdsMeanPerGene <- cdshits[, .(MeanScore = mean(count, na.rm = TRUE)), by = genes] 
CdsMeanPerGene[, `:=` (feature="Cds")]

TrailerMeanPerGene <- trailerhits[, .(MeanScore = mean(count, na.rm = TRUE)), by = genes] 
TrailerMeanPerGene[, `:=` (feature="Trailer")]

MeanScoreFeature <- rbindlist(list(LeaderMeanPerGene,CdsMeanPerGene,TrailerMeanPerGene))
saveRDS(MeanScoreFeature, 'MeanScoreFeature.rds')
ggplot(MeanScoreFeature ,aes(x=MeanScore, fill = feature)) + geom_density(alpha=0.5) + theme(
  axis.text = element_text(size = 18), axis.title  = element_text(size = 20 , face = "bold"),legend.title = element_text(size=20),legend.text= element_text(size=15)) 
