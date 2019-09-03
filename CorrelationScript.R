#Correlation of replicate one and replicate two in vivo (in vitro ) by computing difference percentage of count reads
#Correlation of in vivo and in vitro data
#Correlation of in silico and in vivo as well as in vitro  by computing difference percentage of hits

#Load the aligned reads of in vivo (in vitro)

invivo1 <- GRanges(readGAlignments("NAI-N3_Oblong_S21_L008.bam"))
invivo2 <- GRanges(readGAlignments("NAI-N3_old_Oblong_S24_L008.bam"))


seqlevelsStyle(invivo1) <- "NCBI"
seqlevelsStyle(tx) <- "NCBI"

# Replicate one 
invivo1_new <- ORFik:::convertToOneBasedRanges(invivo1)

txHits <- countOverlaps(tx, invivo1_new)
leaderHits <- countOverlaps(leaders, invivo1_new)/txHits 
cdsHits <- countOverlaps(cds, invivo1_new) /txHits 
trailerHits <- countOverlaps(trailers, invivo1_new) /txHits 
Hits_Invivorep1 <- data.frame(leaderHits, cdsHits, trailerHits) 

# Percentage of hits per leader, cds and trailer in each transcript
Hitst_Invivorep1<- saveRDS(Hits_Invivorep1,'Hits_Invivorep1.rds')

# Replicate two
seqlevelsStyle(invivo2) <- "NCBI"

invivo2_new <- ORFik:::convertToOneBasedRanges(invivo2)

txHits <- countOverlaps(tx,invivo2_new)
leaderHits <- countOverlaps(leaders, invivo2_new) / txHits
cdsHits <- countOverlaps(cds, invivo2_new) / txHits
trailerHits <- countOverlaps(trailers, invivo2_new)/ txHits 
Hits_Invivorep2 <- data.frame(leaderHits, cdsHits, trailerHits) 

# Percentage of hits per leader, cds and trailer
saveRDS(Hits_Invivorep2, 'Hits_Invivorep2.rds')
Hits_Invivorep2 <- readRDS('Hits_Invivorep2.rds')

# Difference percentage between two replicates 
# Replace na value with zero
Hits_Invivorep1[is.na(Hits_Invivorep1)] <- 0
Hits_Invivorep2[is.na(Hits_Invivorep2)] <- 0

# Find total difference percentage for each transcripts between invivo replicate one and replicate two
totalScore <- abs(Hits_Invivorep1$leaderHits - Hits_Invivorep2$leaderHits) + abs(Hits_Invivorep1$cdsHits - Hits_Invivorep2$cdsHits) +
  abs(Hits_Invivorep1$trailerHits - Hits_Invivorep2$trailerHits)
df <- data.frame(totalScore)
df$totalScore <- round (df$totalScore ,2)
df$totalScore <- (df$totalScore) *100

# Plot of percentage difference between in vivo replicates for each transcript
sp <- ggplot(df, aes(x = totalScore)) + geom_freqpoly()
sp +scale_x_continuous(name = "\n% Difference of count reads", limits = c(0,300)) + scale_y_continuous(name = "Frequency of genes\n")+theme(
  axis.text = element_text(size = 18), axis.title  = element_text(size = 20 , face = "bold")
) 


# Pool two replicates in vivo (in vitro)

txHitsInvivo1 <- countOverlaps(tx, invivo1_new)
txHitsInvivo2 <- countOverlaps(tx, invivo2_new)
txHitsInvivo <- (txHitsInvivo1 / sum (htseq_invivo_rep1$reads) + txHitsInvivo2  / sum (htseq_invivo_rep2$reads))* 10 ^ 6

leaderHitsInvivo1 <- countOverlaps(leaders, invivo1_new)
leaderHitsInvivo2 <- countOverlaps(leaders, invivo2_new)
LeaderHitsInvivo <- (leaderHitsInvivo1 / sum (htseq_invivo_rep1$reads) + leaderHitsInvivo2  / sum (htseq_invivo_rep2$reads))* 10 ^ 6
LeaderHitsInvivo <- LeaderHitsInvivo /txHitsInvivo

cdsHitsInvivo1 <- countOverlaps(cds, invivo1_new)
cdsHitsInvivo2 <- countOverlaps(cds, invivo2_new)
cdsHitsInvivo <- (cdsHitsInvivo1 / sum (htseq_invivo_rep1$reads) + cdsHitsInvivo2 / sum (htseq_invivo_rep2$reads))* 10 ^ 6
cdsHitsInvivo <- cdsHitsInvivo / txHitsInvivo

trailerHitsInvivo1 <- countOverlaps(trailers, invivo1_new)
trailerHitsInvivo2 <- countOverlaps(trailers, invivo2_new)
trailerHitsInvivo <- (trailerHitsInvivo1 / sum (htseq_invivo_rep1$reads) + trailerHitsInvivo2 / sum (htseq_invivo_rep2$reads))* 10 ^ 6
trailerHitsInvivo <- trailerHitsInvivo / txHitsInvivo

# percentage of count reads per leader, cds and trailer in each transcript
Hits_Invivo <- data.frame(LeaderHitsInvivo, cdsHitsInvivo, trailerHitsInvivo) 
saveRDS(Hits_Invivo, 'Hits_Invivo.rds')
Hitst_Invivo <- readRDS('Hits_Invivo.rds')

# Do the same script above for in vitro replicates and get in vitro data hits 
# Difference percentage between Invitro and Invivo for each transcript
# Replacing na value with zero

Hits_Invivo[is.na(Hits_Invitro)] <- 0
Hits_Invivo[is.na(Hits_Invivo)] <- 0

# Find total difference percentage for each transcripts between in vivo and in vitro data
totalScore <- abs(Hits_Invivo$LeaderHitsInvivo - Hits_Invitro$leaderHitsInvitro) + abs(Hits_Invivo$cdsHitsInvivo - Hits_Invitro$cdsHitsInvitro) +
  abs(Hits_Invivo$trailerHitsInvivo - Hits_Invitro$trailerHitsInvitro)
df <- data.frame(totalScore)
df$totalScore <- round (df$totalScore ,2)

# Plot of percentage difference between in vivo and in vitro data
sp <- ggplot(df, aes(x = totalScore)) + geom_freqpoly()
sp +scale_x_continuous(name = "Score difference") + scale_y_continuous(name = "Frequency of genes")


# Hits In silico
# Using perGroupLead,perGroupCds and perGroupTrail from ViennaScript.R

perGroupLead <- leaderhits[,.N,by = genes]
perGroupCds <- cdshits[,.N,by = genes]
perGroupTrail <- trailerhits[,.N,by = genes]
# Using cov data from ViennaScript.R
pergroupTx <- cov[, .N,by = genes]

# Find percentage of hits for regions per transcript 

Txlead<- merge(x = pergroupTx, y= perGroupLead, by ="genes", all=TRUE)
TxleadCds <- merge ( x = Txlead , y= perGroupCds , by ="genes", all = TRUE)
names(TxleadCds) <- c("genes","Tx","leader","Cds")

TxleadCdsTrai <- merge (x = TxleadCds , y= perGroupTrail ,by ="genes", all = TRUE)
colnames(TxleadCdsTrai)[colnames(TxleadCdsTrai)=="N"] <- "Trailer"

TxleadCdsTrai[is.na(TxleadCdsTrai)] <- 0

TxleadCdsTrai$leader <-TxleadCdsTrai$leader/TxleadCdsTrai$Tx 
TxleadCdsTrai$Cds <-TxleadCdsTrai$Cds/TxleadCdsTrai$Tx 
TxleadCdsTrai$Trailer <-TxleadCdsTrai$Trailer/TxleadCdsTrai$Tx 
Hits_AllInsilico <- TxleadCdsTrai
saveRDS(Hits_AllInsilico, 'Hits_AllInsilico.rds')
Hitst_AllInsilico <- readRDS('Hits_AllInsilico.rds')

# Difference percentage between Insilico and Invitro or Invivo and Insilico
# Replace na value with zero

totalScore <- abs(Hits_Invivo$LeaderHitsInvivo - Hits_AllInsilico$leader) + abs(Hits_Invivo$cdsHitsInvivo - Hits_AllInsilico$Cds) +
  abs(Hits_Invivo$trailerHitsInvivo - Hits_AllInsilico$Trailer)
df <- data.frame(totalScore)
df$totalScore <- round (df$totalScore ,2)

sp <- ggplot(df, aes(x = totalScore)) + geom_freqpoly()
sp +scale_x_continuous(name = "Score difference") + scale_y_continuous(name = "Frequency of genes")

