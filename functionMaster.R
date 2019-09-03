# Two function are using in Coverage_plot.R
# Get transcriptName from geneName
geneToTxnames <- function(txdb,geneName){
  names <-transcriptLengths(txdb)
  res <- names$gene_id[names$gene_id %in% geneName]
  res <- chmatch(as.character(geneName),as.character(names$gene_id))
  return(names$tx_name[res])
}


# Vienna print top gene
outputFasta <- function(fa,refrence,gene=NULL){
  seqlevelsStyle(refrence) <- seqlevelsStyle(fa)
  if (is.null(gene) ){
    seq <- extractTranscriptSeqs(fa,refrence)
  } else {
    seq <- extractTranscriptSeqs(fa, refrence[gene])
  }
  
  writeXStringSet(seq, filepath = "geneNP_059333.1.fasta")
  return()
}