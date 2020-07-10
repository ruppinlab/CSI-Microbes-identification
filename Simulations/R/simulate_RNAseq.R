library(Rsubread)
library(Biostrings)
set.seed(strtoi(snakemake@wildcards[["sample"]]))

fasta = readDNAStringSet(snakemake@input[[1]])

expr = matrix(1, ncol=1, nrow=length(fasta))


simReads(transcript.file=snakemake@input[[1]], expression.levels=expr,
  output.prefix=snakemake@params[[1]], library.size=100000, simulate.sequencing.error=FALSE)
