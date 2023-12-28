library(dplyr)
library(tidyverse)

source("../../SAHMI/functions/read_kraken_reports.r")

# report = read.delim('output/kraken2/P1/SCAF2963_3_Live.kraken.report.txt', header = F)
report = read.delim(snakemake@input[["report"]], header = F)
report$V8 = trimws(report$V8)
report[report$V8 %in% c('Homo sapiens', 'Bacteria', 'Fungi', 'Viruses'), ]

# sckmer data
# kmer_data = read.table('output/SAHMI/P1/SCAF2963_3_Live.sckmer.txt', header = T)
kmer_data = read.table(snakemake@input[["kmer"]], header = T)
head(kmer_data)

c = kmer_data %>% 
  subset(kmer > 1) %>%
  group_by(taxid) %>%
  mutate(nn = n()) %>%
  subset(nn > 3) %>% 
  group_by(taxid) %>%
  summarize(r = cor.test(kmer, uniq, method = 'spearman')$estimate,
            p = cor.test(kmer, uniq, method = 'spearman')$p.value,
            .groups='keep') %>%
  mutate(p = p.adjust(p))

c$name = report$V8[match(c$taxid, report$V7)] # add taxa names 
c

# k1 <- read_kraken_reports(c("output/kraken2/P1/SCAF2963_3_Live.kraken.report.txt"))
k1 <- read_kraken_reports(c(snakemake@input[["report"]]))
# write.table(c, file="output/SAHMI/P1/SCAF2963_3_Live.barcode.kmer.hits.tsv", sep="\t")
write.table(c, file=snakemake@output[["hits"]], sep="\t")
# write.table(k1, file="output/SAHMI/P1/SCAF2963_3_Live.kraken.report.rpmm.tsv", sep="\t")
write.table(k1, file=snakemake@output[["rpmm"]], sep="\t")

