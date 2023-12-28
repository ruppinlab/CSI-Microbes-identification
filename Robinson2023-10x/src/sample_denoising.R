library(ggplot2)
library(dplyr)

source("../../SAHMI/functions/read_kraken_reports.r")

kr <- read_kraken_reports(c("output/kraken2/P1/SCAF2961_1_Uninfected.kraken.report.txt", "output/kraken2/P1/SCAF2962_2_HK.kraken.report.txt", "output/kraken2/P1/SCAF2963_3_Live.kraken.report.txt", "output/kraken2/P1/SCAF2965_5_Live.kraken.report.txt"))

kr = kr %>%
  group_by(taxid) %>%
  mutate(nn = n()) %>%
  subset(nn > 2) %>% 
  select(-nn)

c2 = kr %>%
  subset(rank %in% c('G', 'S')) %>% 
  group_by(name) %>%
  summarize(r1 = cor(min,uniq,method='spearman'),
            r2 = cor(min,reads,method='spearman'),
            r3 = cor(reads,uniq,method='spearman'),
            p1 = cor.test(min,uniq,method='spearman')$p.value,
            p2 = cor.test(min,reads,method='spearman')$p.value,
            p3 = cor.test(reads,uniq,method='spearman')$p.value
            )

c2 %>% subset(r1>0 & r2>0 & r3>0 & p1<0.05 & p2<0.05 & p3<0.05)
# yields Acinetobacter cumulans and Georgenia faecalis