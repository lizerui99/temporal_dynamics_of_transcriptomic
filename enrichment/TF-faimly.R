library(clusterProfiler)
library(GOplot)
library(tidyverse)

term2gene <- read.delim('TF-class.txt',header = FALSE, sep = '\t')
term2name <- read.delim('TF-anno.txt',header = FALSE,sep = '\t')


gene_data <- read.delim('NEW_N45_padj_phase.txt',header = T,sep = '\t')

gene_select <- filter(gene_data,class == 'up',TF != "#N/A",phase == 5) %>% select(gene_id)

x <- enricher(gene_select$gene,
              TERM2GENE = term2gene,
              TERM2NAME = term2name,
              pvalueCutoff = 0.05,
              pAdjustMethod = 'BH',
              qvalueCutoff = 0.05, 
              maxGSSize = 500)

dotplot(x,showCategory=100,decreasing=T,label_format =50)


write.table(x, 'N45_up_stage4_tf_rich.tab', sep = '\t', row.names = FALSE, quote = FALSE)
