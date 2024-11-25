rm(list=ls())
library(clusterProfiler)
library(GOplot)
library(ggplotify)
library(patchwork)
library(DelayedArray)
library(dbplyr)
library(tidyverse)

##########################################################################################
go_anno <- read.delim('B73_AGPv4_GSMER.annot', header = FALSE, stringsAsFactors = FALSE)
names(go_anno) <- c('gene_id', 'ID')
go_class <- read.delim('go_class.txt', header = FALSE, stringsAsFactors = FALSE)
names(go_class) <- c('ID', 'Description', 'Ontology')
go_anno <- merge(go_anno, go_class, by = 'ID', all.x = TRUE)
##########################################################################################

gene_select <- read.delim('NEW_D45_padj_phase.txt', header = T,stringsAsFactors = FALSE)
##########################################################################################
gene <- gene_select %>%
  filter(class == 'up',phase == "1" )
gene_list <- gene$gene_id

go_rich <- enricher(gene = gene_list,
TERM2GENE = go_anno[c('ID', 'gene_id')],
TERM2NAME = go_anno[c('ID', 'Description')],
pvalueCutoff = 0.05,
pAdjustMethod = 'BH',qvalueCutoff = 0.05, maxGSSize = 500)
write.table(go_rich, 'd45phase_up1.tab', sep = '\t', row.names = FALSE, quote = FALSE)
up1 <- dotplot(go_rich,showCategory=15,decreasing=T,label_format =50,color = "pvalue",title = 'stage1_up')
##########################################################################################
gene <- gene_select %>%
  filter(class == 'up',phase == "2"|phase == "3")
gene_list <- gene$gene_id

go_rich <- enricher(gene = gene_list,
                     TERM2GENE = go_anno[c('ID', 'gene_id')],
                     TERM2NAME = go_anno[c('ID', 'Description')],
                     pvalueCutoff = 0.05,
                     pAdjustMethod = 'BH',qvalueCutoff = 0.05, maxGSSize = 500)
write.table(go_rich, 'd45phase_up2.tab', sep = '\t', row.names = FALSE, quote = FALSE)
up2 <- dotplot(go_rich,showCategory=15,decreasing=T,label_format =50,color = "pvalue",title = 'stage2_up')
##########################################################################################
gene <- gene_select %>%
  filter(class == 'up',phase == "4")
gene_list <- gene$gene_id

go_rich <- enricher(gene = gene_list,
                     TERM2GENE = go_anno[c('ID', 'gene_id')],
                     TERM2NAME = go_anno[c('ID', 'Description')],
                     pvalueCutoff = 0.05,
                     pAdjustMethod = 'BH',qvalueCutoff = 0.05, maxGSSize = 500)
write.table(go_rich, 'd45phase_up3.tab', sep = '\t', row.names = FALSE, quote = FALSE)
up3 <- dotplot(go_rich,showCategory=15,decreasing=T,label_format =50,color = "pvalue",title = 'stage3_up')
##########################################################################################
##########################################################################################
gene <- gene_select %>%
  filter(class == 'up',phase == "5")
gene_list <- gene$gene_id

go_rich <- enricher(gene = gene_list,
                     TERM2GENE = go_anno[c('ID', 'gene_id')],
                     TERM2NAME = go_anno[c('ID', 'Description')],
                     pvalueCutoff = 0.05,
                     pAdjustMethod = 'BH',qvalueCutoff = 0.05, maxGSSize = 500)
write.table(go_rich, 'd45phase_up4.tab', sep = '\t', row.names = FALSE, quote = FALSE)
up4 <- dotplot(go_rich,showCategory=15,decreasing=T,label_format =50,color = "pvalue",title = 'stage4_up')
##########################################################################################

pdf("d45stage_up.pdf",h=10,w=15)
as.ggplot(up1)+as.ggplot(up2)+as.ggplot(up3)+as.ggplot(up4)
dev.off()

