"""
要使用`fimo`获取全部基因启动子的motif，并用R对感兴趣的基因集中的motif进行富集分析，包括计算FDR和enrichment score，你需要遵循以下步骤：
1. 使用`fimo`对基因启动子序列进行motif搜索。
2. 将`fimo`的输出导入R，并进行预处理。
3. 定义感兴趣的基因集。
4. 使用`fisher.test`进行motif富集分析。
5. 计算FDR和enrichment score。
6. 对结果进行排序和可视化。
"""
# 安装和加载所需的R包
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# 加载所需的R包

rm(list=ls())
library(stats)
library(tidyverse)
# 假设你已经有了一个包含基因启动子序列的FASTA文件，以及一个motif数据库文件
# 使用fimo对启动子序列进行motif搜索
# 这一步通常在命令行中完成，这里只是说明
# system("fimo --oc fimo_output --verbosity 1 --thresh 0.001 motifs.meme promoters.fasta")
# 读取fimo的输出文件
fimo_output <- read.delim("fimo.txt", sep="\t", header=TRUE)
# 过滤fimo输出，只保留具有统计显著性的motif
significant_motifs <- fimo_output[fimo_output$p.value < 0.05, ]



# 定义感兴趣的基因集

interested_genes <- read.delim('time_point_gene_list.txt',header=1) %>% 
  filter(trait=='D38',class=='up',phase==1) %>% 
  select('gene_id')# 替换为你的感兴趣基因
interested_genes <- interested_genes$gene_id

# 获取感兴趣的基因集中的motif
interested_motifs <- significant_motifs[significant_motifs$sequence_name %in% interested_genes, ]
# 计算感兴趣基因集中motif的数量
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算所有motif的总数
total_motif_counts <- table(significant_motifs$motif_id,exclude = NULL)
# 计算在sig中出现，但不在insterest中出现的motif
x <- interested_motifs %>% group_by(motif_id) %>% summarise(n=n())
y <- significant_motifs %>% group_by(motif_id) %>% summarise(n=n())
z <- y$motif_id %in% x$motif_id
t <- y$motif_id[!z]
# 将z添加值intersted_motif_counts中
for (non_motif in t) {
  interested_motifs <- interested_motifs %>% 
    add_row(motif_id=non_motif)
}
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算背景基因集中motif的数量（所有motif减去感兴趣的motif）
background_motif_counts <- total_motif_counts - interested_motif_counts

# 过滤掉在感兴趣基因集或背景基因集中计数为0的motif
motif_names <- names(interested_motif_counts)[interested_motif_counts > 0 & background_motif_counts > 0]



# 创建2x2列联表矩阵
# d <- data.frame(gene.not.interest=c(M-k, N-M-n+k), gene.in.interest=c(k, n-k))
# row.names(d) <- c("In_category", "not_in_category")
contingency_tables <- lapply(motif_names, function(motif) {
  matrix(c(interested_motif_counts[[motif]], # 背景，不同motif出现的数目
           background_motif_counts[[motif]]-interested_motif_counts[[motif]], # 感兴趣基因集中motif数目
           nrow(interested_motifs) - interested_motif_counts[[motif]],
           nrow(significant_motifs)  - nrow(interested_motifs) -  background_motif_counts[[motif]]#background_motif_counts[[motif]] + interested_motif_counts[[motif]] # 背景全部motif - 背景中单个motif出现的次数
           # 基因集中的motif - 基因集中单个motif出现的次数
  ), ncol = 2, byrow = TRUE)
})




# 应用Fisher精确检验，并检查每个列联表是否至少是2x2的
fisher_results <- lapply(contingency_tables, function(table) {
  fisher.test(table)$p.value
})

# 将列表转换为向量
fisher_results <- unlist(fisher_results)
# 计算FDR
fdr_results <- p.adjust(fisher_results, method = "BH")
# 计算Enrichment score
enrichment_scores <- lapply(motif_names, function(motif) {
  if (background_motif_counts[[motif]] > 0) {
    log2((interested_motif_counts[[motif]]/nrow(interested_motifs)) / (background_motif_counts[[motif]]/nrow(significant_motifs)))
  } else {
    NA  # 如果背景中的motif数量为0，则无法计算富集度得分，返回NA
  }
})
enrichment_scores <- unlist(enrichment_scores)
# 将结果转换为数据框，并排序
enrichment_results_df <- data.frame(
  Motif = motif_names,
  Enrichment.Score = enrichment_scores,
  Pvalue = fisher_results,
  FDR = fdr_results
)
sorted_enrichment_results <- enrichment_results_df[order(enrichment_results_df$FDR), ]

# 筛选结果，注释motif信息
filter_enrichment_results <- filter(enrichment_results_df,Pvalue<0.05,FDR<0.05,enrichment_scores>0)
motif_anno <- read.delim('meme_info.txt')
anno_enrichment_result <- merge(filter_enrichment_results,motif_anno,by='Motif')
ALL_anno_enrichment_result <- merge(sorted_enrichment_results,motif_anno,by='Motif')


write.table(ALL_anno_enrichment_result, 'D38_up_phase1.txt', sep = '\t', row.names = FALSE, quote = FALSE)


interested_genes <- read.delim('time_point_gene_list.txt',header=1) %>% 
  filter(trait=='D38',class=='up',phase==2) %>% 
  select('gene_id')# 替换为你的感兴趣基因
interested_genes <- interested_genes$gene_id

# 获取感兴趣的基因集中的motif
interested_motifs <- significant_motifs[significant_motifs$sequence_name %in% interested_genes, ]
# 计算感兴趣基因集中motif的数量
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算所有motif的总数
total_motif_counts <- table(significant_motifs$motif_id,exclude = NULL)
# 计算在sig中出现，但不在insterest中出现的motif
x <- interested_motifs %>% group_by(motif_id) %>% summarise(n=n())
y <- significant_motifs %>% group_by(motif_id) %>% summarise(n=n())
z <- y$motif_id %in% x$motif_id
t <- y$motif_id[!z]
# 将z添加值intersted_motif_counts中
for (non_motif in t) {
  interested_motifs <- interested_motifs %>% 
    add_row(motif_id=non_motif)
}
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算背景基因集中motif的数量（所有motif减去感兴趣的motif）
background_motif_counts <- total_motif_counts - interested_motif_counts

# 过滤掉在感兴趣基因集或背景基因集中计数为0的motif
motif_names <- names(interested_motif_counts)[interested_motif_counts > 0 & background_motif_counts > 0]



# 创建2x2列联表矩阵
# d <- data.frame(gene.not.interest=c(M-k, N-M-n+k), gene.in.interest=c(k, n-k))
# row.names(d) <- c("In_category", "not_in_category")
contingency_tables <- lapply(motif_names, function(motif) {
  matrix(c(interested_motif_counts[[motif]], # 背景，不同motif出现的数目
           background_motif_counts[[motif]]-interested_motif_counts[[motif]], # 感兴趣基因集中motif数目
           nrow(interested_motifs) - interested_motif_counts[[motif]],
           nrow(significant_motifs)  - nrow(interested_motifs) -  background_motif_counts[[motif]]#background_motif_counts[[motif]] + interested_motif_counts[[motif]] # 背景全部motif - 背景中单个motif出现的次数
           # 基因集中的motif - 基因集中单个motif出现的次数
  ), ncol = 2, byrow = TRUE)
})




# 应用Fisher精确检验，并检查每个列联表是否至少是2x2的
fisher_results <- lapply(contingency_tables, function(table) {
  fisher.test(table)$p.value
})

# 将列表转换为向量
fisher_results <- unlist(fisher_results)
# 计算FDR
fdr_results <- p.adjust(fisher_results, method = "BH")
# 计算Enrichment score
enrichment_scores <- lapply(motif_names, function(motif) {
  if (background_motif_counts[[motif]] > 0) {
    log2((interested_motif_counts[[motif]]/nrow(interested_motifs)) / (background_motif_counts[[motif]]/nrow(significant_motifs)))
  } else {
    NA  # 如果背景中的motif数量为0，则无法计算富集度得分，返回NA
  }
})
enrichment_scores <- unlist(enrichment_scores)
# 将结果转换为数据框，并排序
enrichment_results_df <- data.frame(
  Motif = motif_names,
  Enrichment.Score = enrichment_scores,
  Pvalue = fisher_results,
  FDR = fdr_results
)
sorted_enrichment_results <- enrichment_results_df[order(enrichment_results_df$FDR), ]

# 筛选结果，注释motif信息
filter_enrichment_results <- filter(enrichment_results_df,Pvalue<0.05,FDR<0.05,enrichment_scores>0)
motif_anno <- read.delim('meme_info.txt')
anno_enrichment_result <- merge(filter_enrichment_results,motif_anno,by='Motif')
ALL_anno_enrichment_result <- merge(sorted_enrichment_results,motif_anno,by='Motif')


write.table(ALL_anno_enrichment_result, 'D38_up_phase2.txt', sep = '\t', row.names = FALSE, quote = FALSE)


interested_genes <- read.delim('time_point_gene_list.txt',header=1) %>% 
  filter(trait=='D38',class=='up',phase==3) %>% 
  select('gene_id')# 替换为你的感兴趣基因
interested_genes <- interested_genes$gene_id

# 获取感兴趣的基因集中的motif
interested_motifs <- significant_motifs[significant_motifs$sequence_name %in% interested_genes, ]
# 计算感兴趣基因集中motif的数量
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算所有motif的总数
total_motif_counts <- table(significant_motifs$motif_id,exclude = NULL)
# 计算在sig中出现，但不在insterest中出现的motif
x <- interested_motifs %>% group_by(motif_id) %>% summarise(n=n())
y <- significant_motifs %>% group_by(motif_id) %>% summarise(n=n())
z <- y$motif_id %in% x$motif_id
t <- y$motif_id[!z]
# 将z添加值intersted_motif_counts中
for (non_motif in t) {
  interested_motifs <- interested_motifs %>% 
    add_row(motif_id=non_motif)
}
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算背景基因集中motif的数量（所有motif减去感兴趣的motif）
background_motif_counts <- total_motif_counts - interested_motif_counts

# 过滤掉在感兴趣基因集或背景基因集中计数为0的motif
motif_names <- names(interested_motif_counts)[interested_motif_counts > 0 & background_motif_counts > 0]



# 创建2x2列联表矩阵
# d <- data.frame(gene.not.interest=c(M-k, N-M-n+k), gene.in.interest=c(k, n-k))
# row.names(d) <- c("In_category", "not_in_category")
contingency_tables <- lapply(motif_names, function(motif) {
  matrix(c(interested_motif_counts[[motif]], # 背景，不同motif出现的数目
           background_motif_counts[[motif]]-interested_motif_counts[[motif]], # 感兴趣基因集中motif数目
           nrow(interested_motifs) - interested_motif_counts[[motif]],
           nrow(significant_motifs)  - nrow(interested_motifs) -  background_motif_counts[[motif]]#background_motif_counts[[motif]] + interested_motif_counts[[motif]] # 背景全部motif - 背景中单个motif出现的次数
           # 基因集中的motif - 基因集中单个motif出现的次数
  ), ncol = 2, byrow = TRUE)
})




# 应用Fisher精确检验，并检查每个列联表是否至少是2x2的
fisher_results <- lapply(contingency_tables, function(table) {
  fisher.test(table)$p.value
})

# 将列表转换为向量
fisher_results <- unlist(fisher_results)
# 计算FDR
fdr_results <- p.adjust(fisher_results, method = "BH")
# 计算Enrichment score
enrichment_scores <- lapply(motif_names, function(motif) {
  if (background_motif_counts[[motif]] > 0) {
    log2((interested_motif_counts[[motif]]/nrow(interested_motifs)) / (background_motif_counts[[motif]]/nrow(significant_motifs)))
  } else {
    NA  # 如果背景中的motif数量为0，则无法计算富集度得分，返回NA
  }
})
enrichment_scores <- unlist(enrichment_scores)
# 将结果转换为数据框，并排序
enrichment_results_df <- data.frame(
  Motif = motif_names,
  Enrichment.Score = enrichment_scores,
  Pvalue = fisher_results,
  FDR = fdr_results
)
sorted_enrichment_results <- enrichment_results_df[order(enrichment_results_df$FDR), ]

# 筛选结果，注释motif信息
filter_enrichment_results <- filter(enrichment_results_df,Pvalue<0.05,FDR<0.05,enrichment_scores>0)
motif_anno <- read.delim('meme_info.txt')
anno_enrichment_result <- merge(filter_enrichment_results,motif_anno,by='Motif')
ALL_anno_enrichment_result <- merge(sorted_enrichment_results,motif_anno,by='Motif')


write.table(ALL_anno_enrichment_result, 'D38_up_phase3.txt', sep = '\t', row.names = FALSE, quote = FALSE)


interested_genes <- read.delim('time_point_gene_list.txt',header=1) %>% 
  filter(trait=='D38',class=='up',phase==4) %>% 
  select('gene_id')# 替换为你的感兴趣基因
interested_genes <- interested_genes$gene_id

# 获取感兴趣的基因集中的motif
interested_motifs <- significant_motifs[significant_motifs$sequence_name %in% interested_genes, ]
# 计算感兴趣基因集中motif的数量
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算所有motif的总数
total_motif_counts <- table(significant_motifs$motif_id,exclude = NULL)
# 计算在sig中出现，但不在insterest中出现的motif
x <- interested_motifs %>% group_by(motif_id) %>% summarise(n=n())
y <- significant_motifs %>% group_by(motif_id) %>% summarise(n=n())
z <- y$motif_id %in% x$motif_id
t <- y$motif_id[!z]
# 将z添加值intersted_motif_counts中
for (non_motif in t) {
  interested_motifs <- interested_motifs %>% 
    add_row(motif_id=non_motif)
}
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算背景基因集中motif的数量（所有motif减去感兴趣的motif）
background_motif_counts <- total_motif_counts - interested_motif_counts

# 过滤掉在感兴趣基因集或背景基因集中计数为0的motif
motif_names <- names(interested_motif_counts)[interested_motif_counts > 0 & background_motif_counts > 0]



# 创建2x2列联表矩阵
# d <- data.frame(gene.not.interest=c(M-k, N-M-n+k), gene.in.interest=c(k, n-k))
# row.names(d) <- c("In_category", "not_in_category")
contingency_tables <- lapply(motif_names, function(motif) {
  matrix(c(interested_motif_counts[[motif]], # 背景，不同motif出现的数目
           background_motif_counts[[motif]]-interested_motif_counts[[motif]], # 感兴趣基因集中motif数目
           nrow(interested_motifs) - interested_motif_counts[[motif]],
           nrow(significant_motifs)  - nrow(interested_motifs) -  background_motif_counts[[motif]]#background_motif_counts[[motif]] + interested_motif_counts[[motif]] # 背景全部motif - 背景中单个motif出现的次数
           # 基因集中的motif - 基因集中单个motif出现的次数
  ), ncol = 2, byrow = TRUE)
})




# 应用Fisher精确检验，并检查每个列联表是否至少是2x2的
fisher_results <- lapply(contingency_tables, function(table) {
  fisher.test(table)$p.value
})

# 将列表转换为向量
fisher_results <- unlist(fisher_results)
# 计算FDR
fdr_results <- p.adjust(fisher_results, method = "BH")
# 计算Enrichment score
enrichment_scores <- lapply(motif_names, function(motif) {
  if (background_motif_counts[[motif]] > 0) {
    log2((interested_motif_counts[[motif]]/nrow(interested_motifs)) / (background_motif_counts[[motif]]/nrow(significant_motifs)))
  } else {
    NA  # 如果背景中的motif数量为0，则无法计算富集度得分，返回NA
  }
})
enrichment_scores <- unlist(enrichment_scores)
# 将结果转换为数据框，并排序
enrichment_results_df <- data.frame(
  Motif = motif_names,
  Enrichment.Score = enrichment_scores,
  Pvalue = fisher_results,
  FDR = fdr_results
)
sorted_enrichment_results <- enrichment_results_df[order(enrichment_results_df$FDR), ]

# 筛选结果，注释motif信息
filter_enrichment_results <- filter(enrichment_results_df,Pvalue<0.05,FDR<0.05,enrichment_scores>0)
motif_anno <- read.delim('meme_info.txt')
anno_enrichment_result <- merge(filter_enrichment_results,motif_anno,by='Motif')
ALL_anno_enrichment_result <- merge(sorted_enrichment_results,motif_anno,by='Motif')


write.table(ALL_anno_enrichment_result, 'D38_up_phase4.txt', sep = '\t', row.names = FALSE, quote = FALSE)


interested_genes <- read.delim('time_point_gene_list.txt',header=1) %>% 
  filter(trait=='D38',class=='up',phase==5) %>% 
  select('gene_id')# 替换为你的感兴趣基因
interested_genes <- interested_genes$gene_id

# 获取感兴趣的基因集中的motif
interested_motifs <- significant_motifs[significant_motifs$sequence_name %in% interested_genes, ]
# 计算感兴趣基因集中motif的数量
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算所有motif的总数
total_motif_counts <- table(significant_motifs$motif_id,exclude = NULL)
# 计算在sig中出现，但不在insterest中出现的motif
x <- interested_motifs %>% group_by(motif_id) %>% summarise(n=n())
y <- significant_motifs %>% group_by(motif_id) %>% summarise(n=n())
z <- y$motif_id %in% x$motif_id
t <- y$motif_id[!z]
# 将z添加值intersted_motif_counts中
for (non_motif in t) {
  interested_motifs <- interested_motifs %>% 
    add_row(motif_id=non_motif)
}
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算背景基因集中motif的数量（所有motif减去感兴趣的motif）
background_motif_counts <- total_motif_counts - interested_motif_counts

# 过滤掉在感兴趣基因集或背景基因集中计数为0的motif
motif_names <- names(interested_motif_counts)[interested_motif_counts > 0 & background_motif_counts > 0]



# 创建2x2列联表矩阵
# d <- data.frame(gene.not.interest=c(M-k, N-M-n+k), gene.in.interest=c(k, n-k))
# row.names(d) <- c("In_category", "not_in_category")
contingency_tables <- lapply(motif_names, function(motif) {
  matrix(c(interested_motif_counts[[motif]], # 背景，不同motif出现的数目
           background_motif_counts[[motif]]-interested_motif_counts[[motif]], # 感兴趣基因集中motif数目
           nrow(interested_motifs) - interested_motif_counts[[motif]],
           nrow(significant_motifs)  - nrow(interested_motifs) -  background_motif_counts[[motif]]#background_motif_counts[[motif]] + interested_motif_counts[[motif]] # 背景全部motif - 背景中单个motif出现的次数
           # 基因集中的motif - 基因集中单个motif出现的次数
  ), ncol = 2, byrow = TRUE)
})




# 应用Fisher精确检验，并检查每个列联表是否至少是2x2的
fisher_results <- lapply(contingency_tables, function(table) {
  fisher.test(table)$p.value
})

# 将列表转换为向量
fisher_results <- unlist(fisher_results)
# 计算FDR
fdr_results <- p.adjust(fisher_results, method = "BH")
# 计算Enrichment score
enrichment_scores <- lapply(motif_names, function(motif) {
  if (background_motif_counts[[motif]] > 0) {
    log2((interested_motif_counts[[motif]]/nrow(interested_motifs)) / (background_motif_counts[[motif]]/nrow(significant_motifs)))
  } else {
    NA  # 如果背景中的motif数量为0，则无法计算富集度得分，返回NA
  }
})
enrichment_scores <- unlist(enrichment_scores)
# 将结果转换为数据框，并排序
enrichment_results_df <- data.frame(
  Motif = motif_names,
  Enrichment.Score = enrichment_scores,
  Pvalue = fisher_results,
  FDR = fdr_results
)
sorted_enrichment_results <- enrichment_results_df[order(enrichment_results_df$FDR), ]

# 筛选结果，注释motif信息
filter_enrichment_results <- filter(enrichment_results_df,Pvalue<0.05,FDR<0.05,enrichment_scores>0)
motif_anno <- read.delim('meme_info.txt')
anno_enrichment_result <- merge(filter_enrichment_results,motif_anno,by='Motif')
ALL_anno_enrichment_result <- merge(sorted_enrichment_results,motif_anno,by='Motif')


write.table(ALL_anno_enrichment_result, 'D38_up_phase5.txt', sep = '\t', row.names = FALSE, quote = FALSE)


interested_genes <- read.delim('time_point_gene_list.txt',header=1) %>% 
  filter(trait=='D38',class=='down',phase==1) %>% 
  select('gene_id')# 替换为你的感兴趣基因
interested_genes <- interested_genes$gene_id

# 获取感兴趣的基因集中的motif
interested_motifs <- significant_motifs[significant_motifs$sequence_name %in% interested_genes, ]
# 计算感兴趣基因集中motif的数量
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算所有motif的总数
total_motif_counts <- table(significant_motifs$motif_id,exclude = NULL)
# 计算在sig中出现，但不在insterest中出现的motif
x <- interested_motifs %>% group_by(motif_id) %>% summarise(n=n())
y <- significant_motifs %>% group_by(motif_id) %>% summarise(n=n())
z <- y$motif_id %in% x$motif_id
t <- y$motif_id[!z]
# 将z添加值intersted_motif_counts中
for (non_motif in t) {
  interested_motifs <- interested_motifs %>% 
    add_row(motif_id=non_motif)
}
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算背景基因集中motif的数量（所有motif减去感兴趣的motif）
background_motif_counts <- total_motif_counts - interested_motif_counts

# 过滤掉在感兴趣基因集或背景基因集中计数为0的motif
motif_names <- names(interested_motif_counts)[interested_motif_counts > 0 & background_motif_counts > 0]



# 创建2x2列联表矩阵
# d <- data.frame(gene.not.interest=c(M-k, N-M-n+k), gene.in.interest=c(k, n-k))
# row.names(d) <- c("In_category", "not_in_category")
contingency_tables <- lapply(motif_names, function(motif) {
  matrix(c(interested_motif_counts[[motif]], # 背景，不同motif出现的数目
           background_motif_counts[[motif]]-interested_motif_counts[[motif]], # 感兴趣基因集中motif数目
           nrow(interested_motifs) - interested_motif_counts[[motif]],
           nrow(significant_motifs)  - nrow(interested_motifs) -  background_motif_counts[[motif]]#background_motif_counts[[motif]] + interested_motif_counts[[motif]] # 背景全部motif - 背景中单个motif出现的次数
           # 基因集中的motif - 基因集中单个motif出现的次数
  ), ncol = 2, byrow = TRUE)
})




# 应用Fisher精确检验，并检查每个列联表是否至少是2x2的
fisher_results <- lapply(contingency_tables, function(table) {
  fisher.test(table)$p.value
})

# 将列表转换为向量
fisher_results <- unlist(fisher_results)
# 计算FDR
fdr_results <- p.adjust(fisher_results, method = "BH")
# 计算Enrichment score
enrichment_scores <- lapply(motif_names, function(motif) {
  if (background_motif_counts[[motif]] > 0) {
    log2((interested_motif_counts[[motif]]/nrow(interested_motifs)) / (background_motif_counts[[motif]]/nrow(significant_motifs)))
  } else {
    NA  # 如果背景中的motif数量为0，则无法计算富集度得分，返回NA
  }
})
enrichment_scores <- unlist(enrichment_scores)
# 将结果转换为数据框，并排序
enrichment_results_df <- data.frame(
  Motif = motif_names,
  Enrichment.Score = enrichment_scores,
  Pvalue = fisher_results,
  FDR = fdr_results
)
sorted_enrichment_results <- enrichment_results_df[order(enrichment_results_df$FDR), ]

# 筛选结果，注释motif信息
filter_enrichment_results <- filter(enrichment_results_df,Pvalue<0.05,FDR<0.05,enrichment_scores>0)
motif_anno <- read.delim('meme_info.txt')
anno_enrichment_result <- merge(filter_enrichment_results,motif_anno,by='Motif')
ALL_anno_enrichment_result <- merge(sorted_enrichment_results,motif_anno,by='Motif')


write.table(ALL_anno_enrichment_result, 'D38_down_phase1.txt', sep = '\t', row.names = FALSE, quote = FALSE)


interested_genes <- read.delim('time_point_gene_list.txt',header=1) %>% 
  filter(trait=='D38',class=='down',phase==2) %>% 
  select('gene_id')# 替换为你的感兴趣基因
interested_genes <- interested_genes$gene_id

# 获取感兴趣的基因集中的motif
interested_motifs <- significant_motifs[significant_motifs$sequence_name %in% interested_genes, ]
# 计算感兴趣基因集中motif的数量
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算所有motif的总数
total_motif_counts <- table(significant_motifs$motif_id,exclude = NULL)
# 计算在sig中出现，但不在insterest中出现的motif
x <- interested_motifs %>% group_by(motif_id) %>% summarise(n=n())
y <- significant_motifs %>% group_by(motif_id) %>% summarise(n=n())
z <- y$motif_id %in% x$motif_id
t <- y$motif_id[!z]
# 将z添加值intersted_motif_counts中
for (non_motif in t) {
  interested_motifs <- interested_motifs %>% 
    add_row(motif_id=non_motif)
}
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算背景基因集中motif的数量（所有motif减去感兴趣的motif）
background_motif_counts <- total_motif_counts - interested_motif_counts

# 过滤掉在感兴趣基因集或背景基因集中计数为0的motif
motif_names <- names(interested_motif_counts)[interested_motif_counts > 0 & background_motif_counts > 0]



# 创建2x2列联表矩阵
# d <- data.frame(gene.not.interest=c(M-k, N-M-n+k), gene.in.interest=c(k, n-k))
# row.names(d) <- c("In_category", "not_in_category")
contingency_tables <- lapply(motif_names, function(motif) {
  matrix(c(interested_motif_counts[[motif]], # 背景，不同motif出现的数目
           background_motif_counts[[motif]]-interested_motif_counts[[motif]], # 感兴趣基因集中motif数目
           nrow(interested_motifs) - interested_motif_counts[[motif]],
           nrow(significant_motifs)  - nrow(interested_motifs) -  background_motif_counts[[motif]]#background_motif_counts[[motif]] + interested_motif_counts[[motif]] # 背景全部motif - 背景中单个motif出现的次数
           # 基因集中的motif - 基因集中单个motif出现的次数
  ), ncol = 2, byrow = TRUE)
})




# 应用Fisher精确检验，并检查每个列联表是否至少是2x2的
fisher_results <- lapply(contingency_tables, function(table) {
  fisher.test(table)$p.value
})

# 将列表转换为向量
fisher_results <- unlist(fisher_results)
# 计算FDR
fdr_results <- p.adjust(fisher_results, method = "BH")
# 计算Enrichment score
enrichment_scores <- lapply(motif_names, function(motif) {
  if (background_motif_counts[[motif]] > 0) {
    log2((interested_motif_counts[[motif]]/nrow(interested_motifs)) / (background_motif_counts[[motif]]/nrow(significant_motifs)))
  } else {
    NA  # 如果背景中的motif数量为0，则无法计算富集度得分，返回NA
  }
})
enrichment_scores <- unlist(enrichment_scores)
# 将结果转换为数据框，并排序
enrichment_results_df <- data.frame(
  Motif = motif_names,
  Enrichment.Score = enrichment_scores,
  Pvalue = fisher_results,
  FDR = fdr_results
)
sorted_enrichment_results <- enrichment_results_df[order(enrichment_results_df$FDR), ]

# 筛选结果，注释motif信息
filter_enrichment_results <- filter(enrichment_results_df,Pvalue<0.05,FDR<0.05,enrichment_scores>0)
motif_anno <- read.delim('meme_info.txt')
anno_enrichment_result <- merge(filter_enrichment_results,motif_anno,by='Motif')
ALL_anno_enrichment_result <- merge(sorted_enrichment_results,motif_anno,by='Motif')


write.table(ALL_anno_enrichment_result, 'D38_down_phase2.txt', sep = '\t', row.names = FALSE, quote = FALSE)


interested_genes <- read.delim('time_point_gene_list.txt',header=1) %>% 
  filter(trait=='D38',class=='down',phase==3) %>% 
  select('gene_id')# 替换为你的感兴趣基因
interested_genes <- interested_genes$gene_id

# 获取感兴趣的基因集中的motif
interested_motifs <- significant_motifs[significant_motifs$sequence_name %in% interested_genes, ]
# 计算感兴趣基因集中motif的数量
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算所有motif的总数
total_motif_counts <- table(significant_motifs$motif_id,exclude = NULL)
# 计算在sig中出现，但不在insterest中出现的motif
x <- interested_motifs %>% group_by(motif_id) %>% summarise(n=n())
y <- significant_motifs %>% group_by(motif_id) %>% summarise(n=n())
z <- y$motif_id %in% x$motif_id
t <- y$motif_id[!z]
# 将z添加值intersted_motif_counts中
for (non_motif in t) {
  interested_motifs <- interested_motifs %>% 
    add_row(motif_id=non_motif)
}
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算背景基因集中motif的数量（所有motif减去感兴趣的motif）
background_motif_counts <- total_motif_counts - interested_motif_counts

# 过滤掉在感兴趣基因集或背景基因集中计数为0的motif
motif_names <- names(interested_motif_counts)[interested_motif_counts > 0 & background_motif_counts > 0]



# 创建2x2列联表矩阵
# d <- data.frame(gene.not.interest=c(M-k, N-M-n+k), gene.in.interest=c(k, n-k))
# row.names(d) <- c("In_category", "not_in_category")
contingency_tables <- lapply(motif_names, function(motif) {
  matrix(c(interested_motif_counts[[motif]], # 背景，不同motif出现的数目
           background_motif_counts[[motif]]-interested_motif_counts[[motif]], # 感兴趣基因集中motif数目
           nrow(interested_motifs) - interested_motif_counts[[motif]],
           nrow(significant_motifs)  - nrow(interested_motifs) -  background_motif_counts[[motif]]#background_motif_counts[[motif]] + interested_motif_counts[[motif]] # 背景全部motif - 背景中单个motif出现的次数
           # 基因集中的motif - 基因集中单个motif出现的次数
  ), ncol = 2, byrow = TRUE)
})




# 应用Fisher精确检验，并检查每个列联表是否至少是2x2的
fisher_results <- lapply(contingency_tables, function(table) {
  fisher.test(table)$p.value
})

# 将列表转换为向量
fisher_results <- unlist(fisher_results)
# 计算FDR
fdr_results <- p.adjust(fisher_results, method = "BH")
# 计算Enrichment score
enrichment_scores <- lapply(motif_names, function(motif) {
  if (background_motif_counts[[motif]] > 0) {
    log2((interested_motif_counts[[motif]]/nrow(interested_motifs)) / (background_motif_counts[[motif]]/nrow(significant_motifs)))
  } else {
    NA  # 如果背景中的motif数量为0，则无法计算富集度得分，返回NA
  }
})
enrichment_scores <- unlist(enrichment_scores)
# 将结果转换为数据框，并排序
enrichment_results_df <- data.frame(
  Motif = motif_names,
  Enrichment.Score = enrichment_scores,
  Pvalue = fisher_results,
  FDR = fdr_results
)
sorted_enrichment_results <- enrichment_results_df[order(enrichment_results_df$FDR), ]

# 筛选结果，注释motif信息
filter_enrichment_results <- filter(enrichment_results_df,Pvalue<0.05,FDR<0.05,enrichment_scores>0)
motif_anno <- read.delim('meme_info.txt')
anno_enrichment_result <- merge(filter_enrichment_results,motif_anno,by='Motif')
ALL_anno_enrichment_result <- merge(sorted_enrichment_results,motif_anno,by='Motif')


write.table(ALL_anno_enrichment_result, 'D38_down_phase3.txt', sep = '\t', row.names = FALSE, quote = FALSE)


interested_genes <- read.delim('time_point_gene_list.txt',header=1) %>% 
  filter(trait=='D38',class=='down',phase==4) %>% 
  select('gene_id')# 替换为你的感兴趣基因
interested_genes <- interested_genes$gene_id

# 获取感兴趣的基因集中的motif
interested_motifs <- significant_motifs[significant_motifs$sequence_name %in% interested_genes, ]
# 计算感兴趣基因集中motif的数量
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算所有motif的总数
total_motif_counts <- table(significant_motifs$motif_id,exclude = NULL)
# 计算在sig中出现，但不在insterest中出现的motif
x <- interested_motifs %>% group_by(motif_id) %>% summarise(n=n())
y <- significant_motifs %>% group_by(motif_id) %>% summarise(n=n())
z <- y$motif_id %in% x$motif_id
t <- y$motif_id[!z]
# 将z添加值intersted_motif_counts中
for (non_motif in t) {
  interested_motifs <- interested_motifs %>% 
    add_row(motif_id=non_motif)
}
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算背景基因集中motif的数量（所有motif减去感兴趣的motif）
background_motif_counts <- total_motif_counts - interested_motif_counts

# 过滤掉在感兴趣基因集或背景基因集中计数为0的motif
motif_names <- names(interested_motif_counts)[interested_motif_counts > 0 & background_motif_counts > 0]



# 创建2x2列联表矩阵
# d <- data.frame(gene.not.interest=c(M-k, N-M-n+k), gene.in.interest=c(k, n-k))
# row.names(d) <- c("In_category", "not_in_category")
contingency_tables <- lapply(motif_names, function(motif) {
  matrix(c(interested_motif_counts[[motif]], # 背景，不同motif出现的数目
           background_motif_counts[[motif]]-interested_motif_counts[[motif]], # 感兴趣基因集中motif数目
           nrow(interested_motifs) - interested_motif_counts[[motif]],
           nrow(significant_motifs)  - nrow(interested_motifs) -  background_motif_counts[[motif]]#background_motif_counts[[motif]] + interested_motif_counts[[motif]] # 背景全部motif - 背景中单个motif出现的次数
           # 基因集中的motif - 基因集中单个motif出现的次数
  ), ncol = 2, byrow = TRUE)
})




# 应用Fisher精确检验，并检查每个列联表是否至少是2x2的
fisher_results <- lapply(contingency_tables, function(table) {
  fisher.test(table)$p.value
})

# 将列表转换为向量
fisher_results <- unlist(fisher_results)
# 计算FDR
fdr_results <- p.adjust(fisher_results, method = "BH")
# 计算Enrichment score
enrichment_scores <- lapply(motif_names, function(motif) {
  if (background_motif_counts[[motif]] > 0) {
    log2((interested_motif_counts[[motif]]/nrow(interested_motifs)) / (background_motif_counts[[motif]]/nrow(significant_motifs)))
  } else {
    NA  # 如果背景中的motif数量为0，则无法计算富集度得分，返回NA
  }
})
enrichment_scores <- unlist(enrichment_scores)
# 将结果转换为数据框，并排序
enrichment_results_df <- data.frame(
  Motif = motif_names,
  Enrichment.Score = enrichment_scores,
  Pvalue = fisher_results,
  FDR = fdr_results
)
sorted_enrichment_results <- enrichment_results_df[order(enrichment_results_df$FDR), ]

# 筛选结果，注释motif信息
filter_enrichment_results <- filter(enrichment_results_df,Pvalue<0.05,FDR<0.05,enrichment_scores>0)
motif_anno <- read.delim('meme_info.txt')
anno_enrichment_result <- merge(filter_enrichment_results,motif_anno,by='Motif')
ALL_anno_enrichment_result <- merge(sorted_enrichment_results,motif_anno,by='Motif')


write.table(ALL_anno_enrichment_result, 'D38_down_phase4.txt', sep = '\t', row.names = FALSE, quote = FALSE)


interested_genes <- read.delim('time_point_gene_list.txt',header=1) %>% 
  filter(trait=='D38',class=='down',phase==5) %>% 
  select('gene_id')# 替换为你的感兴趣基因
interested_genes <- interested_genes$gene_id

# 获取感兴趣的基因集中的motif
interested_motifs <- significant_motifs[significant_motifs$sequence_name %in% interested_genes, ]
# 计算感兴趣基因集中motif的数量
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算所有motif的总数
total_motif_counts <- table(significant_motifs$motif_id,exclude = NULL)
# 计算在sig中出现，但不在insterest中出现的motif
x <- interested_motifs %>% group_by(motif_id) %>% summarise(n=n())
y <- significant_motifs %>% group_by(motif_id) %>% summarise(n=n())
z <- y$motif_id %in% x$motif_id
t <- y$motif_id[!z]
# 将z添加值intersted_motif_counts中
for (non_motif in t) {
  interested_motifs <- interested_motifs %>% 
    add_row(motif_id=non_motif)
}
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算背景基因集中motif的数量（所有motif减去感兴趣的motif）
background_motif_counts <- total_motif_counts - interested_motif_counts

# 过滤掉在感兴趣基因集或背景基因集中计数为0的motif
motif_names <- names(interested_motif_counts)[interested_motif_counts > 0 & background_motif_counts > 0]



# 创建2x2列联表矩阵
# d <- data.frame(gene.not.interest=c(M-k, N-M-n+k), gene.in.interest=c(k, n-k))
# row.names(d) <- c("In_category", "not_in_category")
contingency_tables <- lapply(motif_names, function(motif) {
  matrix(c(interested_motif_counts[[motif]], # 背景，不同motif出现的数目
           background_motif_counts[[motif]]-interested_motif_counts[[motif]], # 感兴趣基因集中motif数目
           nrow(interested_motifs) - interested_motif_counts[[motif]],
           nrow(significant_motifs)  - nrow(interested_motifs) -  background_motif_counts[[motif]]#background_motif_counts[[motif]] + interested_motif_counts[[motif]] # 背景全部motif - 背景中单个motif出现的次数
           # 基因集中的motif - 基因集中单个motif出现的次数
  ), ncol = 2, byrow = TRUE)
})




# 应用Fisher精确检验，并检查每个列联表是否至少是2x2的
fisher_results <- lapply(contingency_tables, function(table) {
  fisher.test(table)$p.value
})

# 将列表转换为向量
fisher_results <- unlist(fisher_results)
# 计算FDR
fdr_results <- p.adjust(fisher_results, method = "BH")
# 计算Enrichment score
enrichment_scores <- lapply(motif_names, function(motif) {
  if (background_motif_counts[[motif]] > 0) {
    log2((interested_motif_counts[[motif]]/nrow(interested_motifs)) / (background_motif_counts[[motif]]/nrow(significant_motifs)))
  } else {
    NA  # 如果背景中的motif数量为0，则无法计算富集度得分，返回NA
  }
})
enrichment_scores <- unlist(enrichment_scores)
# 将结果转换为数据框，并排序
enrichment_results_df <- data.frame(
  Motif = motif_names,
  Enrichment.Score = enrichment_scores,
  Pvalue = fisher_results,
  FDR = fdr_results
)
sorted_enrichment_results <- enrichment_results_df[order(enrichment_results_df$FDR), ]

# 筛选结果，注释motif信息
filter_enrichment_results <- filter(enrichment_results_df,Pvalue<0.05,FDR<0.05,enrichment_scores>0)
motif_anno <- read.delim('meme_info.txt')
anno_enrichment_result <- merge(filter_enrichment_results,motif_anno,by='Motif')
ALL_anno_enrichment_result <- merge(sorted_enrichment_results,motif_anno,by='Motif')


write.table(ALL_anno_enrichment_result, 'D38_down_phase5.txt', sep = '\t', row.names = FALSE, quote = FALSE)


interested_genes <- read.delim('time_point_gene_list.txt',header=1) %>% 
  filter(trait=='D45',class=='up',phase==1) %>% 
  select('gene_id')# 替换为你的感兴趣基因
interested_genes <- interested_genes$gene_id

# 获取感兴趣的基因集中的motif
interested_motifs <- significant_motifs[significant_motifs$sequence_name %in% interested_genes, ]
# 计算感兴趣基因集中motif的数量
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算所有motif的总数
total_motif_counts <- table(significant_motifs$motif_id,exclude = NULL)
# 计算在sig中出现，但不在insterest中出现的motif
x <- interested_motifs %>% group_by(motif_id) %>% summarise(n=n())
y <- significant_motifs %>% group_by(motif_id) %>% summarise(n=n())
z <- y$motif_id %in% x$motif_id
t <- y$motif_id[!z]
# 将z添加值intersted_motif_counts中
for (non_motif in t) {
  interested_motifs <- interested_motifs %>% 
    add_row(motif_id=non_motif)
}
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算背景基因集中motif的数量（所有motif减去感兴趣的motif）
background_motif_counts <- total_motif_counts - interested_motif_counts

# 过滤掉在感兴趣基因集或背景基因集中计数为0的motif
motif_names <- names(interested_motif_counts)[interested_motif_counts > 0 & background_motif_counts > 0]



# 创建2x2列联表矩阵
# d <- data.frame(gene.not.interest=c(M-k, N-M-n+k), gene.in.interest=c(k, n-k))
# row.names(d) <- c("In_category", "not_in_category")
contingency_tables <- lapply(motif_names, function(motif) {
  matrix(c(interested_motif_counts[[motif]], # 背景，不同motif出现的数目
           background_motif_counts[[motif]]-interested_motif_counts[[motif]], # 感兴趣基因集中motif数目
           nrow(interested_motifs) - interested_motif_counts[[motif]],
           nrow(significant_motifs)  - nrow(interested_motifs) -  background_motif_counts[[motif]]#background_motif_counts[[motif]] + interested_motif_counts[[motif]] # 背景全部motif - 背景中单个motif出现的次数
           # 基因集中的motif - 基因集中单个motif出现的次数
  ), ncol = 2, byrow = TRUE)
})




# 应用Fisher精确检验，并检查每个列联表是否至少是2x2的
fisher_results <- lapply(contingency_tables, function(table) {
  fisher.test(table)$p.value
})

# 将列表转换为向量
fisher_results <- unlist(fisher_results)
# 计算FDR
fdr_results <- p.adjust(fisher_results, method = "BH")
# 计算Enrichment score
enrichment_scores <- lapply(motif_names, function(motif) {
  if (background_motif_counts[[motif]] > 0) {
    log2((interested_motif_counts[[motif]]/nrow(interested_motifs)) / (background_motif_counts[[motif]]/nrow(significant_motifs)))
  } else {
    NA  # 如果背景中的motif数量为0，则无法计算富集度得分，返回NA
  }
})
enrichment_scores <- unlist(enrichment_scores)
# 将结果转换为数据框，并排序
enrichment_results_df <- data.frame(
  Motif = motif_names,
  Enrichment.Score = enrichment_scores,
  Pvalue = fisher_results,
  FDR = fdr_results
)
sorted_enrichment_results <- enrichment_results_df[order(enrichment_results_df$FDR), ]

# 筛选结果，注释motif信息
filter_enrichment_results <- filter(enrichment_results_df,Pvalue<0.05,FDR<0.05,enrichment_scores>0)
motif_anno <- read.delim('meme_info.txt')
anno_enrichment_result <- merge(filter_enrichment_results,motif_anno,by='Motif')
ALL_anno_enrichment_result <- merge(sorted_enrichment_results,motif_anno,by='Motif')


write.table(ALL_anno_enrichment_result, 'D45_up_phase1.txt', sep = '\t', row.names = FALSE, quote = FALSE)


interested_genes <- read.delim('time_point_gene_list.txt',header=1) %>% 
  filter(trait=='D45',class=='up',phase==2) %>% 
  select('gene_id')# 替换为你的感兴趣基因
interested_genes <- interested_genes$gene_id

# 获取感兴趣的基因集中的motif
interested_motifs <- significant_motifs[significant_motifs$sequence_name %in% interested_genes, ]
# 计算感兴趣基因集中motif的数量
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算所有motif的总数
total_motif_counts <- table(significant_motifs$motif_id,exclude = NULL)
# 计算在sig中出现，但不在insterest中出现的motif
x <- interested_motifs %>% group_by(motif_id) %>% summarise(n=n())
y <- significant_motifs %>% group_by(motif_id) %>% summarise(n=n())
z <- y$motif_id %in% x$motif_id
t <- y$motif_id[!z]
# 将z添加值intersted_motif_counts中
for (non_motif in t) {
  interested_motifs <- interested_motifs %>% 
    add_row(motif_id=non_motif)
}
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算背景基因集中motif的数量（所有motif减去感兴趣的motif）
background_motif_counts <- total_motif_counts - interested_motif_counts

# 过滤掉在感兴趣基因集或背景基因集中计数为0的motif
motif_names <- names(interested_motif_counts)[interested_motif_counts > 0 & background_motif_counts > 0]



# 创建2x2列联表矩阵
# d <- data.frame(gene.not.interest=c(M-k, N-M-n+k), gene.in.interest=c(k, n-k))
# row.names(d) <- c("In_category", "not_in_category")
contingency_tables <- lapply(motif_names, function(motif) {
  matrix(c(interested_motif_counts[[motif]], # 背景，不同motif出现的数目
           background_motif_counts[[motif]]-interested_motif_counts[[motif]], # 感兴趣基因集中motif数目
           nrow(interested_motifs) - interested_motif_counts[[motif]],
           nrow(significant_motifs)  - nrow(interested_motifs) -  background_motif_counts[[motif]]#background_motif_counts[[motif]] + interested_motif_counts[[motif]] # 背景全部motif - 背景中单个motif出现的次数
           # 基因集中的motif - 基因集中单个motif出现的次数
  ), ncol = 2, byrow = TRUE)
})




# 应用Fisher精确检验，并检查每个列联表是否至少是2x2的
fisher_results <- lapply(contingency_tables, function(table) {
  fisher.test(table)$p.value
})

# 将列表转换为向量
fisher_results <- unlist(fisher_results)
# 计算FDR
fdr_results <- p.adjust(fisher_results, method = "BH")
# 计算Enrichment score
enrichment_scores <- lapply(motif_names, function(motif) {
  if (background_motif_counts[[motif]] > 0) {
    log2((interested_motif_counts[[motif]]/nrow(interested_motifs)) / (background_motif_counts[[motif]]/nrow(significant_motifs)))
  } else {
    NA  # 如果背景中的motif数量为0，则无法计算富集度得分，返回NA
  }
})
enrichment_scores <- unlist(enrichment_scores)
# 将结果转换为数据框，并排序
enrichment_results_df <- data.frame(
  Motif = motif_names,
  Enrichment.Score = enrichment_scores,
  Pvalue = fisher_results,
  FDR = fdr_results
)
sorted_enrichment_results <- enrichment_results_df[order(enrichment_results_df$FDR), ]

# 筛选结果，注释motif信息
filter_enrichment_results <- filter(enrichment_results_df,Pvalue<0.05,FDR<0.05,enrichment_scores>0)
motif_anno <- read.delim('meme_info.txt')
anno_enrichment_result <- merge(filter_enrichment_results,motif_anno,by='Motif')
ALL_anno_enrichment_result <- merge(sorted_enrichment_results,motif_anno,by='Motif')


write.table(ALL_anno_enrichment_result, 'D45_up_phase2.txt', sep = '\t', row.names = FALSE, quote = FALSE)


interested_genes <- read.delim('time_point_gene_list.txt',header=1) %>% 
  filter(trait=='D45',class=='up',phase==3) %>% 
  select('gene_id')# 替换为你的感兴趣基因
interested_genes <- interested_genes$gene_id

# 获取感兴趣的基因集中的motif
interested_motifs <- significant_motifs[significant_motifs$sequence_name %in% interested_genes, ]
# 计算感兴趣基因集中motif的数量
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算所有motif的总数
total_motif_counts <- table(significant_motifs$motif_id,exclude = NULL)
# 计算在sig中出现，但不在insterest中出现的motif
x <- interested_motifs %>% group_by(motif_id) %>% summarise(n=n())
y <- significant_motifs %>% group_by(motif_id) %>% summarise(n=n())
z <- y$motif_id %in% x$motif_id
t <- y$motif_id[!z]
# 将z添加值intersted_motif_counts中
for (non_motif in t) {
  interested_motifs <- interested_motifs %>% 
    add_row(motif_id=non_motif)
}
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算背景基因集中motif的数量（所有motif减去感兴趣的motif）
background_motif_counts <- total_motif_counts - interested_motif_counts

# 过滤掉在感兴趣基因集或背景基因集中计数为0的motif
motif_names <- names(interested_motif_counts)[interested_motif_counts > 0 & background_motif_counts > 0]



# 创建2x2列联表矩阵
# d <- data.frame(gene.not.interest=c(M-k, N-M-n+k), gene.in.interest=c(k, n-k))
# row.names(d) <- c("In_category", "not_in_category")
contingency_tables <- lapply(motif_names, function(motif) {
  matrix(c(interested_motif_counts[[motif]], # 背景，不同motif出现的数目
           background_motif_counts[[motif]]-interested_motif_counts[[motif]], # 感兴趣基因集中motif数目
           nrow(interested_motifs) - interested_motif_counts[[motif]],
           nrow(significant_motifs)  - nrow(interested_motifs) -  background_motif_counts[[motif]]#background_motif_counts[[motif]] + interested_motif_counts[[motif]] # 背景全部motif - 背景中单个motif出现的次数
           # 基因集中的motif - 基因集中单个motif出现的次数
  ), ncol = 2, byrow = TRUE)
})




# 应用Fisher精确检验，并检查每个列联表是否至少是2x2的
fisher_results <- lapply(contingency_tables, function(table) {
  fisher.test(table)$p.value
})

# 将列表转换为向量
fisher_results <- unlist(fisher_results)
# 计算FDR
fdr_results <- p.adjust(fisher_results, method = "BH")
# 计算Enrichment score
enrichment_scores <- lapply(motif_names, function(motif) {
  if (background_motif_counts[[motif]] > 0) {
    log2((interested_motif_counts[[motif]]/nrow(interested_motifs)) / (background_motif_counts[[motif]]/nrow(significant_motifs)))
  } else {
    NA  # 如果背景中的motif数量为0，则无法计算富集度得分，返回NA
  }
})
enrichment_scores <- unlist(enrichment_scores)
# 将结果转换为数据框，并排序
enrichment_results_df <- data.frame(
  Motif = motif_names,
  Enrichment.Score = enrichment_scores,
  Pvalue = fisher_results,
  FDR = fdr_results
)
sorted_enrichment_results <- enrichment_results_df[order(enrichment_results_df$FDR), ]

# 筛选结果，注释motif信息
filter_enrichment_results <- filter(enrichment_results_df,Pvalue<0.05,FDR<0.05,enrichment_scores>0)
motif_anno <- read.delim('meme_info.txt')
anno_enrichment_result <- merge(filter_enrichment_results,motif_anno,by='Motif')
ALL_anno_enrichment_result <- merge(sorted_enrichment_results,motif_anno,by='Motif')


write.table(ALL_anno_enrichment_result, 'D45_up_phase3.txt', sep = '\t', row.names = FALSE, quote = FALSE)


interested_genes <- read.delim('time_point_gene_list.txt',header=1) %>% 
  filter(trait=='D45',class=='up',phase==4) %>% 
  select('gene_id')# 替换为你的感兴趣基因
interested_genes <- interested_genes$gene_id

# 获取感兴趣的基因集中的motif
interested_motifs <- significant_motifs[significant_motifs$sequence_name %in% interested_genes, ]
# 计算感兴趣基因集中motif的数量
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算所有motif的总数
total_motif_counts <- table(significant_motifs$motif_id,exclude = NULL)
# 计算在sig中出现，但不在insterest中出现的motif
x <- interested_motifs %>% group_by(motif_id) %>% summarise(n=n())
y <- significant_motifs %>% group_by(motif_id) %>% summarise(n=n())
z <- y$motif_id %in% x$motif_id
t <- y$motif_id[!z]
# 将z添加值intersted_motif_counts中
for (non_motif in t) {
  interested_motifs <- interested_motifs %>% 
    add_row(motif_id=non_motif)
}
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算背景基因集中motif的数量（所有motif减去感兴趣的motif）
background_motif_counts <- total_motif_counts - interested_motif_counts

# 过滤掉在感兴趣基因集或背景基因集中计数为0的motif
motif_names <- names(interested_motif_counts)[interested_motif_counts > 0 & background_motif_counts > 0]



# 创建2x2列联表矩阵
# d <- data.frame(gene.not.interest=c(M-k, N-M-n+k), gene.in.interest=c(k, n-k))
# row.names(d) <- c("In_category", "not_in_category")
contingency_tables <- lapply(motif_names, function(motif) {
  matrix(c(interested_motif_counts[[motif]], # 背景，不同motif出现的数目
           background_motif_counts[[motif]]-interested_motif_counts[[motif]], # 感兴趣基因集中motif数目
           nrow(interested_motifs) - interested_motif_counts[[motif]],
           nrow(significant_motifs)  - nrow(interested_motifs) -  background_motif_counts[[motif]]#background_motif_counts[[motif]] + interested_motif_counts[[motif]] # 背景全部motif - 背景中单个motif出现的次数
           # 基因集中的motif - 基因集中单个motif出现的次数
  ), ncol = 2, byrow = TRUE)
})




# 应用Fisher精确检验，并检查每个列联表是否至少是2x2的
fisher_results <- lapply(contingency_tables, function(table) {
  fisher.test(table)$p.value
})

# 将列表转换为向量
fisher_results <- unlist(fisher_results)
# 计算FDR
fdr_results <- p.adjust(fisher_results, method = "BH")
# 计算Enrichment score
enrichment_scores <- lapply(motif_names, function(motif) {
  if (background_motif_counts[[motif]] > 0) {
    log2((interested_motif_counts[[motif]]/nrow(interested_motifs)) / (background_motif_counts[[motif]]/nrow(significant_motifs)))
  } else {
    NA  # 如果背景中的motif数量为0，则无法计算富集度得分，返回NA
  }
})
enrichment_scores <- unlist(enrichment_scores)
# 将结果转换为数据框，并排序
enrichment_results_df <- data.frame(
  Motif = motif_names,
  Enrichment.Score = enrichment_scores,
  Pvalue = fisher_results,
  FDR = fdr_results
)
sorted_enrichment_results <- enrichment_results_df[order(enrichment_results_df$FDR), ]

# 筛选结果，注释motif信息
filter_enrichment_results <- filter(enrichment_results_df,Pvalue<0.05,FDR<0.05,enrichment_scores>0)
motif_anno <- read.delim('meme_info.txt')
anno_enrichment_result <- merge(filter_enrichment_results,motif_anno,by='Motif')
ALL_anno_enrichment_result <- merge(sorted_enrichment_results,motif_anno,by='Motif')


write.table(ALL_anno_enrichment_result, 'D45_up_phase4.txt', sep = '\t', row.names = FALSE, quote = FALSE)


interested_genes <- read.delim('time_point_gene_list.txt',header=1) %>% 
  filter(trait=='D45',class=='up',phase==5) %>% 
  select('gene_id')# 替换为你的感兴趣基因
interested_genes <- interested_genes$gene_id

# 获取感兴趣的基因集中的motif
interested_motifs <- significant_motifs[significant_motifs$sequence_name %in% interested_genes, ]
# 计算感兴趣基因集中motif的数量
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算所有motif的总数
total_motif_counts <- table(significant_motifs$motif_id,exclude = NULL)
# 计算在sig中出现，但不在insterest中出现的motif
x <- interested_motifs %>% group_by(motif_id) %>% summarise(n=n())
y <- significant_motifs %>% group_by(motif_id) %>% summarise(n=n())
z <- y$motif_id %in% x$motif_id
t <- y$motif_id[!z]
# 将z添加值intersted_motif_counts中
for (non_motif in t) {
  interested_motifs <- interested_motifs %>% 
    add_row(motif_id=non_motif)
}
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算背景基因集中motif的数量（所有motif减去感兴趣的motif）
background_motif_counts <- total_motif_counts - interested_motif_counts

# 过滤掉在感兴趣基因集或背景基因集中计数为0的motif
motif_names <- names(interested_motif_counts)[interested_motif_counts > 0 & background_motif_counts > 0]



# 创建2x2列联表矩阵
# d <- data.frame(gene.not.interest=c(M-k, N-M-n+k), gene.in.interest=c(k, n-k))
# row.names(d) <- c("In_category", "not_in_category")
contingency_tables <- lapply(motif_names, function(motif) {
  matrix(c(interested_motif_counts[[motif]], # 背景，不同motif出现的数目
           background_motif_counts[[motif]]-interested_motif_counts[[motif]], # 感兴趣基因集中motif数目
           nrow(interested_motifs) - interested_motif_counts[[motif]],
           nrow(significant_motifs)  - nrow(interested_motifs) -  background_motif_counts[[motif]]#background_motif_counts[[motif]] + interested_motif_counts[[motif]] # 背景全部motif - 背景中单个motif出现的次数
           # 基因集中的motif - 基因集中单个motif出现的次数
  ), ncol = 2, byrow = TRUE)
})




# 应用Fisher精确检验，并检查每个列联表是否至少是2x2的
fisher_results <- lapply(contingency_tables, function(table) {
  fisher.test(table)$p.value
})

# 将列表转换为向量
fisher_results <- unlist(fisher_results)
# 计算FDR
fdr_results <- p.adjust(fisher_results, method = "BH")
# 计算Enrichment score
enrichment_scores <- lapply(motif_names, function(motif) {
  if (background_motif_counts[[motif]] > 0) {
    log2((interested_motif_counts[[motif]]/nrow(interested_motifs)) / (background_motif_counts[[motif]]/nrow(significant_motifs)))
  } else {
    NA  # 如果背景中的motif数量为0，则无法计算富集度得分，返回NA
  }
})
enrichment_scores <- unlist(enrichment_scores)
# 将结果转换为数据框，并排序
enrichment_results_df <- data.frame(
  Motif = motif_names,
  Enrichment.Score = enrichment_scores,
  Pvalue = fisher_results,
  FDR = fdr_results
)
sorted_enrichment_results <- enrichment_results_df[order(enrichment_results_df$FDR), ]

# 筛选结果，注释motif信息
filter_enrichment_results <- filter(enrichment_results_df,Pvalue<0.05,FDR<0.05,enrichment_scores>0)
motif_anno <- read.delim('meme_info.txt')
anno_enrichment_result <- merge(filter_enrichment_results,motif_anno,by='Motif')
ALL_anno_enrichment_result <- merge(sorted_enrichment_results,motif_anno,by='Motif')


write.table(ALL_anno_enrichment_result, 'D45_up_phase5.txt', sep = '\t', row.names = FALSE, quote = FALSE)


interested_genes <- read.delim('time_point_gene_list.txt',header=1) %>% 
  filter(trait=='D45',class=='down',phase==1) %>% 
  select('gene_id')# 替换为你的感兴趣基因
interested_genes <- interested_genes$gene_id

# 获取感兴趣的基因集中的motif
interested_motifs <- significant_motifs[significant_motifs$sequence_name %in% interested_genes, ]
# 计算感兴趣基因集中motif的数量
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算所有motif的总数
total_motif_counts <- table(significant_motifs$motif_id,exclude = NULL)
# 计算在sig中出现，但不在insterest中出现的motif
x <- interested_motifs %>% group_by(motif_id) %>% summarise(n=n())
y <- significant_motifs %>% group_by(motif_id) %>% summarise(n=n())
z <- y$motif_id %in% x$motif_id
t <- y$motif_id[!z]
# 将z添加值intersted_motif_counts中
for (non_motif in t) {
  interested_motifs <- interested_motifs %>% 
    add_row(motif_id=non_motif)
}
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算背景基因集中motif的数量（所有motif减去感兴趣的motif）
background_motif_counts <- total_motif_counts - interested_motif_counts

# 过滤掉在感兴趣基因集或背景基因集中计数为0的motif
motif_names <- names(interested_motif_counts)[interested_motif_counts > 0 & background_motif_counts > 0]



# 创建2x2列联表矩阵
# d <- data.frame(gene.not.interest=c(M-k, N-M-n+k), gene.in.interest=c(k, n-k))
# row.names(d) <- c("In_category", "not_in_category")
contingency_tables <- lapply(motif_names, function(motif) {
  matrix(c(interested_motif_counts[[motif]], # 背景，不同motif出现的数目
           background_motif_counts[[motif]]-interested_motif_counts[[motif]], # 感兴趣基因集中motif数目
           nrow(interested_motifs) - interested_motif_counts[[motif]],
           nrow(significant_motifs)  - nrow(interested_motifs) -  background_motif_counts[[motif]]#background_motif_counts[[motif]] + interested_motif_counts[[motif]] # 背景全部motif - 背景中单个motif出现的次数
           # 基因集中的motif - 基因集中单个motif出现的次数
  ), ncol = 2, byrow = TRUE)
})




# 应用Fisher精确检验，并检查每个列联表是否至少是2x2的
fisher_results <- lapply(contingency_tables, function(table) {
  fisher.test(table)$p.value
})

# 将列表转换为向量
fisher_results <- unlist(fisher_results)
# 计算FDR
fdr_results <- p.adjust(fisher_results, method = "BH")
# 计算Enrichment score
enrichment_scores <- lapply(motif_names, function(motif) {
  if (background_motif_counts[[motif]] > 0) {
    log2((interested_motif_counts[[motif]]/nrow(interested_motifs)) / (background_motif_counts[[motif]]/nrow(significant_motifs)))
  } else {
    NA  # 如果背景中的motif数量为0，则无法计算富集度得分，返回NA
  }
})
enrichment_scores <- unlist(enrichment_scores)
# 将结果转换为数据框，并排序
enrichment_results_df <- data.frame(
  Motif = motif_names,
  Enrichment.Score = enrichment_scores,
  Pvalue = fisher_results,
  FDR = fdr_results
)
sorted_enrichment_results <- enrichment_results_df[order(enrichment_results_df$FDR), ]

# 筛选结果，注释motif信息
filter_enrichment_results <- filter(enrichment_results_df,Pvalue<0.05,FDR<0.05,enrichment_scores>0)
motif_anno <- read.delim('meme_info.txt')
anno_enrichment_result <- merge(filter_enrichment_results,motif_anno,by='Motif')
ALL_anno_enrichment_result <- merge(sorted_enrichment_results,motif_anno,by='Motif')


write.table(ALL_anno_enrichment_result, 'D45_down_phase1.txt', sep = '\t', row.names = FALSE, quote = FALSE)


interested_genes <- read.delim('time_point_gene_list.txt',header=1) %>% 
  filter(trait=='D45',class=='down',phase==2) %>% 
  select('gene_id')# 替换为你的感兴趣基因
interested_genes <- interested_genes$gene_id

# 获取感兴趣的基因集中的motif
interested_motifs <- significant_motifs[significant_motifs$sequence_name %in% interested_genes, ]
# 计算感兴趣基因集中motif的数量
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算所有motif的总数
total_motif_counts <- table(significant_motifs$motif_id,exclude = NULL)
# 计算在sig中出现，但不在insterest中出现的motif
x <- interested_motifs %>% group_by(motif_id) %>% summarise(n=n())
y <- significant_motifs %>% group_by(motif_id) %>% summarise(n=n())
z <- y$motif_id %in% x$motif_id
t <- y$motif_id[!z]
# 将z添加值intersted_motif_counts中
for (non_motif in t) {
  interested_motifs <- interested_motifs %>% 
    add_row(motif_id=non_motif)
}
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算背景基因集中motif的数量（所有motif减去感兴趣的motif）
background_motif_counts <- total_motif_counts - interested_motif_counts

# 过滤掉在感兴趣基因集或背景基因集中计数为0的motif
motif_names <- names(interested_motif_counts)[interested_motif_counts > 0 & background_motif_counts > 0]



# 创建2x2列联表矩阵
# d <- data.frame(gene.not.interest=c(M-k, N-M-n+k), gene.in.interest=c(k, n-k))
# row.names(d) <- c("In_category", "not_in_category")
contingency_tables <- lapply(motif_names, function(motif) {
  matrix(c(interested_motif_counts[[motif]], # 背景，不同motif出现的数目
           background_motif_counts[[motif]]-interested_motif_counts[[motif]], # 感兴趣基因集中motif数目
           nrow(interested_motifs) - interested_motif_counts[[motif]],
           nrow(significant_motifs)  - nrow(interested_motifs) -  background_motif_counts[[motif]]#background_motif_counts[[motif]] + interested_motif_counts[[motif]] # 背景全部motif - 背景中单个motif出现的次数
           # 基因集中的motif - 基因集中单个motif出现的次数
  ), ncol = 2, byrow = TRUE)
})




# 应用Fisher精确检验，并检查每个列联表是否至少是2x2的
fisher_results <- lapply(contingency_tables, function(table) {
  fisher.test(table)$p.value
})

# 将列表转换为向量
fisher_results <- unlist(fisher_results)
# 计算FDR
fdr_results <- p.adjust(fisher_results, method = "BH")
# 计算Enrichment score
enrichment_scores <- lapply(motif_names, function(motif) {
  if (background_motif_counts[[motif]] > 0) {
    log2((interested_motif_counts[[motif]]/nrow(interested_motifs)) / (background_motif_counts[[motif]]/nrow(significant_motifs)))
  } else {
    NA  # 如果背景中的motif数量为0，则无法计算富集度得分，返回NA
  }
})
enrichment_scores <- unlist(enrichment_scores)
# 将结果转换为数据框，并排序
enrichment_results_df <- data.frame(
  Motif = motif_names,
  Enrichment.Score = enrichment_scores,
  Pvalue = fisher_results,
  FDR = fdr_results
)
sorted_enrichment_results <- enrichment_results_df[order(enrichment_results_df$FDR), ]

# 筛选结果，注释motif信息
filter_enrichment_results <- filter(enrichment_results_df,Pvalue<0.05,FDR<0.05,enrichment_scores>0)
motif_anno <- read.delim('meme_info.txt')
anno_enrichment_result <- merge(filter_enrichment_results,motif_anno,by='Motif')
ALL_anno_enrichment_result <- merge(sorted_enrichment_results,motif_anno,by='Motif')


write.table(ALL_anno_enrichment_result, 'D45_down_phase2.txt', sep = '\t', row.names = FALSE, quote = FALSE)


interested_genes <- read.delim('time_point_gene_list.txt',header=1) %>% 
  filter(trait=='D45',class=='down',phase==3) %>% 
  select('gene_id')# 替换为你的感兴趣基因
interested_genes <- interested_genes$gene_id

# 获取感兴趣的基因集中的motif
interested_motifs <- significant_motifs[significant_motifs$sequence_name %in% interested_genes, ]
# 计算感兴趣基因集中motif的数量
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算所有motif的总数
total_motif_counts <- table(significant_motifs$motif_id,exclude = NULL)
# 计算在sig中出现，但不在insterest中出现的motif
x <- interested_motifs %>% group_by(motif_id) %>% summarise(n=n())
y <- significant_motifs %>% group_by(motif_id) %>% summarise(n=n())
z <- y$motif_id %in% x$motif_id
t <- y$motif_id[!z]
# 将z添加值intersted_motif_counts中
for (non_motif in t) {
  interested_motifs <- interested_motifs %>% 
    add_row(motif_id=non_motif)
}
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算背景基因集中motif的数量（所有motif减去感兴趣的motif）
background_motif_counts <- total_motif_counts - interested_motif_counts

# 过滤掉在感兴趣基因集或背景基因集中计数为0的motif
motif_names <- names(interested_motif_counts)[interested_motif_counts > 0 & background_motif_counts > 0]



# 创建2x2列联表矩阵
# d <- data.frame(gene.not.interest=c(M-k, N-M-n+k), gene.in.interest=c(k, n-k))
# row.names(d) <- c("In_category", "not_in_category")
contingency_tables <- lapply(motif_names, function(motif) {
  matrix(c(interested_motif_counts[[motif]], # 背景，不同motif出现的数目
           background_motif_counts[[motif]]-interested_motif_counts[[motif]], # 感兴趣基因集中motif数目
           nrow(interested_motifs) - interested_motif_counts[[motif]],
           nrow(significant_motifs)  - nrow(interested_motifs) -  background_motif_counts[[motif]]#background_motif_counts[[motif]] + interested_motif_counts[[motif]] # 背景全部motif - 背景中单个motif出现的次数
           # 基因集中的motif - 基因集中单个motif出现的次数
  ), ncol = 2, byrow = TRUE)
})




# 应用Fisher精确检验，并检查每个列联表是否至少是2x2的
fisher_results <- lapply(contingency_tables, function(table) {
  fisher.test(table)$p.value
})

# 将列表转换为向量
fisher_results <- unlist(fisher_results)
# 计算FDR
fdr_results <- p.adjust(fisher_results, method = "BH")
# 计算Enrichment score
enrichment_scores <- lapply(motif_names, function(motif) {
  if (background_motif_counts[[motif]] > 0) {
    log2((interested_motif_counts[[motif]]/nrow(interested_motifs)) / (background_motif_counts[[motif]]/nrow(significant_motifs)))
  } else {
    NA  # 如果背景中的motif数量为0，则无法计算富集度得分，返回NA
  }
})
enrichment_scores <- unlist(enrichment_scores)
# 将结果转换为数据框，并排序
enrichment_results_df <- data.frame(
  Motif = motif_names,
  Enrichment.Score = enrichment_scores,
  Pvalue = fisher_results,
  FDR = fdr_results
)
sorted_enrichment_results <- enrichment_results_df[order(enrichment_results_df$FDR), ]

# 筛选结果，注释motif信息
filter_enrichment_results <- filter(enrichment_results_df,Pvalue<0.05,FDR<0.05,enrichment_scores>0)
motif_anno <- read.delim('meme_info.txt')
anno_enrichment_result <- merge(filter_enrichment_results,motif_anno,by='Motif')
ALL_anno_enrichment_result <- merge(sorted_enrichment_results,motif_anno,by='Motif')


write.table(ALL_anno_enrichment_result, 'D45_down_phase3.txt', sep = '\t', row.names = FALSE, quote = FALSE)


interested_genes <- read.delim('time_point_gene_list.txt',header=1) %>% 
  filter(trait=='D45',class=='down',phase==4) %>% 
  select('gene_id')# 替换为你的感兴趣基因
interested_genes <- interested_genes$gene_id

# 获取感兴趣的基因集中的motif
interested_motifs <- significant_motifs[significant_motifs$sequence_name %in% interested_genes, ]
# 计算感兴趣基因集中motif的数量
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算所有motif的总数
total_motif_counts <- table(significant_motifs$motif_id,exclude = NULL)
# 计算在sig中出现，但不在insterest中出现的motif
x <- interested_motifs %>% group_by(motif_id) %>% summarise(n=n())
y <- significant_motifs %>% group_by(motif_id) %>% summarise(n=n())
z <- y$motif_id %in% x$motif_id
t <- y$motif_id[!z]
# 将z添加值intersted_motif_counts中
for (non_motif in t) {
  interested_motifs <- interested_motifs %>% 
    add_row(motif_id=non_motif)
}
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算背景基因集中motif的数量（所有motif减去感兴趣的motif）
background_motif_counts <- total_motif_counts - interested_motif_counts

# 过滤掉在感兴趣基因集或背景基因集中计数为0的motif
motif_names <- names(interested_motif_counts)[interested_motif_counts > 0 & background_motif_counts > 0]



# 创建2x2列联表矩阵
# d <- data.frame(gene.not.interest=c(M-k, N-M-n+k), gene.in.interest=c(k, n-k))
# row.names(d) <- c("In_category", "not_in_category")
contingency_tables <- lapply(motif_names, function(motif) {
  matrix(c(interested_motif_counts[[motif]], # 背景，不同motif出现的数目
           background_motif_counts[[motif]]-interested_motif_counts[[motif]], # 感兴趣基因集中motif数目
           nrow(interested_motifs) - interested_motif_counts[[motif]],
           nrow(significant_motifs)  - nrow(interested_motifs) -  background_motif_counts[[motif]]#background_motif_counts[[motif]] + interested_motif_counts[[motif]] # 背景全部motif - 背景中单个motif出现的次数
           # 基因集中的motif - 基因集中单个motif出现的次数
  ), ncol = 2, byrow = TRUE)
})




# 应用Fisher精确检验，并检查每个列联表是否至少是2x2的
fisher_results <- lapply(contingency_tables, function(table) {
  fisher.test(table)$p.value
})

# 将列表转换为向量
fisher_results <- unlist(fisher_results)
# 计算FDR
fdr_results <- p.adjust(fisher_results, method = "BH")
# 计算Enrichment score
enrichment_scores <- lapply(motif_names, function(motif) {
  if (background_motif_counts[[motif]] > 0) {
    log2((interested_motif_counts[[motif]]/nrow(interested_motifs)) / (background_motif_counts[[motif]]/nrow(significant_motifs)))
  } else {
    NA  # 如果背景中的motif数量为0，则无法计算富集度得分，返回NA
  }
})
enrichment_scores <- unlist(enrichment_scores)
# 将结果转换为数据框，并排序
enrichment_results_df <- data.frame(
  Motif = motif_names,
  Enrichment.Score = enrichment_scores,
  Pvalue = fisher_results,
  FDR = fdr_results
)
sorted_enrichment_results <- enrichment_results_df[order(enrichment_results_df$FDR), ]

# 筛选结果，注释motif信息
filter_enrichment_results <- filter(enrichment_results_df,Pvalue<0.05,FDR<0.05,enrichment_scores>0)
motif_anno <- read.delim('meme_info.txt')
anno_enrichment_result <- merge(filter_enrichment_results,motif_anno,by='Motif')
ALL_anno_enrichment_result <- merge(sorted_enrichment_results,motif_anno,by='Motif')


write.table(ALL_anno_enrichment_result, 'D45_down_phase4.txt', sep = '\t', row.names = FALSE, quote = FALSE)


interested_genes <- read.delim('time_point_gene_list.txt',header=1) %>% 
  filter(trait=='D45',class=='down',phase==5) %>% 
  select('gene_id')# 替换为你的感兴趣基因
interested_genes <- interested_genes$gene_id

# 获取感兴趣的基因集中的motif
interested_motifs <- significant_motifs[significant_motifs$sequence_name %in% interested_genes, ]
# 计算感兴趣基因集中motif的数量
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算所有motif的总数
total_motif_counts <- table(significant_motifs$motif_id,exclude = NULL)
# 计算在sig中出现，但不在insterest中出现的motif
x <- interested_motifs %>% group_by(motif_id) %>% summarise(n=n())
y <- significant_motifs %>% group_by(motif_id) %>% summarise(n=n())
z <- y$motif_id %in% x$motif_id
t <- y$motif_id[!z]
# 将z添加值intersted_motif_counts中
for (non_motif in t) {
  interested_motifs <- interested_motifs %>% 
    add_row(motif_id=non_motif)
}
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算背景基因集中motif的数量（所有motif减去感兴趣的motif）
background_motif_counts <- total_motif_counts - interested_motif_counts

# 过滤掉在感兴趣基因集或背景基因集中计数为0的motif
motif_names <- names(interested_motif_counts)[interested_motif_counts > 0 & background_motif_counts > 0]



# 创建2x2列联表矩阵
# d <- data.frame(gene.not.interest=c(M-k, N-M-n+k), gene.in.interest=c(k, n-k))
# row.names(d) <- c("In_category", "not_in_category")
contingency_tables <- lapply(motif_names, function(motif) {
  matrix(c(interested_motif_counts[[motif]], # 背景，不同motif出现的数目
           background_motif_counts[[motif]]-interested_motif_counts[[motif]], # 感兴趣基因集中motif数目
           nrow(interested_motifs) - interested_motif_counts[[motif]],
           nrow(significant_motifs)  - nrow(interested_motifs) -  background_motif_counts[[motif]]#background_motif_counts[[motif]] + interested_motif_counts[[motif]] # 背景全部motif - 背景中单个motif出现的次数
           # 基因集中的motif - 基因集中单个motif出现的次数
  ), ncol = 2, byrow = TRUE)
})




# 应用Fisher精确检验，并检查每个列联表是否至少是2x2的
fisher_results <- lapply(contingency_tables, function(table) {
  fisher.test(table)$p.value
})

# 将列表转换为向量
fisher_results <- unlist(fisher_results)
# 计算FDR
fdr_results <- p.adjust(fisher_results, method = "BH")
# 计算Enrichment score
enrichment_scores <- lapply(motif_names, function(motif) {
  if (background_motif_counts[[motif]] > 0) {
    log2((interested_motif_counts[[motif]]/nrow(interested_motifs)) / (background_motif_counts[[motif]]/nrow(significant_motifs)))
  } else {
    NA  # 如果背景中的motif数量为0，则无法计算富集度得分，返回NA
  }
})
enrichment_scores <- unlist(enrichment_scores)
# 将结果转换为数据框，并排序
enrichment_results_df <- data.frame(
  Motif = motif_names,
  Enrichment.Score = enrichment_scores,
  Pvalue = fisher_results,
  FDR = fdr_results
)
sorted_enrichment_results <- enrichment_results_df[order(enrichment_results_df$FDR), ]

# 筛选结果，注释motif信息
filter_enrichment_results <- filter(enrichment_results_df,Pvalue<0.05,FDR<0.05,enrichment_scores>0)
motif_anno <- read.delim('meme_info.txt')
anno_enrichment_result <- merge(filter_enrichment_results,motif_anno,by='Motif')
ALL_anno_enrichment_result <- merge(sorted_enrichment_results,motif_anno,by='Motif')


write.table(ALL_anno_enrichment_result, 'D45_down_phase5.txt', sep = '\t', row.names = FALSE, quote = FALSE)


interested_genes <- read.delim('time_point_gene_list.txt',header=1) %>% 
  filter(trait=='N38',class=='up',phase==1) %>% 
  select('gene_id')# 替换为你的感兴趣基因
interested_genes <- interested_genes$gene_id

# 获取感兴趣的基因集中的motif
interested_motifs <- significant_motifs[significant_motifs$sequence_name %in% interested_genes, ]
# 计算感兴趣基因集中motif的数量
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算所有motif的总数
total_motif_counts <- table(significant_motifs$motif_id,exclude = NULL)
# 计算在sig中出现，但不在insterest中出现的motif
x <- interested_motifs %>% group_by(motif_id) %>% summarise(n=n())
y <- significant_motifs %>% group_by(motif_id) %>% summarise(n=n())
z <- y$motif_id %in% x$motif_id
t <- y$motif_id[!z]
# 将z添加值intersted_motif_counts中
for (non_motif in t) {
  interested_motifs <- interested_motifs %>% 
    add_row(motif_id=non_motif)
}
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算背景基因集中motif的数量（所有motif减去感兴趣的motif）
background_motif_counts <- total_motif_counts - interested_motif_counts

# 过滤掉在感兴趣基因集或背景基因集中计数为0的motif
motif_names <- names(interested_motif_counts)[interested_motif_counts > 0 & background_motif_counts > 0]



# 创建2x2列联表矩阵
# d <- data.frame(gene.not.interest=c(M-k, N-M-n+k), gene.in.interest=c(k, n-k))
# row.names(d) <- c("In_category", "not_in_category")
contingency_tables <- lapply(motif_names, function(motif) {
  matrix(c(interested_motif_counts[[motif]], # 背景，不同motif出现的数目
           background_motif_counts[[motif]]-interested_motif_counts[[motif]], # 感兴趣基因集中motif数目
           nrow(interested_motifs) - interested_motif_counts[[motif]],
           nrow(significant_motifs)  - nrow(interested_motifs) -  background_motif_counts[[motif]]#background_motif_counts[[motif]] + interested_motif_counts[[motif]] # 背景全部motif - 背景中单个motif出现的次数
           # 基因集中的motif - 基因集中单个motif出现的次数
  ), ncol = 2, byrow = TRUE)
})




# 应用Fisher精确检验，并检查每个列联表是否至少是2x2的
fisher_results <- lapply(contingency_tables, function(table) {
  fisher.test(table)$p.value
})

# 将列表转换为向量
fisher_results <- unlist(fisher_results)
# 计算FDR
fdr_results <- p.adjust(fisher_results, method = "BH")
# 计算Enrichment score
enrichment_scores <- lapply(motif_names, function(motif) {
  if (background_motif_counts[[motif]] > 0) {
    log2((interested_motif_counts[[motif]]/nrow(interested_motifs)) / (background_motif_counts[[motif]]/nrow(significant_motifs)))
  } else {
    NA  # 如果背景中的motif数量为0，则无法计算富集度得分，返回NA
  }
})
enrichment_scores <- unlist(enrichment_scores)
# 将结果转换为数据框，并排序
enrichment_results_df <- data.frame(
  Motif = motif_names,
  Enrichment.Score = enrichment_scores,
  Pvalue = fisher_results,
  FDR = fdr_results
)
sorted_enrichment_results <- enrichment_results_df[order(enrichment_results_df$FDR), ]

# 筛选结果，注释motif信息
filter_enrichment_results <- filter(enrichment_results_df,Pvalue<0.05,FDR<0.05,enrichment_scores>0)
motif_anno <- read.delim('meme_info.txt')
anno_enrichment_result <- merge(filter_enrichment_results,motif_anno,by='Motif')
ALL_anno_enrichment_result <- merge(sorted_enrichment_results,motif_anno,by='Motif')


write.table(ALL_anno_enrichment_result, 'N38_up_phase1.txt', sep = '\t', row.names = FALSE, quote = FALSE)


interested_genes <- read.delim('time_point_gene_list.txt',header=1) %>% 
  filter(trait=='N38',class=='up',phase==2) %>% 
  select('gene_id')# 替换为你的感兴趣基因
interested_genes <- interested_genes$gene_id

# 获取感兴趣的基因集中的motif
interested_motifs <- significant_motifs[significant_motifs$sequence_name %in% interested_genes, ]
# 计算感兴趣基因集中motif的数量
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算所有motif的总数
total_motif_counts <- table(significant_motifs$motif_id,exclude = NULL)
# 计算在sig中出现，但不在insterest中出现的motif
x <- interested_motifs %>% group_by(motif_id) %>% summarise(n=n())
y <- significant_motifs %>% group_by(motif_id) %>% summarise(n=n())
z <- y$motif_id %in% x$motif_id
t <- y$motif_id[!z]
# 将z添加值intersted_motif_counts中
for (non_motif in t) {
  interested_motifs <- interested_motifs %>% 
    add_row(motif_id=non_motif)
}
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算背景基因集中motif的数量（所有motif减去感兴趣的motif）
background_motif_counts <- total_motif_counts - interested_motif_counts

# 过滤掉在感兴趣基因集或背景基因集中计数为0的motif
motif_names <- names(interested_motif_counts)[interested_motif_counts > 0 & background_motif_counts > 0]



# 创建2x2列联表矩阵
# d <- data.frame(gene.not.interest=c(M-k, N-M-n+k), gene.in.interest=c(k, n-k))
# row.names(d) <- c("In_category", "not_in_category")
contingency_tables <- lapply(motif_names, function(motif) {
  matrix(c(interested_motif_counts[[motif]], # 背景，不同motif出现的数目
           background_motif_counts[[motif]]-interested_motif_counts[[motif]], # 感兴趣基因集中motif数目
           nrow(interested_motifs) - interested_motif_counts[[motif]],
           nrow(significant_motifs)  - nrow(interested_motifs) -  background_motif_counts[[motif]]#background_motif_counts[[motif]] + interested_motif_counts[[motif]] # 背景全部motif - 背景中单个motif出现的次数
           # 基因集中的motif - 基因集中单个motif出现的次数
  ), ncol = 2, byrow = TRUE)
})




# 应用Fisher精确检验，并检查每个列联表是否至少是2x2的
fisher_results <- lapply(contingency_tables, function(table) {
  fisher.test(table)$p.value
})

# 将列表转换为向量
fisher_results <- unlist(fisher_results)
# 计算FDR
fdr_results <- p.adjust(fisher_results, method = "BH")
# 计算Enrichment score
enrichment_scores <- lapply(motif_names, function(motif) {
  if (background_motif_counts[[motif]] > 0) {
    log2((interested_motif_counts[[motif]]/nrow(interested_motifs)) / (background_motif_counts[[motif]]/nrow(significant_motifs)))
  } else {
    NA  # 如果背景中的motif数量为0，则无法计算富集度得分，返回NA
  }
})
enrichment_scores <- unlist(enrichment_scores)
# 将结果转换为数据框，并排序
enrichment_results_df <- data.frame(
  Motif = motif_names,
  Enrichment.Score = enrichment_scores,
  Pvalue = fisher_results,
  FDR = fdr_results
)
sorted_enrichment_results <- enrichment_results_df[order(enrichment_results_df$FDR), ]

# 筛选结果，注释motif信息
filter_enrichment_results <- filter(enrichment_results_df,Pvalue<0.05,FDR<0.05,enrichment_scores>0)
motif_anno <- read.delim('meme_info.txt')
anno_enrichment_result <- merge(filter_enrichment_results,motif_anno,by='Motif')
ALL_anno_enrichment_result <- merge(sorted_enrichment_results,motif_anno,by='Motif')


write.table(ALL_anno_enrichment_result, 'N38_up_phase2.txt', sep = '\t', row.names = FALSE, quote = FALSE)


interested_genes <- read.delim('time_point_gene_list.txt',header=1) %>% 
  filter(trait=='N38',class=='up',phase==3) %>% 
  select('gene_id')# 替换为你的感兴趣基因
interested_genes <- interested_genes$gene_id

# 获取感兴趣的基因集中的motif
interested_motifs <- significant_motifs[significant_motifs$sequence_name %in% interested_genes, ]
# 计算感兴趣基因集中motif的数量
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算所有motif的总数
total_motif_counts <- table(significant_motifs$motif_id,exclude = NULL)
# 计算在sig中出现，但不在insterest中出现的motif
x <- interested_motifs %>% group_by(motif_id) %>% summarise(n=n())
y <- significant_motifs %>% group_by(motif_id) %>% summarise(n=n())
z <- y$motif_id %in% x$motif_id
t <- y$motif_id[!z]
# 将z添加值intersted_motif_counts中
for (non_motif in t) {
  interested_motifs <- interested_motifs %>% 
    add_row(motif_id=non_motif)
}
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算背景基因集中motif的数量（所有motif减去感兴趣的motif）
background_motif_counts <- total_motif_counts - interested_motif_counts

# 过滤掉在感兴趣基因集或背景基因集中计数为0的motif
motif_names <- names(interested_motif_counts)[interested_motif_counts > 0 & background_motif_counts > 0]



# 创建2x2列联表矩阵
# d <- data.frame(gene.not.interest=c(M-k, N-M-n+k), gene.in.interest=c(k, n-k))
# row.names(d) <- c("In_category", "not_in_category")
contingency_tables <- lapply(motif_names, function(motif) {
  matrix(c(interested_motif_counts[[motif]], # 背景，不同motif出现的数目
           background_motif_counts[[motif]]-interested_motif_counts[[motif]], # 感兴趣基因集中motif数目
           nrow(interested_motifs) - interested_motif_counts[[motif]],
           nrow(significant_motifs)  - nrow(interested_motifs) -  background_motif_counts[[motif]]#background_motif_counts[[motif]] + interested_motif_counts[[motif]] # 背景全部motif - 背景中单个motif出现的次数
           # 基因集中的motif - 基因集中单个motif出现的次数
  ), ncol = 2, byrow = TRUE)
})




# 应用Fisher精确检验，并检查每个列联表是否至少是2x2的
fisher_results <- lapply(contingency_tables, function(table) {
  fisher.test(table)$p.value
})

# 将列表转换为向量
fisher_results <- unlist(fisher_results)
# 计算FDR
fdr_results <- p.adjust(fisher_results, method = "BH")
# 计算Enrichment score
enrichment_scores <- lapply(motif_names, function(motif) {
  if (background_motif_counts[[motif]] > 0) {
    log2((interested_motif_counts[[motif]]/nrow(interested_motifs)) / (background_motif_counts[[motif]]/nrow(significant_motifs)))
  } else {
    NA  # 如果背景中的motif数量为0，则无法计算富集度得分，返回NA
  }
})
enrichment_scores <- unlist(enrichment_scores)
# 将结果转换为数据框，并排序
enrichment_results_df <- data.frame(
  Motif = motif_names,
  Enrichment.Score = enrichment_scores,
  Pvalue = fisher_results,
  FDR = fdr_results
)
sorted_enrichment_results <- enrichment_results_df[order(enrichment_results_df$FDR), ]

# 筛选结果，注释motif信息
filter_enrichment_results <- filter(enrichment_results_df,Pvalue<0.05,FDR<0.05,enrichment_scores>0)
motif_anno <- read.delim('meme_info.txt')
anno_enrichment_result <- merge(filter_enrichment_results,motif_anno,by='Motif')
ALL_anno_enrichment_result <- merge(sorted_enrichment_results,motif_anno,by='Motif')


write.table(ALL_anno_enrichment_result, 'N38_up_phase3.txt', sep = '\t', row.names = FALSE, quote = FALSE)


interested_genes <- read.delim('time_point_gene_list.txt',header=1) %>% 
  filter(trait=='N38',class=='up',phase==4) %>% 
  select('gene_id')# 替换为你的感兴趣基因
interested_genes <- interested_genes$gene_id

# 获取感兴趣的基因集中的motif
interested_motifs <- significant_motifs[significant_motifs$sequence_name %in% interested_genes, ]
# 计算感兴趣基因集中motif的数量
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算所有motif的总数
total_motif_counts <- table(significant_motifs$motif_id,exclude = NULL)
# 计算在sig中出现，但不在insterest中出现的motif
x <- interested_motifs %>% group_by(motif_id) %>% summarise(n=n())
y <- significant_motifs %>% group_by(motif_id) %>% summarise(n=n())
z <- y$motif_id %in% x$motif_id
t <- y$motif_id[!z]
# 将z添加值intersted_motif_counts中
for (non_motif in t) {
  interested_motifs <- interested_motifs %>% 
    add_row(motif_id=non_motif)
}
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算背景基因集中motif的数量（所有motif减去感兴趣的motif）
background_motif_counts <- total_motif_counts - interested_motif_counts

# 过滤掉在感兴趣基因集或背景基因集中计数为0的motif
motif_names <- names(interested_motif_counts)[interested_motif_counts > 0 & background_motif_counts > 0]



# 创建2x2列联表矩阵
# d <- data.frame(gene.not.interest=c(M-k, N-M-n+k), gene.in.interest=c(k, n-k))
# row.names(d) <- c("In_category", "not_in_category")
contingency_tables <- lapply(motif_names, function(motif) {
  matrix(c(interested_motif_counts[[motif]], # 背景，不同motif出现的数目
           background_motif_counts[[motif]]-interested_motif_counts[[motif]], # 感兴趣基因集中motif数目
           nrow(interested_motifs) - interested_motif_counts[[motif]],
           nrow(significant_motifs)  - nrow(interested_motifs) -  background_motif_counts[[motif]]#background_motif_counts[[motif]] + interested_motif_counts[[motif]] # 背景全部motif - 背景中单个motif出现的次数
           # 基因集中的motif - 基因集中单个motif出现的次数
  ), ncol = 2, byrow = TRUE)
})




# 应用Fisher精确检验，并检查每个列联表是否至少是2x2的
fisher_results <- lapply(contingency_tables, function(table) {
  fisher.test(table)$p.value
})

# 将列表转换为向量
fisher_results <- unlist(fisher_results)
# 计算FDR
fdr_results <- p.adjust(fisher_results, method = "BH")
# 计算Enrichment score
enrichment_scores <- lapply(motif_names, function(motif) {
  if (background_motif_counts[[motif]] > 0) {
    log2((interested_motif_counts[[motif]]/nrow(interested_motifs)) / (background_motif_counts[[motif]]/nrow(significant_motifs)))
  } else {
    NA  # 如果背景中的motif数量为0，则无法计算富集度得分，返回NA
  }
})
enrichment_scores <- unlist(enrichment_scores)
# 将结果转换为数据框，并排序
enrichment_results_df <- data.frame(
  Motif = motif_names,
  Enrichment.Score = enrichment_scores,
  Pvalue = fisher_results,
  FDR = fdr_results
)
sorted_enrichment_results <- enrichment_results_df[order(enrichment_results_df$FDR), ]

# 筛选结果，注释motif信息
filter_enrichment_results <- filter(enrichment_results_df,Pvalue<0.05,FDR<0.05,enrichment_scores>0)
motif_anno <- read.delim('meme_info.txt')
anno_enrichment_result <- merge(filter_enrichment_results,motif_anno,by='Motif')
ALL_anno_enrichment_result <- merge(sorted_enrichment_results,motif_anno,by='Motif')


write.table(ALL_anno_enrichment_result, 'N38_up_phase4.txt', sep = '\t', row.names = FALSE, quote = FALSE)


interested_genes <- read.delim('time_point_gene_list.txt',header=1) %>% 
  filter(trait=='N38',class=='up',phase==5) %>% 
  select('gene_id')# 替换为你的感兴趣基因
interested_genes <- interested_genes$gene_id

# 获取感兴趣的基因集中的motif
interested_motifs <- significant_motifs[significant_motifs$sequence_name %in% interested_genes, ]
# 计算感兴趣基因集中motif的数量
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算所有motif的总数
total_motif_counts <- table(significant_motifs$motif_id,exclude = NULL)
# 计算在sig中出现，但不在insterest中出现的motif
x <- interested_motifs %>% group_by(motif_id) %>% summarise(n=n())
y <- significant_motifs %>% group_by(motif_id) %>% summarise(n=n())
z <- y$motif_id %in% x$motif_id
t <- y$motif_id[!z]
# 将z添加值intersted_motif_counts中
for (non_motif in t) {
  interested_motifs <- interested_motifs %>% 
    add_row(motif_id=non_motif)
}
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算背景基因集中motif的数量（所有motif减去感兴趣的motif）
background_motif_counts <- total_motif_counts - interested_motif_counts

# 过滤掉在感兴趣基因集或背景基因集中计数为0的motif
motif_names <- names(interested_motif_counts)[interested_motif_counts > 0 & background_motif_counts > 0]



# 创建2x2列联表矩阵
# d <- data.frame(gene.not.interest=c(M-k, N-M-n+k), gene.in.interest=c(k, n-k))
# row.names(d) <- c("In_category", "not_in_category")
contingency_tables <- lapply(motif_names, function(motif) {
  matrix(c(interested_motif_counts[[motif]], # 背景，不同motif出现的数目
           background_motif_counts[[motif]]-interested_motif_counts[[motif]], # 感兴趣基因集中motif数目
           nrow(interested_motifs) - interested_motif_counts[[motif]],
           nrow(significant_motifs)  - nrow(interested_motifs) -  background_motif_counts[[motif]]#background_motif_counts[[motif]] + interested_motif_counts[[motif]] # 背景全部motif - 背景中单个motif出现的次数
           # 基因集中的motif - 基因集中单个motif出现的次数
  ), ncol = 2, byrow = TRUE)
})




# 应用Fisher精确检验，并检查每个列联表是否至少是2x2的
fisher_results <- lapply(contingency_tables, function(table) {
  fisher.test(table)$p.value
})

# 将列表转换为向量
fisher_results <- unlist(fisher_results)
# 计算FDR
fdr_results <- p.adjust(fisher_results, method = "BH")
# 计算Enrichment score
enrichment_scores <- lapply(motif_names, function(motif) {
  if (background_motif_counts[[motif]] > 0) {
    log2((interested_motif_counts[[motif]]/nrow(interested_motifs)) / (background_motif_counts[[motif]]/nrow(significant_motifs)))
  } else {
    NA  # 如果背景中的motif数量为0，则无法计算富集度得分，返回NA
  }
})
enrichment_scores <- unlist(enrichment_scores)
# 将结果转换为数据框，并排序
enrichment_results_df <- data.frame(
  Motif = motif_names,
  Enrichment.Score = enrichment_scores,
  Pvalue = fisher_results,
  FDR = fdr_results
)
sorted_enrichment_results <- enrichment_results_df[order(enrichment_results_df$FDR), ]

# 筛选结果，注释motif信息
filter_enrichment_results <- filter(enrichment_results_df,Pvalue<0.05,FDR<0.05,enrichment_scores>0)
motif_anno <- read.delim('meme_info.txt')
anno_enrichment_result <- merge(filter_enrichment_results,motif_anno,by='Motif')
ALL_anno_enrichment_result <- merge(sorted_enrichment_results,motif_anno,by='Motif')


write.table(ALL_anno_enrichment_result, 'N38_up_phase5.txt', sep = '\t', row.names = FALSE, quote = FALSE)


interested_genes <- read.delim('time_point_gene_list.txt',header=1) %>% 
  filter(trait=='N38',class=='down',phase==1) %>% 
  select('gene_id')# 替换为你的感兴趣基因
interested_genes <- interested_genes$gene_id

# 获取感兴趣的基因集中的motif
interested_motifs <- significant_motifs[significant_motifs$sequence_name %in% interested_genes, ]
# 计算感兴趣基因集中motif的数量
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算所有motif的总数
total_motif_counts <- table(significant_motifs$motif_id,exclude = NULL)
# 计算在sig中出现，但不在insterest中出现的motif
x <- interested_motifs %>% group_by(motif_id) %>% summarise(n=n())
y <- significant_motifs %>% group_by(motif_id) %>% summarise(n=n())
z <- y$motif_id %in% x$motif_id
t <- y$motif_id[!z]
# 将z添加值intersted_motif_counts中
for (non_motif in t) {
  interested_motifs <- interested_motifs %>% 
    add_row(motif_id=non_motif)
}
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算背景基因集中motif的数量（所有motif减去感兴趣的motif）
background_motif_counts <- total_motif_counts - interested_motif_counts

# 过滤掉在感兴趣基因集或背景基因集中计数为0的motif
motif_names <- names(interested_motif_counts)[interested_motif_counts > 0 & background_motif_counts > 0]



# 创建2x2列联表矩阵
# d <- data.frame(gene.not.interest=c(M-k, N-M-n+k), gene.in.interest=c(k, n-k))
# row.names(d) <- c("In_category", "not_in_category")
contingency_tables <- lapply(motif_names, function(motif) {
  matrix(c(interested_motif_counts[[motif]], # 背景，不同motif出现的数目
           background_motif_counts[[motif]]-interested_motif_counts[[motif]], # 感兴趣基因集中motif数目
           nrow(interested_motifs) - interested_motif_counts[[motif]],
           nrow(significant_motifs)  - nrow(interested_motifs) -  background_motif_counts[[motif]]#background_motif_counts[[motif]] + interested_motif_counts[[motif]] # 背景全部motif - 背景中单个motif出现的次数
           # 基因集中的motif - 基因集中单个motif出现的次数
  ), ncol = 2, byrow = TRUE)
})




# 应用Fisher精确检验，并检查每个列联表是否至少是2x2的
fisher_results <- lapply(contingency_tables, function(table) {
  fisher.test(table)$p.value
})

# 将列表转换为向量
fisher_results <- unlist(fisher_results)
# 计算FDR
fdr_results <- p.adjust(fisher_results, method = "BH")
# 计算Enrichment score
enrichment_scores <- lapply(motif_names, function(motif) {
  if (background_motif_counts[[motif]] > 0) {
    log2((interested_motif_counts[[motif]]/nrow(interested_motifs)) / (background_motif_counts[[motif]]/nrow(significant_motifs)))
  } else {
    NA  # 如果背景中的motif数量为0，则无法计算富集度得分，返回NA
  }
})
enrichment_scores <- unlist(enrichment_scores)
# 将结果转换为数据框，并排序
enrichment_results_df <- data.frame(
  Motif = motif_names,
  Enrichment.Score = enrichment_scores,
  Pvalue = fisher_results,
  FDR = fdr_results
)
sorted_enrichment_results <- enrichment_results_df[order(enrichment_results_df$FDR), ]

# 筛选结果，注释motif信息
filter_enrichment_results <- filter(enrichment_results_df,Pvalue<0.05,FDR<0.05,enrichment_scores>0)
motif_anno <- read.delim('meme_info.txt')
anno_enrichment_result <- merge(filter_enrichment_results,motif_anno,by='Motif')
ALL_anno_enrichment_result <- merge(sorted_enrichment_results,motif_anno,by='Motif')


write.table(ALL_anno_enrichment_result, 'N38_down_phase1.txt', sep = '\t', row.names = FALSE, quote = FALSE)


interested_genes <- read.delim('time_point_gene_list.txt',header=1) %>% 
  filter(trait=='N38',class=='down',phase==2) %>% 
  select('gene_id')# 替换为你的感兴趣基因
interested_genes <- interested_genes$gene_id

# 获取感兴趣的基因集中的motif
interested_motifs <- significant_motifs[significant_motifs$sequence_name %in% interested_genes, ]
# 计算感兴趣基因集中motif的数量
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算所有motif的总数
total_motif_counts <- table(significant_motifs$motif_id,exclude = NULL)
# 计算在sig中出现，但不在insterest中出现的motif
x <- interested_motifs %>% group_by(motif_id) %>% summarise(n=n())
y <- significant_motifs %>% group_by(motif_id) %>% summarise(n=n())
z <- y$motif_id %in% x$motif_id
t <- y$motif_id[!z]
# 将z添加值intersted_motif_counts中
for (non_motif in t) {
  interested_motifs <- interested_motifs %>% 
    add_row(motif_id=non_motif)
}
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算背景基因集中motif的数量（所有motif减去感兴趣的motif）
background_motif_counts <- total_motif_counts - interested_motif_counts

# 过滤掉在感兴趣基因集或背景基因集中计数为0的motif
motif_names <- names(interested_motif_counts)[interested_motif_counts > 0 & background_motif_counts > 0]



# 创建2x2列联表矩阵
# d <- data.frame(gene.not.interest=c(M-k, N-M-n+k), gene.in.interest=c(k, n-k))
# row.names(d) <- c("In_category", "not_in_category")
contingency_tables <- lapply(motif_names, function(motif) {
  matrix(c(interested_motif_counts[[motif]], # 背景，不同motif出现的数目
           background_motif_counts[[motif]]-interested_motif_counts[[motif]], # 感兴趣基因集中motif数目
           nrow(interested_motifs) - interested_motif_counts[[motif]],
           nrow(significant_motifs)  - nrow(interested_motifs) -  background_motif_counts[[motif]]#background_motif_counts[[motif]] + interested_motif_counts[[motif]] # 背景全部motif - 背景中单个motif出现的次数
           # 基因集中的motif - 基因集中单个motif出现的次数
  ), ncol = 2, byrow = TRUE)
})




# 应用Fisher精确检验，并检查每个列联表是否至少是2x2的
fisher_results <- lapply(contingency_tables, function(table) {
  fisher.test(table)$p.value
})

# 将列表转换为向量
fisher_results <- unlist(fisher_results)
# 计算FDR
fdr_results <- p.adjust(fisher_results, method = "BH")
# 计算Enrichment score
enrichment_scores <- lapply(motif_names, function(motif) {
  if (background_motif_counts[[motif]] > 0) {
    log2((interested_motif_counts[[motif]]/nrow(interested_motifs)) / (background_motif_counts[[motif]]/nrow(significant_motifs)))
  } else {
    NA  # 如果背景中的motif数量为0，则无法计算富集度得分，返回NA
  }
})
enrichment_scores <- unlist(enrichment_scores)
# 将结果转换为数据框，并排序
enrichment_results_df <- data.frame(
  Motif = motif_names,
  Enrichment.Score = enrichment_scores,
  Pvalue = fisher_results,
  FDR = fdr_results
)
sorted_enrichment_results <- enrichment_results_df[order(enrichment_results_df$FDR), ]

# 筛选结果，注释motif信息
filter_enrichment_results <- filter(enrichment_results_df,Pvalue<0.05,FDR<0.05,enrichment_scores>0)
motif_anno <- read.delim('meme_info.txt')
anno_enrichment_result <- merge(filter_enrichment_results,motif_anno,by='Motif')
ALL_anno_enrichment_result <- merge(sorted_enrichment_results,motif_anno,by='Motif')


write.table(ALL_anno_enrichment_result, 'N38_down_phase2.txt', sep = '\t', row.names = FALSE, quote = FALSE)


interested_genes <- read.delim('time_point_gene_list.txt',header=1) %>% 
  filter(trait=='N38',class=='down',phase==3) %>% 
  select('gene_id')# 替换为你的感兴趣基因
interested_genes <- interested_genes$gene_id

# 获取感兴趣的基因集中的motif
interested_motifs <- significant_motifs[significant_motifs$sequence_name %in% interested_genes, ]
# 计算感兴趣基因集中motif的数量
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算所有motif的总数
total_motif_counts <- table(significant_motifs$motif_id,exclude = NULL)
# 计算在sig中出现，但不在insterest中出现的motif
x <- interested_motifs %>% group_by(motif_id) %>% summarise(n=n())
y <- significant_motifs %>% group_by(motif_id) %>% summarise(n=n())
z <- y$motif_id %in% x$motif_id
t <- y$motif_id[!z]
# 将z添加值intersted_motif_counts中
for (non_motif in t) {
  interested_motifs <- interested_motifs %>% 
    add_row(motif_id=non_motif)
}
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算背景基因集中motif的数量（所有motif减去感兴趣的motif）
background_motif_counts <- total_motif_counts - interested_motif_counts

# 过滤掉在感兴趣基因集或背景基因集中计数为0的motif
motif_names <- names(interested_motif_counts)[interested_motif_counts > 0 & background_motif_counts > 0]



# 创建2x2列联表矩阵
# d <- data.frame(gene.not.interest=c(M-k, N-M-n+k), gene.in.interest=c(k, n-k))
# row.names(d) <- c("In_category", "not_in_category")
contingency_tables <- lapply(motif_names, function(motif) {
  matrix(c(interested_motif_counts[[motif]], # 背景，不同motif出现的数目
           background_motif_counts[[motif]]-interested_motif_counts[[motif]], # 感兴趣基因集中motif数目
           nrow(interested_motifs) - interested_motif_counts[[motif]],
           nrow(significant_motifs)  - nrow(interested_motifs) -  background_motif_counts[[motif]]#background_motif_counts[[motif]] + interested_motif_counts[[motif]] # 背景全部motif - 背景中单个motif出现的次数
           # 基因集中的motif - 基因集中单个motif出现的次数
  ), ncol = 2, byrow = TRUE)
})




# 应用Fisher精确检验，并检查每个列联表是否至少是2x2的
fisher_results <- lapply(contingency_tables, function(table) {
  fisher.test(table)$p.value
})

# 将列表转换为向量
fisher_results <- unlist(fisher_results)
# 计算FDR
fdr_results <- p.adjust(fisher_results, method = "BH")
# 计算Enrichment score
enrichment_scores <- lapply(motif_names, function(motif) {
  if (background_motif_counts[[motif]] > 0) {
    log2((interested_motif_counts[[motif]]/nrow(interested_motifs)) / (background_motif_counts[[motif]]/nrow(significant_motifs)))
  } else {
    NA  # 如果背景中的motif数量为0，则无法计算富集度得分，返回NA
  }
})
enrichment_scores <- unlist(enrichment_scores)
# 将结果转换为数据框，并排序
enrichment_results_df <- data.frame(
  Motif = motif_names,
  Enrichment.Score = enrichment_scores,
  Pvalue = fisher_results,
  FDR = fdr_results
)
sorted_enrichment_results <- enrichment_results_df[order(enrichment_results_df$FDR), ]

# 筛选结果，注释motif信息
filter_enrichment_results <- filter(enrichment_results_df,Pvalue<0.05,FDR<0.05,enrichment_scores>0)
motif_anno <- read.delim('meme_info.txt')
anno_enrichment_result <- merge(filter_enrichment_results,motif_anno,by='Motif')
ALL_anno_enrichment_result <- merge(sorted_enrichment_results,motif_anno,by='Motif')


write.table(ALL_anno_enrichment_result, 'N38_down_phase3.txt', sep = '\t', row.names = FALSE, quote = FALSE)


interested_genes <- read.delim('time_point_gene_list.txt',header=1) %>% 
  filter(trait=='N38',class=='down',phase==4) %>% 
  select('gene_id')# 替换为你的感兴趣基因
interested_genes <- interested_genes$gene_id

# 获取感兴趣的基因集中的motif
interested_motifs <- significant_motifs[significant_motifs$sequence_name %in% interested_genes, ]
# 计算感兴趣基因集中motif的数量
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算所有motif的总数
total_motif_counts <- table(significant_motifs$motif_id,exclude = NULL)
# 计算在sig中出现，但不在insterest中出现的motif
x <- interested_motifs %>% group_by(motif_id) %>% summarise(n=n())
y <- significant_motifs %>% group_by(motif_id) %>% summarise(n=n())
z <- y$motif_id %in% x$motif_id
t <- y$motif_id[!z]
# 将z添加值intersted_motif_counts中
for (non_motif in t) {
  interested_motifs <- interested_motifs %>% 
    add_row(motif_id=non_motif)
}
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算背景基因集中motif的数量（所有motif减去感兴趣的motif）
background_motif_counts <- total_motif_counts - interested_motif_counts

# 过滤掉在感兴趣基因集或背景基因集中计数为0的motif
motif_names <- names(interested_motif_counts)[interested_motif_counts > 0 & background_motif_counts > 0]



# 创建2x2列联表矩阵
# d <- data.frame(gene.not.interest=c(M-k, N-M-n+k), gene.in.interest=c(k, n-k))
# row.names(d) <- c("In_category", "not_in_category")
contingency_tables <- lapply(motif_names, function(motif) {
  matrix(c(interested_motif_counts[[motif]], # 背景，不同motif出现的数目
           background_motif_counts[[motif]]-interested_motif_counts[[motif]], # 感兴趣基因集中motif数目
           nrow(interested_motifs) - interested_motif_counts[[motif]],
           nrow(significant_motifs)  - nrow(interested_motifs) -  background_motif_counts[[motif]]#background_motif_counts[[motif]] + interested_motif_counts[[motif]] # 背景全部motif - 背景中单个motif出现的次数
           # 基因集中的motif - 基因集中单个motif出现的次数
  ), ncol = 2, byrow = TRUE)
})




# 应用Fisher精确检验，并检查每个列联表是否至少是2x2的
fisher_results <- lapply(contingency_tables, function(table) {
  fisher.test(table)$p.value
})

# 将列表转换为向量
fisher_results <- unlist(fisher_results)
# 计算FDR
fdr_results <- p.adjust(fisher_results, method = "BH")
# 计算Enrichment score
enrichment_scores <- lapply(motif_names, function(motif) {
  if (background_motif_counts[[motif]] > 0) {
    log2((interested_motif_counts[[motif]]/nrow(interested_motifs)) / (background_motif_counts[[motif]]/nrow(significant_motifs)))
  } else {
    NA  # 如果背景中的motif数量为0，则无法计算富集度得分，返回NA
  }
})
enrichment_scores <- unlist(enrichment_scores)
# 将结果转换为数据框，并排序
enrichment_results_df <- data.frame(
  Motif = motif_names,
  Enrichment.Score = enrichment_scores,
  Pvalue = fisher_results,
  FDR = fdr_results
)
sorted_enrichment_results <- enrichment_results_df[order(enrichment_results_df$FDR), ]

# 筛选结果，注释motif信息
filter_enrichment_results <- filter(enrichment_results_df,Pvalue<0.05,FDR<0.05,enrichment_scores>0)
motif_anno <- read.delim('meme_info.txt')
anno_enrichment_result <- merge(filter_enrichment_results,motif_anno,by='Motif')
ALL_anno_enrichment_result <- merge(sorted_enrichment_results,motif_anno,by='Motif')


write.table(ALL_anno_enrichment_result, 'N38_down_phase4.txt', sep = '\t', row.names = FALSE, quote = FALSE)


interested_genes <- read.delim('time_point_gene_list.txt',header=1) %>% 
  filter(trait=='N38',class=='down',phase==5) %>% 
  select('gene_id')# 替换为你的感兴趣基因
interested_genes <- interested_genes$gene_id

# 获取感兴趣的基因集中的motif
interested_motifs <- significant_motifs[significant_motifs$sequence_name %in% interested_genes, ]
# 计算感兴趣基因集中motif的数量
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算所有motif的总数
total_motif_counts <- table(significant_motifs$motif_id,exclude = NULL)
# 计算在sig中出现，但不在insterest中出现的motif
x <- interested_motifs %>% group_by(motif_id) %>% summarise(n=n())
y <- significant_motifs %>% group_by(motif_id) %>% summarise(n=n())
z <- y$motif_id %in% x$motif_id
t <- y$motif_id[!z]
# 将z添加值intersted_motif_counts中
for (non_motif in t) {
  interested_motifs <- interested_motifs %>% 
    add_row(motif_id=non_motif)
}
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算背景基因集中motif的数量（所有motif减去感兴趣的motif）
background_motif_counts <- total_motif_counts - interested_motif_counts

# 过滤掉在感兴趣基因集或背景基因集中计数为0的motif
motif_names <- names(interested_motif_counts)[interested_motif_counts > 0 & background_motif_counts > 0]



# 创建2x2列联表矩阵
# d <- data.frame(gene.not.interest=c(M-k, N-M-n+k), gene.in.interest=c(k, n-k))
# row.names(d) <- c("In_category", "not_in_category")
contingency_tables <- lapply(motif_names, function(motif) {
  matrix(c(interested_motif_counts[[motif]], # 背景，不同motif出现的数目
           background_motif_counts[[motif]]-interested_motif_counts[[motif]], # 感兴趣基因集中motif数目
           nrow(interested_motifs) - interested_motif_counts[[motif]],
           nrow(significant_motifs)  - nrow(interested_motifs) -  background_motif_counts[[motif]]#background_motif_counts[[motif]] + interested_motif_counts[[motif]] # 背景全部motif - 背景中单个motif出现的次数
           # 基因集中的motif - 基因集中单个motif出现的次数
  ), ncol = 2, byrow = TRUE)
})




# 应用Fisher精确检验，并检查每个列联表是否至少是2x2的
fisher_results <- lapply(contingency_tables, function(table) {
  fisher.test(table)$p.value
})

# 将列表转换为向量
fisher_results <- unlist(fisher_results)
# 计算FDR
fdr_results <- p.adjust(fisher_results, method = "BH")
# 计算Enrichment score
enrichment_scores <- lapply(motif_names, function(motif) {
  if (background_motif_counts[[motif]] > 0) {
    log2((interested_motif_counts[[motif]]/nrow(interested_motifs)) / (background_motif_counts[[motif]]/nrow(significant_motifs)))
  } else {
    NA  # 如果背景中的motif数量为0，则无法计算富集度得分，返回NA
  }
})
enrichment_scores <- unlist(enrichment_scores)
# 将结果转换为数据框，并排序
enrichment_results_df <- data.frame(
  Motif = motif_names,
  Enrichment.Score = enrichment_scores,
  Pvalue = fisher_results,
  FDR = fdr_results
)
sorted_enrichment_results <- enrichment_results_df[order(enrichment_results_df$FDR), ]

# 筛选结果，注释motif信息
filter_enrichment_results <- filter(enrichment_results_df,Pvalue<0.05,FDR<0.05,enrichment_scores>0)
motif_anno <- read.delim('meme_info.txt')
anno_enrichment_result <- merge(filter_enrichment_results,motif_anno,by='Motif')
ALL_anno_enrichment_result <- merge(sorted_enrichment_results,motif_anno,by='Motif')


write.table(ALL_anno_enrichment_result, 'N38_down_phase5.txt', sep = '\t', row.names = FALSE, quote = FALSE)


interested_genes <- read.delim('time_point_gene_list.txt',header=1) %>% 
  filter(trait=='N45',class=='up',phase==1) %>% 
  select('gene_id')# 替换为你的感兴趣基因
interested_genes <- interested_genes$gene_id

# 获取感兴趣的基因集中的motif
interested_motifs <- significant_motifs[significant_motifs$sequence_name %in% interested_genes, ]
# 计算感兴趣基因集中motif的数量
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算所有motif的总数
total_motif_counts <- table(significant_motifs$motif_id,exclude = NULL)
# 计算在sig中出现，但不在insterest中出现的motif
x <- interested_motifs %>% group_by(motif_id) %>% summarise(n=n())
y <- significant_motifs %>% group_by(motif_id) %>% summarise(n=n())
z <- y$motif_id %in% x$motif_id
t <- y$motif_id[!z]
# 将z添加值intersted_motif_counts中
for (non_motif in t) {
  interested_motifs <- interested_motifs %>% 
    add_row(motif_id=non_motif)
}
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算背景基因集中motif的数量（所有motif减去感兴趣的motif）
background_motif_counts <- total_motif_counts - interested_motif_counts

# 过滤掉在感兴趣基因集或背景基因集中计数为0的motif
motif_names <- names(interested_motif_counts)[interested_motif_counts > 0 & background_motif_counts > 0]



# 创建2x2列联表矩阵
# d <- data.frame(gene.not.interest=c(M-k, N-M-n+k), gene.in.interest=c(k, n-k))
# row.names(d) <- c("In_category", "not_in_category")
contingency_tables <- lapply(motif_names, function(motif) {
  matrix(c(interested_motif_counts[[motif]], # 背景，不同motif出现的数目
           background_motif_counts[[motif]]-interested_motif_counts[[motif]], # 感兴趣基因集中motif数目
           nrow(interested_motifs) - interested_motif_counts[[motif]],
           nrow(significant_motifs)  - nrow(interested_motifs) -  background_motif_counts[[motif]]#background_motif_counts[[motif]] + interested_motif_counts[[motif]] # 背景全部motif - 背景中单个motif出现的次数
           # 基因集中的motif - 基因集中单个motif出现的次数
  ), ncol = 2, byrow = TRUE)
})




# 应用Fisher精确检验，并检查每个列联表是否至少是2x2的
fisher_results <- lapply(contingency_tables, function(table) {
  fisher.test(table)$p.value
})

# 将列表转换为向量
fisher_results <- unlist(fisher_results)
# 计算FDR
fdr_results <- p.adjust(fisher_results, method = "BH")
# 计算Enrichment score
enrichment_scores <- lapply(motif_names, function(motif) {
  if (background_motif_counts[[motif]] > 0) {
    log2((interested_motif_counts[[motif]]/nrow(interested_motifs)) / (background_motif_counts[[motif]]/nrow(significant_motifs)))
  } else {
    NA  # 如果背景中的motif数量为0，则无法计算富集度得分，返回NA
  }
})
enrichment_scores <- unlist(enrichment_scores)
# 将结果转换为数据框，并排序
enrichment_results_df <- data.frame(
  Motif = motif_names,
  Enrichment.Score = enrichment_scores,
  Pvalue = fisher_results,
  FDR = fdr_results
)
sorted_enrichment_results <- enrichment_results_df[order(enrichment_results_df$FDR), ]

# 筛选结果，注释motif信息
filter_enrichment_results <- filter(enrichment_results_df,Pvalue<0.05,FDR<0.05,enrichment_scores>0)
motif_anno <- read.delim('meme_info.txt')
anno_enrichment_result <- merge(filter_enrichment_results,motif_anno,by='Motif')
ALL_anno_enrichment_result <- merge(sorted_enrichment_results,motif_anno,by='Motif')


write.table(ALL_anno_enrichment_result, 'N45_up_phase1.txt', sep = '\t', row.names = FALSE, quote = FALSE)


interested_genes <- read.delim('time_point_gene_list.txt',header=1) %>% 
  filter(trait=='N45',class=='up',phase==2) %>% 
  select('gene_id')# 替换为你的感兴趣基因
interested_genes <- interested_genes$gene_id

# 获取感兴趣的基因集中的motif
interested_motifs <- significant_motifs[significant_motifs$sequence_name %in% interested_genes, ]
# 计算感兴趣基因集中motif的数量
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算所有motif的总数
total_motif_counts <- table(significant_motifs$motif_id,exclude = NULL)
# 计算在sig中出现，但不在insterest中出现的motif
x <- interested_motifs %>% group_by(motif_id) %>% summarise(n=n())
y <- significant_motifs %>% group_by(motif_id) %>% summarise(n=n())
z <- y$motif_id %in% x$motif_id
t <- y$motif_id[!z]
# 将z添加值intersted_motif_counts中
for (non_motif in t) {
  interested_motifs <- interested_motifs %>% 
    add_row(motif_id=non_motif)
}
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算背景基因集中motif的数量（所有motif减去感兴趣的motif）
background_motif_counts <- total_motif_counts - interested_motif_counts

# 过滤掉在感兴趣基因集或背景基因集中计数为0的motif
motif_names <- names(interested_motif_counts)[interested_motif_counts > 0 & background_motif_counts > 0]



# 创建2x2列联表矩阵
# d <- data.frame(gene.not.interest=c(M-k, N-M-n+k), gene.in.interest=c(k, n-k))
# row.names(d) <- c("In_category", "not_in_category")
contingency_tables <- lapply(motif_names, function(motif) {
  matrix(c(interested_motif_counts[[motif]], # 背景，不同motif出现的数目
           background_motif_counts[[motif]]-interested_motif_counts[[motif]], # 感兴趣基因集中motif数目
           nrow(interested_motifs) - interested_motif_counts[[motif]],
           nrow(significant_motifs)  - nrow(interested_motifs) -  background_motif_counts[[motif]]#background_motif_counts[[motif]] + interested_motif_counts[[motif]] # 背景全部motif - 背景中单个motif出现的次数
           # 基因集中的motif - 基因集中单个motif出现的次数
  ), ncol = 2, byrow = TRUE)
})




# 应用Fisher精确检验，并检查每个列联表是否至少是2x2的
fisher_results <- lapply(contingency_tables, function(table) {
  fisher.test(table)$p.value
})

# 将列表转换为向量
fisher_results <- unlist(fisher_results)
# 计算FDR
fdr_results <- p.adjust(fisher_results, method = "BH")
# 计算Enrichment score
enrichment_scores <- lapply(motif_names, function(motif) {
  if (background_motif_counts[[motif]] > 0) {
    log2((interested_motif_counts[[motif]]/nrow(interested_motifs)) / (background_motif_counts[[motif]]/nrow(significant_motifs)))
  } else {
    NA  # 如果背景中的motif数量为0，则无法计算富集度得分，返回NA
  }
})
enrichment_scores <- unlist(enrichment_scores)
# 将结果转换为数据框，并排序
enrichment_results_df <- data.frame(
  Motif = motif_names,
  Enrichment.Score = enrichment_scores,
  Pvalue = fisher_results,
  FDR = fdr_results
)
sorted_enrichment_results <- enrichment_results_df[order(enrichment_results_df$FDR), ]

# 筛选结果，注释motif信息
filter_enrichment_results <- filter(enrichment_results_df,Pvalue<0.05,FDR<0.05,enrichment_scores>0)
motif_anno <- read.delim('meme_info.txt')
anno_enrichment_result <- merge(filter_enrichment_results,motif_anno,by='Motif')
ALL_anno_enrichment_result <- merge(sorted_enrichment_results,motif_anno,by='Motif')


write.table(ALL_anno_enrichment_result, 'N45_up_phase2.txt', sep = '\t', row.names = FALSE, quote = FALSE)


interested_genes <- read.delim('time_point_gene_list.txt',header=1) %>% 
  filter(trait=='N45',class=='up',phase==3) %>% 
  select('gene_id')# 替换为你的感兴趣基因
interested_genes <- interested_genes$gene_id

# 获取感兴趣的基因集中的motif
interested_motifs <- significant_motifs[significant_motifs$sequence_name %in% interested_genes, ]
# 计算感兴趣基因集中motif的数量
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算所有motif的总数
total_motif_counts <- table(significant_motifs$motif_id,exclude = NULL)
# 计算在sig中出现，但不在insterest中出现的motif
x <- interested_motifs %>% group_by(motif_id) %>% summarise(n=n())
y <- significant_motifs %>% group_by(motif_id) %>% summarise(n=n())
z <- y$motif_id %in% x$motif_id
t <- y$motif_id[!z]
# 将z添加值intersted_motif_counts中
for (non_motif in t) {
  interested_motifs <- interested_motifs %>% 
    add_row(motif_id=non_motif)
}
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算背景基因集中motif的数量（所有motif减去感兴趣的motif）
background_motif_counts <- total_motif_counts - interested_motif_counts

# 过滤掉在感兴趣基因集或背景基因集中计数为0的motif
motif_names <- names(interested_motif_counts)[interested_motif_counts > 0 & background_motif_counts > 0]



# 创建2x2列联表矩阵
# d <- data.frame(gene.not.interest=c(M-k, N-M-n+k), gene.in.interest=c(k, n-k))
# row.names(d) <- c("In_category", "not_in_category")
contingency_tables <- lapply(motif_names, function(motif) {
  matrix(c(interested_motif_counts[[motif]], # 背景，不同motif出现的数目
           background_motif_counts[[motif]]-interested_motif_counts[[motif]], # 感兴趣基因集中motif数目
           nrow(interested_motifs) - interested_motif_counts[[motif]],
           nrow(significant_motifs)  - nrow(interested_motifs) -  background_motif_counts[[motif]]#background_motif_counts[[motif]] + interested_motif_counts[[motif]] # 背景全部motif - 背景中单个motif出现的次数
           # 基因集中的motif - 基因集中单个motif出现的次数
  ), ncol = 2, byrow = TRUE)
})




# 应用Fisher精确检验，并检查每个列联表是否至少是2x2的
fisher_results <- lapply(contingency_tables, function(table) {
  fisher.test(table)$p.value
})

# 将列表转换为向量
fisher_results <- unlist(fisher_results)
# 计算FDR
fdr_results <- p.adjust(fisher_results, method = "BH")
# 计算Enrichment score
enrichment_scores <- lapply(motif_names, function(motif) {
  if (background_motif_counts[[motif]] > 0) {
    log2((interested_motif_counts[[motif]]/nrow(interested_motifs)) / (background_motif_counts[[motif]]/nrow(significant_motifs)))
  } else {
    NA  # 如果背景中的motif数量为0，则无法计算富集度得分，返回NA
  }
})
enrichment_scores <- unlist(enrichment_scores)
# 将结果转换为数据框，并排序
enrichment_results_df <- data.frame(
  Motif = motif_names,
  Enrichment.Score = enrichment_scores,
  Pvalue = fisher_results,
  FDR = fdr_results
)
sorted_enrichment_results <- enrichment_results_df[order(enrichment_results_df$FDR), ]

# 筛选结果，注释motif信息
filter_enrichment_results <- filter(enrichment_results_df,Pvalue<0.05,FDR<0.05,enrichment_scores>0)
motif_anno <- read.delim('meme_info.txt')
anno_enrichment_result <- merge(filter_enrichment_results,motif_anno,by='Motif')
ALL_anno_enrichment_result <- merge(sorted_enrichment_results,motif_anno,by='Motif')


write.table(ALL_anno_enrichment_result, 'N45_up_phase3.txt', sep = '\t', row.names = FALSE, quote = FALSE)


interested_genes <- read.delim('time_point_gene_list.txt',header=1) %>% 
  filter(trait=='N45',class=='up',phase==4) %>% 
  select('gene_id')# 替换为你的感兴趣基因
interested_genes <- interested_genes$gene_id

# 获取感兴趣的基因集中的motif
interested_motifs <- significant_motifs[significant_motifs$sequence_name %in% interested_genes, ]
# 计算感兴趣基因集中motif的数量
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算所有motif的总数
total_motif_counts <- table(significant_motifs$motif_id,exclude = NULL)
# 计算在sig中出现，但不在insterest中出现的motif
x <- interested_motifs %>% group_by(motif_id) %>% summarise(n=n())
y <- significant_motifs %>% group_by(motif_id) %>% summarise(n=n())
z <- y$motif_id %in% x$motif_id
t <- y$motif_id[!z]
# 将z添加值intersted_motif_counts中
for (non_motif in t) {
  interested_motifs <- interested_motifs %>% 
    add_row(motif_id=non_motif)
}
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算背景基因集中motif的数量（所有motif减去感兴趣的motif）
background_motif_counts <- total_motif_counts - interested_motif_counts

# 过滤掉在感兴趣基因集或背景基因集中计数为0的motif
motif_names <- names(interested_motif_counts)[interested_motif_counts > 0 & background_motif_counts > 0]



# 创建2x2列联表矩阵
# d <- data.frame(gene.not.interest=c(M-k, N-M-n+k), gene.in.interest=c(k, n-k))
# row.names(d) <- c("In_category", "not_in_category")
contingency_tables <- lapply(motif_names, function(motif) {
  matrix(c(interested_motif_counts[[motif]], # 背景，不同motif出现的数目
           background_motif_counts[[motif]]-interested_motif_counts[[motif]], # 感兴趣基因集中motif数目
           nrow(interested_motifs) - interested_motif_counts[[motif]],
           nrow(significant_motifs)  - nrow(interested_motifs) -  background_motif_counts[[motif]]#background_motif_counts[[motif]] + interested_motif_counts[[motif]] # 背景全部motif - 背景中单个motif出现的次数
           # 基因集中的motif - 基因集中单个motif出现的次数
  ), ncol = 2, byrow = TRUE)
})




# 应用Fisher精确检验，并检查每个列联表是否至少是2x2的
fisher_results <- lapply(contingency_tables, function(table) {
  fisher.test(table)$p.value
})

# 将列表转换为向量
fisher_results <- unlist(fisher_results)
# 计算FDR
fdr_results <- p.adjust(fisher_results, method = "BH")
# 计算Enrichment score
enrichment_scores <- lapply(motif_names, function(motif) {
  if (background_motif_counts[[motif]] > 0) {
    log2((interested_motif_counts[[motif]]/nrow(interested_motifs)) / (background_motif_counts[[motif]]/nrow(significant_motifs)))
  } else {
    NA  # 如果背景中的motif数量为0，则无法计算富集度得分，返回NA
  }
})
enrichment_scores <- unlist(enrichment_scores)
# 将结果转换为数据框，并排序
enrichment_results_df <- data.frame(
  Motif = motif_names,
  Enrichment.Score = enrichment_scores,
  Pvalue = fisher_results,
  FDR = fdr_results
)
sorted_enrichment_results <- enrichment_results_df[order(enrichment_results_df$FDR), ]

# 筛选结果，注释motif信息
filter_enrichment_results <- filter(enrichment_results_df,Pvalue<0.05,FDR<0.05,enrichment_scores>0)
motif_anno <- read.delim('meme_info.txt')
anno_enrichment_result <- merge(filter_enrichment_results,motif_anno,by='Motif')
ALL_anno_enrichment_result <- merge(sorted_enrichment_results,motif_anno,by='Motif')


write.table(ALL_anno_enrichment_result, 'N45_up_phase4.txt', sep = '\t', row.names = FALSE, quote = FALSE)


interested_genes <- read.delim('time_point_gene_list.txt',header=1) %>% 
  filter(trait=='N45',class=='up',phase==5) %>% 
  select('gene_id')# 替换为你的感兴趣基因
interested_genes <- interested_genes$gene_id

# 获取感兴趣的基因集中的motif
interested_motifs <- significant_motifs[significant_motifs$sequence_name %in% interested_genes, ]
# 计算感兴趣基因集中motif的数量
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算所有motif的总数
total_motif_counts <- table(significant_motifs$motif_id,exclude = NULL)
# 计算在sig中出现，但不在insterest中出现的motif
x <- interested_motifs %>% group_by(motif_id) %>% summarise(n=n())
y <- significant_motifs %>% group_by(motif_id) %>% summarise(n=n())
z <- y$motif_id %in% x$motif_id
t <- y$motif_id[!z]
# 将z添加值intersted_motif_counts中
for (non_motif in t) {
  interested_motifs <- interested_motifs %>% 
    add_row(motif_id=non_motif)
}
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算背景基因集中motif的数量（所有motif减去感兴趣的motif）
background_motif_counts <- total_motif_counts - interested_motif_counts

# 过滤掉在感兴趣基因集或背景基因集中计数为0的motif
motif_names <- names(interested_motif_counts)[interested_motif_counts > 0 & background_motif_counts > 0]



# 创建2x2列联表矩阵
# d <- data.frame(gene.not.interest=c(M-k, N-M-n+k), gene.in.interest=c(k, n-k))
# row.names(d) <- c("In_category", "not_in_category")
contingency_tables <- lapply(motif_names, function(motif) {
  matrix(c(interested_motif_counts[[motif]], # 背景，不同motif出现的数目
           background_motif_counts[[motif]]-interested_motif_counts[[motif]], # 感兴趣基因集中motif数目
           nrow(interested_motifs) - interested_motif_counts[[motif]],
           nrow(significant_motifs)  - nrow(interested_motifs) -  background_motif_counts[[motif]]#background_motif_counts[[motif]] + interested_motif_counts[[motif]] # 背景全部motif - 背景中单个motif出现的次数
           # 基因集中的motif - 基因集中单个motif出现的次数
  ), ncol = 2, byrow = TRUE)
})




# 应用Fisher精确检验，并检查每个列联表是否至少是2x2的
fisher_results <- lapply(contingency_tables, function(table) {
  fisher.test(table)$p.value
})

# 将列表转换为向量
fisher_results <- unlist(fisher_results)
# 计算FDR
fdr_results <- p.adjust(fisher_results, method = "BH")
# 计算Enrichment score
enrichment_scores <- lapply(motif_names, function(motif) {
  if (background_motif_counts[[motif]] > 0) {
    log2((interested_motif_counts[[motif]]/nrow(interested_motifs)) / (background_motif_counts[[motif]]/nrow(significant_motifs)))
  } else {
    NA  # 如果背景中的motif数量为0，则无法计算富集度得分，返回NA
  }
})
enrichment_scores <- unlist(enrichment_scores)
# 将结果转换为数据框，并排序
enrichment_results_df <- data.frame(
  Motif = motif_names,
  Enrichment.Score = enrichment_scores,
  Pvalue = fisher_results,
  FDR = fdr_results
)
sorted_enrichment_results <- enrichment_results_df[order(enrichment_results_df$FDR), ]

# 筛选结果，注释motif信息
filter_enrichment_results <- filter(enrichment_results_df,Pvalue<0.05,FDR<0.05,enrichment_scores>0)
motif_anno <- read.delim('meme_info.txt')
anno_enrichment_result <- merge(filter_enrichment_results,motif_anno,by='Motif')
ALL_anno_enrichment_result <- merge(sorted_enrichment_results,motif_anno,by='Motif')


write.table(ALL_anno_enrichment_result, 'N45_up_phase5.txt', sep = '\t', row.names = FALSE, quote = FALSE)


interested_genes <- read.delim('time_point_gene_list.txt',header=1) %>% 
  filter(trait=='N45',class=='down',phase==1) %>% 
  select('gene_id')# 替换为你的感兴趣基因
interested_genes <- interested_genes$gene_id

# 获取感兴趣的基因集中的motif
interested_motifs <- significant_motifs[significant_motifs$sequence_name %in% interested_genes, ]
# 计算感兴趣基因集中motif的数量
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算所有motif的总数
total_motif_counts <- table(significant_motifs$motif_id,exclude = NULL)
# 计算在sig中出现，但不在insterest中出现的motif
x <- interested_motifs %>% group_by(motif_id) %>% summarise(n=n())
y <- significant_motifs %>% group_by(motif_id) %>% summarise(n=n())
z <- y$motif_id %in% x$motif_id
t <- y$motif_id[!z]
# 将z添加值intersted_motif_counts中
for (non_motif in t) {
  interested_motifs <- interested_motifs %>% 
    add_row(motif_id=non_motif)
}
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算背景基因集中motif的数量（所有motif减去感兴趣的motif）
background_motif_counts <- total_motif_counts - interested_motif_counts

# 过滤掉在感兴趣基因集或背景基因集中计数为0的motif
motif_names <- names(interested_motif_counts)[interested_motif_counts > 0 & background_motif_counts > 0]



# 创建2x2列联表矩阵
# d <- data.frame(gene.not.interest=c(M-k, N-M-n+k), gene.in.interest=c(k, n-k))
# row.names(d) <- c("In_category", "not_in_category")
contingency_tables <- lapply(motif_names, function(motif) {
  matrix(c(interested_motif_counts[[motif]], # 背景，不同motif出现的数目
           background_motif_counts[[motif]]-interested_motif_counts[[motif]], # 感兴趣基因集中motif数目
           nrow(interested_motifs) - interested_motif_counts[[motif]],
           nrow(significant_motifs)  - nrow(interested_motifs) -  background_motif_counts[[motif]]#background_motif_counts[[motif]] + interested_motif_counts[[motif]] # 背景全部motif - 背景中单个motif出现的次数
           # 基因集中的motif - 基因集中单个motif出现的次数
  ), ncol = 2, byrow = TRUE)
})




# 应用Fisher精确检验，并检查每个列联表是否至少是2x2的
fisher_results <- lapply(contingency_tables, function(table) {
  fisher.test(table)$p.value
})

# 将列表转换为向量
fisher_results <- unlist(fisher_results)
# 计算FDR
fdr_results <- p.adjust(fisher_results, method = "BH")
# 计算Enrichment score
enrichment_scores <- lapply(motif_names, function(motif) {
  if (background_motif_counts[[motif]] > 0) {
    log2((interested_motif_counts[[motif]]/nrow(interested_motifs)) / (background_motif_counts[[motif]]/nrow(significant_motifs)))
  } else {
    NA  # 如果背景中的motif数量为0，则无法计算富集度得分，返回NA
  }
})
enrichment_scores <- unlist(enrichment_scores)
# 将结果转换为数据框，并排序
enrichment_results_df <- data.frame(
  Motif = motif_names,
  Enrichment.Score = enrichment_scores,
  Pvalue = fisher_results,
  FDR = fdr_results
)
sorted_enrichment_results <- enrichment_results_df[order(enrichment_results_df$FDR), ]

# 筛选结果，注释motif信息
filter_enrichment_results <- filter(enrichment_results_df,Pvalue<0.05,FDR<0.05,enrichment_scores>0)
motif_anno <- read.delim('meme_info.txt')
anno_enrichment_result <- merge(filter_enrichment_results,motif_anno,by='Motif')
ALL_anno_enrichment_result <- merge(sorted_enrichment_results,motif_anno,by='Motif')


write.table(ALL_anno_enrichment_result, 'N45_down_phase1.txt', sep = '\t', row.names = FALSE, quote = FALSE)


interested_genes <- read.delim('time_point_gene_list.txt',header=1) %>% 
  filter(trait=='N45',class=='down',phase==2) %>% 
  select('gene_id')# 替换为你的感兴趣基因
interested_genes <- interested_genes$gene_id

# 获取感兴趣的基因集中的motif
interested_motifs <- significant_motifs[significant_motifs$sequence_name %in% interested_genes, ]
# 计算感兴趣基因集中motif的数量
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算所有motif的总数
total_motif_counts <- table(significant_motifs$motif_id,exclude = NULL)
# 计算在sig中出现，但不在insterest中出现的motif
x <- interested_motifs %>% group_by(motif_id) %>% summarise(n=n())
y <- significant_motifs %>% group_by(motif_id) %>% summarise(n=n())
z <- y$motif_id %in% x$motif_id
t <- y$motif_id[!z]
# 将z添加值intersted_motif_counts中
for (non_motif in t) {
  interested_motifs <- interested_motifs %>% 
    add_row(motif_id=non_motif)
}
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算背景基因集中motif的数量（所有motif减去感兴趣的motif）
background_motif_counts <- total_motif_counts - interested_motif_counts

# 过滤掉在感兴趣基因集或背景基因集中计数为0的motif
motif_names <- names(interested_motif_counts)[interested_motif_counts > 0 & background_motif_counts > 0]



# 创建2x2列联表矩阵
# d <- data.frame(gene.not.interest=c(M-k, N-M-n+k), gene.in.interest=c(k, n-k))
# row.names(d) <- c("In_category", "not_in_category")
contingency_tables <- lapply(motif_names, function(motif) {
  matrix(c(interested_motif_counts[[motif]], # 背景，不同motif出现的数目
           background_motif_counts[[motif]]-interested_motif_counts[[motif]], # 感兴趣基因集中motif数目
           nrow(interested_motifs) - interested_motif_counts[[motif]],
           nrow(significant_motifs)  - nrow(interested_motifs) -  background_motif_counts[[motif]]#background_motif_counts[[motif]] + interested_motif_counts[[motif]] # 背景全部motif - 背景中单个motif出现的次数
           # 基因集中的motif - 基因集中单个motif出现的次数
  ), ncol = 2, byrow = TRUE)
})




# 应用Fisher精确检验，并检查每个列联表是否至少是2x2的
fisher_results <- lapply(contingency_tables, function(table) {
  fisher.test(table)$p.value
})

# 将列表转换为向量
fisher_results <- unlist(fisher_results)
# 计算FDR
fdr_results <- p.adjust(fisher_results, method = "BH")
# 计算Enrichment score
enrichment_scores <- lapply(motif_names, function(motif) {
  if (background_motif_counts[[motif]] > 0) {
    log2((interested_motif_counts[[motif]]/nrow(interested_motifs)) / (background_motif_counts[[motif]]/nrow(significant_motifs)))
  } else {
    NA  # 如果背景中的motif数量为0，则无法计算富集度得分，返回NA
  }
})
enrichment_scores <- unlist(enrichment_scores)
# 将结果转换为数据框，并排序
enrichment_results_df <- data.frame(
  Motif = motif_names,
  Enrichment.Score = enrichment_scores,
  Pvalue = fisher_results,
  FDR = fdr_results
)
sorted_enrichment_results <- enrichment_results_df[order(enrichment_results_df$FDR), ]

# 筛选结果，注释motif信息
filter_enrichment_results <- filter(enrichment_results_df,Pvalue<0.05,FDR<0.05,enrichment_scores>0)
motif_anno <- read.delim('meme_info.txt')
anno_enrichment_result <- merge(filter_enrichment_results,motif_anno,by='Motif')
ALL_anno_enrichment_result <- merge(sorted_enrichment_results,motif_anno,by='Motif')


write.table(ALL_anno_enrichment_result, 'N45_down_phase2.txt', sep = '\t', row.names = FALSE, quote = FALSE)


interested_genes <- read.delim('time_point_gene_list.txt',header=1) %>% 
  filter(trait=='N45',class=='down',phase==3) %>% 
  select('gene_id')# 替换为你的感兴趣基因
interested_genes <- interested_genes$gene_id

# 获取感兴趣的基因集中的motif
interested_motifs <- significant_motifs[significant_motifs$sequence_name %in% interested_genes, ]
# 计算感兴趣基因集中motif的数量
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算所有motif的总数
total_motif_counts <- table(significant_motifs$motif_id,exclude = NULL)
# 计算在sig中出现，但不在insterest中出现的motif
x <- interested_motifs %>% group_by(motif_id) %>% summarise(n=n())
y <- significant_motifs %>% group_by(motif_id) %>% summarise(n=n())
z <- y$motif_id %in% x$motif_id
t <- y$motif_id[!z]
# 将z添加值intersted_motif_counts中
for (non_motif in t) {
  interested_motifs <- interested_motifs %>% 
    add_row(motif_id=non_motif)
}
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算背景基因集中motif的数量（所有motif减去感兴趣的motif）
background_motif_counts <- total_motif_counts - interested_motif_counts

# 过滤掉在感兴趣基因集或背景基因集中计数为0的motif
motif_names <- names(interested_motif_counts)[interested_motif_counts > 0 & background_motif_counts > 0]



# 创建2x2列联表矩阵
# d <- data.frame(gene.not.interest=c(M-k, N-M-n+k), gene.in.interest=c(k, n-k))
# row.names(d) <- c("In_category", "not_in_category")
contingency_tables <- lapply(motif_names, function(motif) {
  matrix(c(interested_motif_counts[[motif]], # 背景，不同motif出现的数目
           background_motif_counts[[motif]]-interested_motif_counts[[motif]], # 感兴趣基因集中motif数目
           nrow(interested_motifs) - interested_motif_counts[[motif]],
           nrow(significant_motifs)  - nrow(interested_motifs) -  background_motif_counts[[motif]]#background_motif_counts[[motif]] + interested_motif_counts[[motif]] # 背景全部motif - 背景中单个motif出现的次数
           # 基因集中的motif - 基因集中单个motif出现的次数
  ), ncol = 2, byrow = TRUE)
})




# 应用Fisher精确检验，并检查每个列联表是否至少是2x2的
fisher_results <- lapply(contingency_tables, function(table) {
  fisher.test(table)$p.value
})

# 将列表转换为向量
fisher_results <- unlist(fisher_results)
# 计算FDR
fdr_results <- p.adjust(fisher_results, method = "BH")
# 计算Enrichment score
enrichment_scores <- lapply(motif_names, function(motif) {
  if (background_motif_counts[[motif]] > 0) {
    log2((interested_motif_counts[[motif]]/nrow(interested_motifs)) / (background_motif_counts[[motif]]/nrow(significant_motifs)))
  } else {
    NA  # 如果背景中的motif数量为0，则无法计算富集度得分，返回NA
  }
})
enrichment_scores <- unlist(enrichment_scores)
# 将结果转换为数据框，并排序
enrichment_results_df <- data.frame(
  Motif = motif_names,
  Enrichment.Score = enrichment_scores,
  Pvalue = fisher_results,
  FDR = fdr_results
)
sorted_enrichment_results <- enrichment_results_df[order(enrichment_results_df$FDR), ]

# 筛选结果，注释motif信息
filter_enrichment_results <- filter(enrichment_results_df,Pvalue<0.05,FDR<0.05,enrichment_scores>0)
motif_anno <- read.delim('meme_info.txt')
anno_enrichment_result <- merge(filter_enrichment_results,motif_anno,by='Motif')
ALL_anno_enrichment_result <- merge(sorted_enrichment_results,motif_anno,by='Motif')


write.table(ALL_anno_enrichment_result, 'N45_down_phase3.txt', sep = '\t', row.names = FALSE, quote = FALSE)


interested_genes <- read.delim('time_point_gene_list.txt',header=1) %>% 
  filter(trait=='N45',class=='down',phase==4) %>% 
  select('gene_id')# 替换为你的感兴趣基因
interested_genes <- interested_genes$gene_id

# 获取感兴趣的基因集中的motif
interested_motifs <- significant_motifs[significant_motifs$sequence_name %in% interested_genes, ]
# 计算感兴趣基因集中motif的数量
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算所有motif的总数
total_motif_counts <- table(significant_motifs$motif_id,exclude = NULL)
# 计算在sig中出现，但不在insterest中出现的motif
x <- interested_motifs %>% group_by(motif_id) %>% summarise(n=n())
y <- significant_motifs %>% group_by(motif_id) %>% summarise(n=n())
z <- y$motif_id %in% x$motif_id
t <- y$motif_id[!z]
# 将z添加值intersted_motif_counts中
for (non_motif in t) {
  interested_motifs <- interested_motifs %>% 
    add_row(motif_id=non_motif)
}
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算背景基因集中motif的数量（所有motif减去感兴趣的motif）
background_motif_counts <- total_motif_counts - interested_motif_counts

# 过滤掉在感兴趣基因集或背景基因集中计数为0的motif
motif_names <- names(interested_motif_counts)[interested_motif_counts > 0 & background_motif_counts > 0]



# 创建2x2列联表矩阵
# d <- data.frame(gene.not.interest=c(M-k, N-M-n+k), gene.in.interest=c(k, n-k))
# row.names(d) <- c("In_category", "not_in_category")
contingency_tables <- lapply(motif_names, function(motif) {
  matrix(c(interested_motif_counts[[motif]], # 背景，不同motif出现的数目
           background_motif_counts[[motif]]-interested_motif_counts[[motif]], # 感兴趣基因集中motif数目
           nrow(interested_motifs) - interested_motif_counts[[motif]],
           nrow(significant_motifs)  - nrow(interested_motifs) -  background_motif_counts[[motif]]#background_motif_counts[[motif]] + interested_motif_counts[[motif]] # 背景全部motif - 背景中单个motif出现的次数
           # 基因集中的motif - 基因集中单个motif出现的次数
  ), ncol = 2, byrow = TRUE)
})




# 应用Fisher精确检验，并检查每个列联表是否至少是2x2的
fisher_results <- lapply(contingency_tables, function(table) {
  fisher.test(table)$p.value
})

# 将列表转换为向量
fisher_results <- unlist(fisher_results)
# 计算FDR
fdr_results <- p.adjust(fisher_results, method = "BH")
# 计算Enrichment score
enrichment_scores <- lapply(motif_names, function(motif) {
  if (background_motif_counts[[motif]] > 0) {
    log2((interested_motif_counts[[motif]]/nrow(interested_motifs)) / (background_motif_counts[[motif]]/nrow(significant_motifs)))
  } else {
    NA  # 如果背景中的motif数量为0，则无法计算富集度得分，返回NA
  }
})
enrichment_scores <- unlist(enrichment_scores)
# 将结果转换为数据框，并排序
enrichment_results_df <- data.frame(
  Motif = motif_names,
  Enrichment.Score = enrichment_scores,
  Pvalue = fisher_results,
  FDR = fdr_results
)
sorted_enrichment_results <- enrichment_results_df[order(enrichment_results_df$FDR), ]

# 筛选结果，注释motif信息
filter_enrichment_results <- filter(enrichment_results_df,Pvalue<0.05,FDR<0.05,enrichment_scores>0)
motif_anno <- read.delim('meme_info.txt')
anno_enrichment_result <- merge(filter_enrichment_results,motif_anno,by='Motif')
ALL_anno_enrichment_result <- merge(sorted_enrichment_results,motif_anno,by='Motif')


write.table(ALL_anno_enrichment_result, 'N45_down_phase4.txt', sep = '\t', row.names = FALSE, quote = FALSE)


interested_genes <- read.delim('time_point_gene_list.txt',header=1) %>% 
  filter(trait=='N45',class=='down',phase==5) %>% 
  select('gene_id')# 替换为你的感兴趣基因
interested_genes <- interested_genes$gene_id

# 获取感兴趣的基因集中的motif
interested_motifs <- significant_motifs[significant_motifs$sequence_name %in% interested_genes, ]
# 计算感兴趣基因集中motif的数量
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算所有motif的总数
total_motif_counts <- table(significant_motifs$motif_id,exclude = NULL)
# 计算在sig中出现，但不在insterest中出现的motif
x <- interested_motifs %>% group_by(motif_id) %>% summarise(n=n())
y <- significant_motifs %>% group_by(motif_id) %>% summarise(n=n())
z <- y$motif_id %in% x$motif_id
t <- y$motif_id[!z]
# 将z添加值intersted_motif_counts中
for (non_motif in t) {
  interested_motifs <- interested_motifs %>% 
    add_row(motif_id=non_motif)
}
interested_motif_counts <- table(interested_motifs$motif_id,exclude = NULL)
# 计算背景基因集中motif的数量（所有motif减去感兴趣的motif）
background_motif_counts <- total_motif_counts - interested_motif_counts

# 过滤掉在感兴趣基因集或背景基因集中计数为0的motif
motif_names <- names(interested_motif_counts)[interested_motif_counts > 0 & background_motif_counts > 0]



# 创建2x2列联表矩阵
# d <- data.frame(gene.not.interest=c(M-k, N-M-n+k), gene.in.interest=c(k, n-k))
# row.names(d) <- c("In_category", "not_in_category")
contingency_tables <- lapply(motif_names, function(motif) {
  matrix(c(interested_motif_counts[[motif]], # 背景，不同motif出现的数目
           background_motif_counts[[motif]]-interested_motif_counts[[motif]], # 感兴趣基因集中motif数目
           nrow(interested_motifs) - interested_motif_counts[[motif]],
           nrow(significant_motifs)  - nrow(interested_motifs) -  background_motif_counts[[motif]]#background_motif_counts[[motif]] + interested_motif_counts[[motif]] # 背景全部motif - 背景中单个motif出现的次数
           # 基因集中的motif - 基因集中单个motif出现的次数
  ), ncol = 2, byrow = TRUE)
})




# 应用Fisher精确检验，并检查每个列联表是否至少是2x2的
fisher_results <- lapply(contingency_tables, function(table) {
  fisher.test(table)$p.value
})

# 将列表转换为向量
fisher_results <- unlist(fisher_results)
# 计算FDR
fdr_results <- p.adjust(fisher_results, method = "BH")
# 计算Enrichment score
enrichment_scores <- lapply(motif_names, function(motif) {
  if (background_motif_counts[[motif]] > 0) {
    log2((interested_motif_counts[[motif]]/nrow(interested_motifs)) / (background_motif_counts[[motif]]/nrow(significant_motifs)))
  } else {
    NA  # 如果背景中的motif数量为0，则无法计算富集度得分，返回NA
  }
})
enrichment_scores <- unlist(enrichment_scores)
# 将结果转换为数据框，并排序
enrichment_results_df <- data.frame(
  Motif = motif_names,
  Enrichment.Score = enrichment_scores,
  Pvalue = fisher_results,
  FDR = fdr_results
)
sorted_enrichment_results <- enrichment_results_df[order(enrichment_results_df$FDR), ]

# 筛选结果，注释motif信息
filter_enrichment_results <- filter(enrichment_results_df,Pvalue<0.05,FDR<0.05,enrichment_scores>0)
motif_anno <- read.delim('meme_info.txt')
anno_enrichment_result <- merge(filter_enrichment_results,motif_anno,by='Motif')
ALL_anno_enrichment_result <- merge(sorted_enrichment_results,motif_anno,by='Motif')


write.table(ALL_anno_enrichment_result, 'N45_down_phase5.txt', sep = '\t', row.names = FALSE, quote = FALSE)
