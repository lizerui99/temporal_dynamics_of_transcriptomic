rm(list=ls())

library(WGCNA)
# setwd('C:/Users/Public/Desktop/Training/wgcna')
options(stringsAsFactors = FALSE) # 允许使用多线程
enableWGCNAThreads()

RNAseq=read.delim('n45wgcna.tab',header=T,sep="\t",check.names=F) # 第一行为表头，第一列是行的名称
RNAseq <- RNAseq[!duplicated(RNAseq$Gene_ID),]
row.names(RNAseq) <- RNAseq$Gene_ID
RNAseq_voom <- RNAseq[-1]
# 4-8g，10000 ，32 - 20000多
#RNAseq_voom=RNAseq_voom[seq(1,2000,by=1),] # 为了计算的快，使用2000个基因


datExprOri=t(RNAseq_voom)    #转置，每行为一个样品，每列为一个基因

#======================= 过滤样本 过滤缺失基因数目多于10%的样本：===========
nGenes = ncol(datExprOri)      #基因数目
NumberMissingBySample=apply(is.na(data.frame(datExprOri)),1, sum)   #计算每行NA的数目
KeepSample= NumberMissingBySample<0.1*nGenes #判断哪些样品的缺失率<10%，
table(KeepSample)  #过滤的样品数目统计

datExpr=datExprOri[KeepSample,]  # 表达矩阵过滤样品
SampleName=rownames(datExprOri)[KeepSample] #样品名称变量过滤样品

#=====================过滤方差为0的基因，以及缺失样本超过10%的基因=====
nSamples = nrow(datExpr)  # 统计样品数目
variancedatExpr=as.vector(apply(as.matrix(datExpr),2,var, na.rm=T))  #按列（基因）取方差
no.missingdatExpr=as.vector(apply(is.na(as.matrix(datExpr)),2, sum) )#按列（基因）统计缺失数目
KeepGenes= variancedatExpr>0 & no.missingdatExpr<0.1*nSamples  # 保留方差不等于0，且缺失低于10%的基因
table(KeepGenes)# 过滤统计

datExpr=datExpr[, KeepGenes] #过滤基因
GeneName=colnames(datExpr) #过滤前的基因名称


#根据层次聚类绘制样本树
tree=hclust(dist(datExpr),method = 'average')
pdf(file='sample clustering.pdf',w=10,h=7)
plot(tree,xlab="", sub="", cex = 0.7,main="Sample Clustering") 
dev.off()


#================设置网络构建参数选择范围======
powers = c(1:20)  #设置beta值的取值范围，临界举证β--软阈值
cor <- WGCNA::cor 
sft = pickSoftThreshold(datExpr, RsquaredCut = 0.9,powerVector = powers, verbose = 5, networkType = 'unsigned')  #选择软阈值,powervector,verbose输出
sft$powerEstimate # 选取的β值

#每个候选软阈值的计算R2绘图
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# 绘制一条设定的阈值线
abline(h=0.90,col="red")
abline(h=0.85,col="blue")

#绘制每个候选beta的平均连接度
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
    main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


#===============自动构建WGCNA模型==================
#========表达矩阵转换成邻接矩阵，然后再将邻接矩阵转换成拓扑矩阵,识别模块=============
net = blockwiseModules(
  datExpr,
  power = 7,#sft$powerEstimate,             #软阈值，前面计算出来的，β值
  maxBlockSize = 10000,                  #最大block大小，将所有基因放在一个block中，3w 需要32g内存
  TOMType = "unsigned",                  #选择unsigned，使用标准TOM矩阵
  deepSplit = 2, minModuleSize = 50,     #剪切树参数，deepSplit取值0-4
  mergeCutHeight = 0.25,                 #模块合并参数，越大模块越少
  numericLabels = TRUE,                             # T返回数字，F返回颜色
  pamRespectsDendro = FALSE,  
  saveTOMs = TRUE,                       #计算后可以存储，然后载入 
  saveTOMFileBase = "FPKM-TOM",
  loadTOMs = FALSE,
  verbose = 3                            # 
)

#获取模块颜色
moduleColors = labels2colors(net$colors)
#获取每个模块的特征值
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes # 计算模块特征值
MEs = orderMEs(MEs0); #重新排序，相似颜色的挨在一起

# 输出每个基因所在的模块，以及与该模块的KME值
#file.remove('All_Gene_KME.txt')
for(module in substring(colnames(MEs),3)){
  if(module == "grey") next
  ME=as.data.frame(MEs[,paste("ME",module,sep="")])
  colnames(ME)=module
  datModExpr=datExpr[,moduleColors==module]
  datKME = signedKME(datModExpr, ME)
  datKME=cbind(datKME,rep(module,length(datKME)))
  write.table(datKME,quote = F,row.names = T,append = T,file = "All_Gene_KME.txt",col.names = F)
}

#===============模块绘图==================
#将表达矩阵转换为一个颜色矩阵，使用log10（FPKM+1）
expColor=t(numbers2colors(log10(datExpr+1),colors=blueWhiteRed(100),naColor="grey"))
colnames(expColor)=rownames(datExpr)
#绘制基因的树形图，模块图，以及每个样品的表达图
pdf("wgcna.dendroColors.pdf",height = 7,width = 9)

# 绘制
plotDendroAndColors(net$dendrograms[[1]], 
                    colors=cbind(moduleColors[net$blockGenes[[1]]],expColor),
                    c("Module",colnames(expColor)),
                    dendroLabels = F, hang = 0.03,
                    addGuide = T, guideHang = 0.05,
                    cex.rowText=0.5)
dev.off()

#绘制两两模块间的邻接矩阵
pdf("wgcna.adjacency.heatmap.pdf",height = 10,width = 9)
plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap",plotDendrograms = F,
                      marDendro = c(4,4,2,4))
dev.off()


#绘制所有模块的表达值热图与特征值条形图
for(module in substring(colnames(MEs),3)){
  if(module == "grey") next
  ME=MEs[,paste("ME",module,sep="")]
  pdf(paste("wgcna.", module, ".express.barplot.pdf", sep=""),height = 7,width = 9)
  par(mfrow=c(2,1),mar=c(0.3,5.5,3,2))
  plotMat(t(scale(datExpr[,moduleColors==module])),
          rlabels=T,main=module,cex.main=2,clabels=T)
  
  par(mar=c(5,4.2,0,0.7))
  barplot(ME,col=module,main="",cex.main=2,ylab="eigengene expression",xlab="sample")
  dev.off()
}
#将结果保存成cytoscape的输入文件格式

#读入构建时保存的TOM矩阵
load('FPKM-TOM-block.1.RData') 
TOM=as.matrix(TOM)

#或重新计算TOM矩阵
#TOM = TOMsimilarityFromExpr(datExpr, power =sft$powerEstimate,TOMType = "unsigned"); 

for(module in substring(colnames(MEs),3)){
  if(module == "grey") next
  probes = colnames(datExpr)
  inModule = is.finite(match(moduleColors, module))
  modProbes = probes[inModule]
  modTOM = TOM[inModule, inModule]
  dimnames(modTOM) = list(modProbes, modProbes)
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = paste("CytoscapeInput-edges-", module, ".txt", sep=""),
                                 nodeFile = paste("CytoscapeInput-nodes-", module, ".txt", sep=""),
                                 weighted = TRUE,
                                 threshold = 0.15,
                                 nodeNames = modProbes,
                                 nodeAttr = moduleColors[inModule])
}

#与表型的关联
allTraits<-read.delim('pheno.txt',header=T,sep="\t",check.names = F)  #读取表型值

#去除删除的样本
traitRows = match(SampleName, allTraits[,1]);
datTraits = allTraits[traitRows, -1];    
rownames(datTraits) = allTraits[traitRows, 1];

#计算ME与表型的相关系数，并计算p值
moduleTraitCor = cor(MEs, datTraits, use = "p") # use=p代表去掉缺失计算
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)#

#整理要显示在图中的数字,第一行为相关系数，第二行为p值
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)

#绘制关联热图
pdf("wgcna.Module-trait.heatmap.pdf", width = 10, height = 15)
par(mar = c(6, 8.8, 3, 2.2))  
labeledHeatmap(Matrix = moduleTraitCor,   #相关系数
               xLabels = colnames(datTraits), #x轴为表型
               yLabels = names(MEs),#y轴为模块
               ySymbols = names(MEs),
               colorLabels = FALSE, 
               colors = blueWhiteRed(50),
               textMatrix = textMatrix, #每个单元格的内容
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-1,1),
               main = "Module-trait relationships")
dev.off()

adjacency = adjacency(datExpr,power=9)
kIM <- intramodularConnectivity(adjacency, moduleColors)
kIM$module <- moduleColors
write.csv(kIM,'kim.csv')

