rm(list=ls())# "清除环境"
getwd()
setwd("D:\\various cancer genes\\seventeen types cancer\\Acute myeloid leukemia")
source("https://bioconductor.org/biocLite.R")
library(DOSE)
library(org.Hs.eg.db)
library(topGO)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")
library(clusterProfiler) #安装富集分析包，若已安装，无需再次安装
library(pathview)
library(GOSim)
datarty <-as.matrix(read.table("leukemiagenes.txt",header = FALSE)) #读取数据
datartym <-as.matrix(read.table("malacard genes.txt",header = FALSE))
kmintergene <-as.matrix(intersect(datarty,datartym))
test1 = bitr(kmintergene, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")#将SYMBOL格式转为ENSEMBL和ENTERZID格式 
    allgenes <-read.table("D:\\my postgraduate project\\hcc\\HiSeqV2",header = TRUE)
    allgenes <-as.matrix(allgenes[,1]) #富集分析的背景基因集
    background <-bitr(allgenes, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")
    ego_BP <- enrichGO(gene = test1$ENTREZID, 
                       universe = background$ENTREZID , #背景基因集
                       OrgDb = org.Hs.eg.db, #没有organism="human"，改为OrgDb=org.Hs.eg.db
                       #keytype = 'ENSEMBL',
                       ont = "BP", #也可以是 CC  BP  MF中的一种
                       pAdjustMethod = "BH", #矫正方式 holm”, “hochberg”, “hommel”, “bonferroni”, “BH”, “BY”, “fdr”, “none”中的一种
                       pvalueCutoff = 0.05, #P值会过滤掉很多，可以全部输出
                       qvalueCutoff = 0.05,
                       readable = TRUE) #Gene ID 转成gene Symbol ，易读
    dotplot(ego_BP,title="EnrichmentGO_BP_dot")#点图，按富集的数从大到小的
    barplot(ego_BP,showCategory=20,title="EnrichmentGO_BP")#条状图，按p从小到大排，绘制前20个Term
    ##可视化--
    plotGOgraph(ego_BP)
    dev.off() #关闭图形界面
    