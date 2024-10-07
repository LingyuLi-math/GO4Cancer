rm(list=ls())
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
library(data.table)
library(DESeq2)
library(readr)
library(caret)
library(pheatmap) 
library(doMC)
library(tidyr)
library(stringr)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggplot2)
registerDoMC(cores = 4)
getwd()
setwd("D:/my postgraduate project")
rawda <- read.table("D:\\my postgraduate project\\hcc\\LIHC_raw_counts.txt",sep="\t",header = TRUE)
#找出正常组织样本矩阵
x <- colnames(rawda)  
geneID <- rawda[, 1]
colnames(rawda)[1] <- "geneID" 
x <- x[-1]
x <- substr(x, 14, 15) 
x <- as.integer(x)
x <- x %in% (10:19)
y <- colnames(rawda)[-1]
y <- y[x]
tumor_matched_normal_counts <- cbind(geneID, subset(rawda, select = y))

#筛选对应癌症组织样本矩阵
x <- colnames(tumor_matched_normal_counts)[-1]
y <- substr(x, 9, 12)
z <- subset(rawda, select = grep(paste(y, collapse = "|"), colnames(rawda), value = TRUE))
colnames(z)
tumor_counts <- subset(z, select = setdiff(colnames(z), colnames(tumor_matched_normal_counts)))
tumor_counts <- cbind(geneID, tumor_counts)
colnames(tumor_counts)
cts_LIHC <- cbind(tumor_matched_normal_counts, tumor_counts[-1]) ##表达矩阵
colnames(cts_LIHC)
rm(x, y, z)
gc()

coldata_LIHC <- data.frame(c(rep("tumor_matched_normal",50), rep("tumor",50)))
row.names(coldata_LIHC) <- colnames(cts_LIHC)[-1]
colnames(coldata_LIHC) <- "condition"
head(coldata_LIHC)
rownames(cts_LIHC) <- cts_LIHC$geneID
cts_LIHC <- cts_LIHC[, -1]
cts_LIHC <- round(cts_LIHC)  ##对数据取整

#countData-表达矩阵；colData-分组矩阵；design-分组矩阵中包含组别信息的列名。
dds_LIHC <- DESeqDataSetFromMatrix(countData = cts_LIHC, colData = coldata_LIHC, design = ~condition)
dds_LIHC 
dds_LIHC$condition <- relevel(dds_LIHC$condition, ref = "tumor_matched_normal")##将对照组放在首位
dds_LIHC <- DESeq(dds_LIHC)

results_LIHC <- results(dds_LIHC, alpha = 0.01)##设置FDR cutoff值为0.05(5 %)
table(dds_LIHC$padj<0.05)
diff_gene_deseq2 <-subset(results_LIHC, padj < 0.05 & abs(log2FoldChange) > 1)
dim(diff_gene_deseq2)
write.csv(diff_gene_deseq2,file= "LIHC DEG_treat_vs_control.csv")
#####以上为基因差异性表达



#相关系数分析
myfile<- read.csv("LIHC DEG_treat_vs_control.csv",header = TRUE)
myfile1 <-as.matrix(myfile[,1])
myfile2 <-myfile1[1:5]

mysubf <-str_sub(myfile2,3,8)
gene.df <- bitr(mysubf, fromType = "ENTREZID", #fromType是指你的数据ID类型是属于哪一类的
               toType = "SYMBOL", #toType是指你要转换成哪种ID类型，可以写多种，也可以只写一种
               OrgDb = org.Hs.eg.db)#Orgdb是指对应的注释包是哪个
myfile3 <-as.matrix(myfile1[6:4191])
colnames(myfile3) <- 'lm'
MYFILEEE <- myfile[-c(1:5),]
hhh <-sub("\\|.*", "", myfile3)
hhh1 <-as.matrix(gene.df[,2])
hhh2 <-as.matrix(rbind(hhh,hhh1))
sample11 <-read.table("D:\\my postgraduate project\\hcc\\bq.txt",sep="\t",header = TRUE)
mysample <-t(read.table("D:\\my postgraduate project\\hcc\\bq.txt",sep="\t",header = TRUE))
myfile <- read.table("D:/my postgraduate project/hcc/HiSeqV2",header = TRUE,sep ="\t")
head(myfile)
myfile1 <-myfile[,-1]
rtitle <-myfile[,1]
colnames(myfile1) <- mysample
rownames(myfile1) <- rtitle
dim(myfile1)
head(myfile1)
#overl <-intersect(hhh2,hccp1)
myfile1 <- myfile1[(rownames(myfile1) %in% hccp1),]
head(myfile1)
cfile <- myfile1[, (colnames(myfile1) %in% 0)]
nfile <- myfile1[, (colnames(myfile1) %in% 1)]
#head(nfile)
#计算疾病基因间的相关系数
nfilet <-t(nfile)
cfilet <-t(cfile)
matrixN=cor(nfilet)
matrix=cor(cfilet)
diag(matrixN)=0
#matrix2 <- data.frame(matrix)
hist(matrixN, col="skyblue", xlab="coefficient", main="correlation histogram",plot = TRUE,freq = FALSE)
#rug(jitter(matrixN))
lines(density(matrixN),col ="red")
curve(dnorm(matrix, mean=mean(matrix), sd=sd(ww)), add=TRUE, col="darkblue", lwd=2)
#qqnorm(matrixN)
length(matrixN)
mean(matrixN) 
ww <-as.vector(matrixN)
le <-mean(matrix)-1.96*sd(ww)
re <-mean(matrix)+1.96*sd(ww)#,μ+1.96σ）内的面积为95.449974%。
rcn <-qnorm(0.975, mean = mean(matrixN), sd = sd(ww), lower.tail = TRUE, log.p = FALSE) #的返回值是给定概率p后的下分位点
lcn <-qnorm(0.025, mean = mean(matrixN), sd = sd(ww), lower.tail = TRUE, log.p = FALSE)
diag(matrix)=0
#rrr<-matrix(which[!is.na(matrix)])
#write.table(matrix,"coefficient_matrix.txt",sep="\t")   
DF = matrixN
value = lcn
which.names <- function(DF, value){
  ind <- which(DF<=value, arr.ind=TRUE)
  print(paste(rownames(DF)[ind[,"row"]],  colnames(DF)[ind[,"col"]], sep=' '))
}
namereN <-as.matrix(which.names(DF,value))
namere1N <-as.matrix(which.names(DF,value))
namere2N<- rbind(namereN,namere1N)
write.table(namere2N,"Cedn.txt",sep="\t")   
REGENES <- read.table("Cedn.txt",header =TRUE)
#namere1 <-as.matrix(namere)
Cedge <-separate(REGENES, col = V1, into = c("V1","V2"), sep = " ")
write.table(Cedge,"cedn.txt",sep="\t")   

fff <-read.table("cedn.txt")
write.table(fff,"cancer pathway genes/Cedgen.txt",sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE)
