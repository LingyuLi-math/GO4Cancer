###kegg and malacard genes enriched go terms
rm(list=ls())
getwd()
setwd("D:\\various cancer genes\\Acute myeloid leukemia")
source("https://bioconductor.org/biocLite.R")
library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)
library(GOSim)
datarty <-as.matrix(read.table("leukemiagenes.txt",header = FALSE))
datartym <-as.matrix(read.table("malacard genes.txt",header = FALSE))
kmintergene <-as.matrix(intersect(datarty,datartym))
HGenriched <-function(datarty,datartym){
  a=1
  if(a==1){
    datarty$V1 <- as.character(datarty$V1) #需要character格式，然后进行ID转化
    #将SYMBOL格式转为ENSEMBL和ENTERZID格式 
    test1 = bitr(kmintergene, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")
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
    pdf("D:\\My master thesis\\picture in sci\\KEGG_BP.pdf")
    ego_BPP <-as.matrix(ego_BP)
    barplot(ego_BP,showCategory=3,title="EnrichmentGO_BP")#条状图，按p从小到大排，绘制前20个Term
    dev.off()
    dev.off()
    barplot(ego_BP[,2])
    write.csv(ego_BP,"chong zuo\\kmintergenes enriched go terms.csv",row.names =FALSE)
  }
  a=2
  if(a==2){
    datartym$V1 <- as.character(datartym$V1) #需要character格式，然后进行ID转化
    #将SYMBOL格式转为ENSEMBL和ENTERZID格式 
    test1 = bitr(datartym$V1, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")
    ego_BP1 <- enrichGO(gene = test1$ENTREZID, 
                       universe = background$ENTREZID , #背景基因集
                       OrgDb = org.Hs.eg.db, #没有organism="human"，改为OrgDb=org.Hs.eg.db
                       #keytype = 'ENSEMBL',
                       ont = "BP", #也可以是 CC  BP  MF中的一种
                       pAdjustMethod = "BH", #矫正方式 holm”, “hochberg”, “hommel”, “bonferroni”, “BH”, “BY”, “fdr”, “none”中的一种
                       pvalueCutoff = 0.0001, #P值会过滤掉很多，可以全部输出
                       qvalueCutoff = 0.05,
                       readable = TRUE) #Gene ID 转成gene Symbol ，易读
    write.csv(ego_BP1,"chong zuo\\mala enriched go terms.csv",row.names =FALSE)
  }
  
  interkm <-as.matrix(intersect(as.matrix(ego_BP[,1]),as.matrix(ego_BP1[,1])))
  kinter <-as.matrix(ego_BP[ego_BP[,1] %in% interkm,])
  minter <-as.matrix(ego_BP1[ego_BP1[,1] %in% interkm,])
  kinter <-as.matrix(kinter[order(kinter[,1]),])
  minter <-as.matrix(minter[order(minter[,1]),])
  kinter <-kinter[,c(1:2,6)]
  minter <-minter[,c(1:2,6)]
  interkmpp <-as.matrix(cbind(kinter,minter[,3]))
  rownames(interkmpp) <-NULL
  colnames(interkmpp) <-c("ID","Description","k_p.adjust","m_p.adjust")
  write.csv(interkmpp,"chong zuo\\kegg and mala inter enriched go terms.csv",row.names =FALSE)
  return(interkmpp)
}
a1 <-HGenriched(datarty,datartym)
###去除go terms冗余
a01 <-as.matrix(read.csv("chong zuo\\kegg and mala inter enriched go terms.csv"))
goandminp <-function(a01){
  a0134mean <-a01[,c(3,4)]
  a0134mean1 <-apply(a0134mean, 2, as.numeric)
  a0134mean2 <-as.matrix(apply(a0134mean1,1,min))
  a01 <-as.matrix(cbind(a01,a0134mean2))
  colnames(a01) <-c("ID","Description","k_p.adjust","m_p.adjust","min_p.adjust")
  write.csv(a01,"chong zuo\\kegg and mala inter enriched go terms.csv",row.names =FALSE)
  a0115 <-a01[,c(1,5)]
  write.csv(a0115,"chong zuo\\kegg and mala inter enriched go terms and min p.csv",row.names =FALSE)
  return(a0115)
}
jg <-goandminp(a01)
###生成聚类矩阵
##用REVIGO去完冗余，对应原来的
setwd("D:\\various cancer genes")
B1 <-as.matrix(read.csv("chong zuo\\REVIGO_0.7.csv"))
By <-as.matrix(read.csv("chong zuo\\kegg and mala inter enriched go terms and min p.csv"))
notf <-as.matrix(c("GO:0104004","GO:0062012","GO:0062013")) #,"GO:0110110","GO:0150063","GO:0099177","GO:0062009"
#xuyaogo <-as.matrix(B1[,1])
xuyaogo <-as.matrix(rbind(notf,as.matrix(B1[,1])))
by <-as.matrix(By[By[,1] %in% xuyaogo,])
write.csv(by,"chong zuo\\kegg and mala shao enriched go terms and min p.csv",row.names =FALSE)
##生成矩阵
data1 <-as.matrix(read.csv("Acute myeloid leukemia\\chong zuo\\kegg and mala inter enriched go terms.csv"))
data1 <-as.matrix(data1[,c(1,5)])
colnames(data1) <-NULL
data2 <-as.matrix(read.csv("Basal cell carcinoma\\chong zuo\\kegg and mala shao enriched go terms and min p.csv"))
data3 <-as.matrix(read.csv("Bladder cancer\\chong zuo\\kegg and mala inter enriched go terms and min p.csv"))
data4 <-as.matrix(read.csv("Breast cancer\\chong zuo\\kegg and mala shao enriched go terms and min p.csv"))
data5 <-as.matrix(read.csv("Chronic myeloid leukemia\\chong zuo\\kegg and mala shao enriched go terms and min p.csv"))
data6 <-as.matrix(read.csv("Colorectal cancer\\chong zuo\\kegg and mala elite inter enriched go terms and min p.csv"))
data7 <-as.matrix(read.csv("Endometrial cancer\\chong zuo\\kegg and mala shao enriched go terms and min p.csv"))
data8 <-as.matrix(read.csv("Gastric cancer\\chong zuo\\kegg and mala shao enriched go terms and min p.csv"))
data9 <-as.matrix(read.csv("Glioma\\chong zuo\\kegg and mala shao enriched go terms and min p.csv"))
data10 <-as.matrix(read.csv("Hepatocellular carcinoma\\chong zuo\\kegg and mala shao enriched go terms and min p.csv"))
data11 <-as.matrix(read.csv("Non-small cell lung cancer\\chong zuo\\kegg and mala shao enriched go terms and min p.csv"))
data12 <-as.matrix(read.csv("Pancreatic cancer\\chong zuo\\kegg and mala shao enriched go terms and min p.csv"))
data13 <-as.matrix(read.csv("Prostate cancer\\chong zuo\\kegg and mala shao enriched go terms and min p.csv"))
data14 <-as.matrix(read.csv("Renal cell carcinoma\\chong zuo\\kegg and mala inter enriched go terms and min p.csv"))
data15 <-as.matrix(read.csv("Thyroid cancer\\chong zuo\\kegg and mala inter enriched go terms and min p.csv"))
allgo <-as.matrix(rbind(as.matrix(data1[,1]),as.matrix(data2[,1]),as.matrix(data3[,1]),as.matrix(data4[,1]),as.matrix(data5[,1]),as.matrix(data6[,1]),as.matrix(data7[,1]),as.matrix(data8[,1]),as.matrix(data9[,1]),as.matrix(data10[,1]),as.matrix(data11[,1]),as.matrix(data12[,1]),as.matrix(data13[,1]),as.matrix(data14[,1]),as.matrix(data15[,1])))
allgo <-as.matrix(unique(allgo))
#对应
hmatrix <-c("AML","BCC","BLC","BC","CML","CRC","ENDMC","GASC","Glioma","HCC","NSCLC","PNCA","PC","RCC","THCA")
m1 <- matrix(1,nrow=567,ncol=15,dimnames=list(allgo,hmatrix))
m1 <-apply(m1,2,as.character )
rownames(m1) <-allgo
duiy <- function(data1,e,m1) {
  for (i in 1:length(allgo)) {
    for (j in 1:nrow(data1)) {
      if(data1[j,1]==rownames(m1)[i]){
        m1[i,e]=data1[j,2]
        }
    }
  }
  return(m1)
}
c1 <-duiy(data15,15,c1)
###画聚类图
rrr<-rownames(c1)
c1 <-apply(c1,2,as.numeric)
X <-log10(c1)
X <-abs(X)
rownames(X) <-rrr 
write.csv(X,"picture of my paper/log hou go and p.csv")
rc <- rainbow(nrow(X), start=0, end=.3)
cc <- rainbow(ncol(X), start=0, end=.3)
library("gplots")
col=bluered
dev.off()  ##消除热图过大现象
my_palette <- colorRampPalette(c("aliceblue","red" ,"blue","green","yellow"))(n = 567)
resu <-heatmap.2(X,col=my_palette,scale="none",trace="none",density.info = "none",
                 key.title = NA,key.xlab =NA,dendrogram = "none",labRow=FALSE,
                 lhei = c(0.7,7))   #lhei = c(0.7,7),,lwid=c(1,3)
x_reorder <- X[rev(resu$rowInd), resu$colInd]
# x_reorder 为重新排序后的数据
library(pRoloc)
goidx_reorder <-as.matrix(row.names(x_reorder))
gotote <- as.matrix(goIdToTerm(goidx_reorder))
colnames(gotote) <-"go terms"
x_reordert <-as.matrix(cbind(gotote,x_reorder))
write.csv(x_reordert,"picture of my paper/cancer15 cluster go terms and p.csv")
x_consis <-as.matrix(x_reordert[c(484:567),])
write.csv(x_consis,"picture of my paper/cancer15 consistancy go terms and p.csv")
x_consisim <-as.matrix(rownames(x_consis))
###计算这些的相似性
library(GOSemSim)
library(GOSim)
d <- godata('org.Hs.eg.db', ont="BP", computeIC=FALSE)
smiallc <-as.matrix(mgoSim(rownames(hccgo),cj,semData=d,measure = "Wang",combine=NULL))
write.csv(smiall,"picture of my paper/go terms84 similarity each.csv")
cell10go <-as.matrix(read.table("D:\\my postgraduate project\\cancer pathway genes\\sixteen go terms.txt",header = FALSE,sep = ","))
cell10go <-as.matrix(cell10go[,1])
smiall_ten <-as.matrix(mgoSim(x_consisim,cell10go,semData=d,measure = "Wang",combine=NULL))
write.csv(smiall_ten,"picture of my paper/go terms84 and cell 10.csv")

###画84个go terms聚类散点图
library(factoextra)
df <- scale(smiall) 
head(df, n = 5)
fviz_nbclust(df, kmeans, method = "wss") + geom_vline(xintercept = 9, linetype = 2)
set.seed(123)
km_result <- kmeans(df, 40, nstart = 24)
print(km_result)
dd <- cbind(smiall, cluster = km_result$cluster)
head(dd)
table(km_result$cluster)
fviz_cluster(km_result, data = df,
             palette = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07","#FF0000","#00FF00","#FF7F00","#A67D3D","#DB70DB"),
             ellipse.type = "euclid",
             star.plot = FALSE, 
             repel = FALSE,
             ggtheme = theme_minimal()
)
smiall1 <-as.data.frame(smiall)
ggplot(smiall1)
