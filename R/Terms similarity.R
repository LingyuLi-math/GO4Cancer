getwd()
setwd("D:/various cancer genes")
library(GOSemSim)
library(GOSim)
ce <- as.matrix(read.csv("lihc edge terms.csv",header = TRUE))
cn <-as.matrix(read.csv("lihc node terms.csv",header = TRUE))
ne <- as.matrix(read.csv("lihcn edge terms.csv",header = TRUE))
nn <-as.matrix(read.csv("lihcn node terms.csv",header = TRUE))
ce <-as.matrix(ce[,2])
cn <-as.matrix(cn[,2])
ne <-as.matrix(ne[,2])
nn <-as.matrix(nn[,2])
allterms<-unique(rbind(ce,cn,ne,nn))

d <- godata('org.Hs.eg.db', ont="BP", computeIC=TRUE)
simmmmmm <-goSim("GO:0001709","GO:0001708",semData=d,measure = "Lin")
smie <-mgoSim(ne,nn,semData=d,measure = "Wang",combine="BMA")
smin <-mgoSim(cn,nn,semData=d,measure = "Wang",combine="BMA")
smiall <-mgoSim(allterms,allterms,semData=d,measure = "Wang",combine=NULL)
unen <-unique(rbind(ce,cn))
sixtgo <-read.table("sixteen go terms.txt",header = FALSE,sep = ",")
sixgoid <-as.matrix(sixtgo[,1])
rowNAME <-read.csv("initial cluster matrix.csv")
rowNAME <-as.matrix(rowNAME[,1])
smiall1 <-mgoSim(rowNAME,rowNAME,semData=d,measure = "Wang",combine=NULL)
goSim("GO:0001953","GO:0001952",semData=d,measure = "Wang")
pac=c("GO:0007265","GO:0048011","GO:0016772","GO:0016303","GO:0004713")
smiPPP <-mgoSim(pac,sixgoid,semData=d,measure = "Wang",combine="BMA")
hccg <-read.table("hcc.txt",header = FALSE,sep = ";")
hccg <-as.matrix(hccg[,1])
head(hccg)
#ddd <-sub("^\\S+\\s+", '', hccg)
head(ddd)
write.table(hccg,"station gene.txt",sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE)
hccg <-read.table("station gene.txt",header = FALSE,sep = " ")
hccg <-as.matrix(hccg[,9])
write.table(hccg,"hcc 168genes.txt",quote = FALSE,row.names = FALSE,col.names = FALSE)
enrichnetRE <-read.table("enrichnet_ranking_table.txt",header = TRUE)
enrichnetRE <-as.matrix(enrichnetRE[,1])
write.table(enrichnetRE,"enrichnet go.txt",quote = FALSE,row.names = FALSE,col.names = FALSE)
enrichgo <-read.table("enrichnet go.txt",header = FALSE,sep = "~")
enrichgo <-enrichgo[c(1:20),]
enrichgo <-as.matrix(enrichgo[,1])
smiall1 <-mgoSim(ce,enrichgo,semData=d,measure = "Wang",combine=NULL)
smiPPP <-mgoSim(enrichgo,sixgoid,semData=d,measure = "Wang",combine="BMA")
cancersim <-mgoSim(allgo,allgo,semData=d,measure = "Wang",combine=NULL)
uni <-unique(allgo)
cancersimuni <-mgoSim(uni,uni,semData=d,measure = "Wang",combine=NULL)
upp  <-as.matrix(cancersimuni[upper.tri(cancersimuni)])
#allanc <-as.matrix(largec[!allanc])
hist(upp, col="skyblue", xlab="sem coefficient", main="semsim histogram",plot = TRUE,freq = FALSE)
#rug(jitter(matrixN))
lines(density(upp),col ="red")
Ucn <-qnorm(0.05, mean = mean(upp), sd = sd(upp), lower.tail = FALSE, log.p = FALSE)
cancersimuni[lower.tri(cancersimuni)] <- 0  
diag(cancersimuni) <-0
DF = cancersimuni
value = lcn
which.names <- function(DF, value){
  ind <- which(DF>=value, arr.ind=TRUE)
  print(paste(rownames(DF)[ind[,"row"]],  colnames(DF)[ind[,"col"]], sep=' '))
}
namesim <-as.matrix(which.names(DF,value))



###compare enriched go terms of kegg and malacard
KEGGgo17 <-as.matrix(read.table("Thyroid cancer\\go ID.txt",header = FALSE))
malacardgo17 <-as.matrix(read.table("Thyroid cancer\\malacard ID.txt",header = FALSE))
BLC <-mgoSim(KEGGgo17,malacardgo17,semData=d,measure = "Wang",combine=NULL)
write.csv(BLC,"Thyroid cancer\\kegg and malacard similarity.csv")
###求kegg和malacard的go term 与six go terms 的相似性
A1 <-as.matrix(read.table("Thyroid cancer\\go ID.txt",header = FALSE))
A1m <-as.matrix(read.table("Thyroid cancer\\malacard ID.txt",header = FALSE))
A1me <-as.matrix(read.table("Thyroid cancer\\malacard elite ID.txt",header = FALSE))
allthreeid <-rbind(A1,A1m,A1me)
A1GOS <-mgoSim(A1A,sixchpa,semData=d,measure = "Wang",combine=NULL)
write.csv(smiall1,"218 terms similarity.csv")
###整体的相似性

#allthreeid <-as.matrix(read.table("Basal cell carcinoma\\all three ID.txt",header = FALSE))
allthreeid <-unique(allthreeid)
malacardgo <-as.matrix(read.table("Bladder cancer\\malacard ID.txt",header = FALSE))
malaelitego <-as.matrix(read.table("Bladder cancer\\malacard elite ID.txt",header = FALSE))
allthsim <-as.matrix(mgoSim(allthreeid,allthreeid,semData=d,measure = "Wang",combine=NULL))
write.csv(allthsim,"Thyroid cancer\\all three kinds similarity.csv")
kmesimilar <-as.matrix(mgoSim(KEGGgo,malaelitego,semData=d,measure = "Wang",combine="BMA"))
kmallgsimilar <-rbind(kmsimilar,kmesimilar)
colnames(kmallgsimilar) <-"kegg"
write.csv(kmallgsimilar,"Thyroid cancer\\kmallsimilar.csv",row.names = c("mala","malae"))

##聚类分析
#install.packages("factoextra")
library(factoextra)
data <- as.matrix(read.csv("D:/various cancer genes/malacard terms similarity.csv",header = TRUE,row.names = 1))
df <- scale(data)
head(df, n = 5)
fviz_nbclust(df, kmeans, method = "wss") + geom_vline(xintercept = 4, linetype = 2)
set.seed(123)   #设置随机数种子，保证实验的可重复进行
km_result <- kmeans(df, 8, nstart = 25)  #利用k-means是进行聚类
print(km_result)
dd <- as.data.frame(cbind(data, cluster = km_result$cluster))
head(dd)
table(dd$cluster)
#进行可视化展示
fviz_cluster(km_result, data = df,
             palette = c("#2E9FDF", "#00AFBB", "#E7B800","#FC4E07","#E7B888","#5C4B07","#10FC40","#A10E48","#B357FB"), #"#FC4E07","#E7B888","#FC4B07"
             ellipse.type = "euclid",
             star.plot = TRUE, 
             repel = TRUE,
             ggtheme = theme_minimal()
)

##层次聚类
#先求样本之间两两相似性
result <- dist(df, method = "euclidean")
#产生层次结构
result_hc <- hclust(d = result, method = "ward.D2")
#进行初步展示
fviz_dend(result_hc, cex = 0.6)
fviz_dend(result_hc, k = 4, 
          cex = 0.5, 
          k_colors = c("#2E9FDF", "#00AFBB", "#E7B800","#FC4E07","#E7B888","#5C4B07","#10FC40","#A10E48","#B357FB"),
          color_labels_by_k = TRUE, 
          rect = TRUE          
)
###
library(gplots)
x <- as.matrix(mtcars)
rc <- rainbow(nrow(x), start=0, end=.3)
cc <- rainbow(ncol(x), start=0, end=.3)
###每个疾病去冗余
data <-read.csv("Thyroid cancer/Thyroid-enrich.csv",header = TRUE)  
data1 <-read.csv("Thyroid cancer/malacard enriched.txt.csv",header = TRUE)
data2 <-read.csv("Thyroid cancer/malacard elite enrich.csv",header = TRUE)
DATA <-rbind(data,data1,data2)
DATA <-DATA[,c(1,5)]
write.csv(DATA,"Thyroid cancer\\three go and p.csv")
data <-as.matrix(data[,1])
data1 <-as.matrix(data1[,1])
data2 <-as.matrix(data2[,1])
inter2 <-as.matrix(intersect(data1,data2))

