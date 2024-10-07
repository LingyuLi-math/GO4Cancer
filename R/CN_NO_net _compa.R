library(GOSemSim)
library(GOSim)
library(ggplot2)
library(carData)
library(corrplot)
# Install
if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/ggpubr")
library(ggpubr)
library(plyr)
library(gridExtra)

setOntology(ont = "BP", loadIC=TRUE, DIR=NULL)
#H <-getAncestors()
HH <-getOffsprings()
HHH <-unlist(HH)
d <- godata('org.Hs.eg.db', ont="BP", computeIC=FALSE)
a <-as.matrix(read.csv("D:\\My master thesis\\picture in sci\\spe and con\\ba dao shi unique.csv",header = F))
aa <-as.matrix(read.table("D:\\my postgraduate project\\cancer pathway genes\\sixteen go terms.txt",sep = ","))
amlsim <-as.matrix(mgoSim(a,aa[,1],semData=d,measure = "Wang",combine=NULL))
write.csv(amlsim,"D:\\My master thesis\\picture in sci\\go4cancer\\consimi.csv")
amlsim <-amlsim[-34,]
corrplot(amlsim)  #[c(1:10),]
dasim <- as.matrix(apply(amlsim,1,max))
dasimc <-as.matrix(apply(amlsim,2,max))
congoterm <-as.matrix(rbind(dasim,dasimc))
type1 <-as.matrix(rep("consistent GO terms",50))
congoterm1 <-data.frame(cbind(congoterm,type1))
colnames(congoterm1) <-c("sim","type")
congoterm1$sim <-as.numeric(as.vector(congoterm1$sim))
congoterm1hou <-as.matrix(congoterm1)

dasimean <-mean(dasim)
dasimcmean <-mean(dasimc)

###随机go terms 1次
set.seed(1)
rand <-as.matrix(sample(HHH, size=41, replace = FALSE))
row.names(rand) <-NULL
randsim <-as.matrix(mgoSim(rand,aa[,1],semData=d,measure = "Wang",combine=NULL))
randsimr <- as.matrix(apply(randsim,1,max))
randsimrc <-as.matrix(apply(randsim,2,max))



randomterm <-as.matrix(rbind(randsimr,randsimrc))
type2 <-as.matrix(rep("random GO terms",51))
randomterm1 <-data.frame(cbind(randomterm,type2))
randomterm1$X1 <-as.numeric(as.vector(randomterm1$X1))
randomterm1hou <-as.matrix(randomterm1)
###一致的go term和随机的整合
boxdata <-data.frame(rbind(congoterm1hou,randomterm1hou))
boxdata$sim <-as.numeric(as.vector(boxdata$sim))


#randsimrmean <-mean(randsimr)
#randsimrcmean <-mean(randsimrc)
###画箱线图
pdf("D:\\My master thesis\\picture in sci\\go4cancer\\xiangshiyitu32230ci.pdf",width=6,height=5) #

compare_means(sim~type, data = boxdata,method = "t.test")
# Change method
my_comparisons <- list( c("consistent GO terms", "random GO terms") )
ggboxplot(boxdata, x = "type", y = "sim",
          color = "type", palette = "jco",
          add = "jitter", fill = c("#B5D2BB", "#59B6CA"))+   #"jco", "ucscgb", "uchicago", "simpsons" and "rickandmorty".
  stat_compare_means(comparisons = my_comparisons,method = "t.test")+ # Add pairwise comparisons p-value
  stat_summary(fun = mean,geom = "point", shape = 18, size = 4, color = "#E57373")# Add global p-value

dev.off()  



###normalcancer specific sin 箱线图
a <-as.matrix(read.csv("D:\\various cancer genes\\seventeen types cancer\\Endometrial cancer\\kminter\\czigong0.8 GO.csv",header = F))
b <-as.matrix(read.csv("D:\\various cancer genes\\seventeen types cancer\\Endometrial cancer\\kminter\\nzigong0.8 GO.csv",header = F))
endspe <-as.matrix(read.csv("D:\\My master thesis\\picture in sci\\spe and con\\ENDMC\\ENDMC.csv",header = F))

amlsim <-as.matrix(mgoSim(a,endspe,semData=d,measure = "Wang",combine=NULL))
amlsimN <-as.matrix(mgoSim(b,endspe,semData=d,measure = "Wang",combine=NULL))
amlsim <-amlsim[,-2]
amlsimN <-amlsimN[,-2]
dasim <- as.matrix(apply(amlsim,1,max))
dasimc <-as.matrix(apply(amlsim,2,max))
congoterm <-as.matrix(rbind(dasim,dasimc))
type1 <-as.matrix(rep("GO terms of disease network",31))
congoterm1 <-data.frame(cbind(congoterm,type1))
colnames(congoterm1) <-c("similarity","type")
congoterm1$similarity <-as.numeric(as.vector(congoterm1$similarity))
congoterm1hou <-as.matrix(congoterm1)

###正常网络

randsimr <- as.matrix(apply(amlsimN,1,max))
randsimrc <-as.matrix(apply(amlsimN,2,max))
randomterm <-as.matrix(rbind(randsimr,randsimrc))
type2 <-as.matrix(rep("GO terms of normal network",18))
randomterm1 <-data.frame(cbind(randomterm,type2))
randomterm1$X1 <-as.numeric(as.vector(randomterm1$X1))
randomterm1hou <-as.matrix(randomterm1)

###疾病和正常组合
boxdata <-data.frame(rbind(congoterm1hou,randomterm1hou))
boxdata$similarity <-as.numeric(as.vector(boxdata$similarity))

###画箱线图
pdf("D:\\various cancer genes\\seventeen types cancer\\Endometrial cancer\\kminter\\cnxiangxiantuhhh.pdf",width=6,height=5) #

compare_means(similarity~type, data = boxdata,method = "t.test")
# Change method
my_comparisons <- list( c("GO terms of disease network", "GO terms of normal network") )
ggboxplot(boxdata, x = "type", y = "similarity",
          color = "type", palette = "jco",
          add = "jitter", fill = c("#FF6666", "#009999"))+   #"jco", "ucscgb", "uchicago", "simpsons" and "rickandmorty".
  stat_compare_means(comparisons = my_comparisons,method = "t.test")+ # Add pairwise comparisons p-value
  stat_summary(fun = mean,geom = "point", shape = 18, size = 4, color = "purple")# Add global p-value

dev.off()  


pdf("D:\\various cancer genes\\seventeen types cancer\\Endometrial cancer\\kminter\\cnxiangxiantu.pdf",width=6,height=4) #
p1<-ggplot(data=boxdata,aes(x=type,y=sim,fill=type))+geom_boxplot()+
  labs(x = NULL, y = "similarity")+
  theme(axis.text.x = element_text(size=10,family="serif",face = "bold",colour = "black"),
        axis.text.y = element_text(family="serif",face = "bold",colour = "black",size=10))+
  scale_fill_manual(values=c("#FF6666", "#009999"))

p1
dev.off()  
