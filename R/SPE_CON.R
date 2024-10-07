##批量读kmala deeper go terms
getwd()
setwd("D:/various cancer genes/seventeen types cancer")
doc.names <- dir("D:/various cancer genes/seventeen types cancer")
doc.path <- sapply(doc.names,function(names) paste(names,"chong zuo/kmala deeper go terms.csv",sep='/'))
doc <- sapply(doc.path, function(doc) as.matrix(read.csv(doc)))

doc.path1 <-sapply(doc.names,function(names) paste(names,"chong zuo/kegg and mala inter enriched go terms.csv",sep='/'))
doc1 <- sapply(doc.path1, function(doc1) as.matrix(read.csv(doc1)))
c1<-as.matrix(read.csv("D:\\various cancer genes\\seventeen types cancer\\Acute myeloid leukemia\\chong zuo\\kegg and mala inter enriched go terms.csv",header = TRUE))
c2<-as.matrix(read.csv("D:\\various cancer genes\\seventeen types cancer\\Basal cell carcinoma\\chong zuo\\kegg and mala shao enriched go terms and min p.csv",header = TRUE))
c3<-as.matrix(read.csv("D:\\various cancer genes\\seventeen types cancer\\Bladder cancer\\chong zuo\\kegg and mala inter enriched go terms.csv",header = TRUE))
c4<-as.matrix(read.csv("D:\\various cancer genes\\seventeen types cancer\\Breast cancer\\chong zuo\\kegg and mala shao enriched go terms and min p.csv",header = TRUE))
c5<-as.matrix(read.csv("D:\\various cancer genes\\seventeen types cancer\\Chronic myeloid leukemia\\chong zuo\\kegg and mala shao enriched go terms and min p.csv",header = TRUE))
c6<-as.matrix(read.csv("D:\\various cancer genes\\seventeen types cancer\\Colorectal cancer\\chong zuo\\kegg and mala inter enriched go terms.csv",header = TRUE))
c7<-as.matrix(read.csv("D:\\various cancer genes\\seventeen types cancer\\Endometrial cancer\\chong zuo\\kegg and mala shao enriched go terms and min p.csv",header = TRUE))
c8<-as.matrix(read.csv("D:\\various cancer genes\\seventeen types cancer\\Gastric cancer\\chong zuo\\kegg and mala shao enriched go terms and min p.csv",header = TRUE))
c9<-as.matrix(read.csv("D:\\various cancer genes\\seventeen types cancer\\Glioma\\chong zuo\\kegg and mala shao enriched go terms and min p.csv",header = TRUE))
c10<-as.matrix(read.csv("D:\\various cancer genes\\seventeen types cancer\\Hepatocellular carcinoma\\chong zuo\\kegg and mala shao enriched go terms and min p.csv",header = TRUE))
c11<-as.matrix(read.csv("D:\\various cancer genes\\seventeen types cancer\\Melanoma\\chong zuo\\kegg enriched go terms.csv",header = TRUE))
c12<-as.matrix(read.csv("D:\\various cancer genes\\seventeen types cancer\\Non-small cell lung cancer\\chong zuo\\kegg and mala shao enriched go terms and min p.csv",header = TRUE))
c13<-as.matrix(read.csv("D:\\various cancer genes\\seventeen types cancer\\Pancreatic cancer\\chong zuo\\kegg and mala shao enriched go terms and min p.csv",header = TRUE))
c14<-as.matrix(read.csv("D:\\various cancer genes\\seventeen types cancer\\Prostate cancer\\chong zuo\\kegg and mala shao enriched go terms and min p.csv",header = TRUE))
c15<-as.matrix(read.csv("D:\\various cancer genes\\seventeen types cancer\\Renal cell carcinoma\\chong zuo\\kegg and mala inter enriched go terms.csv",header = TRUE))
c16<-as.matrix(read.csv("D:\\various cancer genes\\seventeen types cancer\\Small cell lung cancer\\chong zuo\\kegg enriched go terms.csv",header = TRUE))
c17<-as.matrix(read.csv("D:\\various cancer genes\\seventeen types cancer\\Thyroid cancer\\chong zuo\\kegg and mala inter enriched go terms.csv",header = TRUE))
#####
a1 <-as.matrix(doc$`Acute myeloid leukemia`)
a1q<-as.matrix(c1[c1[,1] %in% a1,c(1,5)])
a1z <-as.matrix(cbind(a1[,2],a1q))
a2 <-as.matrix(doc$`Basal cell carcinoma`)
a2q<-as.matrix(c2[c2[,1] %in% a2,])
a2z <-as.matrix(cbind(a2[,2],a2q))

a3 <-as.matrix(doc$`Bladder cancer`)
a3q<-as.matrix(c3[c3[,1] %in% a3,c(1,5)])
a3z <-as.matrix(cbind(a3[,2],a3q))

a4 <-as.matrix(doc$`Breast cancer`)
a4q<-as.matrix(c4[c4[,1] %in% a4,])
a4z <-as.matrix(cbind(a4[,2],a4q))

a5 <-as.matrix(doc$`Chronic myeloid leukemia`)
a5q<-as.matrix(c5[c5[,1] %in% a5,])
a5z <-as.matrix(cbind(a5[,2],a5q))

a6 <-as.matrix(doc$`Colorectal cancer`)
a6q<-as.matrix(c6[c6[,1] %in% a6,c(1,5)])
a6z <-as.matrix(cbind(a6[,2],a6q))

a7 <-as.matrix(doc$`Endometrial cancer`)
a7q<-as.matrix(c7[c7[,1] %in% a7,])
a7z <-as.matrix(cbind(a7[,2],a7q))

a8 <-as.matrix(doc$`Gastric cancer`)
a8q<-as.matrix(c8[c8[,1] %in% a8,])
a8z <-as.matrix(cbind(a8[,2],a8q))

a9 <-as.matrix(doc$Glioma)
a9q<-as.matrix(c9[c9[,1] %in% a9,])
a9z <-as.matrix(cbind(a9[,2],a9q))

a10 <-as.matrix(doc$`Hepatocellular carcinoma`)
a10q<-as.matrix(c10[c10[,1] %in% a10,])
a10z <-as.matrix(cbind(a10[,2],a10q))

a11 <-as.matrix(doc$Melanoma)
a11q<-as.matrix(c11[c11[,1] %in% a11,c(1,6)])
a11z <-as.matrix(cbind(a11[,2],a11q))

a12 <-as.matrix(doc$`Non-small cell lung cancer`)  
a12q<-as.matrix(c12[c12[,1] %in% a12,])###
a12z <-as.matrix(cbind(a12[,2],a12q))

a13 <-as.matrix(doc$`Pancreatic cancer`)
a13q<-as.matrix(c13[c13[,1] %in% a13,])
a13z <-as.matrix(cbind(a13[,2],a13q))

a14 <-as.matrix(doc$`Prostate cancer`)
a14q<-as.matrix(c14[c14[,1] %in% a14,])
a14z <-as.matrix(cbind(a14[,2],a14q))

a15 <-as.matrix(doc$`Renal cell carcinoma`)
a15q<-as.matrix(c15[c15[,1] %in% a15,c(1,5)])
a15z <-as.matrix(cbind(a15[,2],a15q))

a16 <-as.matrix(doc$`Small cell lung cancer`)
a16q<-as.matrix(c16[c16[,1] %in% a16,c(1,6)])
a16z <-as.matrix(cbind(a16[,2],a16q))

a17 <-as.matrix(doc$`Thyroid cancer`)
a17q<-as.matrix(c17[c17[,1] %in% a17,c(1,5)])
a17z <-as.matrix(cbind(a17[,2],a17q))
a <-as.matrix(rbind(a1z,a2z,a3z,a4z,a5z,a6z,a7z,a8z,a9z,a10z,a11z,a12z,a13z,a14z,a15z,a16z,a17z))
b  <-as.matrix(cbind(as.matrix(a[,2]),as.matrix(a[,1])))
write.csv(a,"D:\\My master thesis\\picture in sci\\17 go terms specificity consistencyzzz.csv",row.names =FALSE,quote=FALSE)


####每个go terms对应几个疾病
amlgo <-as.matrix(read.csv("D:\\My master thesis\\picture in sci\\17 go terms specificity consistencyzzz.csv",header = TRUE))
www <-as.matrix(unique(amlgo[,2]))
qq1 <-c()
for (i in 1:566) {
  w=as.matrix(apply(amlgo,1,function(x) any(x==www[i,1])))
  qq <-as.matrix(amlgo[w,])
  if (nrow(qq)==3&ncol(qq)==1) {
    qq <-t(qq)
    qq1 <-as.matrix(rbind(qq,qq1))
  }
}
write.csv(qq1,"D:\\My master thesis\\picture in sci\\spe and con\\one you.csv",row.names =FALSE)

a11111 <-as.matrix(unique(a[,2]))
