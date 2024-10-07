rm(list=ls())
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GO.db")
library(stringr)
library(GOSemSim)
library(GOSim)
d <- godata('org.Hs.eg.db', ont="BP", computeIC=FALSE)
setOntology(ont = "BP", loadIC=TRUE, DIR=NULL)
H <-getAncestors()
HH <-getOffsprings()
amlgo <-as.matrix(read.csv("D:\\various cancer genes\\Thyroid cancer\\chong zuo\\kegg and mala inter enriched go terms.csv",header = TRUE))
amlgo <-as.matrix(amlgo[,1])
#去冗余函数
deredun <-function(amlgo){
  amlsim <-as.matrix(mgoSim(amlgo,amlgo,semData=d,measure = "Wang",combine=NULL))
  amlsim[upper.tri(amlsim)] <-0
  diag(amlsim)<-0
  DF = amlsim
  value = 0.7
  which.names <- function(DF, value){
    ind <- which(DF>=value, arr.ind=TRUE)
    print(paste(rownames(DF)[ind[,"row"]],  colnames(DF)[ind[,"col"]], sep=' '))
  }
  amlgo0.7 <-as.matrix(which.names(DF,value))
  amlgo0.7 <-str_split_fixed(amlgo0.7, " ", 2)
  #将相似性大于0.7的go term中specific的go terms保留
  a <-amlgo0.7
  yiyang <- c()
  for (i in 1:nrow(a)) {
    if(a[i,2]%in%H[a[i,1]][[1]]) {you <-a[i,1]}
    if(a[i,2]%in%HH[a[i,1]][[1]]){you <-a[i,2]}
    if(!(a[i,2]%in%HH[a[i,1]][[1]])&!(a[i,2]%in%H[a[i,1]][[1]])){you <-rbind(a[i,1],a[i,2])}
    yiyang <-as.matrix(rbind(you,yiyang))
  }
  yiyang1 <-unique(yiyang)
  b <-unique(as.matrix(rbind(as.matrix(amlgo0.7[,1]),as.matrix(amlgo0.7[,2]))))
  meiyouxiang <-as.matrix(setdiff(amlgo,b))
  labc <- list(yiyang1=yiyang1,b=b,meiyouxiang=meiyouxiang)
}

###
for (j in 1:10) {
  ww[j] <-deredun(amlgo)
  you[j] <-as.matrix(ww$meiyouxiang)
  xiayilun[j] <-as.matrix(ww$yiyang1)
  if(xiayilun[j]==xiayilun[j-1]){
    
    AMLGO <-unique(as.matrix(rbind(xiayilun[j-1],you1,you2,you3)))
  }
  
}
###
ww <-deredun(aaaa1)
you1 <-as.matrix(ww$meiyouxiang)
xiayilun <-as.matrix(ww$yiyang1)

ww1 <-deredun(xiayilun)
you2 <-as.matrix(ww1$meiyouxiang)
xiayilun2 <-as.matrix(ww1$yiyang1)

ww2 <-deredun(xiayilun2)
you3 <-as.matrix(ww2$meiyouxiang)
xiayilun3 <-as.matrix(ww2$yiyang1)

ww3 <-deredun(xiayilun3)
you4 <-as.matrix(ww3$meiyouxiang)
xiayilun4 <-as.matrix(ww3$yiyang1)

ww4 <-deredun(xiayilun4)
you5 <-as.matrix(ww4$meiyouxiang)
xiayilun5 <-as.matrix(ww4$yiyang1)

ww5 <-deredun(xiayilun5)
you6 <-as.matrix(ww5$meiyouxiang)
xiayilun6 <-as.matrix(ww5$yiyang1)

ww6 <-deredun(xiayilun6)
you7 <-as.matrix(ww6$meiyouxiang)
xiayilun7 <-as.matrix(ww6$yiyang1)


AMLGO <-unique(as.matrix(rbind(xiayilun4,you1,you2,you3,you4)))
colnames(AMLGO) <-NULL
write.csv(AMLGO,"D:\\various cancer genes\\Small cell lung cancer\\chong zuo\\kmala deeper go terms.csv",row.names =FALSE)
