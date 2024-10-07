###随机go terms ，取30次，取均值
randr <-matrix(data=0,nrow=41,ncol=1)
randc <-matrix(data=0,nrow=10,ncol=1)
set.seed(1)
for (i in 1:60) {
  bianl <- paste("rand",i,sep = "")
  bianl <-as.matrix(sample(HHH, size=41, replace = FALSE))
  if (sum(is.na(bianl))==0) {
    row.names(bianl) <-NULL
    rands <- paste("randsim",i,sep = "")
    rands <-as.matrix(mgoSim(bianl,aa[,1],semData=d,measure = "Wang",combine=NULL))
    randsimr <- as.matrix(apply(rands,1,max))
    randsimrc<-as.matrix(apply(rands,2,max))
    randr <-as.matrix(cbind(randr,randsimr))
    randc <-as.matrix(cbind(randc,randsimrc))
    
  }
}
randr <-randr[,-1]
randc <-randc[,-1]
randsimr1 <-as.matrix(apply(randr,1,mean))
randsimrc1 <-as.matrix(apply(randc,1,mean))

randomterm <-as.matrix(rbind(randsimr1,randsimrc1))
type2 <-as.matrix(rep("random GO terms",51))
randomterm1 <-data.frame(cbind(randomterm,type2))
randomterm1$X1 <-as.numeric(as.vector(randomterm1$X1))
randomterm1hou <-as.matrix(randomterm1)


###一致的go term和随机的整合
boxdata <-data.frame(rbind(congoterm1hou,randomterm1hou))
boxdata$sim <-as.numeric(as.vector(boxdata$sim))


randsimrmean <-mean(randsimr1)
randsimrcmean <-mean(randsimrc1)
