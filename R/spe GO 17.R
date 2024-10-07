
# 画特异性结果AML ---------------------------------------------------------------


a <-as.matrix(read.csv("D:\\My master thesis\\picture in sci\\spe and con\\Glioma\\Glioma zhi.csv",header=F))
gotote <- as.matrix(goIdToTerm(a[,2]))
a <-as.matrix(a[,-1])
aml <-as.matrix(cbind(gotote,a))
aml <-as.matrix(aml[,-2])
write.csv(aml,"D:\\My master thesis\\picture in sci\\spe and con\\Glioma\\YOU go term.csv")
A5 <-as.matrix(aml[,3])
A5=apply(A5,2,as.numeric)
A5 <--log10(A5)
AML <-as.matrix(cbind(aml[,c(1:2)],A5))
colnames(AML) <- c("Description","ID","P.adj")
AML1 <-data.frame(AML)
library(dplyr)
AML1 <-unite(AML1,"Description",Description, ID)
AML1$P.adj <-as.numeric(as.vector(AML1$P.adj))

df1<-AML1%>%
  arrange(desc(P.adj))
library(tidyverse)
pdf("D:\\My master thesis\\picture in sci\\spe and con\\PC\\tiaoxingtu.pdf",width=12,height=7)
df1 %>% 
  ggplot(aes(reorder(Description, P.adj), P.adj,width = 0.6)) + 
  geom_col(aes(fill = P.adj)) + 
  scale_fill_gradient2(low = "red", 
                       high = "#F701C6" ) +    #midpoint = median(df1$P.adj)由一个色到另一个色
  coord_flip() + 
  labs(x = "Description",y = "-log10(P.adj)")
dev.off() 
 
 





