##AML»­Í¼
rm(list=ls())
install.packages("ggplot2")
library(ggplot2)
amlgo1111 <-as.matrix(read.csv("D:\\various cancer genes\\seventeen types cancer\\Endometrial cancer\\chong zuo\\kegg and mala inter enriched go terms.csv",header = TRUE))
AMLGO <-as.matrix(read.csv("D:\\various cancer genes\\seventeen types cancer\\Endometrial cancer\\chong zuo\\kmala deeper go terms.csv"))
A1 <- amlgo1111[amlgo1111[,1] %in% AMLGO[,1],]
A5 <-as.matrix(A1[,c(3:5)])
A5=apply(A5,2,as.numeric)
A5 <--log10(A5)
zz <-as.matrix(seq(1:19))
ZZ1 <-as.matrix(rep(1,94))
z <-as.matrix(rbind(zz,ZZ1))
A2 <-as.matrix(cbind(A1[,c(1:2)],A5))
A2 <-as.matrix(cbind(A2,z))
colnames(A2) <- c("ID","Description","K_P.adj","M_P.adj","min","-log10(P.adj)")
A11 <-data.frame(A2)
#A11»­go termµÄµãÍ¼
library(ggplot2)
A11$K_P.adj <-as.numeric(as.vector(A11$K_P.adj))
A11$M_P.adj <-as.numeric(as.vector(A11$M_P.adj))
pdf("D:\\My master thesis\\picture in sci\\shishiENDMC.pdf",width=6,height=7)
ggplot(A11,mapping = aes(x=K_P.adj,y=Description))+
  geom_point(aes(x=K_P.adj,y=Description,colour = "#61AC79"),size=1)+     #angle=90, 
  geom_point(aes(x=M_P.adj,y=Description,colour = "#F69790"),size=1)+
  theme(axis.text.x = element_text(hjust=1, vjust=.6,size=6,family="serif",face = "bold",colour = "black"),
        axis.text.y = element_text(family="serif",face = "bold",colour = "black",size=5),
        panel.background = element_rect(fill = "white", colour = "white"))+
  labs(x =  "-log10(P.adj)", y =NULL)+
  scale_x_continuous(breaks=seq(3, 19, 1))+
  theme(panel.grid.major.x = element_blank(),panel.grid.major.y = element_line(size = 0.3, linetype = 'solid',colour = "#6D9EC1"),
        axis.ticks = element_blank(),axis.title.y = element_text(colour="black",face="bold",family="serif",size=7),
        legend.text = element_text(family="serif",size = 6,face = "bold"),legend.title = element_text(size = 6,family="serif",face = "bold"))+
  scale_color_identity(name = "Database",
                       breaks = c("#61AC79", "#F69790"),
                       labels = c("KEGG", "MalaCards"),
                       guide = "legend")

dev.off()  



library(GO.db)
library(org.Hs.eg.db) ## human annotation package from bioconductor

getAllBPChildren <- function(goids)
{
  ans <- unique(unlist(mget(goids, GOBPCHILDREN), use.names=FALSE))
  ans <- ans[!is.na(ans)]
}

level3_terms <- getAllBPChildren("GO:0002376")
level4_terms <- getAllBPChildren(level3_terms)

level4_genes <- mget(intersect(level4_terms, keys(org.Hs.eg.db)), org.Hs.eg.db)
length(intersect(level3_terms, level4_terms))
