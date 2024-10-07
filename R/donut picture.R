library(reshape2)
library(ggplot2)
a <-data.frame(read.csv("D:\\My master thesis\\picture in sci\\seventeen genes.csv"))
a$kzhanbi = a$KEGG / sum(a$KEGG)
a$mzhanbi = a$MalaCards/ sum(a$MalaCards)
a$ezhanbi = a$Elite.genes.in.MalaCards / sum(a$Elite.genes.in.MalaCards)
mdat = melt(a,id.vars="Name",measure.vars=c("ezhanbi","mzhanbi","kzhanbi"))
mdat$value <-as.numeric(as.vector(mdat$value))
mdat$width <-c(rep(.8,17),rep(2,17),rep(4,17))
mdat$width <-as.numeric(as.vector(mdat$width))
pdf("D:\\My master thesis\\picture in sci\\donut picture.pdf")
p <-ggplot(mdat, aes(x=variable, fill=Name, y=value,width=1)) +
  geom_bar(aes(width=width), position="stack", stat="identity")+
  scale_x_discrete(limits = c(" "," ","ezhanbi","mzhanbi"," ","kzhanbi"))+
  scale_fill_manual(values=c("#C7B797","#A1ADCD","#8E76A7","#D4D9E9","#8E8886","#E3EDE4","#755591",
                             "#F7A5A0","#C6B696","#647FAD","#9FD2DE","#FAC1B8","#F48680","#B5A9C8",
                             "#33B375","#D2C5AC","#A5C261"))+
  geom_segment(x=3, y=0, xend=3,yend=1,colour="#524B48",size=0.3)+
  geom_segment(x=2.6, y=0, xend=2.6,yend=1,colour="#524B48",size=0.3)+
  coord_polar("y")+
  theme_void()
B = seq(0.9541029,1,0.001835883)
for (i in B) {
  p <- p + geom_segment(x=4, y=i, xend=8,yend=i,colour="#B8A47C",size=0.02)
}  
B2 =seq(0.9102921,0.9541029,0.00292072)
for (i in B2) {
  p <- p + geom_segment(x=4, y=i, xend=8,yend=i,colour="#43649B",size=0.02)
} 
B3 =seq(0.8817803,0.9102921,0.001055993)
for (i in B3) {
  p <- p + geom_segment(x=4, y=i, xend=8,yend=i,colour="#552D77",size=0.01)
} 
B4 =seq(0.779555,0.8817803,0.00159727)
for (i in B4) {
  p <- p + geom_segment(x=4, y=i, xend=8,yend=i,colour="#8C9BC2",size=0.01)
} 
B5 =seq(0.7267038,0.779555,0.003303199)
for (i in B5) {
  p <- p + geom_segment(x=4, y=i, xend=8,yend=i,colour="#544D4A",size=0.01)
} 
B6 =seq(0.6668985,0.7267038,0.001128402)
for (i in B6) {
  p <- p + geom_segment(x=4, y=i, xend=8,yend=i,colour="#9EC6A9",size=0.01)
} 
B7 =seq(0.6265647,0.6668985,0.001833355)
for (i in B7) {
  p <- p + geom_segment(x=4, y=i, xend=8,yend=i,colour="#050306",size=0.01)
} 
B8 =seq(0.5229486,0.6265647,0.002656823)
for (i in B8) {
  p <- p + geom_segment(x=4, y=i, xend=8,yend=i,colour="#F26D6D",size=0.01)
} 
B9 =seq(0.4707928,0.5229486,0.002745041)
for (i in B9) {
  p <- p + geom_segment(x=4, y=i, xend=8,yend=i,colour="#756D6A",size=0.01)
}
B10 =seq(0.3539639,0.4707928,0.00208623)
for (i in B10) {
  p <- p + geom_segment(x=4, y=i, xend=8,yend=i,colour="#0C4887",size=0.01)
}
B11 =seq(0.3038944,0.3539639,0.006258693)
for (i in B11) {
  p <- p + geom_segment(x=4, y=i, xend=8,yend=i,colour="#06A6BE",size=0.01)
}
B12 =seq(0.2579973,0.3038944,0.001240462)
for (i in B12) {
  p <- p + geom_segment(x=4, y=i, xend=8,yend=i,colour="#F5908A",size=0.01)
}
B13 =seq(0.2058415,0.2579973,0.001448771)
for (i in B13) {
  p <- p + geom_segment(x=4, y=i, xend=8,yend=i,colour="#9F536B",size=0.01)
}
B14 =seq(0.1383867,0.2058415,0.001163014)
for (i in B14) {
  p <- p + geom_segment(x=4, y=i, xend=8,yend=i,colour="#755793",size=0.01)
}
B15 =seq(0.09040339,0.1383867,0.003198887)
for (i in B15) {
  p <- p + geom_segment(x=4, y=i, xend=8,yend=i,colour="#029047",size=0.01)
}
B16 =seq(0.02573023,0.09040339,0.01293463)
for (i in B16) {
  p <- p + geom_segment(x=4, y=i, xend=8,yend=i,colour="#C7B898",size=0.01)
}
B17 =seq(0,0.02573023,0.001513543)
for (i in B17) {
  p <- p + geom_segment(x=4, y=i, xend=8,yend=i,colour="#90C4A2",size=0.01)
}
p

dev.off()





