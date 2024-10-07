getwd()
library(reshape2)
a <-read.csv("D:\\My master thesis\\enrich go seven.csv")
a <-a[c(1:17),]
#将数据转换为长形式reshape2包。
mdat = melt(a,id.vars="Type",measure.vars=c("inter","MalaCards","KEGG"))
pdf("D:\\My master thesis\\picture in sci\\enrich go seventheen1.pdf",width=2.5,height=1)
plot_1 = ggplot(mdat, aes(x=Type, y=value, fill=variable)) +
  theme_bw() +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values=c("#D3C5AC","#8898C1","#3D5F99"))+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.6,size=2,family="serif",face = "bold",colour = "black",margin = margin(t = 0.4)),
       axis.text.y = element_text(family="serif",face = "bold",colour = "black",size=2,margin = margin(r = 0.4)),
       panel.background = element_rect(fill = "white", colour = "white"))+
  theme( plot.background = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.border = element_blank(),panel.grid.major.x = element_blank(),panel.grid.major.y = element_blank(),
        axis.title.y = element_text(colour="black",face="bold",family="serif",size=2),
        legend.text = element_text(family="serif",size = 6,face = "bold"),legend.title = element_text(size = 6,family="serif",face = "bold"))+
  theme(axis.line = element_line(color = 'black',size =0.08),axis.ticks.x = element_line(size = 0.08),
        axis.ticks.y = element_line(size = 0.08),axis.ticks.length =unit(.025, "cm"))+
  scale_x_discrete(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0),breaks=seq(0,3000,200))+
  labs(x = NULL, y = "Number of enriched GO terms")+
  scale_color_identity(labels = c("Enriched GO terms in both KEGG and MalaCards","Enriched GO terms in MalaCards","Enriched GO terms in KEGG"),
                       guide = "legend")
#"#3D5F99","#8898C1","#D3C5AC"
plot_1
dev.off()

