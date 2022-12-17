
########## Fig 2 #############



library(pheatmap)
library(RColorBrewer)
library('openxlsx')


rm(list=ls()) 


pool_meta2<-read.xlsx("heatmap.xlsx",sheet = 1,startRow = 1)
rock_meta2<-read.xlsx("heatmap.xlsx",sheet = 2,startRow = 1)


##top25
pre25<-read.xlsx("heatmap.xlsx",sheet = 3,startRow = 1)

pool_meta3<-merge(pool_meta2,pre25,by="id_order_pool")
pool_meta4<-pool_meta3[,c(-ncol(pool_meta3))]

rock_meta3<-merge(rock_meta2,pre25,by="id_order_rock")
rock_meta4<-rock_meta3[,c(-ncol(rock_meta3))]

pool<-pool_meta4

row.names(pool)<-pool[,1]
pool<-pool[,c(-1)]

class(pool)

colormap<-colorRampPalette(rev(brewer.pal(n=7,name="RdYlBu")))(100)
breaks=seq(min(unlist(c(pool))),max(unlist(c(pool))),length.out = 100)
breaks2=seq((-1),2,length.out = 100)
pool_heatmap<-pheatmap(pool,scale="row",color = colormap,breaks = breaks2,border_color = "grey",
            cluster_cols = F,
            
            main = "Pool",
            cutree_rows = 4,
            angle_col =0
)
pool_heatmap

rock<-rock_meta4

row.names(rock)<-rock[,1]
rock<-rock[,c(-1)]

class(rock)

colormap<-colorRampPalette(rev(brewer.pal(n=7,name="RdYlBu")))(100)
breaks=seq(min(unlist(c(rock))),max(unlist(c(rock))),length.out = 100)
breaks2=seq((-1),2,length.out = 100)
pool_heatmap<-pheatmap(rock,scale="row",color = colormap,breaks = breaks2,border_color = "grey",
                       cluster_cols = F,
                       
                       main = "Rock",
                       cutree_rows = 4,
                       angle_col =0
)
pool_heatmap


############# Fig 4  ################



library(ggplot2)
library(ggsci)
library(ggsignif)
library(ggpubr)

rm(list=ls()) 

pool_kegg_pathway<-read.xlsx("boxplot.xlsx",sheet =1,startRow = 1)
pool_kegg_pathway2<-pool_kegg_pathway

pool_rock_meta<-read.xlsx("boxplot.xlsx",sheet = 2,startRow = 1)


pool_rock_meta2<-scale(pool_rock_meta[,4:ncol(pool_rock_meta)],center = T,scale = T)
pool_rock_meta3<-cbind(pool_rock_meta[1:3],pool_rock_meta2)
colnames(pool_rock_meta3)<-colnames(pool_rock_meta)


key_ggplot<-pool_rock_meta3[,c(1:3,4)]
names(key_ggplot)[4]<-"meta"
key_ggplot$Compound_ID<-names(pool_rock_meta3)[4]

for (i in 4:ncol(pool_rock_meta3)){
  onedat<-pool_rock_meta3[,c(1:3,i)]
  names(onedat)[4]<-"meta"
  onedat$Compound_ID<-names(pool_rock_meta3)[i]
  key_ggplot<-rbind(key_ggplot,onedat)
}


pool_ggplot_pathway2<-merge(pool_kegg_pathway2,key_ggplot,by="Compound_ID")



mytheme <- theme_bw() + theme(axis.text.x = element_text(size=rel(1),colour="black")) + theme(axis.text.y = element_text(size=rel(1),colour="black")) 
mytheme <- mytheme + theme(axis.title.x = element_text(size=rel(1))) + theme(axis.title.y = element_text(size=rel(1)))
mytheme <- mytheme + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())



b<-pool_ggplot_pathway2


### Cellular_stress_response
b1<-subset(b,b$pathway=="cellular_stress_response")
b1$names <- factor(b1$names, levels=c("Pyridoxamine","Pyridoxal","4-Pyridoxic acid","Glutathione","Cys-Gly"), ordered=TRUE)


###  energy_metabolism
b1<-subset(b,b$pathway=="energy_metabolism")
b1$names <- factor(b1$names, levels=c("Pantothenate","Pantetheine","Phenylalanine","Adrenaline",
                                      "Citrate","Isocitrate","2-Oxoglutarate","Malate"), ordered=TRUE)

###  osmoregulation
b1<-subset(b,b$pathway=="osmoregulation")
b1$names <- factor(b1$names, levels=c("G6P","F6P","Threonine","Serine","Methionine"), ordered=TRUE)


###  pum
b1<-subset(b,b$pathway=="pum")
b1$names <- factor(b1$names, levels=c("GDP","Inosine","AMP","Hypoxanthine","Xanthine"), ordered=TRUE)


a<-ggplot(data=b1,aes(x=group,y=meta,fill=group,shape=microhabitat,color=group))+
  geom_boxplot(alpha=.8,outlier.size =0.8,lwd=0.4)+
  scale_fill_manual(values=c("#D73128","#E36F4F","#E08B52","#6997C5", "#5C8ABE","#4575B3"))+ 
  scale_color_manual(values=c("#D73128","#E36F4F","#E08B52","#6997C5", "#5C8ABE","#4575B3"))+ 
  facet_wrap(.~names,ncol = 5)+
  scale_shape_manual(values=c(1,15))+
  scale_y_continuous(limits=c(-2,4),breaks = c(-2,0,2,4))+
  xlab("Sampling times")+ylab("Normalized concentration")+
  scale_x_discrete(labels=c("B0"="D0","B3"="D3","B5"="D5","W0"="N0","W3"="N3","W5"="N5"))+
  mytheme+
  theme(legend.position = "none")
a



################# Fig 5 ################


rm(list=ls())


df<-read.xlsx("kegg.xlsx",
              sheet = 1,startRow = 1)



mytheme <- theme_bw() + theme(axis.text.x = element_text(face="bold",size=rel(1.6),colour="black")) + theme(axis.text.y = element_text(size=rel(1.6),colour="black")) 
mytheme <- mytheme + theme(axis.title.x = element_text(face="bold",size=rel(1.4))) + theme(axis.title.y = element_text(face="bold",size=rel(1.4)))
mytheme <- mytheme + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())


a<-ggplot() +
  geom_point(data=df,aes(x=kegg_pw,y=-log(p), fill=-log(p),size = Enrichment_Ratio),shape = 21, alpha=.8,color="black")+
  scale_shape_manual(values = c(21))+
  scale_size_area(max_size = 12)+
  facet_grid(rows = vars(microhabitat))+
  scale_fill_distiller(palette = "Spectral")+
  labs(size="Pathway impact")+

  mytheme
a


#############  Fig 6 ################
library(ggrepel)
rm(list=ls())


rock3<-read.xlsx("rda.xlsx",sheet = 1,startRow = 1)

env2<-read.xlsx("rda.xlsx",sheet = 2,startRow = 1)

rda_tb<-cca(spec~.,env2,na.action = na.omit)

ef=envfit(rda_tb,env2,permu=999)
ef
vif.cca(rda_tb)

rda_tb_forward_p <- ordistep(rda(spec~1, env2, scale = T), scope = formula(rda_tb), direction = 'forward', permutations = 999,na.action = na.omit)

mod2 <- rda(spec ~ dn + time + microhabitat  +  temperature  ,env2,scale = T)      
vif.cca(mod2) 
anova(mod2)   
summary(mod2) 


ef=envfit(mod2,env2,permu=999)
ef


library(rdacca.hp)
rdacca.hp(spec,env2,scale = T )



r2 <- RsquareAdj(mod2)
rda_noadj <- r2$r.squared  
rda_adj <- r2$adj.r.squared 
rda_adj
sp<- scores(mod2,choices = 1:2, display = 'sp')
lab5<-c(
  'pum','pum' ,'pum','pum','pum',
  'energy','energy' ,'energy','energy',
  'energy','energy' ,'energy','energy',
  'osmoregulation','osmoregulation' ,'osmoregulation','osmoregulation','osmoregulation',
  'Cellular_stress_response','Cellular_stress_response' ,'Cellular_stress_response','Cellular_stress_response','Cellular_stress_response'
)
sp<-as.data.frame(sp)
sp$lab<-lab5
cn<-mod2$CCA$biplot[,1:2]
cn<-as.data.frame(cn)
row.names(cn)
row.names(sp)

0.7073 * rda_adj 
0.1511 * rda_adj 
myCol <- c("darkolivegreen", "blue2", "yellow2","coral2","plum3")
windowsFonts(TNM= windowsFont("Arial"))

b<-ggplot(data=sp) +
  geom_hline(yintercept=0, linetype=2,color='black') + 
  geom_vline(xintercept=0, linetype=2,color='black') + 
  geom_point(mapping=aes(x=RDA1,y=RDA2,fill = lab,shape=lab),size =5,alpha=.7)+ 
  geom_segment(data = cn,aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(length =unit(0.3,"cm"),type = "closed",angle=22.5),linetype=1,
               colour ="black",size=0.6)+ 
  geom_text_repel(data =cn,aes(x=RDA1/2,y=RDA2/2,label=row.names(cn)),size=3,
                  color='black')+ 
  theme_bw()+labs( x= "RDA1 (14.15%)", y ="RDA2 (3.02%)",color ='')+
  scale_fill_brewer(palette = "Set3")+
  scale_shape_manual(values = c(21,22,23,24,25))+
  theme(axis.title = element_text(color='black',size=10.5),
        axis.text = element_text(color='black',size=10.5),
        panel.grid = element_blank())+
  theme(legend.position = "right")
b



#########  PERMANOVA   #############
library('openxlsx')
library('vegan')


rm(list=ls())

df<-read.xlsx("permanova.xlsx",sheet = 1,startRow = 1)

df.env<-read.xlsx("permanova.xlsx" ,sheet = 2,startRow = 1)

df.dist <- vegdist(df, method="bray")

dune.div <- adonis2(df ~ time*microhabitat, data = df.env, permutations = 999, method="bray")
dune.div