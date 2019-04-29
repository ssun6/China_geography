library(vegan)
library(ggpubr)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(tidyverse)
library(rstatix)   

setwd("/users/ssun5/china_git/China_geography/PCoA_PERMANOVA/")
gen_tab=read.table(file="/Users/shansun/Google\ Drive/China/otu_table_12500_wtax_L6.txt",sep="\t",header=T,row.names=1)
gen_ab=gen_tab*12500
tab=log10(gen_tab*12500+1)
map=read.csv(file="/Users/shansun/Google\ Drive/China/Metadata_01082019/metadata_subset_01082019.csv",row.names=1, header=T)
map1=map[match(colnames(gen_tab),rownames(map)),]
map1$t11=factor(c(1:15)[factor(map1$t1)])
map1$t11t2=paste(map1$t11,map1$t2,sep="_")

col11=c("red","blue","green","orange","turquoise1","deeppink","black","royalblue","darkgreen","hotpink","darkcyan","goldenrod1","brown","grey","purple")
col1=col11[factor(map1$t11)]
col2=c("red","blue")[factor(map1$t2)]
pro_names=c("11Beijing","21Liaoning","23Heilongjiang","31Shanghai","32Jiangsu","33Zhejiang","37Shandong","41Henan","42Hubei","43Hunan","45Guangxi","52Guizhou","53Yunnan","55Chongqin","61Shaanxi")
pro_names1=c("1Beijing","2Liaoning","3Heilongjiang","4Shanghai","5Jiangsu","6Zhejiang","7Shandong","8Henan","9Hubei","10Hunan","11Guangxi","12Guizhou","13Yunnan","14Chongqin","15Shaanxi")

#Leave obe province out
#random forest categorical metadata

rf2=vector()
pvalue1=vector()
pvalue_within=vector()
pvalue_betw=vector()
m=1
for (i in c(8,77,9,79,80,83)){
  print (m)
  fname=paste("/Users/shansun/Google\ Drive/China/figures/RF/l1o_RF_",colnames(map)[i],".csv",sep="")
  rse_tab=read.csv(file=fname,header=F)
  rse_tab=apply(rse_tab,1, as.character)
  rse_tab=apply(rse_tab,1, as.numeric)
  fname1=paste("/Users/shansun/Google\ Drive/China/figures/RF/l1o_RF_shuffle_",colnames(map)[i],".csv",sep="")
  rse_tab1=read.csv(file=fname1,header=F)
  rse_tab1=apply(rse_tab1,1, as.character)
  rse_tab1=apply(rse_tab1,1, as.numeric)
  pvalue1[m]=t.test(rse_tab[1:15,1:6],rse_tab1[1:15,7])$p.value
  pvalue_within[m]=t.test(rse_tab[1:15,1:6],rse_tab1[1:150,1:6])$p.value
  pvalue_betw[m]=t.test(rse_tab[1:15,7],rse_tab1[1:150,7])$p.value
  rf4=cbind(rep(colnames(map)[i],1155),c(as.numeric(rse_tab[1:15,1:6]),as.numeric(rse_tab[1:15,7]),as.numeric(rse_tab1[1:150,1:6]),as.numeric(rse_tab1[1:150,7])),c(rep("within",90),rep("between",15),rep("within_control",900),rep("between_control",150)))
  rf2=rbind(rf2,rf4)
  m=m+1
}

colnames(rf2)=c("Metadata","tpr","cat")
rf2=data.frame(rf2)
rf2$cat=factor(rf2$cat,levels=levels(factor(rf2$cat))[c(3,4,1,2)])
rf2[,2]=as.numeric(as.character(rf2[,2]))


stat.test=compare_means(tpr~cat, rf2, method = "t.test", group.by = "Metadata")
stat.test1=stat.test
stat.test1$p=p.adjust(stat.test$p,method="fdr")
stat.test1$p.format=format.pval(stat.test1$p,digits=1)
stat.test2=stat.test1[sort(c(seq(1,36,6),seq(2,36,6),seq(6,36,6))),]
fdr_perc=unlist(stat.test1[sort(seq(1,36,6)),5])
length(fdr_perc) #6
length(which(fdr_perc<0.05))#4

myggplot=list()
for (n in c("t2","EDUC","excreta","male","occup","toilet" )){
  stat.test3=stat.test2[stat.test2$Metadata==n,]
  a=rf2[rf2$Metadata==n,]
  a_max=max(a[,2])
  a$title=n
  if (all((stat.test3[,5][1]<0.05)[1:2])){
    myggplot[n][[1]]=ggboxplot(a, x = "cat", y = "tpr", color = "cat", palette = "jco", add = "jitter", ggtheme = theme_bw() ,ylim=c(0,a_max+0.2))+
      font("title", size = 14, color = "black", face = "bold")+
      stat_pvalue_manual(stat.test3, label = "p.format", y.position = c(a_max,a_max+0.05,a_max+0.1)) +facet_grid(. ~ title)+
      theme(strip.text.x = element_text(size = 20, angle = 0),
            legend.position="none",
            axis.text.y   = element_text(size=14),
            axis.text.x   = element_text(size=14),
            axis.title.y  = element_text(size=18),
            axis.title.x  = element_text(size=18),
            panel.background = element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.line = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=0.5))
  }
}

tiff("/Users/shansun/Documents/GitHub/RandomForest/l1o_RF_tpr_cat.tiff", width=24,height=5,units="in",res=200)
par(mfrow=c(2,2),mar=c(5,5,5,5))
grid.arrange(myggplot[[1]],myggplot[[2]],myggplot[[3]],myggplot[[4]],ncol=4)  
dev.off() 

#RF for continuous metadata
mean1=vector()
mean_inc=vector()
rf2_1=vector()
m=1
for (i in c(6,10:12,15:67,70,71,84,85)){
  fname=paste("/Users/shansun/Google\ Drive/China/figures/RF/l1o_RF_",colnames(map)[i],".csv",sep="")
  rse_tab=read.csv(file=fname,header=F)
  rse_tab=apply(rse_tab,1, as.character)
  rse_tab=apply(rse_tab,1, as.numeric)
  fname1=paste("/Users/shansun/Google\ Drive/China/figures/RF/l1o_RF_shuffle_",colnames(map)[i],".csv",sep="")
  rse_tab1=read.csv(file=fname1,header=F)
  rse_tab1=apply(rse_tab1,1, as.character)
  rse_tab1=apply(rse_tab1,1, as.numeric)
  rf4=cbind(rep(colnames(map)[i],1155),c(as.numeric(rse_tab[1:15,1:6]),as.numeric(rse_tab[1:15,7]),as.numeric(rse_tab1[1:150,1:6]),as.numeric(rse_tab1[1:150,7])),c(rep("within",90),rep("between",15),rep("within_control",900),rep("between_control",150)))
  rf2_1=rbind(rf2_1,rf4)
  mean1[m]=mean(rse_tab)
  mean_inc[m]=mean(rse_tab[,7])/mean(rse_tab[,1:6])
  m=m+1
}

mean_increase=mean_inc[which(fdr_perc<0.05)][order(mean1[which(fdr_perc<0.05)])]
names(mean_increase)=colnames(map)[c(6,10:12,15:67,70,71,84,85)][which(fdr_perc<0.05)][order(mean1[which(fdr_perc<0.05)])]

colnames(rf2_1)=c("Metadata","rmse","cat")
rf2_1=data.frame(rf2_1)
rf2_1$cat=factor(rf2_1$cat,levels=levels(factor(rf2_1$cat))[c(3,4,1,2)])
rf2_1[,2]=as.numeric(as.character(rf2_1[,2]))

stat.test=compare_means(rmse~cat, rf2_1, method = "t.test", group.by = "Metadata")
stat.test1=stat.test
stat.test1$p=p.adjust(stat.test$p,method="fdr")
stat.test1$p.format=format.pval(stat.test1$p,digits=1)
stat.test2=stat.test1[sort(c(seq(1,366,6),seq(2,366,6),seq(6,366,6))),]
fdr_perc=unlist(stat.test1[sort(seq(1,366,6)),5])
length(fdr_perc) #61
length(which(fdr_perc<0.05))#27

myggplot=list()
for (n in colnames(map)[c(6,10:12,15:67,70,71,84,85)][order(mean1)]){
  print(n)
  stat.test3=stat.test2[stat.test2$Metadata==n,]
  a=rf2_1[rf2_1$Metadata==n,]
  a_max=max(a[,2])
  a_gra=max(a[,2])/10
  a$title=n
  if (all((stat.test3[,5][1]<0.05)[1])){
    myggplot[n][[1]]=ggboxplot(a, x = "cat", y = "rmse", color = "cat", palette = "jco", add = "jitter", ggtheme = theme_bw() ,ylim=c(0,a_max+a_gra*4))+
      stat_pvalue_manual(stat.test3, label = "p.format", y.position = c(a_max,a_max+a_gra,a_max+a_gra*2)) +facet_grid(. ~ title)+
      theme(strip.text.x = element_text(size = 20, angle = 0),
            legend.position="top",
            axis.text.y   = element_text(size=14),
            axis.text.x   = element_text(size=14),
            axis.title.y  = element_text(size=18),
            axis.title.x  = element_text(size=18),
            panel.background = element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.line = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=0.5))
  }
}


for (i in 0:2){
  fname=paste("/Users/shansun/Documents/GitHub/RandomForest/l1o_RF_rmse_cat_filt_",i,".tiff",sep="")
  tiff(fname, width=24,height=5,units="in",res=200)
  par(mfrow=c(1,4),mar=c(5,5,5,5))
  if(i==2){
    grid.arrange(myggplot[[i*4+1]],myggplot[[i*4+2]],ncol=4)
  }else{
    grid.arrange(myggplot[[i*4+1]],myggplot[[i*4+2]],myggplot[[i*4+3]],myggplot[[i*4+4]],ncol=4)
  }
  dev.off() 
}
