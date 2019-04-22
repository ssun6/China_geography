library(vegan)
setwd("/Users/shansun/Documents/GitHub/China_geography/tax_anova/")
#ANOVA test of all taxa at all levels
gen_tab=read.table(file="otu_table_12500_wtax_L6.txt",sep="\t",header=T,row.names=1)
gen_ab=gen_tab*12500
tab=log10(gen_tab*12500+1)
map=read.csv(file="metadata_subset_01082019.csv",row.names=1, header=T)
map=map[match(colnames(gen_tab),rownames(map)),]
col11=c("red","blue","green","orange","turquoise1","deeppink","black","royalblue","darkgreen","hotpink","darkcyan","goldenrod1","brown","grey","purple")
col1=col11[factor(map$t1)]
col2=c("red","blue")[factor(map$t2)]
pro_names=c("1Beijing","2Liaoning","3Heilongjiang","4Shanghai","5Jiangsu","6Zhejiang","7Shandong","8Henan","9Hubei","10Hunan","11Guangxi","12Guizhou","13Yunnan","14Chongqin","15Shaanxi")
map$t11=c(1:15)[factor(map$t1)]

tab_ano=list()
tab_fdr=list()
for (i in 2:6){
  fname=paste("/Users/shansun/Google\ Drive/China/china1_tax/otu_table_12500_wtax_L",i,".txt",sep="")
  tab1=read.table(file=fname,sep="\t",header=T,row.names=1)
  tab2=log10(tab1*12500+1)
  tab=tab2[which(apply(tab1,1,function(i){length(which(i!=0))})/2164>0.1),]#present in more than 10% samples
  tab_ano1=list()
  for (n in 1:dim(tab)[1]){
    tab_ano1[n]=summary(aov(lm(as.numeric(tab[n,])~factor(map$t1))))[[1]][1,5]
  }
  names(tab_ano1)=rownames(tab)
  tab_ano[[i-1]]=unlist(tab_ano1)
  tab_fdr[[i-1]]=p.adjust(tab_ano[[i-1]],method="fdr")
  print(length(tab_fdr[[i-1]]))
  print(length(tab_fdr[[i-1]][tab_fdr[[i-1]]<0.05]))
  tab3=apply(tab1[names(which(tab_fdr[[i-1]]<0.05)),],1,function(i){aggregate(i ~ map$t1, FUN=mean)})
  tab4=t(do.call(cbind,tab3)[,seq(2,length(tab3)*2,2)])
  tab_names1=sapply(strsplit(rownames(tab4),"\\.i"), "[", 1)
  tab_names=sapply(strsplit(tab_names1,";.__"), "[", i)
  tab5=tab4[!is.na(tab_names)&tab_names!="Other",]
  rownames(tab5)=na.omit(tab_names[tab_names!="Other"])
  colnames(tab5)=c(1:15)
  tab6=tab5[which(rowMeans(tab5)>0.0001),] #relative abundance >0.01%
  tab7=tab6[order(rowMeans(tab6),decreasing = T),]
  tab8=t(scale(t(tab7)))
  print(dim(tab5)[1])
  print(dim(tab6)[1])
  fname1=paste("heatmap_level",i,".png",sep="")
  png(fname1,width=1000,height=dim(tab8)[1]^0.75*60)
  pheatmap(tab8,cluster_rows = T,cluster_cols = T,border_color=NA,cex=1.3,cellheight=15)
  dev.off()
}

