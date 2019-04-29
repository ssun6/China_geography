library(vegan)
setwd("/Users/shansun/Documents/GitHub/China_geography/gen_meta_association/")
setwd("/users/ssun5/china_git/China_geography/PCoA_PERMANOVA/")
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
pro_n=sort(unique(map$t1))
#metadata
#continuous variable #70 c(6,10:68,70:71,84:91)
#categorical variable #8: 2 levels:8, 77
#                     >2 levels:7 (t1),9,78,79,80,83

#genus metadata associations across all provinces
#kendall correlation for continuous metadata
gen_names=vector()
meta_names=vector()
70*1343 #94010
pvalues=vector()
coef=vector()
n=1
for (i in 1:dim(tab)[1]){
  for (j in c(6,10:68,70:71,84:91)){
      tax_sub=tab[i,]
      map_sub=map[,j]
      if (all(is.na(tax_sub))|all(is.na(map_sub))){
        gen_names[n]=rownames(tab)[i]
        meta_names[n]=colnames(map)[j]
        pvalues[n]=NA
        coef[n]=NA
        next
      }
      t=cor.test(as.numeric(tax_sub),as.numeric(map_sub),method = "kendall")
      if (!is.na(t$p.value)){
        gen_names[n]=rownames(tab)[i]
        meta_names[n]=colnames(map)[j]
        pvalues[n]=t$p.value
        coef[n]=t$estimate
      }else{
        gen_names[n]=rownames(tab)[i]
        meta_names[n]=colnames(map)[j]
        pvalues[n]=NA
        coef[n]=NA
      }
    n=n+1
    print(n)
  } 
}
fdr=p.adjust(pvalues,methods="fdr")
gen_meta_kendall=cbind(gen_names,meta_names,coef,pvalues,fdr)
write.csv(gen_meta_kendall,file="China_kendall_all.csv")

#find the taxa overrepresented in associations with metadata
gn_fdr0.05=gen_meta_kendall[gen_meta_kendall$fdr<0.05,]
dim(gn_fdr0.05) #1832
-sort(-table(gn_fdr0.05[,2]))
"
pct_kani         pct_fani        mcarotene            mvite         mretinol           mfiber 
128              103               87               85               75               73 
sani_15               U9           WEIGHT          index15            mvita          mniacin 
71               59               46               40               39               29 
mcu        count_fg2              mmg              mzn      mp_enercarb       hhincgross 
28               26               26               26               25               23 
mvitc       mp_enerpro           HEIGHT       mp_enerfat        mriboflav             mash 
23               22               21               21               21               20 
mprotein              BMI              U10        dbpavg_v2              mca        mgramscon 
19               18               18               17               17               17 
mthiamin           mwater               TG            TG_mg               TP          insulin 
17               17               17               17               17               14 
mfat             CHOl          CHOl_mg              mna      totalPA_est              age 
14               13               13               13               13               12 
mse              mfe           CRP_hs     mderive_kcal        sbpavg_v2            Apo_A 
12               10                9                9                9                8 
mcarbo          menergy               mk              mmn               mp      change_kcal 
8                8                8                8                7                6 
HDL_C           HDL_mg            LDL_C           LDL_mg              h1c       change_fat 
6                6                6                6                5                4 
Apo_B     count_fg2_t1        glu_field           P500_2           P507_2     change_count 
3                3                2                2                2                1 
count_fg2_commid            LP_a_ 
0                0 "

order0.05=data.frame(table(as.character(data.frame(strsplit(as.character(gn_fdr0.05[,1]),";"),stringsAsFactors = F)[4,])))
order0.05$Perc=order0.05$Freq/sum(order0.05$Freq)
order1=data.frame(table(as.character(data.frame(strsplit(as.character(gn_fdr_all[,1]),";"),stringsAsFactors = F)[4,])))
order1$Perc=order1$Freq/sum(order1$Freq)
order_merge=data.frame(as.matrix(merge(order0.05,order1,by="Var1",all=TRUE)),stringsAsFactors = F)

#get the top 10 metadata with most associated taxa
metanames=names(sort(table(gn_fdr0.05[,2]),decreasing=T))[1:10]
for (mn in metanames){
  order_n=data.frame(table(as.character(data.frame(strsplit(as.character(gn_fdr0.05[gn_fdr0.05$meta_names==mn,1]),";"),stringsAsFactors = F)[4,])))
  order_n$Perc=order_n$Freq/sum(order_n$Freq)
  order_merge=data.frame(as.matrix(merge(order_merge,order_n,by=1,all=T)),stringsAsFactors = F)
}
order_merge[is.na(order_merge)]=0

order_merge1=data.frame(order_merge[order_merge$Var1!="o__"&order_merge$Var1!="Other",seq(1,25,2)])
order_merge2=data.frame(order_merge1[order_merge1$Perc.x>0.015,])
order_merge3=matrix(as.numeric(unlist(order_merge2[,-1])),ncol=dim(order_merge2)[2]-1)
order_merge4=as.matrix(rbind(order_merge3,1-colSums(order_merge3)))
rownames(order_merge4)=c(order_merge2[,1],"Other")

order_merge_1=matrix(as.numeric(unlist(order_merge[,seq(2,25,2)])),nrow=dim(order_merge)[1])
order_merge1_1=data.frame(order_merge[order_merge$Var1!="o__"&order_merge$Var1!="Other",seq(2,dim(order_merge)[2],2)])
order_merge2_1=data.frame(order_merge1_1[order_merge1$Perc.x>0.015,])
order_merge3_1=matrix(as.numeric(unlist(order_merge2_1)),ncol=dim(order_merge2_1)[2])
rownames(order_merge3_1)=order_merge2[,1]

#Fisher test for the overrepresented and underrepresented taxa and averaged taxa
fdr_gen_meta1=matrix(nrow=9,ncol=10)
all_meta=colSums(order_merge_1)
for (i in 1:10){
  p_gen_meta=vector()
  for (n in 1:9){
    order_overa=order_merge3_1[n,2]
    p_gen_meta[n]=fisher.test(matrix(c(order_merge3_1[n,i+2],order_overa-order_merge3_1[n,i+2],all_meta[i+2]-order_merge3_1[n,i+2],all_meta[2]-order_overa-(all_meta[i+2]-order_merge3_1[n,i+2])),nrow=2))$p.value
  }
  fdr_gen_meta1[,i]=p.adjust(p_gen_meta,method="fdr")
}
rownames(fdr_gen_meta1)=order_merge2[,1]
colnames(fdr_gen_meta1)=metanames
write.csv(fdr_gen_meta1,file="fisher_fdr_gen_meta_kendall.csv")

#barplots for the taxa associated with metadata at order levels
order_merge5=order_merge4[,c(3:12,1:2)]
write.csv(order_merge5,file="gen_meta_kendall_for_barplots.csv")

order_merge5=read.csv(file="gen_meta_kendall_for_barplots.csv",row.names = 1)
order_merge5=as.matrix(order_merge5[,-11])
col11=c("hotpink","blue","beige","orange","grey","turquoise1","yellow","purple","red","royalblue","green","brown","darkcyan","goldenrod1","lavender","black","darkgreen")

meta_lab=c("animal kcal","animal fat","carotene intake","vitamin E intake","retinol intake","fiber intake", "sanitation index","hip circum","weight","urbanization index")
png("barplot_gen_meta_kendall.png",width=3900,height=800)
par(mfrow=c(1,2),mar=c(5,5,5,5),las=1)
barplot(order_merge5,col=rev(col11[1:dim(order_merge5)[1]]),width=1,names.arg=c(meta_lab,"all"),ylab="Percentage",cex.lab=2.5,cex.axis=1.5, cex.names=1.5)
plot.new()
legend("left",c("Other",rev(as.character(order_merge2[,1]))),col=col11[1:10],pch=16,cex=2)
dev.off()

-sort(-table(gn_fdr0.05[,3]))
"index15       sani_15   totalPA_est        mfiber     mcarotene 
181           163            83            77            76 
mniacin      mretinol     mvite           mcu        mcarbo 
64            63            60            59            56 
mvita    hhincgross        mmg    mp_enercarb            U9 
55            53            53            53            46 
mvitc    mp_enerfat        WEIGHT        mzn        HEIGHT 
45            43            42            38            36 "

#wilcox test for 2 level metadata (t2 and gender) to determine the directions
gen_names=vector()
meta_names=vector()
pvalues=vector()
coef=vector()
conf1=vector()
conf2=vector()
n=1
for (i in 1:dim(tab)[1]){
  for (j in c(8,77)){
    tax_sub1=as.numeric(tab[i,map[,j]==levels(factor(map[,j]))[1]])
    tax_sub2=as.numeric(tab[i,map[,j]==levels(factor(map[,j]))[2]])
    t=wilcox.test(tax_sub1,tax_sub2, conf.int = TRUE,exact=TRUE)
    if (!is.na(t$p.value)){
      gen_names[n]=rownames(tab)[i]
      meta_names[n]=colnames(map)[j]
      pvalues[n]=t$p.value
      coef[n]=t$estimate
      conf1[n]=t$conf.int[1]
      conf2[n]=t$conf.int[2]
    }else{
      gen_names[n]=rownames(tab)[i]
      meta_names[n]=colnames(map)[j]
      pvalues[n]=NA
      coef[n]=NA
      conf1[n]=NA
      conf2[n]=NA
    }
    print(n)
    n=n+1
  } 
}
wilcox_fdr=p.adjust(pvalues,method="fdr")
gen_meta_wilcox=cbind(gen_names,meta_names,coef,conf,pvalues,wilcox_fdr)
write.csv(gen_meta_wilcox,file="China_wilcox_all.csv")

#find the order level difference
gen_meta_wilcox=read.csv(file="China_wilcox_all.csv",row.names=1)
gn_fdr0.05=gen_meta_wilcox[gen_meta_wilcox$fdr<0.05,] 
order0.05=data.frame(table(as.character(data.frame(strsplit(as.character(gn_fdr0.05[,1]),";"),stringsAsFactors = F)[4,])))
order0.05$Perc=order0.05$Freq/sum(order0.05$Freq)
order1=data.frame(table(as.character(data.frame(strsplit(as.character(gn_fdr_all[,1]),";"),stringsAsFactors = F)[4,])))
order1$Perc=order1$Freq/sum(order1$Freq)
order_merge=data.frame(as.matrix(merge(order0.05,order1,by="Var1",all=TRUE)),stringsAsFactors = F)

metanames=names(sort(table(gn_fdr_sig_wil[,2]),decreasing=T))
for (mn in metanames){
  order_n=data.frame(table(as.character(data.frame(strsplit(as.character(gn_fdr_sig_wil[gn_fdr_sig_wil$meta_names==mn,1]),";"),stringsAsFactors = F)[4,])))
  order_n$Perc=order_n$Freq/sum(order_n$Freq)
  order_merge=data.frame(as.matrix(merge(order_merge,order_n,by="Var1",all=T)),stringsAsFactors = F)
}
order_merge[is.na(order_merge)]=0

order_merge1=data.frame(order_merge[order_merge$Var1!="o__"&order_merge$Var1!="Other",seq(1,dim(order_merge)[2],2)])
order_merge2=data.frame(order_merge1[order_merge1$Perc.x>0.015,])
order_merge3=matrix(as.numeric(unlist(order_merge2[,-1])),ncol=length(metanames)+2)
order_merge4_wil=as.matrix(rbind(order_merge3,1-colSums(order_merge3)))
rownames(order_merge4_wil)=c(order_merge2[,1],"Other")

order_merge_1=matrix(as.numeric(unlist(order_merge[,seq(2,dim(order_merge)[2],2)])),nrow=dim(order_merge)[1])
order_merge1_1=data.frame(order_merge[order_merge$Var1!="o__"&order_merge$Var1!="Other",seq(2,dim(order_merge)[2],2)])
order_merge2_1=data.frame(order_merge1_1[order_merge1$Perc.x>0.015,])
order_merge3_1=matrix(as.numeric(unlist(order_merge2_1)),ncol=length(metanames)+2)
rownames(order_merge3_1)=order_merge2[,1]

fdr_gen_meta1=matrix(nrow=dim(order_merge3_1)[1],ncol=length(metanames))
all_meta=colSums(order_merge_1)
for (i in 1:length(metanames)){
  p_gen_meta=vector()
  for (n in 1:dim(order_merge3_1)[1]){
    order_overa=order_merge3_1[n,2]
    p_gen_meta[n]=fisher.test(matrix(c(order_merge3_1[n,i+2],order_overa-order_merge3_1[n,i+2],all_meta[i+2]-order_merge3_1[n,i+2],all_meta[2]-order_overa-(all_meta[i+2]-order_merge3_1[n,i+2])),nrow=2))$p.value
  }
  fdr_gen_meta1[,i]=p.adjust(p_gen_meta,method="fdr")
}
rownames(fdr_gen_meta1)=order_merge2[,1]
fdr_gen_meta_wil=fdr_gen_meta1
colnames(fdr_gen_meta_wil)=c("t2","gender")

#kw for more than 2 level categorical factors
gen_names_1=vector()
meta_names_1=vector()
pvalues_1=vector()
coef_1=vector()
n=1
for (i in 1:dim(tab)[1]){
  for (j in c(7,9,78:80,83)){
    tax_sub=tab[i,]
    map_sub=map[,j]
    if (all(is.na(tax_sub))|all(is.na(map_sub))){
      gen_names_1[n]=rownames(tab)[i]
      meta_names_1[n]=colnames(map)[j]
      pvalues_1[n]=NA
      coef_1[n]=NA
      next
    }
    if (length(levels(factor(map_sub)))==1){
      gen_names_1[n]=rownames(tab)[i]
      meta_names_1[n]=colnames(map)[j]
      pvalues_1[n]=NA
      coef_1[n]=NA
      next
    }
    t=kruskal.test(as.numeric(tax_sub)~factor(map_sub))
    if (!is.na(t$p.value)){
      gen_names_1[n]=rownames(tab)[i]
      meta_names_1[n]=colnames(map)[j]
      pvalues_1[n]=t$p.value
      coef_1[n]=t$statistic
    }else{
      gen_names_1[n]=rownames(tab)[i]
      meta_names_1[n]=colnames(map)[j]
      pvalues_1[n]=NA
      coef_1[n]=NA
    }
    print(n)
    n=n+1
  } 
}
gn_fdr_1=data.frame(as.matrix(cbind(gen_names_1,meta_names_1,coef_1,pvalues_1)),stringsAsFactors=FALSE)
gn_fdr_1$fdr=p.adjust(gn_fdr_1$pvalues_1,method="fdr")
gn_fdr_sig=na.omit(gn_fdr_1[gn_fdr_1$fdr<0.05,])
dim(gn_fdr_sig)#240
gn_fdr_sig=gn_fdr_sig[order(gn_fdr_sig$fdr),]
table(gn_fdr_sig[,2])
"drk_water      EDUC   excreta     occup    toilet 
68        24        18        19       111 "
write.csv(gn_fdr_1,file="genus_metadata_kw.csv")

gn_fdr_1=read.csv(file="genus_metadata_kw.csv",row.names=1)
gn_fdr_sig=na.omit(gn_fdr_1[gn_fdr_1$fdr<0.05,])

order0.05=data.frame(table(as.character(data.frame(strsplit(as.character(gn_fdr_sig[,1]),";"),stringsAsFactors = F)[4,])))
order0.05$Perc=order0.05$Freq/sum(order0.05$Freq)
order1=data.frame(table(as.character(data.frame(strsplit(as.character(gn_fdr_1[,1]),";"),stringsAsFactors = F)[4,])))
order1$Perc=order1$Freq/sum(order1$Freq)
order_merge=data.frame(as.matrix(merge(order0.05,order1,by="Var1",all=TRUE)),stringsAsFactors = F)

metanames=names(sort(table(gn_fdr_sig[,2]),decreasing=T))
for (mn in metanames){
  order_n=data.frame(table(as.character(data.frame(strsplit(as.character(gn_fdr_sig[gn_fdr_sig$meta_names==mn,1]),";"),stringsAsFactors = F)[4,])))
  order_n$Perc=order_n$Freq/sum(order_n$Freq)
  order_merge=data.frame(as.matrix(merge(order_merge,order_n,by="Var1",all=T)),stringsAsFactors = F)
}
order_merge[is.na(order_merge)]=0

order_merge1=data.frame(order_merge[order_merge$Var1!="o__"&order_merge$Var1!="Other",seq(1,dim(order_merge)[2],2)])
order_merge2=data.frame(order_merge1[order_merge1$Perc.x>0.015,])
order_merge3=matrix(as.numeric(unlist(order_merge2[,-1])),ncol=length(metanames)+2)
order_merge4_kw=as.matrix(rbind(order_merge3,1-colSums(order_merge3)))
rownames(order_merge4_kw)=c(order_merge2[,1],"Other")

order_merge_1=matrix(as.numeric(unlist(order_merge[,seq(2,dim(order_merge)[2],2)])),nrow=dim(order_merge)[1])
order_merge1_1=data.frame(order_merge[order_merge$Var1!="o__"&order_merge$Var1!="Other",seq(2,dim(order_merge)[2],2)])
order_merge2_1=data.frame(order_merge1_1[order_merge1$Perc.x>0.015,])
order_merge3_1=matrix(as.numeric(unlist(order_merge2_1)),ncol=length(metanames)+2)
rownames(order_merge3_1)=order_merge2[,1]

fdr_gen_meta1=matrix(nrow=dim(order_merge3_1)[1],ncol=length(metanames))
all_meta=colSums(order_merge_1)
for (i in 1:length(metanames)){
  p_gen_meta=vector()
  for (n in 1:dim(order_merge3_1)[1]){
    order_overa=order_merge3_1[n,2]
    p_gen_meta[n]=fisher.test(matrix(c(order_merge3_1[n,i+2],order_overa-order_merge3_1[n,i+2],all_meta[i+2]-order_merge3_1[n,i+2],all_meta[2]-order_overa-(all_meta[i+2]-order_merge3_1[n,i+2])),nrow=2))$p.value
  }
  fdr_gen_meta1[,i]=p.adjust(p_gen_meta,method="fdr")
}
rownames(fdr_gen_meta1)=order_merge2[,1]
colnames(fdr_gen_meta1)=metanames
fdr_gen_meta_kw=fdr_gen_meta1

fdr_merge4_both=merge(fdr_gen_meta_wil,fdr_gen_meta_kw,by="row.names",all=T)
write.csv(fdr_merge4_both,file="figuresfisher_fdr_metadata_categorical.csv")

col11=c("hotpink","blue","darkgreen","orange","grey","red","yellow","purple","turquoise1","royalblue","green","brown","darkcyan","goldenrod1","lavender","beige","black")

metanames_lab1=c("urbanization","gender")
png("barplot_order_categorical_meta_wilcoxon.png",width=1500,height=800)
par(mfrow=c(1,2),mar=c(5,5,5,5),las=1)
barplot(order_merge4_wil[,c(3,4,2)],width=1,col=rev(col11[1:dim(order_merge4_wil)[1]]),names.arg=c(metanames_lab1,"all"),ylab="Percentage",cex.lab=2.5,cex.axis=1.5, cex.names=2)
plot.new()
legend("left",rev(rownames(order_merge4_wil)),col=col11[1:dim(order_merge4_both)[1]],pch=16,cex=2)
dev.off()

metanames_lab2=c("toilet","drinking water","education","occupation","excreta")
png("barplot_order_categorical_meta_kw.png",width=3000,height=800)
par(mfrow=c(1,2),mar=c(5,5,5,5),las=1)
barplot(order_merge4_kw[,c(3:7,2)],width=1,col=rev(col11[1:dim(order_merge4_kw)[1]]),names.arg=c(metanames_lab2,"all"),ylab="Percentage",cex.lab=2.5,cex.axis=1.5, cex.names=2)
plot.new()
legend("left",rev(rownames(order_merge4_kw)),col=col11[1:dim(order_merge4_kw)[1]],pch=16,cex=2)
dev.off()