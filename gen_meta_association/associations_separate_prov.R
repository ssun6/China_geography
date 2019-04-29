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

#metadata and genera correlation in separate provinces
#kendall correlation for continuous metadata

1343*77 #103411

gen_names=vector()
meta_names=vector()
pvalues=matrix(nrow=103411,ncol=15)
coef=matrix(nrow=103411,ncol=15 )
n=1
for (i in 1:dim(tab)[1]){
  for (j in c(6,10:68,70:71,84:91)){
    for (m in 1:15){
      tax_sub=tab[i,map$t1==pro_n[m]]
      map_sub=map[map$t1==pro_n[m],j]
      if (all(is.na(tax_sub))|all(is.na(map_sub))){
        gen_names[n]=rownames(tab)[i]
        meta_names[n]=colnames(map)[j]
        pvalues[n,m]=NA
        coef[n,m]=NA
        next
      }
      t=cor.test(as.numeric(tax_sub),as.numeric(map_sub),method = "kendall")
      if (!is.na(t$p.value)){
        gen_names[n]=rownames(tab)[i]
        meta_names[n]=colnames(map)[j]
        pvalues[n,m]=t$p.value
        coef[n,m]=t$estimate
      }else{
        gen_names[n]=rownames(tab)[i]
        meta_names[n]=colnames(map)[j]
        pvalues[n,m]=NA
        coef[n,m]=NA
      }
    }
    n=n+1
    print(n)
  } 
}
for (i in 1:dim(tab)[1]){
  for (j in c(9,78:80,83)){
    for (m in 1:15){
      tax_sub1=as.numeric(tab[i,map[,j]==levels(factor(map[,j]))[1]])
      tax_sub2=as.numeric(tab[i,map[,j]==levels(factor(map[,j]))[2]])
      t=wilcox.test(tax_sub1,tax_sub2, conf.int = TRUE,exact=TRUE)
      if (!is.na(t$p.value)){
        gen_names[n]=rownames(tab)[i]
        meta_names[n]=colnames(map)[j]
        pvalues[n,m]=t$p.value
        coef[n,m]=t$estimate
      }else{
        gen_names[n]=rownames(tab)[i]
        meta_names[n]=colnames(map)[j]
        pvalues[n,m]=NA
        coef[n,m]=NA
      }
    }
    print(n)
    n=n+1
  } 
}
for (i in 1:dim(tab)[1]){
  for (j in c(8,77))){
    for (m in 1:15){
      tax_sub1=as.numeric(tab[i,map[,j]==levels(factor(map[,j]))[1]])
      tax_sub2=as.numeric(tab[i,map[,j]==levels(factor(map[,j]))[2]])
      t=wilcox.test(tax_sub1,tax_sub2, conf.int = TRUE,exact=TRUE)
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
    }
    n=n+1
    print(n)
  } 
}
p_adjs=matrix(p.adjust(pvalues,method="fdr"),nrow=103411)
rownames(p_adjs)=paste(gen_names,meta_names,sep="~")
rownames(coef)=paste(gen_names,meta_names,sep="~")
write.csv(p_adjs,file="China_association_prov_fdr.csv")
write.csv(coef,file="China_association_prov_estimate.csv")

#permutation test to see whether the number of significant ones happen by chance or not.
library(picante)
#for kendall test
1343*70 #94010
p_adjs1=p_adjs[1:94010,]
num_sig=apply(p_adjs1,1,function(i){length(which(i<0.05))})
kend_tab_ken=data.frame(table(num_sig))
colnames(kend_tab_ken)[1]="num_sig_1"
p_adjs_1=randomizeMatrix(p_adjs1,null.model = "frequency",iterations = 100000000)
for (j in 1:1000){
  p_adjs_1=randomizeMatrix(p_adjs_1,null.model = "frequency",iterations = 1000000000)
  num_sig_1=apply(p_adjs_1,1,function(i){length(which(i<0.05))})
  kend_tab1=data.frame(table(num_sig_1))
  kend_tab_ken=merge(kend_tab_ken, kend_tab1,by="num_sig_1",all=T)
  print(j)
}
pvals_ken=apply(kend_tab_ken[,3:1002],1,function(y) sum(length(which(!is.na(y)))))/1000
names(pvals_ken)=kend_tab_ken[,1]
print(pvals_ken)

#for kruskal test
1343*5 #6715
94010+1343*5 #100725
p_adjs2=p_adjs[94011:100725,]
num_sig=apply(p_adjs2,1,function(i){length(which(i<0.05))})
kend_tab_kw=data.frame(table(num_sig))
colnames(kend_tab_kw)[1]="num_sig_1"
p_adjs_2=randomizeMatrix(p_adjs2,null.model = "frequency",iterations = 100000000)
for (j in 1:1000){
  p_adjs_2=randomizeMatrix(p_adjs_2,null.model = "frequency",iterations = 1000000000)
  num_sig_1=apply(p_adjs_2,1,function(i){length(which(i<0.05))})
  kend_tab1=data.frame(table(num_sig_1))
  kend_tab_kw=merge(kend_tab_kw, kend_tab1,by="num_sig_1",all=T)
  print(j)
}

pvals_kw=apply(kend_tab_kw[,3:1002],1,function(y) sum(length(which(!is.na(y)))))/1000
names(pvals_kw)=kend_tab_kw[,1]
print(pvals_kw)

#for wilcoxin test considering the direction of the test
1343*5 #6715
94010+1343*5 #100725
p_adjs2=p_adjs[94011:100725,]
num_sig=apply(p_adjs2,1,function(i){length(which(i<0.05))})
kend_tab_wil=data.frame(table(num_sig))
colnames(kend_tab_wil)[1]="num_sig_1"
p_adjs_2=randomizeMatrix(p_adjs2,null.model = "frequency",iterations = 100000000)
for (j in 1:1000){
  p_adjs_2=randomizeMatrix(p_adjs_2,null.model = "frequency",iterations = 1000000000)
  num_sig_1=apply(p_adjs_2,1,function(i){length(which(i<0.05))})
  kend_tab1=data.frame(table(num_sig_1))
  kend_tab_wil=merge(kend_tab_wil, kend_tab1,by="num_sig_1",all=T)
  print(j)
}
pvals_wil=apply(kend_tab_wil[,3:1002],1,function(y) sum(length(which(!is.na(y)))))/1000
names(pvals_wil)=kend_tab_wil[,1]
print(pvals_wil)




