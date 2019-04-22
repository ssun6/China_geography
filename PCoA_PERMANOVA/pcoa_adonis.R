library(vegan)
#setwd("/Users/shansun/Documents/GitHub/China_geography/PCoA_PERMANOVA/")
setwd("/users/ssun5/china_git/China_geography/PCoA_PERMANOVA/")
#PCoA analysis based on genus abundance
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

gen_pcoa=capscale(t(tab)~1,distance="bray")
(gen_pcoa$CA$eig/sum(gen_pcoa$CA$eig))[1:6]

png("pcoa12_t1.png",width=900,height=750)
par(mar=c(5,5,5,5),mfrow=c(1,1))
pcoa12=ordiplot(gen_pcoa,choices=c(1,2),type="none",cex.lab=2,xlab="PCoA1 (14.2%)",ylab="PCoA2 (5.5%)")
points(pcoa12,"sites",col=adjustcolor(col1, alpha.f = 0.05),pch=16,cex=2.5)
ordiellipse(pcoa12, map$t11, kind="se", conf=0.95, lwd=4, draw = "lines", col=col11[1],show.groups=1,label=T,font=2,cex=1) #beijing
ordiellipse(pcoa12, map$t11, kind="se", conf=0.95, lwd=4, draw = "lines", col=col11[2],show.groups=2,label=T,font=2,cex=1) #liaoning
ordiellipse(pcoa12, map$t11, kind="se", conf=0.95, lwd=4, draw = "lines", col=col11[3],show.groups=3,label=T,font=2,cex=1) #heilongjiang
ordiellipse(pcoa12, map$t11, kind="se", conf=0.95, lwd=4, draw = "lines", col=col11[4],show.groups=4,label=T,font=2,cex=1) #shanghai
ordiellipse(pcoa12, map$t11, kind="se", conf=0.95, lwd=4, draw = "lines", col=col11[5],show.groups=5,label=T,font=2,cex=1) #jiangsu
ordiellipse(pcoa12, map$t11, kind="se", conf=0.95, lwd=4, draw = "lines", col=col11[6],show.groups=6,label=T,font=2,cex=1) #zhejiang
ordiellipse(pcoa12, map$t11, kind="se", conf=0.95, lwd=4, draw = "lines", col=col11[7],show.groups=7,label=T,font=2,cex=1) #shandong
ordiellipse(pcoa12, map$t11, kind="se", conf=0.95, lwd=4, draw = "lines", col=col11[8],show.groups=8,label=T,font=2,cex=1) #henan
ordiellipse(pcoa12, map$t11, kind="se", conf=0.95, lwd=4, draw = "lines", col=col11[9],show.groups=9,label=T,font=2,cex=1) #hubei
ordiellipse(pcoa12, map$t11, kind="se", conf=0.95, lwd=4, draw = "lines", col=col11[10],show.groups=10,label=T,font=2,cex=1) #hunan
ordiellipse(pcoa12, map$t11, kind="se", conf=0.95, lwd=4, draw = "lines", col=col11[11],show.groups=11,label=T,font=2,cex=1) #guangxi
ordiellipse(pcoa12, map$t11, kind="se", conf=0.95, lwd=4, draw = "lines", col=col11[12],show.groups=12,label=T,font=2,cex=1) #guizhou
ordiellipse(pcoa12, map$t11, kind="se", conf=0.95, lwd=4, draw = "lines", col=col11[13],show.groups=13,label=T,font=2,cex=1) #yunnan
ordiellipse(pcoa12, map$t11, kind="se", conf=0.95, lwd=4, draw = "lines", col=col11[14],show.groups=14,label=T,font=2,cex=1) #chongqin
ordiellipse(pcoa12, map$t11, kind="se", conf=0.95, lwd=4, draw = "lines", col=col11[15],show.groups=15,label=T,font=2,cex=1) # Shaanxi
legend("topleft",pro_names,col=col11,cex=1.5,pch=16,bty = "n")
dev.off()

pro_n=sort(unique(map$t1))

ado_R2_China=matrix(nrow=78,ncol=2)
n=1
for (i in c(6:68,70,71,77:80,83:91)){
  tab1=na.omit(data.frame(cbind(t(tab),map1[,i])))#remove na rows
  dim(tab1) #2164 1344
  if(length(table(tab1[,1344]))>20){
    fit=adonis(tab1[,1:1343] ~ tab1[,1344], data=tab1, permutations=999)
  }else{
    fit=adonis(tab1[,1:1343] ~ factor(tab1[,1344]), data=tab1, permutations=999)
  }
  ado_R2_China[n,1]=fit$aov.tab[1,5]
  ado_R2_China[n,2]=fit$aov.tab[1,6]
  n=n+1
}
rownames(ado_R2_China)=colnames(map)[c(6:68,70,71,77:80,83:91)]
colnames(ado_R2_China)=c("R2","P")
ado_R2_China=data.frame(ado_R2_China)
ado_R2_China$FDR=p.adjust(ado_R2_China$P,method="fdr")
write.csv(ado_R2_China,file="ado_R2_China.csv")

ado_R2_China1=ado_R2_China[order(ado_R2_China[,1],decreasing=T),]
ado_R2_China2=ado_R2_China1[1:25,]
bar_lab=c("Province","Occupation","Toilet" ,"Urban/rural","Animal source kcal","Animal source fat","Drink water source","Carotene intake","Vitamin E intake", "Sanitation index","Education level","Retinol intake","Urbanization index","Weight","Animal source calories change","Niacin intake","Food diversity","%energy from carbs","Fiber intake","Height","Vitamin A intake","Hip circumference","Zinc intake","Copper intake","Total protein intake")
c("%energy from fat", "Total food intake","DBP","%energy from protein","Water intake","Vitamin C intake","Calcium intake","Protein intake","Magnesium","BMI","Age","Waist circumference","Carbs intake")
col_bar=c("pink","yellow","lightgreen","lightblue")[factor(rev(c(1,1,2,2,3,3,2,3,3,2,1,3,2,4,3,3,3,3,3,4,3,4,rep(3,3))))]
png("pcoa12_t1.png",width=900,height=750)
par(mar=c(5,15,5,5),mfrow=c(1,2))
barplot(rev(as.numeric(ado_R2_China2[,1])),col=col_bar,xlim=c(0,0.2),names.arg=rev(bar_lab),las=2,main="Variance explained",horiz=T)
plot.new()
legend("topleft",c("Demography","Lifestyle","Diet","Physiology"),col="black",pt.bg=c("pink","yellow","lightgreen","lightblue"),cex=1.5,pch=22,bty = "n")
dev.off()

#interactions between province (t1) and other factors
ado_R2_inter=matrix(nrow=77,ncol=6)
for (i in c(6,8:68,70,71,77:80,83:91)){
  print(i)
  print(colnames(map1)[i])
  print(length(table(map1[,i])))
  #remove na rows
  if(any(is.na(map1[,i]))){
    tab_m=t(tab)[-which(is.na(map1[,i])),]
    map_m=na.omit(map1[,i])
  }else{
    tab_m=t(tab)
    map_m=map1[,i]
  }
  if(length(table(map1[,i]))>20){
    fit=adonis(tab_m ~ as.numeric(as.character(map1[,i]))*factor(map1$t1), permutations=999)
  }else{
    fit=adonis(tab_m ~ factor(map1[,i])*factor(map1$t1), permutations=999)
  }
  ado_R2_inter[i-5,1]=fit$aov.tab[1,5]
  ado_R2_inter[i-5,2]=fit$aov.tab[1,6]
  ado_R2_inter[i-5,3]=fit$aov.tab[2,5]
  ado_R2_inter[i-5,4]=fit$aov.tab[2,6]
  ado_R2_inter[i-5,5]=fit$aov.tab[3,5]
  ado_R2_inter[i-5,6]=fit$aov.tab[3,6]
}
rownames(ado_R2_inter)=colnames(map1)[c(6,8:68,70,71,77:80,83:91)]
colnames(ado_R2_inter)=c("meta_R2","meta_p","t1_R2","t1_p","inter_R2","inter_p")
write.csv(ado_R2_inter,file="china1_t1_interaction_effect_size.csv")

#PERMANOVA test in each province
ado_R2_t1_all=matrix(nrow=77,ncol=30)
for (j in 1:15){
  print(j)
  tab_m1=t(tab)[map$t1==pro_n[j],]
  tab_m=tab_m1[,-which(colSums(tab_m1)==0)]
  map_m=map[map$t1==pro_n[j],]
  m=1
  for (i in c(6,8:68,70,71,77:80,83:91)){
    print(i)
    #remove na rows
    if(any(is.na(map_m[,i]))){
      tab_n=tab_m[-which(is.na(map_m[,i])),]
      map_n=na.omit(map_m[,i])
    }else{
      tab_n=tab_m
      map_n=map_m[,i]
    }
    if (length(table(map_n))<2){
      next
    }
    if(length(table(map_n))>20){
      fit=adonis(tab_n ~ as.numeric(as.character(map_n)), permutations=999)
    }else{
      fit=adonis(tab_n ~ factor(map_n), permutations=999)
    }
    ado_R2_t1_all[m,j*2-1]=fit$aov.tab[1,5]
    ado_R2_t1_all[m,j*2]=fit$aov.tab[1,6]
    m=m+1
  }
}
rownames(ado_R2_t1_all)=colnames(map)[c(6,8:12,15:71,77:91)]
colnames(ado_R2_t1_all)=paste(sort(rep(c(1:15),2)),c("effect","P"),sep="_")
write.csv(ado_R2_t1_all,file="china1_t1_effect_size_each_prov.csv")

#American gut data
amg_gen=read.table(file="amg_db_L6.txt",sep="\t",row.names=1,header=T)
map1=read.csv(file="us_meta_filt.csv",row.names=1)
map1=map1[map1$body_habitat=="UBERON:feces",]
sample_inter=intersect(colnames(amg_gen),rownames(map1))
map2=map1[sample_inter,]
map2$state=as.character(map2$state)

amg_gen1=t(t(amg_gen)/colSums(amg_gen)*mean(colSums(amg_gen)))
amg_gen2=amg_gen1[,sample_inter]
amg_gen3=amg_gen2[apply(amg_gen2,1,function(i){sum(i!=0)})!=0,]
amg_gen4=log10(amg_gen3+1)

gen1_pcoa=capscale(t(amg_gen4)~1,distance="bray")
percentVariance<- eigenvals(gen1_pcoa)/sum(eigenvals(gen1_pcoa))
percentVariance[1:3]

map2[map2=="Not provided"]=NA
map2$state[map2$state%in%names(table(map2$state))[table(map2$state)<20]]="Other"
map2$state=factor(map2$state)
col1=col11[map2$state]
png("amg_pcoa12_state.png",width=700,height=500)
par(mar=c(5,5,5,5))
pcoa12=ordiplot(gen1_pcoa,choices=c(1,2),type="none",cex.lab=2,xlab="PCoA1 (8.61%)",ylab="PCoA2 (5.12%)")
points(pcoa12,"sites",col=adjustcolor(col1, alpha.f = 0.1),pch=16,cex=1.5)
for (n in c(1:11,13:16)){
  ordiellipse(pcoa12, map2$state, kind="se", conf=0.95, lwd=4, draw = "lines", col=col11[n],show.groups=levels(factor(map2$state))[n],label=T,font=2,cex=1)
}
legend("topleft",levels(factor(map2$state))[c(1:11,13:16)],col=col11,cex=1.5,pch=16,bty = "n")
dev.off()

meta_list=match(c("age_cat","alcohol_frequency","antibiotic_history","bmi_cat","census_region","consume_animal_products_abx","country_of_birth","deodorant_use","diet_type","dog","drinking_water_source","exercise_frequency","flossing_frequency","flu_vaccine_date","frozen_dessert_frequency",
                  "fruit_frequency","height_cm","high_fat_red_meat_frequency","homecooked_meals_frequency","level_of_education","meat_eggs_frequency","milk_cheese_frequency","multivitamin","nail_biter","race","sex","sugary_sweets_frequency","types_of_plants","vegetable_frequency"),colnames(map2))
meta_list1=match(c("alcohol_frequency","antibiotic_history","census_region","consume_animal_products_abx","country_of_birth","deodorant_use","diet_type","dog","drinking_water_source","exercise_frequency","flossing_frequency","flu_vaccine_date","frozen_dessert_frequency",
                   "fruit_frequency","high_fat_red_meat_frequency","homecooked_meals_frequency","level_of_education","meat_eggs_frequency","milk_cheese_frequency","multivitamin","nail_biter","race","sex","sugary_sweets_frequency","types_of_plants","vegetable_frequency","state","age_corrected","bmi_corrected","height_cm"),colnames(map2))
map2$state[map2$state%in%names(table(map2$state))[table(map2$state)<5]]="Other"
map2$state=factor(map2$state)

fit_m=matrix(nrow=30,ncol=2)
m=1
for (n in meta_list1){
  print(colnames(map2)[n])
  print(length(levels(factor(map2[,n]))))
  if(any(is.na(map2[,n]))){
    tab_m=t(amg_gen4)[-which(is.na(map2[,n])),]
    map_m=na.omit(map2[,n])
  }else{
    tab_m=t(amg_gen4)
    map_m=map2[,n]
  }
  if (length(levels(factor(map2[,n])))>70){
    fit1=adonis(tab_m ~ as.numeric(as.character(map_m)),permutations = 999)
  }else{
    fit1=adonis(tab_m ~ factor(map_m),permutations = 999)
  }
  fit_m[m,1]=fit1$aov.tab[1,5]
  fit_m[m,2]=fit1$aov.tab[1,6]
  print(m)
  m=m+1
}

rownames(fit_m)=c("alcohol_frequency","antibiotic_history","census_region","consume_animal_products_abx","country_of_birth","deodorant_use","diet_type","dog","drinking_water_source","exercise_frequency","flossing_frequency","flu_vaccine_date","frozen_dessert_frequency",
                  "fruit_frequency","high_fat_red_meat_frequency","homecooked_meals_frequency","level_of_education","meat_eggs_frequency","milk_cheese_frequency","multivitamin","nail_biter","race","sex","sugary_sweets_frequency","types_of_plants","vegetable_frequency","state","age_corrected","bmi_corrected","height_cm")
fit_m1=fit_m[fit_m[,2]<0.05,]
fit_m1=fit_m1[order(fit_m1[,1],decreasing = T),]

bar_lab=c("State","Race","Antibiotic history","Age","Vegetable frequency","Education","Types of plants","Alcohol frequency", "Fruit frequency","Exercise frequency","High fat red meat frequency","Homecooked meals frequency","Diet type","Sex","Sugary sweets frequency")
col_bar=c("pink","yellow","lightgreen","lightblue")[factor(rev(c(1,1,4,1,3,1,3,2,3,2,3,3,3,1,3)))]
tiff("Amg_var_adonis.tiff",res=300, width =6, height = 4,units="in")
par(mar=c(5,15,5,5),mfrow=c(1,1))
barplot(rev(as.numeric(fit_m1[c(1:14,15),1])),col=col_bar,xlim=c(0,0.2),names.arg=rev(bar_lab),las=2,main="AMG variance explained",horiz=T,cex.axis=0.8,cex.names=0.8,cex.main=1)
dev.off()




