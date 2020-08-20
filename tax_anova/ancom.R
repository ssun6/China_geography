library(exactRankTests)
library(nlme)
library(dplyr)
library(ggplot2)
library(compositions)
library(readr)
library(tidyverse)

setwd("/Users/shansun/ANCOM/")
source("scripts/ancom_v2.1.R")

t_run=vector()
for (i in 2:6){
  gen_tab=read.table(file=paste0("/Users/shansun/Google\ Drive/China/revision/tax_anova2/otu_table_nc_wtax_L",i,".txt"),sep="\t",header=T,row.names=1)
  map=read.csv(file="/Users/shansun/Google\ Drive/China/Metadata_01082019/metadata_subset_01082019.csv",row.names=1, header=T)
  map=map[match(colnames(gen_tab),rownames(map)),]
  map$SampleID=rownames(map)
  pro_names=c("1Beijing","2Liaoning","3Heilongjiang","4Shanghai","5Jiangsu","6Zhejiang","7Shandong","8Henan","9Hubei","10Hunan","11Guangxi","12Guizhou","13Yunnan","14Chongqin","15Shaanxi")
  map$t11=factor(c(1:15)[factor(map$t1)])
  
  # Step 1: Data preprocessing

  feature_table = gen_tab; meta_data=map;sample_var = "SampleID"; group_var = NULL
  out_cut = 0.05; zero_cut = 0.90; lib_cut = 1000; neg_lb = FALSE
  prepro = feature_table_pre_process(feature_table, meta_data, sample_var, group_var,
                                     out_cut, zero_cut, lib_cut, neg_lb)
  feature_table = prepro$feature_table # Preprocessed feature table
  meta_data = prepro$meta_data # Preprocessed metadata
  struc_zero = prepro$structure_zeros # Structural zero info

  # Step 2: ANCOM

  main_var = "t11"; p_adj_method = "BH"; alpha = 0.05
  adj_formula = NULL; rand_formula = NULL
  t_start = Sys.time()
  res = ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method,
              alpha, adj_formula, rand_formula)
  t_end = Sys.time()
  t_run[i-1] = t_end - t_start # around 30s

  write_csv(res$out, paste0("/Users/shansun/Google\ Drive/China/revision/tax_anova2/ancom_l",i,".csv"))
  
  tab_fdr=read.csv(paste0("/Users/shansun/Google\ Drive/China/revision/tax_anova2/ancom_l",i,".csv"),row.names=1,header=T)
  print(length(which(!tab_fdr[,2])))
 
  tab3=sapply(by(t(gen_tab[rownames(tab_fdr)[which(tab_fdr[,2])],]),map$t11, colMeans),identity)
  tab_names1=sapply(strsplit(rownames(tab3),"\\.i"), "[", 1)
  tab_names=sapply(strsplit(tab_names1,";.__"), "[", i)
  tab5=tab3[!is.na(tab_names)&tab_names!="Other",]
  rownames(tab5)=na.omit(tab_names[tab_names!="Other"])
  colnames(tab5)=c(1:15)
  tab6=tab5[which(rowMeans(tab5)/sum(rowMeans(tab5))>0.0002),] #relative abundance >0.01%
  tab7=tab6[order(rowMeans(tab6),decreasing = T),]
  tab8=t(scale(t(tab7)))
  fname1=paste("/Users/shansun/Google\ Drive/China/revision/tax_anova2/ancom_heatmap_level",i,".pdf",sep="")
  pdf(fname1,width=10,height=dim(tab8)[1]^0.75*0.8)
  pheatmap(tab8,cluster_rows = T,cluster_cols = T,border_color=NA,cex=1.3,cellheight=15)
  dev.off()
}

library(grid)
draw_colnames_45 <- function (coln, ...) {
  m = length(coln)
  x = (1:m)/m - 1/2/m
  grid.text(coln, x = x, y = unit(0.96, "npc"), vjust = 1.5, 
            hjust = 0.5, rot = 0, gp = gpar(...)) ## Was 'hjust=0' and 'rot=270'
}

## 'Overwrite' default draw_colnames with your own version 
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))

meta_data = meta_data %>% rownames_to_column("sampleID")
sample_id = "sampleID"; formula = NULL
p_adj_method = "BH"; zero_cut = 0.90; lib_cut = 1000; group = "t11"
struc_zero = TRUE; neg_lb = TRUE; tol = 1e-05; max_iter = 100
conserve = TRUE; alpha = 0.05; global = TRUE

out = ancombc(feature_table, meta_data, sample_id, formula,
              p_adj_method, zero_cut, lib_cut, group, struc_zero,
              neg_lb, tol, max_iter, conserve, alpha, global)
