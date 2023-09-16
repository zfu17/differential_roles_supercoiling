## Computes loci local contact sum within each group (PP.NPNP.)
## Author: Ziqi Fu
## Email: zfu17@jhu.edu
## This script uses multi-core processing (2 cores). If undesired, change mc.cores to 1. 

##### packages #####
{
library(ggplot2)
library(tidyquant)
library(tidyverse)
library(data.table)
library(reshape2)
library(parallel)
library(tidyquant)
}

##### functions do-not-modify #####
## special version, raw
toMat = function(f){
  temp = fread(f)%>%as.data.frame()
  mat = matrix(0,nrow = max(temp[,1:2]),ncol=max(temp[,1:2]))
  for (i in seq_along(temp[,1])){
    mat[temp[i,2],temp[i,1]] = temp[i,3]
    mat[temp[i,1],temp[i,2]] = temp[i,3]
  }
  colnames(mat) = 1:dim(mat)[2]
  rownames(mat) = 1:dim(mat)[1]
  return((mat)) ## spatial distance
} 
dis = function(x,y){
  temp = rep(NA,length(x))
  for (i in 1:length(x)){
    if (abs(x[i]-y[i])<ceiling(929/2)){
      temp[i] = (abs(x[i]-y[i]))
    } else{
      temp[i] = (929-abs(x[i]-y[i]))
    }
  }
  return(temp)
}
dis_500 = function(x,y){
  temp = rep(NA,length(x))
  for (i in 1:length(x)){
    if (abs(x[i]-y[i])<ceiling(9284/2)){
      temp[i] = (abs(x[i]-y[i]))
    } else{
      temp[i] = (9284-abs(x[i]-y[i]))
    }
  }
  return(temp)
}

mytheme = theme_bw(base_size = 14)+
  theme(
    #plot.background = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = 'black'),
    legend.text=element_text(size=12)
) 


##### Inputting HiC matrices, res=5k or .5k, ICED #####
setwd('/Users/ziqi_fu/Documents/R_scripts/Ecoli/hic_matrices/')
#load('hic_matDF_hires_list.Rdata') 
load('hic_matDF_list.Rdata')

##### Inputting reproducible peaks (no filter applied)
setwd('/Users/ziqi_fu/Documents/R_scripts/Ecoli/peaks')
gapr_df = fread('gapr_peaks.csv')
topo_df = fread('topo_peaks.csv')
gapr_peak = gapr_df$`5k_mdpt_bin` %>% unique() 
topo_peak = topo_df$`5k_mdpt_bin` %>% unique() 
overlap = intersect(gapr_df$`5k_mdpt_bin`,topo_df$`5k_mdpt_bin`)
nopeak = setdiff(1:929,union(gapr_df$`5k_mdpt_bin`,topo_df$`5k_mdpt_bin`))
gapR_only = setdiff(gapr_df$`5k_mdpt_bin`,topo_df$`5k_mdpt_bin`)
topo_only = setdiff(topo_df$`5k_mdpt_bin`,gapr_df$`5k_mdpt_bin`)
no_gapR = setdiff(1:929,gapr_peak)
no_topo = setdiff(1:929,topo_peak)

##### create ggplot data frame with peak information ##### 
get_extended_5k = function(mat_list) {
  mat=lapply(mat_list, function(x){
    temp = as.data.frame(melt(x))
    colnames(temp) = c('loci1','loci2','contact')
    temp = mutate(temp,distance=dis(loci1,loci2))
    temp = temp %>% mutate(
      l1_gapr = loci1 %in% gapr_peak,
      l1_topo = loci1 %in% topo_peak,
      l1_overlap = loci1 %in% overlap,
      l1_gapr_only = loci1 %in% gapR_only,
      l1_topo_only = loci1 %in% topo_only,
      l1_nopeak = loci1 %in% nopeak,
      l1_no_gapr = loci1 %in% no_gapR,
      l1_no_topo = loci1 %in% no_topo,
      l2_gapr = loci2 %in% gapr_peak,
      l2_topo = loci2 %in% topo_peak,
      l2_overlap = loci2 %in% overlap,
      l2_gapr_only = loci2 %in% gapR_only,
      l2_topo_only = loci2 %in% topo_only,
      l2_nopeak = loci2 %in% nopeak,
      l2_no_gapr = loci2 %in% no_gapR,
      l2_no_topo = loci2 %in% no_topo
    )
    return(temp%>%as.data.frame())})
  return(mat)
}
hic_mat_extended = get_extended_5k(hic_matDF_list)

##### summing up for 5k #####
setwd('/Users/ziqi_fu/Documents/R_scripts/Ecoli/figures/local_group_sum/')
N = 37 #breaking the chromosome into N pieces
 
#### pp vs npnp for gapR
i=1 #matrix index
gpgp = mclapply(1:928,function(pos){
  m = filter(hic_mat_extended[[i]],loci1==pos)
  if(prod(m$l1_gapr)==1){
    m=filter(m,l2_gapr==1 & distance>1)
  }else{
    m=filter(m,l2_gapr==0 & distance>1)
  } ##p-p only
  return(log2(1/m$contact) %>% mean())
},mc.cores = 5) %>% unlist()

df = data.frame(pos=1:928*5/1000,gpgp,Peak = ifelse(1:928 %in% gapr_peak,'PS','aPS')) %>% filter(gpgp>0)
df %>% ggplot(aes(x=pos,y=gpgp,color=Peak)) +
  geom_point(size=1,alpha=0.3)+
  geom_ma(ma_fun = SMA,n=20,linetype='solid',size=4) +
  mytheme+
  labs(y='Average contacts',x='Genomic Coordinate (Mb)')

#### pall vs npall for gapR
i=1 #matrix index
gpgp = mclapply(1:928,function(pos){
  m = filter(hic_mat_extended[[i]],loci1==pos & distance > 20)
  return(sum(log2(1/m$contact)))
},mc.cores = 5) %>% unlist()

df = data.frame(pos=1:928,gpgp,Peak = ifelse(1:928 %in% gapr_peak,'PS','aPS')) %>% filter(gpgp>0)
df %>% ggplot(aes(x=pos*5/1000,y=gpgp,color=Peak)) +
  geom_point(size=1,alpha=0.3)+
  geom_ma(ma_fun = SMA,n=20,linetype='solid',size=4) +
  mytheme+
  labs(y='Average contacts',x='Genomic Coordinate (Mb)')


## pick a window, then we can perform a statistical test
win = 40
lapply(seq(0,850,80),function(x){
  w = x:(x+win)
  ps = filter(df,pos%in%w & Peak == 'PS')$gpgp 
  aps = filter(df,pos%in%w & Peak == 'aPS')$gpgp 
  if(length(ps)*length(aps)>0){
    test = t.test(ps,aps,alternative = 'greater')
    return(test$p.value)
  } else(
    return(NA)
  )
}) %>% unlist() %>% as.data.frame() %>% View()



#### pp vs npnp for topo
i=1 #matrix index
gpgp = mclapply(1:928,function(pos){
  m = filter(hic_mat_extended[[i]],loci1==pos)
  if(prod(m$l1_topo)==1){
    m=filter(m,l2_topo==1 & distance>20)
  }else{
    m=filter(m,l2_topo==0 & distance>20)
  } ##p-p only
  return(mean(log2(1/m$contact)))
},mc.cores = 5) %>% unlist()

df = data.frame(pos=1:928*5/1000,gpgp,Peak = ifelse(1:928 %in% topo_peak,'NS','aNS')) %>% filter(gpgp>0)
df %>% ggplot(aes(x=pos,y=gpgp,color=Peak)) +
  geom_point(size=1,alpha=0.3)+
  geom_ma(ma_fun = SMA,n=20,linetype='solid',size=4) +
  mytheme+
  labs(y='Average contacts',x='Genomic Coordinate (Mb)')

#### pall vs npall for topo
i=1 #matrix index
gpgp = mclapply(1:928,function(pos){
  m = filter(hic_mat_extended[[i]],loci1==pos & distance > 20)
  return(sum(log2(1/m$contact)))
},mc.cores = 5) %>% unlist()

df = data.frame(pos=1:928*5/1000,gpgp,Peak = ifelse(1:928 %in% topo_peak,'NS','aNS')) %>% filter(gpgp>0)
df %>% ggplot(aes(x=pos,y=gpgp,color=Peak)) +
  geom_point(size=1,alpha=0.3)+
  geom_ma(ma_fun = SMA,n=20,linetype='solid',size=4) +
  mytheme+
  labs(y='Average contacts',x='Genomic Coordinate (Mb)')










