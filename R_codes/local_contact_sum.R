## Computes loci local contact sum, Fig 1F-H
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
  return(1/(1+mat))
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

mytheme = theme_bw(base_size = 12)+
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
# filenames = c('lioy_LB1.matrix','lioy_LB2.matrix',
#               'monica_LB1.matrix','monica_LB2.matrix')
# filenames = c('HUab_1.matrix','HUab_2.matrix')
# filenames = c('CC_HpaII_rep1_500.matrix','CC_HpaII_rep2_500.matrix','CC_MluCI_500.matrix','CC_HpaII_MluCI_500.matrix')
# filenames_hires = c('CC_HpaII_rep1_500.matrix','CC_HpaII_rep2_500.matrix')
# filenames = c('lioy_LB1.matrix','monica_LB1.matrix')
# hic_matDF_list = mclapply(filenames,function(x) {toMat(x)},mc.cores=4) # toMat returns 1/(1+count)
# hic_matDF_hires_list = mclapply(filenames_hires,function(x) {toMat(x)},mc.cores=2) # can be loaded directly
# names(hic_matDF_list) = filenames
# names(hic_matDF_hires_list) = filenames_hires
# save(hic_matDF_list,file='hic_matDF_list.Rdata')
# save(hic_matDF_hires_list,file='hic_matDF_hires_list.Rdata')

load('hic_matDF_hires_list.Rdata') 
load('hic_matDF_list.Rdata')

##### summing up for 5k #####
setwd('/Users/ziqi_fu/Documents/R_scripts/Ecoli/figures/window_local_contact_sum/')
win_size = 1000000 / 5000 #resolution=5k
gap = 100000 /5000
 
loc_sum = lapply(hic_matDF_list,function(mat){
  mat = cbind(mat,mat)
  return(lapply(1:dim(mat)[1],function(row){
    f = row:(row+win_size) + gap  
    return(1/(mat[row,f]%>%sum()))
  })%>%unlist())
})

loc_plot_5k=lapply(1:length(loc_sum),function(i){
  p=data.frame(bin=1:length(loc_sum[[i]]),sum=loc_sum[[i]]/max(loc_sum[[i]])) %>% ggplot() +
    geom_point(aes(x=5*bin,y=sum))+
    geom_smooth(aes(x=5*bin,y=sum),method='gam')+
    mytheme+
    theme(axis.title.y = element_text(size=11))+
    labs(x='Genomic coordinate (Kb)',y='Norm. contact sum in \n window 100Kb - 1Mb')
  ggsave(paste0('local_sum_res5k_',i,'.pdf'),p,width = 4,height = 2)
})

##### summing up for .5K #####
win_size = 1000000 / 500 #resolution=.5k
gap = 100000 /500

loc_sum_hires = lapply(hic_matDF_hires_list,function(mat){
  mat = cbind(mat,mat)
  return(mclapply(1:dim(mat)[1],function(row){
    f = row:(row+win_size) + gap  
    return(1/(mat[row,f]%>%sum()))
  },mc.cores = 5)%>%unlist())
})

loc_plot_0.5k=lapply(1:length(loc_sum_hires),function(i){
  p=data.frame(bin=1:length(loc_sum_hires[[i]]),sum=loc_sum_hires[[i]]/max(loc_sum_hires[[i]])) %>% ggplot() +
    geom_point(aes(x=.5*bin,y=sum),size=.8)+
    geom_smooth(aes(x=.5*bin,y=sum),method='gam')+
    mytheme+
    theme(axis.title.y = element_text(size=11))+
    labs(x='Genomic coordinate (Kb)',y='Norm. contact sum in \n window 100Kb - 1Mb')
  ggsave(paste0('local_sum_res0.5k_',i,'.pdf'),p,width = 4,height = 2)
})


##### chromosome structure all-all #####
hic_mat_extended = lapply(hic_matDF_list, function(x){
  temp = as.data.frame(melt(x))
  colnames(temp) = c('loci1','loci2','contact')
  temp = mutate(temp,distance=dis(loci1,loci2))
  return(temp)
  }
)
plot_df = lapply(names(hic_mat_extended), function(x){
  all = hic_mat_extended[[x]] %>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
  df = data.frame(contact=all,dataset=x)
  return (df)
})
all_5k_df=do.call(rbind.data.frame,plot_df) %>% mutate(contact.distance=5*contact.distance)

all_5k_df%>% ggplot(aes(x=contact.distance,y=contact.avg,color=dataset))+
  geom_point()+
  mytheme

hic_mat_extended_hires = mclapply(hic_matDF_hires_list, function(x){
  temp = as.data.frame(melt(x))
  colnames(temp) = c('loci1','loci2','contact')
  temp = mutate(temp,distance=dis_500(loci1,loci2))
  return(temp)
},mc.cores=2)
plot_df_hires = mclapply(names(hic_mat_extended_hires), function(x){
  all = hic_mat_extended_hires[[x]] %>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
  df = data.frame(contact=all,dataset=x)
  return (df)
},mc.cores=2)

all_500_df=do.call(rbind.data.frame,plot_df_hires) %>% mutate(contact.distance=.5*contact.distance)

all_plot=rbind(all_500_df,all_5k_df) %>% ggplot(aes(x=contact.distance,y=contact.avg,color=dataset))+
  geom_ma(ma_fun = SMA, n = 25,lwd=2,linetype='solid')+
  mytheme+
  labs(color='',x='Loci linear separation distance (Kb)',y='1 / (1+contact)')+
  theme(legend.position = c(.8,.2))+
  theme(legend.background=element_rect(fill = alpha("white", 0.3))) #7,5
ggsave('/Users/ziqi_fu/Documents/R_scripts/Ecoli/figures/all_all_structure.pdf',all_plot,width=7,height=5)






