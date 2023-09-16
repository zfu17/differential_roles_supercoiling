## This script plots all-all contacts for all Hi-C data sets Supplement S1D
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
{
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
}

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

load('~/Documents/R_scripts/Ecoli/hic_matrices/hic_matDF_hires_list.Rdata') ## 500
load('~/Documents/R_scripts/Ecoli/hic_matrices/hic_matDF_list.Rdata') ## 5k

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
  theme(legend.position = c(.8,.25),legend.key.width = unit(35,"pt"))+
  theme(legend.background=element_rect(fill = alpha("white", 0.1))) +#7,5 
  scale_color_discrete(labels=c('Cockram 1','Cockram 2','Lioy 1','Lioy 2','Guo 1','Guo 2'))
all_plot
ggsave('/Users/ziqi_fu/Documents/R_scripts/Ecoli/figures/all_all_structure.pdf',all_plot,width=6,height=4)






