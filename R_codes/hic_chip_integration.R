## HiC and ChIP-seq integration
## Author: Ziqi Fu
## Email: zfu17@jhu.edu
## This script uses multi-core processing (4 cores). Change mc.cores to 1 to avoid parallel processing. 

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

mytheme = theme_bw(base_size = 14)+
  theme(
    #plot.background = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = 'black'),
    legend.text=element_text(size=12)+
    theme(legend.background=element_rect(fill = alpha("white", 0.3)))
) 
}

##### Inputting HiC matrices, res=5k, ICED #####
setwd('/Users/ziqi_fu/Documents/R_scripts/Ecoli/hic_matrices/')
# filenames = c('lioy_LB1.matrix','lioy_LB2.matrix',
              # 'monica_LB1.matrix','monica_LB2.matrix')
# filenames_hires = c('CC_HpaII_rep1_500.matrix','CC_HpaII_rep2_500.matrix')
# filenames = c('HUab_1.matrix','HUab_2.matrix')
# filenames = c('CC_HpaII_rep1_500.matrix','CC_HpaII_rep2_500.matrix','CC_MluCI_500.matrix','CC_HpaII_MluCI_500.matrix')
# hic_matDF_list = mclapply(filenames,function(x) {toMat(x)},mc.cores = 4) # toMat returns 1/(1+count)
load('hic_matDF_hires_list.Rdata') #0.5 Kb resolution
load('hic_matDF_list.Rdata')
# ter_reg = 91:555 ## define ter region bins
# ori_reg = setdiff(1:928,ter_reg)
# hic_matDF_list_ter = lapply(hic_matDF_list,function(x) x[ter_reg,ter_reg])
# hic_matDF_list_ori = lapply(hic_matDF_list,function(x) x[ori_reg,ori_reg])
gc()

##### Inputting reproducible peaks (no filter applied)
setwd('/Users/ziqi_fu/Documents/R_scripts/Ecoli/peaks')
gapr_df = fread('gapr_peaks.csv')
# gapr_df = fread('gapr76_peaks.csv') ## make gapr76 figures
topo_df = fread('topo_peaks.csv')
# gapr76_df = fread('gapr76_peaks.csv')
# gapr_df = mutate(gapr_df,reg = ifelse(`5k_mdpt_bin`%in% ter_reg,'ter','ori'))
# topo_df = mutate(topo_df,reg = ifelse(`5k_mdpt_bin`%in% ter_reg,'ter','ori'))
# gapr76_df = mutate(gapr76_df,reg = ifelse(`5k_mdpt_bin`%in% ter_reg,'ter','ori'))
# see reproducible peak processing script written by Ken

##### creating peak vectors #####
gapr_peak = gapr_df$`5k_mdpt_bin` %>% unique() 
topo_peak = topo_df$`5k_mdpt_bin` %>% unique() 
overlap = intersect(gapr_df$`5k_mdpt_bin`,topo_df$`5k_mdpt_bin`)
nopeak = setdiff(1:929,union(gapr_df$`5k_mdpt_bin`,topo_df$`5k_mdpt_bin`))
gapR_only = setdiff(gapr_df$`5k_mdpt_bin`,topo_df$`5k_mdpt_bin`)
topo_only = setdiff(topo_df$`5k_mdpt_bin`,gapr_df$`5k_mdpt_bin`)
no_gapR = setdiff(1:929,gapr_peak)
no_topo = setdiff(1:929,topo_peak)

#######################################################  
##### gapr and topo peaks distribution statistics #####
####################################################### 
peak_width_dis_plot = data.frame(width=c(gapr_df$width,topo_df$width,gapr76_df$width),
           peak=c(rep('GapR (N=609)',length(gapr_df$width)),rep('Topo (N=446)',
                  length(topo_df$width)),rep('GapR1-76 (N=482)',length(gapr76_df$width)))) %>% ggplot() +
  geom_density(aes(x=log2(width),fill=peak),alpha=.3) +
  mytheme+
  theme(legend.position = c(.8,.8)) +
  theme(legend.background=element_rect(fill = alpha("white", .3)))+
  labs(fill='',x='Log2 normalized peak width')
ggsave('/Users/ziqi_fu/Documents/R_scripts/Ecoli/figures/peak_width_dis_plot.pdf',
       peak_width_dis_plot,width=5.5,height=3.5)

## gapr peak distribution in ter and ori
gapr_width_ter_ori = gapr_df %>% ggplot(aes(x=log2(width),after_stat(count),fill=reg))+
  geom_density(alpha=.5,position = "stack")+
  mytheme+
  labs(fill='',x='Log2 normalized GapR peak width')+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 8))+
  theme(legend.position =c(.8,.8))+
  theme(legend.background=element_rect(fill = alpha("white", .3)))
gapr_width_ter_ori
ggsave('/Users/ziqi_fu/Documents/R_scripts/Ecoli/figures/gapr_width_ter_ori.pdf',
       gapr_width_ter_ori,width=4.5,height=3)  

## topo peak distribution in ter and ori
topo_width_ter_ori = topo_df %>% ggplot(aes(x=log2(width),after_stat(count),fill=reg))+
  geom_density(alpha=.5,position = "stack")+
  mytheme+
  labs(fill='',x='Log2 normalized topo peak width')+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 8))+
  theme(legend.position =c(.8,.8))+
  theme(legend.background=element_rect(fill = alpha("white", .3)))
topo_width_ter_ori
ggsave('/Users/ziqi_fu/Documents/R_scripts/Ecoli/figures/topo_width_ter_ori.pdf',
       topo_width_ter_ori,width=4.5,height=3) 

## gapr1-76 peak distribution in ter and ori
gapr76_width_ter_ori = gapr76_df %>% ggplot(aes(x=log2(width),after_stat(count),fill=reg))+
  geom_density(alpha=.5,position = "stack")+
  mytheme+
  labs(fill='',x='Log2 normalized GapR peak width')+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 8))+
  theme(legend.position =c(.8,.8))+
  theme(legend.background=element_rect(fill = alpha("white", .3)))
gapr76_width_ter_ori
ggsave('/Users/ziqi_fu/Documents/R_scripts/Ecoli/figures/gapr76_width_ter_ori.pdf',
       gapr76_width_ter_ori,width=4.5,height=3) 

########################################################## 
##### create ggplot data frame with peak information #####
########################################################## 
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
# hic_mat_extended_ter = get_extended_5k(hic_matDF_list_ter)
# hic_mat_extended_ori = get_extended_5k(hic_matDF_list_ori)


##### generate plots for res = 5k ##### 
setwd('/Users/ziqi_fu/Documents/R_scripts/Ecoli/figures/res5k_new/')  

get_plot_df_5k = function(df_list,peak) {
  if (peak=='gapr'){
    res=mclapply(df_list,function(x){
      pp = filter(x,l1_gapr==T & l2_gapr==T)%>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
      p_all = filter(x,l1_gapr==T)%>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
      npnp = filter(x,l1_gapr==F & l2_gapr==F)%>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
      np_all = filter(x,l1_gapr==F)%>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
      all = x %>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
      df = rbind(pp,npnp,all)
      s = c(rep('PS-PS',dim(pp)[1]),rep('aPS-aPS',dim(npnp)[1]),rep('all-all',dim(all)[1]))
      df = cbind(df,status=s)
      return(df)
    },mc.cores=4) }
  if (peak == 'topo') {
    res=mclapply(df_list,function(x){
      pp = filter(x,l1_topo==T & l2_topo==T)%>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
      p_all = filter(x,l1_topo==T)%>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
      npnp = filter(x,l1_topo==F & l2_gapr==F)%>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
      np_all = filter(x,l1_topo==F)%>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
      all = x %>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
      df = rbind(pp,npnp,all)
      s = c(rep('NS-NS',dim(pp)[1]),rep('aNS-aNS',dim(npnp)[1]),rep('all-all',dim(all)[1]))
      df = cbind(df,status=s)
      return(df)
    },mc.cores=4)
  }
  return(res)
}
# plot_df_ter = get_plot_df_5k(hic_mat_extended_ter,'gapr')
# plot_df_ori = get_plot_df_5k(hic_mat_extended_ori,'gapr')
plot_df = get_plot_df_5k(hic_mat_extended,'gapr')

plot_pp_5k = function(plot_df_ls){
  plots = lapply(plot_df_ls,function(x){
    p = ggplot(data=x,aes(x=5*distance,y=avg,col=status))+
      geom_point(size=0.05,alpha=0.3)+
      geom_ma(ma_fun = SMA, n=15,linetype='solid',size=0.8)+
      mytheme+
      theme(legend.position = c(.8,.35))+
      theme(legend.background=element_rect(fill = alpha("white", 0.2)))+
      labs(y="1 / (1 + contact)",x="Loci linear separation distance (Kb)",color='')
    return (p)
  })
  return(plots)
}
gapr_pp_plots = plot_pp_5k(plot_df)
# gapr_pp_plots_ter = plot_pp_5k(plot_df_ter)
# gapr_pp_plots_ori = plot_pp_5k(plot_df_ori)

for(i in 1:4){
  print(gapr_pp_plots[[i]])
  ggsave(paste0('gapr_500_',toString(i),'.pdf'),device="pdf",height=2.5,width=4.5)
}

##### plotting gapR difference plots ##### 
plot_df = lapply(hic_mat_extended,function(x){
  pp = filter(x,l1_gapr==T & l2_gapr==T)%>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
  p_all = filter(x,l1_gapr==T)%>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
  npnp = filter(x,l1_gapr==F & l2_gapr==F)%>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
  np_all = filter(x,l1_gapr==F)%>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
  aa = x %>% group_by(distance) %>% summarise_at(vars(contact),list(avg=mean))
  df = data.frame((npnp-pp)/max(npnp-pp),dis=1:dim(npnp-pp)[1])
  return(df)
})
gapr_diff_plots = lapply(plot_df,function(x){
  p = ggplot(data=x,aes(x=5*dis,y=avg))+
    geom_ma(ma_fun = SMA, n=20,linetype='solid',size=1.2)+
    geom_point(size=0.3,alpha=0.3)+
    geom_hline(yintercept = 0,color='gray2',lwd=0.8,linetype='dashed')+
    mytheme+
    # coord_cartesian(ylim=c(0,.6))+ # for gapR_diff
    labs(y="Diff. in 1 / (1+contact)",x="Loci linear separation distance (Kb)",color='')
    # theme(legend.position = c(0.75,0.8))
    # coord_cartesian(ylim=c(-0.05,0.25))
  return (p)
})

for(i in 1:4){
  print(gapr_diff_plots[[i]])
  ggsave(paste0('gapR_diff_',toString(i),'.pdf'),units='in',
         height=2.5,width=4.5)
}

##### plotting topo plots ##### 
# plot_df_topo = lapply(hic_mat_extended,function(x){
#   pp = filter(x,l1_topo==T & l2_topo==T)%>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
#   p_all = filter(x,l1_topo==T)%>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
#   npnp = filter(x,l1_topo==F & l2_gapr==F)%>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
#   np_all = filter(x,l1_topo==F)%>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
#   all = x %>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
#   
#   df = rbind(pp,npnp,all)
#   s = c(rep('NS-NS',dim(pp)[1]),rep('aNS-aNS',dim(npnp)[1]),rep('all-all',dim(all)[1]))
#   df = cbind(df,status=s)
#   return(df)
# })

plot_df_topo = get_plot_df_5k(hic_mat_extended,'topo')

topo_plots = lapply(plot_df_topo,function(x){
  p = ggplot(data=x,aes(x=5*distance,y=avg,col=status))+
    geom_point(size=0.2,alpha=0.3)+
    geom_ma(ma_fun = SMA, n=15,linetype='solid',size=0.8)+
    mytheme+
    theme(legend.position = c(.78,.35))+
    labs(y="1 / (1 + contact)",x="Loci linear separation distance (Kb)",color='')+
    theme(legend.background=element_rect(fill = alpha("white", .3)))
  return (p)
})
# setwd('/Users/ziqi_fu/Documents/R_scripts/Ecoli/figures/')
for(i in 1:4){
  print(topo_plots[[i]])
  ggsave(paste0('topo',toString(i),'.pdf'),units='in',
         height=2.5,width=4.5)
}

##### plotting topo difference #####
plot_df_topo_diff = lapply(hic_mat_extended,function(x){
  pp = filter(x,l1_topo==T & l2_topo==T)%>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
  p_all = filter(x,l1_topo==T)%>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
  npnp = filter(x,l1_topo==F & l2_gapr==F)%>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
  np_all = filter(x,l1_topo==F)%>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
  #all = x %>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
  
  # df = rbind(npnp-pp,np_all-p_all)
  # s = c(rep('tNP-tNP vs tP-tP',dim(npnp/pp)[1]),rep('tNP-all vs tP-all',dim(np_all/p_all)[1]))
  # df = cbind(df,status=s)
  aa = x %>% group_by(distance) %>% summarise_at(vars(contact),list(avg=mean))
  df = data.frame((npnp-pp)/max(npnp-pp),dis=rep(1:465,2))
  # df = mutate(df,distance=rep(1:4643,2)) ## for 500bp resolution
  return(df)
})
topo_diff_plots = lapply(plot_df_topo_diff,function(x){
  p = ggplot(data=x,aes(x=5*dis,y=avg))+
    geom_ma(ma_fun = SMA, n=30,linetype='solid',size=1.4)+
    geom_point(size=0.3,alpha=0.2)+
    geom_hline(yintercept = 0,color='gray2',lwd=0.8,linetype='dashed')+
    mytheme+
    labs(y="Diff. in 1 / (1+contact)",x="Loci linear separation distance (Kb)",color='')+
    theme(legend.position = c(0.75,0.8))
    # coord_cartesian(ylim=c(-0.05,0.15))
  return (p)
})
# setwd('/Users/ziqi_fu/Documents/R_scripts/Ecoli/figures/')
for(i in 1:4){
  print(topo_diff_plots[[i]])
  ggsave(paste0('topo_diff',toString(i),'.pdf'),units='in',
         height=2.5,width=4.5)
}


##### plotting gapR vs topo plots #####
plot_df_gt = lapply(hic_mat_extended,function(x){
  gg = filter(x,l1_gapr_only==T & l2_gapr_only==T)%>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
  tt = filter(x,l1_topo_only==T & l2_topo_only==T)%>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
  npnp = filter(x,l1_nopeak==T & l2_nopeak==T)%>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
  #all = x %>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
  
  df = rbind(gg,tt,npnp)
  s = c(rep('PSonly-PSonly',dim(gg)[1]),rep('NSonly-NSonly',dim(tt)[1]),
        rep('aa-aa',dim(npnp)[1]))
  df = cbind(df,status=s)
  return(df)
})
gapRtopo_plots = lapply(plot_df_gt,function(x){
  p = ggplot(data=x,aes(x=5*distance,y=avg,col=status))+
    geom_point(size=0.1,alpha=0.3)+
    geom_ma(ma_fun = SMA, n=15,linetype='solid',size=0.6)+
    mytheme+
    theme(legend.position = c(0.75,0.35))+
    labs(y="1 / (1 + contact)",x="Loci linear separation distance (Kb)",color='')+
    theme(legend.background=element_rect(fill = alpha("white", .3)))
  return (p)
})
# setwd('/Users/ziqi_fu/Documents/R_scripts/Ecoli/figures/')
for(i in 1:4){
  print(gapRtopo_plots[[i]])
  ggsave(paste0('gapR_vs_topo_',toString(i),'.pdf'),units='in',
         height=2.5,width=4.5)
}

## test statistical significant

temp = filter(plot_df[[1]],distance>300) 
lm(data=temp,avg ~ distance + status) %>% summary()
temp = filter(plot_df[[1]],distance<350)
lm(data=temp,avg ~ distance + status) %>% summary()

############################################
##### plotting gapR vs topo difference #####
############################################
plot_df = lapply(hic_mat_extended,function(x){
  gg = filter(x,l1_gapr_only==T & l2_gapr_only==T)%>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
  tt = filter(x,l1_topo_only==T & l2_topo_only==T)%>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
  npnp = filter(x,l1_nopeak==T & l2_nopeak==T)%>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
  #all = x %>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
  
  df = rbind(npnp-tt,npnp-gg)
  s = c(rep('NN vs tOnly',dim(gg)[1]),rep('NN vs gOnly',dim(tt)[1]))
  df = cbind(df,status=s)
  df = mutate(df,distance=rep(1:465,2))
  # df = mutate(df,distance=rep(1:4643,2))
  return(df)
})
gapRtopo_diff_plots = lapply(plot_df,function(x){
  p = ggplot(data=x,aes(x=5*distance,y=avg,col=status))+
    geom_point(size=0.3,alpha=0.7)+
    geom_ma(ma_fun = SMA, n=15,linetype='solid',lwd=0.8)+
    geom_hline(yintercept = 0,color='gray2',lwd=0.8,linetype='dashed')+
    mytheme+
    theme(legend.position = c(0.8,0.8))+
    coord_cartesian(ylim=c(-0.1,0.3))+
    #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    labs(y="Difference in 1 / (1+ICED contact counts)",
         x="Loci linear separation distance (Kb)",color='')
  return (p)
})

setwd('/Users/ziqi_fu/Documents/R_scripts/Ecoli/figures/')
for(i in 1:4){
  print(gapRtopo_diff_plots[[i]])
  ggsave(paste0('gapr_topo_diff_',toString(i),'.png'),device="png",units='in',
         height=3,width=4.5)
}

### test
filter(plot_df[[1]],distance>400 & status == 'NN vs tOnly')$avg %>% t.test(
  y=filter(plot_df[[1]],distance>400 & status == 'NN vs gOnly')$avg,paired=T
)


##### plotting gapR vs topo difference 2 #####
plot_df_gt_diff = lapply(hic_mat_extended,function(x){
  gg = filter(x,l1_gapr_only==T & l2_gapr_only==T)%>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
  tt = filter(x,l1_topo_only==T & l2_topo_only==T)%>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
  npnp = filter(x,l1_nopeak==T & l2_nopeak==T)%>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
  #all = x %>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
  
  df = data.frame((tt-gg)/max(tt-gg),status=rep('tOnly vs gOnly',dim(tt)[1]))
  df= mutate(df,distance=1:465)
  # df = mutate(df,distance=rep(1:4643,1))
  return(df)
})

gapRtopo_diff_plots = lapply(plot_df_gt_diff,function(x){
  p = ggplot(data=x,aes(x=5*distance,y=avg))+
    geom_point(size=0.3,alpha=0.25)+
    geom_ma(ma_fun = SMA, n=20,linetype='solid',size=1.2)+
    geom_hline(yintercept = 0,color='gray2',lwd=0.8,linetype='dashed')+
    mytheme+
    # theme(legend.position = c(.8,.85))+
    # coord_cartesian(ylim=c(-0.1,0.3))+
    #scale_color_manual(values=c('black'))+
    #scale_color_discrete(name = "tOnly : gOnly",labels = "tOnly : gOnly")+
    #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    labs(y="Diff. in 1 / (1+contact)",
         x="Loci linear separation distance (Kb)",color='')
  return (p)
})
# setwd('/Users/ziqi_fu/Documents/R_scripts/Ecoli/figures/')
for(i in 1:4){
  print(gapRtopo_diff_plots[[i]])
  ggsave(paste0('gapr_topo_diff_',toString(i),'.pdf'),units='in',
         height=2.5,width=4.5)
}

#############################
##### resolution = 500 #####
#############################
gapr_df$start_bin_500 = ceiling(gapr_df$start/500)
gapr_df$end_bin_500 = ceiling (gapr_df$end/500)
gapr_int_ls_500= list()
for (i in 1:dim(gapr_df)[1]){
  gapr_int_ls_500[[i]] = gapr_df$start_bin_500[i]:gapr_df$end_bin_500[i]
}
gapr_peak_500 = unlist(gapr_int_ls_500) %>% unique()

topo_df$start_bin_500 = ceiling(topo_df$start/500)
topo_df$end_bin_500 = ceiling (topo_df$end/500)
topo_int_ls_500= list()
for (i in 1:dim(topo_df)[1]){
  topo_int_ls_500[[i]] = topo_df$start_bin_500[i]:topo_df$end_bin_500[i]
}
topo_peak_500 = unlist(topo_int_ls_500) %>% unique()

overlap_500 = intersect(gapr_peak_500,topo_peak_500)
nopeak_500 = setdiff(1:9284,union(gapr_peak_500,topo_peak_500))
gapR_only_500 = setdiff(gapr_peak_500,topo_peak_500)
topo_only_500 = setdiff(topo_peak_500,gapr_peak_500)
no_gapR_500 = setdiff(1:9284,gapr_peak_500)
no_topo_500 = setdiff(1:9284,topo_peak_500)

# ter_reg = 898:5540
# ori_reg = setdiff(1:9284,ter_reg)
# hic_matDF_hires_list_ter = mclapply(hic_matDF_hires_list,function(x) {x[ter_reg,ter_reg]},mc.cores=2)
# hic_matDF_hires_list_ori = mclapply(hic_matDF_hires_list,function(x) {x[ori_reg,ori_reg]},mc.cores = 2)

###### 500 bp plots #####
get_extended_500 = function(mat_list) {
  mat=lapply(mat_list, function(x){
    temp = as.data.frame(melt(x))
    colnames(temp) = c('loci1','loci2','contact')
    temp = mutate(temp,distance=dis_500(loci1,loci2))
    temp = temp %>% mutate(
      l1_gapr = loci1 %in% gapr_peak_500,
      l1_topo = loci1 %in% topo_peak_500,
      l1_overlap = loci1 %in% overlap_500,
      l1_gapr_only = loci1 %in% gapR_only_500,
      l1_topo_only = loci1 %in% topo_only_500,
      l1_nopeak = loci1 %in% nopeak_500,
      l1_no_gapr = loci1 %in% no_gapR_500,
      l1_no_topo = loci1 %in% no_topo_500,
      l2_gapr = loci2 %in% gapr_peak_500,
      l2_topo = loci2 %in% topo_peak_500,
      l2_overlap = loci2 %in% overlap_500,
      l2_gapr_only = loci2 %in% gapR_only_500,
      l2_topo_only = loci2 %in% topo_only_500,
      l2_nopeak = loci2 %in% nopeak_500,
      l2_no_gapr = loci2 %in% no_gapR_500,
      l2_no_topo = loci2 %in% no_topo_500
    )
    return(temp)}) 
  return(mat)
}
hic_mat_extended_500 = get_extended_500(hic_matDF_hires_list)

# hic_mat_extended_500_ter = get_extended_500(hic_matDF_hires_list_ter)
# hic_mat_extended_500_ori = get_extended_500(hic_matDF_hires_list_ori)

##### plotting gapR plots ##### 
setwd('/Users/ziqi_fu/Documents/R_scripts/Ecoli/figures/res500_new/')  
gapr_plot_df = mclapply(hic_mat_extended_500,function(x){
  pp = filter(x,l1_gapr==T & l2_gapr==T)%>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
  p_all = filter(x,l1_gapr==T)%>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
  npnp = filter(x,l1_gapr==F & l2_gapr==F)%>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
  np_all = filter(x,l1_gapr==F)%>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
  all = x %>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
  
  df = rbind(pp,npnp,all)
  s = c(rep('PS-PS',dim(pp)[1]),rep('aPS-aPS',dim(npnp)[1]),rep('all-all',dim(all)[1]))
  df = cbind(df,status=s)
  return(df)
},mc.cores=2)

gapr_plots_500 = lapply(gapr_plot_df,function(x){
  p = ggplot(data=x,aes(x=.5*distance,y=avg,col=status))+
    geom_point(size=0.05,alpha=0.3)+
    geom_ma(ma_fun = SMA, n=30,linetype='solid',size=0.8)+
    mytheme+
    labs(y="1 / (1 + contact)",x="Loci linear separation distance (Kb)",color='')+
    theme(legend.position = c(.8,.35))+
    theme(legend.background=element_rect(fill = alpha("white", .3)))
  return (p)
})
for(i in 1:2){
  print(gapr_plots_500[[i]])
  ggsave(paste0('gapr_500_',toString(i),'.pdf'),device="pdf",units='in',
         height=2.5,width=4.5)
}
##### plotting gapR difference plots ##### 
gapr_diff_plot_df = mclapply(hic_mat_extended_500,function(x){
  pp = filter(x,l1_gapr==T & l2_gapr==T)%>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
  p_all = filter(x,l1_gapr==T)%>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
  npnp = filter(x,l1_gapr==F & l2_gapr==F)%>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
  np_all = filter(x,l1_gapr==F)%>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
  #all = x %>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
  
  # df = data.frame(npnp-pp)
  # s = c(rep('gNP-gNP vs gP-gP',dim(npnp-pp)[1]),rep('gNP-all vs gP-all',dim(np_all-p_all)[1]))
  # df = cbind(df,status=s)
  if (dim(npnp)[1]<dim(pp)[1]) {pp=filter(pp,distance %in% npnp$distance)}
  if (dim(npnp)[1]>dim(pp)[1]) {npnp=filter(npnp,distance %in% pp$distance)}
  df = data.frame((npnp-pp)/max(npnp-pp),dis=1:dim(npnp-pp)[1])
  # df = mutate(df,distance=rep(1:4643,2)) ## for 500 resolution
  return(df)
},mc.cores=2)

gapr_diff_plots = lapply(gapr_diff_plot_df,function(x){
  p = ggplot(data=x,aes(x=.5*dis,y=avg))+
    geom_ma(ma_fun = SMA, n=100,linetype='solid',size=2)+
    geom_point(size=0.1,alpha=0.08)+
    geom_hline(yintercept = 0,color='gray2',lwd=0.8,linetype='dashed')+
    mytheme+
    # coord_cartesian(ylim=c(-.2,.2))+
    labs(y="Diff. in 1 / (1+contact)",x="Loci linear separation distance (Kb)",color='')
  # theme(legend.position = c(0.75,0.8))
  # coord_cartesian(ylim=c(-0.05,0.25))
  return (p)
})
# setwd('/Users/ziqi_fu/Documents/R_scripts/Ecoli/figures/')
for(i in 1:2){
  print(gapr_diff_plots[[i]])
  ggsave(paste0('gapR_diff_',toString(i),'.pdf'),units='in',
         height=2.5,width=4.5)
}

##### plotting topo plots ##### 
plot_df_topo = mclapply(hic_mat_extended_500,function(x){
  pp = filter(x,l1_topo==T & l2_topo==T)%>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
  p_all = filter(x,l1_topo==T)%>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
  npnp = filter(x,l1_topo==F & l2_gapr==F)%>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
  np_all = filter(x,l1_topo==F)%>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
  all = x %>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
  
  df = rbind(pp,npnp,all)
  s = c(rep('NS-NS',dim(pp)[1]),rep('aNS-aNS',dim(npnp)[1]),rep('all-all',dim(all)[1]))
  df = cbind(df,status=s)
  return(df)
},mc.cores=2)
topo_plots = lapply(plot_df_topo,function(x){
  p = ggplot(data=x,aes(x=.5*distance,y=avg,col=status))+
    geom_point(size=0.2,alpha=0.3)+
    geom_ma(ma_fun = SMA, n=15,linetype='solid',size=0.8)+
    mytheme+
    labs(y="1 / (1 + contact)",x="Loci linear separation distance (Kb)",color='')+
    theme(legend.position = c(.8,.35))+
    theme(legend.background=element_rect(fill = alpha("white", .3)))
  return (p)
})
# setwd('/Users/ziqi_fu/Documents/R_scripts/Ecoli/figures/')
for(i in 1:2){
  print(topo_plots[[i]])
  ggsave(paste0('topo_500_',toString(i),'.pdf'),units='in',
         height=2.5,width=4.5)
}

##### plotting topo difference #####
plot_df_topo_diff = mclapply(hic_mat_extended_500,function(x){
  pp = filter(x,l1_topo==T & l2_topo==T)%>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
  p_all = filter(x,l1_topo==T)%>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
  npnp = filter(x,l1_topo==F & l2_gapr==F)%>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
  np_all = filter(x,l1_topo==F)%>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
  #all = x %>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
  
  # df = rbind(npnp-pp,np_all-p_all)
  # s = c(rep('tNP-tNP vs tP-tP',dim(npnp/pp)[1]),rep('tNP-all vs tP-all',dim(np_all/p_all)[1]))
  # df = cbind(df,status=s)
  df = data.frame((npnp-pp)/max(npnp-pp),dis=1:dim(npnp-pp)[1])
  # df = mutate(df,distance=rep(1:4643,2)) ## for 500bp resolution
  return(df)
},mc.cores = 2)
topo_diff_plots = lapply(plot_df_topo_diff,function(x){
  p = ggplot(data=x,aes(x=.5*dis,y=avg))+
    geom_ma(ma_fun = SMA, n=100,linetype='solid',size=1.2)+
    geom_point(size=0.1,alpha=0.08)+
    geom_hline(yintercept = 0,color='gray2',lwd=0.8,linetype='dashed')+
    mytheme+
    labs(y="Diff. in 1 / (1+contact)",x="Loci linear separation distance (Kb)",color='')+
    theme(legend.position = c(0.75,0.8))+
    coord_cartesian(ylim=c(-0.4,1))
  return (p)
})
# setwd('/Users/ziqi_fu/Documents/R_scripts/Ecoli/figures/')
for(i in 1:2){
  print(topo_diff_plots[[i]])
  ggsave(paste0('topo_diff',toString(i),'.pdf'),units='in',
         height=2.5,width=4.5)
}


##### plotting gapR vs topo plots #####
plot_df_gt = lapply(hic_mat_extended_500,function(x){
  gg = filter(x,l1_gapr_only==T & l2_gapr_only==T)%>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
  tt = filter(x,l1_topo_only==T & l2_topo_only==T)%>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
  npnp = filter(x,l1_nopeak==T & l2_nopeak==T)%>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
  #all = x %>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
  
  df = rbind(gg,tt,npnp)
  s = c(rep('PSonly-PSonly',dim(gg)[1]),rep('NSonly-NSonly',dim(tt)[1]),
        rep('aa-aa',dim(npnp)[1]))
  df = cbind(df,status=s)
  return(df)
})
gapRtopo_plots = lapply(plot_df_gt,function(x){
  p = ggplot(data=x,aes(x=.5*distance,y=avg,col=status))+
    geom_point(size=0.1,alpha=0.1)+
    geom_ma(ma_fun = SMA, n=60,linetype='solid',size=0.6)+
    mytheme+
    # coord_cartesian(ylim=c(.4,.95))+ ## fig1 (.2,.9), fig2 (.4,.95)
    theme(legend.position = c(0.7,0.35))+
    labs(y="1 / (1+contact)",x="Loci linear separation distance (Kb)",color='')+
    theme(legend.background=element_rect(fill = alpha("white", .3)))
  return (p)
})

for(i in 1:2){
  print(gapRtopo_plots[[i]])
  ggsave(paste0('gapR_vs_topo_',toString(i),'.pdf'),units='in',
         height=2.5,width=4.5)
}

### gapr vs topo difference
plot_df_gt_diff = lapply(hic_mat_extended_500,function(x){
  gg = filter(x,l1_gapr_only==T & l2_gapr_only==T)%>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
  tt = filter(x,l1_topo_only==T & l2_topo_only==T)%>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
  npnp = filter(x,l1_nopeak==T & l2_nopeak==T)%>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
  #all = x %>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
  
  df = data.frame((tt-gg)/max(tt-gg),status=rep('tOnly vs gOnly',dim(tt)[1]))
  # df= mutate(df,distance=1:465)
  df = mutate(df,distance=rep(1:4643,1))
  return(df)
})
gapRtopo_diff_plots = lapply(plot_df_gt_diff,function(x){
  p = ggplot(data=x,aes(x=.5*distance,y=avg))+
    geom_point(size=0.3,alpha=0.08)+
    geom_ma(ma_fun = SMA, n=100,linetype='solid',size=1.2)+
    geom_hline(yintercept = 0,color='gray2',lwd=0.8,linetype='dashed')+
    mytheme+
    labs(y="Diff. in 1 / (1+contact)",
         x="Loci linear separation distance (Kb)",color='')+
    theme(legend.background=element_rect(fill = alpha("white", .3)))
  return (p)
})
# setwd('/Users/ziqi_fu/Documents/R_scripts/Ecoli/figures/')
for(i in 1:2){
  print(gapRtopo_diff_plots[[i]])
  ggsave(paste0('gapr_topo_diff_',toString(i),'.pdf'),units='in',
         height=2.5,width=4.5)
}

##### 0.5 Kb all-all comparison #####
all_ter_ori_plot = lapply(1:4,function(x){
  df_plot = rbind(filter(plot_df_ter[[x]],status=='all-all')%>%mutate(reg='ter'),
                  filter(plot_df_ori[[x]],status=='all-all')%>%mutate(reg='ori'))
  plot = df_plot %>% ggplot(aes(x=.5*distance,y=avg,color=reg)) +
    geom_point(size=.3,alpha=.5)+
    geom_ma(ma_fun = SMA, n=30,linetype='solid',lwd=2)+
    mytheme+
    theme(legend.position = c(.8,.2))+
    labs(color='Domain',y="1 / (1+ICED contact counts)",x="Loci linear separation distance (Kb)")
  return(plot)
})
grid.arrange(grobs=all_ter_ori_plot)
lapply(1:4,function(x){
  ggsave(paste0('/Users/ziqi_fu/Documents/R_scripts/Ecoli/figures/compaction_ter_ori/',filenames[x],'_ter_ori_compact.pdf'),
         all_ter_ori_plot[[x]],height=3,width=4.5)
})










