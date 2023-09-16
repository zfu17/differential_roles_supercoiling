## HiC and ChIP-seq integration for machine learning section
## Author: Ziqi Fu
## Email: zfu17@jhu.edu
## This script uses multi-core processing (4 cores). Set mc.cores to 1 to avoid parallel processing. 

##### packages #####
{
  library(ggplot2)
  library(tidyquant)
  library(tidyverse)
  library(data.table)
  library(reshape2)
  library(parallel)
  library(glue)
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
      legend.text=element_text(size=12)
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
load('hic_matDF_hires_list.Rdata') # precomputed frequency matrices
load('hic_matDF_list.Rdata') # precomputed frequency matrices
# ter_reg = 91:555 ## define ter region bins
# ori_reg = setdiff(1:928,ter_reg)
# hic_matDF_list_ter = lapply(hic_matDF_list,function(x) x[ter_reg,ter_reg])
# hic_matDF_list_ori = lapply(hic_matDF_list,function(x) x[ori_reg,ori_reg])
gc()

##### Inputting reproducible peaks (no filter applied)
setwd('/Users/ziqi_fu/Documents/R_scripts/Ecoli/peaks')
gapr_df = fread('gapr_peaks.csv')
topo_df = fread('topo_peaks.csv')
gapr76_df = fread('gapr76_peaks.csv')
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


##### plots ######

get_extended_5k = function(mat_list) {
  mat=lapply(mat_list, function(x){
    temp = as.data.frame(melt(x))
    colnames(temp) = c('loci1','loci2','contact')
    # temp = mutate(temp,distance=dis(loci1,loci2))
    temp = mutate(temp,distance=loci1-loci2)
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
get_plot_df_5k = function(df_list,peak) {
  if (peak=='gapr'){
    res=mclapply(df_list,function(x){
      pp = filter(x,l1_gapr==T & l2_gapr==T)%>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
      p_all = filter(x,l1_gapr==T)%>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
      npnp = filter(x,l1_gapr==F & l2_gapr==F)%>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
      np_all = filter(x,l1_gapr==F)%>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
      all = x %>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
      df = rbind(pp,p_all,npnp,np_all,all)
      s = c(rep('PS-PS',dim(pp)[1]),
            rep('PS-all',dim(p_all)[1]),
            rep('aPS-aPS',dim(npnp)[1]),
            rep('aPS-all',dim(np_all)[1]),
            rep('all-all',dim(all)[1]))
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
      df = rbind(pp,p_all,npnp,np_all,all)
      s = c(rep('NS-NS',dim(pp)[1]),
            rep('NS-all',dim(p_all)[1]),
            rep('aNS-aNS',dim(npnp)[1]),
            rep('aNS-all',dim(np_all)[1]),
            rep('all-all',dim(all)[1]))
      df = cbind(df,status=s)
      return(df)
    },mc.cores=4)}
  return(res)
}
plot_pp_5k = function(plot_df_ls){
  plots = lapply(plot_df_ls,function(x){
    p = ggplot(data=x,aes(x=5*distance,y=avg,col=status))+
      geom_point(size=0.05,alpha=0.15)+
      geom_ma(ma_fun = SMA, n=30,linetype='solid',size=0.8)+
      mytheme+
      theme(legend.position = 'right')+
      theme(legend.background=element_rect(fill = alpha("white", 0.2)))+
      labs(y="1 / (1 + contact)",x="Loci linear separation distance (Kb)",color='')+
      scale_x_continuous(breaks=seq(-5000,5000,1000))
    return (p)
  })
  return(plots)
}

hic_mat_extended = get_extended_5k(hic_matDF_list)

plot_df_gapr = get_plot_df_5k(hic_mat_extended,'gapr')
gapr_pp_plots = plot_pp_5k(plot_df_gapr)

plot_df_topo = get_plot_df_5k(hic_mat_extended,'topo')
topo_pp_plots = plot_pp_5k(plot_df_topo)

setwd('~/Documents/R_scripts/Ecoli/figures/rf_peak_full/')
lapply(1:4,function(p){
  ggsave(glue('gapr_full_{p}.pdf'),gapr_pp_plots[[p]],width=7,height=3)
  ggsave(glue('topo_full_{p}.pdf'),topo_pp_plots[[p]],width=7,height=3)
})

##### generate peak to all only plots ##### 
get_plot_df_5k_peak_all = function(df_list,peak) {
  if (peak=='gapr'){
    res=mclapply(df_list,function(x){
      pp = filter(x,l1_gapr==T & l2_gapr==T)%>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
      p_all = filter(x,l1_gapr==T)%>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
      npnp = filter(x,l1_gapr==F & l2_gapr==F)%>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
      np_all = filter(x,l1_gapr==F)%>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
      all = x %>%group_by(distance) %>%summarise_at(vars(contact),list(avg=mean))
      df = rbind(p_all,np_all,all)
      s = c(rep('PS-all',dim(p_all)[1]),
            rep('aPS-all',dim(np_all)[1]),
            rep('all-all',dim(all)[1]))
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
      df = rbind(p_all,np_all,all)
      s = c(rep('NS-all',dim(p_all)[1]),
            rep('aNS-all',dim(np_all)[1]),
            rep('all-all',dim(all)[1]))
      df = cbind(df,status=s)
      return(df)
    },mc.cores=4)}
  return(res)
}
hic_mat_extended = get_extended_5k(hic_matDF_list)
plot_pp_5k_peak_all = function(plot_df_ls){
  plots = lapply(plot_df_ls,function(x){
    p = ggplot(data=x,aes(x=5*distance,y=avg,col=status))+
      geom_point(size=0.05,alpha=0.15)+
      geom_ma(ma_fun = SMA, n=30,linetype='solid',size=0.8)+
      mytheme+
      theme(legend.position = c(.7,.3),
            legend.background=element_rect(fill = alpha("white", 0.2)))+
      labs(y="1 / (1 + contact)",x="Loci linear separation distance (Kb)",color='')+
      scale_x_continuous(breaks=seq(-5000,5000,1000))
    return (p)
  })
  return(plots)
}

plot_df_gapr = get_plot_df_5k_peak_all(hic_mat_extended,'gapr')
gapr_pp_plots = plot_pp_5k_peak_all(plot_df_gapr)

plot_df_topo = get_plot_df_5k_peak_all(hic_mat_extended,'topo')
topo_pp_plots = plot_pp_5k_peak_all(plot_df_topo)

setwd('~/Documents/R_scripts/Ecoli/figures/rf_peak_full/')
lapply(1:4,function(p){
  ggsave(glue('gapr_peak_all_{p}.pdf'),gapr_pp_plots[[p]],width=5.5,height=3)
  ggsave(glue('topo_peak_all_{p}.pdf'),topo_pp_plots[[p]],width=5.5,height=3)
})





