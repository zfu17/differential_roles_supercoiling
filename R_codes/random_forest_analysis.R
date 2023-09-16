## This script loads and analyses trained random forests.
## Author: Ziqi Fu
## Email: zfu17@jhu.edu
## This script uses multi-core processing (5 cores). If undesired, change mc.cores to 1. 

##### packages #####
{
  library(tidyverse)
  library(tidyquant)
  library(data.table)
  library(reshape2)
  library(dplyr)
  library(data.table)
  library(ggplot2)
  library(randomForest)
  library(ROCR)
  library(pROC)
  library(caret)
  library(parallel)
  library(PRROC)
  library(stringr)
  library(glue)
  library(gmodels)
}

{
  toMat = function(f){
    temp = fread(f)%>%as.data.frame()
    mat = matrix(0,nrow = max(temp[,1:2]),ncol=max(temp[,1:2]))
    for (i in seq_along(temp[,1])){
      mat[temp[i,2],temp[i,1]] = temp[i,3]
      mat[temp[i,1],temp[i,2]] = temp[i,3]
    }
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
  make_nb_mat = function(temp_mat,LEN){
    temp_mat = 1/(1+temp_mat)[1:928,1:928] #log2 transformed, 928x928 size
    temp = matrix(0,nrow=928,ncol=2*LEN+1) #initialization
    for (i in 1:928){
      for (j in -LEN:LEN){
        if ((i+j)%%928!=0) {temp[i,j+LEN+1] = temp_mat[i,(i+j)%%928]}
        else {temp[i,j+LEN+1] = temp_mat[i,928]}
      }
    }
    return (temp)
  } #function that returns the mat
  
  make_nb_mat_500 = function(temp_mat,LEN){
    temp_mat = 1/(1+temp_mat)[1:9284,1:9284] #log2 transformed, 928x928 size
    temp = matrix(0,nrow=9284,ncol=2*LEN+1) #initialization
    for (i in 1:9284){
      for (j in -LEN:LEN){
        if ((i+j)%%9284!=0) {temp[i,j+LEN+1] = temp_mat[i,(i+j)%%9284]}
        else {temp[i,j+LEN+1] = temp_mat[i,9284]}
      }
    }
    return (temp)
  }
  
  make_df = function(temp_mat,peaks){
    temp_df = temp_mat %>% as.data.frame() %>%
      mutate(index = 1:928) %>% 
      mutate(RESPONSE = ifelse(index %in% peaks,"peak","nonpeak")) %>%
      #mutate(RESPONSE = peaks[1:928]) %>%
      transform(RESPONSE=as.factor(RESPONSE)) %>%
      select(-index)
    return (temp_df)
  }
  
  make_df_500 = function(temp_mat,peaks){
    temp_df = temp_mat %>% as.data.frame() %>%
      mutate(index = 1:9284) %>% 
      mutate(RESPONSE = ifelse(index %in% peaks,"peak","nonpeak")) %>%
      #mutate(RESPONSE = peaks[1:928]) %>%
      transform(RESPONSE=as.factor(RESPONSE)) %>%
      select(-index)
    return (temp_df)
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
  
  pr_theme = theme_bw(base_size = 12)+
    theme(
      #plot.background = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      #panel.border = element_blank(),
      #axis.line = element_line(color = 'black'),
      legend.text=element_text(size=12),
      legend.position=c(.8,.8)
    ) 
  
  roc_theme = theme_bw(base_size = 12)+
    theme(
      #plot.background = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      #panel.border = element_blank(),
      #axis.line = element_line(color = 'black'),
      legend.text=element_text(size=12),
      legend.position=c(.8,.2)
    ) 
  
  ## precision-recall curve and roc curve
  get.pr.curve = function(pred,testing){
    pred_prob=pred$peak
    fg <- pred_prob[testing$RESPONSE == 'peak']
    bg <- pred_prob[testing$RESPONSE == 'nonpeak']
    return (pr.curve(fg,bg,curve = T))
  }
  
  get.roc.curve = function(pred,testing){
    pred_prob=pred$peak
    fg <- pred_prob[testing$RESPONSE == 'peak']
    bg <- pred_prob[testing$RESPONSE == 'nonpeak']
    return (roc.curve(fg,bg,curve = T))
  }
}

##### loading #####
setwd('~/Documents/R_scripts/Ecoli/')
# folds_4_5k = readRDS('~/Documents/R_scripts/Ecoli/folds_4_5k.rds')
load('rf_result_new/testing_training_index.Rdata')
##### Inputting HiC matrices, res=5k, ICED #####
setwd('~/Documents/R_scripts/Ecoli/')
lioy_LB1 = toMat('hic_matrices/lioy_LB1.matrix')
lioy_LB2 = toMat('hic_matrices/lioy_LB2.matrix')
monica_LB1 = toMat('hic_matrices/monica_LB1.matrix')
# monica_LB2 = toMat('hic_matrices/monica_LB2.matrix')
##### Inputting reproducible peaks (no filter applied)
setwd('~/Documents/R_scripts/Ecoli/peaks')
gapr_df = fread('gapr_peaks.csv')
topo_df = fread('topo_peaks.csv')
gapr_peaks = gapr_df$`5k_mdpt_bin` %>% unique()  
topo_peaks = topo_df$`5k_mdpt_bin` %>% unique()  

##### creating data matrix #####
lioy1_gapr = make_nb_mat(lioy_LB1,464) %>% make_df(peaks=gapr_peaks)
lioy1_topo = make_nb_mat(lioy_LB1,464) %>% make_df(peaks=topo_peaks)
lioy1_gapr_half = make_nb_mat(lioy_LB1,232) %>% make_df(peaks=gapr_peaks)

lioy2_gapr = make_nb_mat(lioy_LB2,464) %>% make_df(peaks=gapr_peaks)
lioy2_topo = make_nb_mat(lioy_LB2,464) %>% make_df(peaks=topo_peaks)

monica1_gapr = make_nb_mat(monica_LB1,464) %>% make_df(peaks=gapr_peaks)
monica1_topo = make_nb_mat(monica_LB1,464) %>% make_df(peaks=topo_peaks)


##### gapr trained samples, 25% #####
setwd('~/Documents/R_scripts/Ecoli/rf_result_new/')
gapr_models_75 = list.files(pattern='gapr_lioy1_75')
gapr_models_50 = list.files(pattern='gapr_lioy1_50') #14
gapr_models_25 = list.files(pattern='^gapr_lioy1_25') #all
gapr_models_25_half = list.files(pattern='gapr_lioy1_half') #all

topo_models_75 = list.files(pattern='topo_lioy1_75')
topo_models_50 = list.files(pattern='topo_lioy1_50') #12
topo_models_25 = list.files(pattern='^topo_lioy1_25') #all

rand_gapr = list.files(pattern='rand_gapr')
rand_topo = list.files(pattern='rand_topo')

##### function loads all models #####
load_model = function(mod_list) {tmp = lapply(mod_list,load,.GlobalEnv)} ## loading the models}

##### function returns average AUC across all folds (1-4) #####
get_avg_auroc = function(mod_list) {
  if (grepl('topo',as.name(mod_list))) {dat = lioy1_topo}
  if (grepl('gapr',as.name(mod_list))) {dat = lioy1_gapr}
  if (grepl('half',as.name(mod_list))) {dat = lioy1_gapr_half}
  m_list = lapply(stringr::str_sub(mod_list,1,-7),function(x) return(get(x))) ## make a list
  res = lapply(1:length(m_list),function(test_id){
    current_fold = m_list[[test_id]]
    pred = mclapply(1:10,function(rf){
      model = lapply(1:20,function(m){
      if (is.null(current_fold[[rf]][[m]])==F){
      return(mutate(current_fold[[rf]][[m]]$results,ntree = m))}
    })
      results = do.call(rbind.data.frame,model) %>% arrange(desc(ROC))
      best_model = current_fold[[rf]][[results[1,'ntree']]]
      test = dat[testing_idx[,test_id],]
      pred = predict(best_model,newdata=test,type='prob')
      return(pred)
    },mc.cores=5)
    auc = lapply(pred,function(mod){
      x = get.roc.curve(mod,testing = dat[testing_idx[,test_id],])
      x$auc %>% round(4) 
    }) %>% unlist() %>% mean()
    return(auc)
  }) %>% unlist() %>% mean() %>% round(3)
  return(res)
}

get_avg_auprc = function(mod_list) {
  if (grepl('topo',as.name(mod_list))) {dat = lioy1_topo}
  if (grepl('gapr',as.name(mod_list))) {dat = lioy1_gapr}
  if (grepl('half',as.name(mod_list))) {dat = lioy1_gapr_half}
  m_list = lapply(stringr::str_sub(mod_list,1,-7),function(x) return(get(x))) ## make a list
  res = lapply(1:length(m_list),function(test_id){
    current_fold = m_list[[test_id]]
    pred = mclapply(1:10,function(rf){
      model = lapply(1:20,function(m){
        if(is.null(current_fold[[rf]][[m]])==F){
        return(mutate(current_fold[[rf]][[m]]$results,ntree = m))}
      })
      results = do.call(rbind.data.frame,model) %>% arrange(desc(ROC))
      best_model = current_fold[[rf]][[results[1,'ntree']]]
      test = dat[testing_idx[,test_id],]
      pred = predict(best_model,newdata=test,type='prob')
      return(pred)
    },mc.cores=5)
    auc = lapply(pred,function(mod){
      x = get.pr.curve(mod,testing = dat[testing_idx[,test_id],])
      x$auc.integral %>% round(4) 
    }) %>% unlist() %>% mean()
    return(auc)
  }) %>% unlist() %>% mean() %>% round(3)
  return(res)
}

##### load models and compute the average AUCs #####
set.seed(19980415)

load_model(gapr_models_25)
load_model(gapr_models_50)
load_model(rand_gapr)
load_model(gapr_models_75)
# load_model(gapr_models_25_half)
load_model(topo_models_25)
load_model(topo_models_50)
load_model(topo_models_75)
load_model(rand_topo)

get_avg_auroc(gapr_models_25)
get_avg_auroc(gapr_models_50)
get_avg_auroc(gapr_models_75)
# get_avg_auroc(gapr_models_25_half)
get_avg_auroc(topo_models_25)
get_avg_auroc(topo_models_50)
get_avg_auroc(topo_models_75)
# get_avg_auroc(gapr_models_75)
get_avg_auroc(rand_gapr)
get_avg_auroc(rand_topo)

get_avg_auprc(gapr_models_25)
get_avg_auprc(gapr_models_50)
get_avg_auprc(gapr_models_75)
# get_avg_auprc(gapr_models_25_half)
get_avg_auprc(topo_models_25)
get_avg_auprc(topo_models_50)
get_avg_auprc(topo_models_75)
# get_avg_auprc(gapr_models_75)
get_avg_auprc(rand_gapr)
get_avg_auprc(rand_topo)

##### feature importance  #####
get_mda = function(mod_list) {
  if (grepl('topo',as.name(mod_list))) {dat = lioy1_topo}
  if (grepl('gapr',as.name(mod_list))) {dat = lioy1_gapr}
  if (grepl('half',as.name(mod_list))) {dat = lioy1_gapr_half}
  m_list = lapply(stringr::str_sub(mod_list,1,-7),function(x) return(get(x))) ## make a list
  res = lapply(1:length(m_list),function(test_id){
    current_fold = m_list[[test_id]]
    pred = mclapply(1:10,function(rf){
      model = lapply(1:20,function(m){
        if(is.null(current_fold[[rf]][[m]])==F) {
        return(mutate(current_fold[[rf]][[m]]$results,ntree = m))}
      })
      results = do.call(rbind.data.frame,model) %>% arrange(desc(ROC))
      best_model = current_fold[[rf]][[results[1,'ntree']]]
      return(best_model$finalModel)
    },mc.cores=5)
    mda = lapply(pred,function(mod){importance(mod,type=1) %>% return()}) 
    return(do.call(cbind.data.frame,mda) %>% rowMeans())
  }) 
  return(do.call(cbind.data.frame,res) %>% rowMeans())
}
setwd('~/Documents/R_scripts/Ecoli/rf_plots/')

gapr_impt = data.frame(GapR_25 = get_mda(gapr_models_25),
                       GapR_50 = get_mda(gapr_models_50),
                       GapR_75 = get_mda(gapr_models_75),
           dis = -464:464 * 5) %>%
  melt(id.vars=c('dis'),value.name = 'MDA',variable.name='groups') %>% 
  ggplot(aes(x=dis,y=MDA,color=groups)) +
  geom_point(size=.8,alpha=.4) +
  geom_ma(aes(x=dis,y=MDA,color=groups),ma_fun = SMA, n = 50,lwd=1,linetype = "solid") +
  mytheme+
  labs(x = 'Loci linear separation distance (Kb)',y='Mean Decrease Accuracy',color='') +
  theme(legend.position = c(.85,.8), legend.background = element_rect(fill = alpha("white", .3)))
gapr_impt
ggsave('gapr_impt.pdf',gapr_impt,width=4.5,height=3)

topo_impt = data.frame(Topo_25 = get_mda(topo_models_25),
                       Topo_50 = get_mda(topo_models_50),
                       Topo_75 = get_mda(topo_models_75),
                       dis = -464:464 * 5) %>%
  melt(id.vars=c('dis'),value.name = 'MDA',variable.name='groups') %>% 
  ggplot(aes(x=dis,y=MDA,color=groups)) +
  geom_point(size=.8,alpha=.4) +
  geom_ma(aes(x=dis,y=MDA,color=groups),ma_fun = SMA, n = 30,lwd=1,linetype = "solid") +
  mytheme+
  labs(x = 'Loci linear separation distance (Kb)',y='Mean Decrease Accuracy',color='') +
  theme(legend.position = c(.85,.8), legend.background = element_rect(fill = alpha("white", .3)))
topo_impt
ggsave('topo_impt.pdf',topo_impt,width=4.5,height=3)

##### model performance AUC  #####
get_all_roc = function(mod_list) {
  if (grepl('topo',as.name(mod_list))) {dat = lioy1_topo}
  if (grepl('gapr',as.name(mod_list))) {dat = lioy1_gapr}
  if (grepl('half',as.name(mod_list))) {dat = lioy1_gapr_half}
  m_list = lapply(stringr::str_sub(mod_list,1,-7),function(x) return(get(x))) ## make a list
  res = lapply(1:length(m_list),function(test_id){
    current_fold = m_list[[test_id]]
    pred = mclapply(1:10,function(rf){
      model = lapply(1:20,function(m){
        if(is.null(current_fold[[rf]][[m]])==F) {
          return(mutate(current_fold[[rf]][[m]]$results,ntree = m))}
      })
      results = do.call(rbind.data.frame,model) %>% arrange(desc(ROC))
      best_model = current_fold[[rf]][[results[1,'ntree']]]
      test = dat[testing_idx[,test_id],]
      pred = predict(best_model,newdata=test,type='prob')
      return(pred)
    },mc.cores=5)
    roc_ls = lapply(pred,function(mod){
      x = get.roc.curve(mod,testing = dat[testing_idx[,test_id],])
      return(x)
    })  
    return(roc_ls)
  })
  return(do.call(c,res))
}

plot_all_auc = function(peak_type){
  options(warn=-1)
  if (peak_type=='gapr') {mod_list = list(gapr_models_25,gapr_models_50,gapr_models_75,rand_gapr)}
  if (peak_type=='topo') {mod_list = list(topo_models_25,topo_models_50,topo_models_75,rand_topo)}
  
  get_all_list_25 = get_all_roc(mod_list[[1]])
  df_25 = do.call(rbind,mclapply(1:length(get_all_list_25),function(i){
    tmp = data.frame(get_all_list_25[[i]]$curve[,c(1,2)])
    return(mutate(tmp,model=glue('25_{i}')))
  },mc.cores=5))
  auc_25 = sapply(get_all_list_25,function(x) x$auc)
  ci_25 = ci(auc_25)[1:3] %>% round(3)
  mid_auc_25 = order(auc_25)[20]
  mid_auc_25 = filter(df_25,model==glue('25_{mid_auc_25}'))
  
  get_all_list_50 = get_all_roc(mod_list[[2]])
  df_50 = do.call(rbind,lapply(1:length(get_all_list_50),function(i){
    tmp = data.frame(get_all_list_50[[i]]$curve[,c(1,2)])
    return(mutate(tmp,model=glue('50_{i}')))
  }))
  auc_50 = sapply(get_all_list_50,function(x) x$auc)
  ci_50 = ci(auc_50)[1:3] %>% round(3)
  mid_auc_50 = order(auc_50)[20]
  mid_auc_50 = filter(df_50,model==glue('50_{mid_auc_50}'))
  
  get_all_list_75 = get_all_roc(mod_list[[3]])
  df_75 = do.call(rbind,mclapply(1:length(get_all_list_75),function(i){
    tmp = data.frame(get_all_list_75[[i]]$curve[,c(1,2)])
    return(mutate(tmp,model=glue('75_{i}')))
  },mc.cores=5))
  auc_75 = sapply(get_all_list_75,function(x) x$auc)
  ci_75 = ci(auc_75)[1:3] %>% round(3)
  mid_auc_75 = order(auc_75)[20]
  mid_auc_75 = filter(df_75,model==glue('75_{mid_auc_75}'))
  
  get_all_list_rand = get_all_roc(mod_list[[4]])
  df_rand = do.call(rbind,lapply(1:length(get_all_list_rand),function(i){
    tmp = data.frame(get_all_list_rand[[i]]$curve[,c(1,2)])
    return(mutate(tmp,model=glue('rand_{i}')))
  }))
  auc_rand = sapply(get_all_list_rand,function(x) x$auc)
  ci_rand = ci(auc_rand)[1:3] %>% round(3)
  mid_auc_rand = order(auc_rand)[20]
  mid_auc_rand = filter(df_rand,model==glue('rand_{mid_auc_rand}'))
  df = mutate(rbind(df_25,df_50,df_75,df_rand),perc=as.factor(c(rep(paste0(peak_type,"_25"),dim(df_25)[1]),
                                                rep(paste0(peak_type,"_50"),dim(df_50)[1]),
                                                rep(paste0(peak_type,"_75"),dim(df_75)[1]),
                                                rep('baseline',dim(df_rand)[1]))))
  # levels(df$perc) = c("gapr_25" ,"gapr_50", "baseline")
  return(list(df,mid_auc_25,mid_auc_50,mid_auc_75,mid_auc_rand))
}

setwd('~/Documents/R_scripts/Ecoli/rf_plots/')
plot_info_gapr = plot_all_auc('gapr')
gapr_roc_plot = ggplot() + 
  geom_line(data=plot_info_gapr[[1]],aes(X1,X2,color=perc,group=model),
            alpha=ifelse(plot_info_gapr[[1]]$perc=='baseline',.1,.3),lwd=.15) +
  geom_line(data=plot_info_gapr[[2]],aes(X1,X2),alpha=.7,lwd=1.5,color='#F8766D') +
  geom_line(data=plot_info_gapr[[3]],aes(X1,X2),alpha=.7,lwd=1.5,color='#00BA38') +
  geom_line(data=plot_info_gapr[[4]],aes(X1,X2),alpha=.7,lwd=1.5,color='#00BFC4') +
  geom_line(data=plot_info_gapr[[5]],aes(X1,X2),alpha=.7,lwd=1.5,color='#444444') +
  scale_color_manual(breaks = unique(plot_info_gapr[[1]]$perc),
                     values = c('#F8766D','#00BA38','#00BFC4','#444444'),
                     labels = c('GapR, 25% (0.74)','GapR, 50% (0.74)','GapR, 75% (0.74)','Baseline (0.43)'),
                     guide=guide_legend(override.aes = list(size = c(1.5,1.5,1.5,1.5),alpha=c(2,2,2,2)))) +
  geom_abline(slope = 1,intercept = 0,color='black',linetype='dashed',lwd=.5)+ ## replaced with baseline model
  roc_theme +
  labs(x='1 - Specificity',y='Sensitivity',color='')+
  theme(legend.position = c(.78,.2),legend.key.width = unit(25,"pt"),
        legend.background = element_rect(fill = alpha("white", .3)),
        legend.text = element_text(size=10)) 
gapr_roc_plot
ggsave('gapr_roc_plot.pdf',gapr_roc_plot,width=4.5,height=4)

plot_info_topo = plot_all_auc('topo')
topo_roc_plot = ggplot() + 
  geom_line(data=plot_info_topo[[1]],aes(X1,X2,color=perc,group=model),
            alpha=ifelse(plot_info_topo[[1]]$perc=='baseline',.1,.3),lwd=.15) +
  geom_line(data=plot_info_topo[[2]],aes(X1,X2),alpha=.7,lwd=1.5,color='#F8766D') +
  geom_line(data=plot_info_topo[[3]],aes(X1,X2),alpha=.7,lwd=1.5,color='#00BA38') +
  geom_line(data=plot_info_topo[[4]],aes(X1,X2),alpha=.7,lwd=1.5,color='#00BFC4') +
  geom_line(data=plot_info_topo[[5]],aes(X1,X2),alpha=.7,lwd=1.5,color='#444444') +
  scale_color_manual(breaks = unique(plot_info_topo[[1]]$perc),
                     values = c('#F8766D','#00BA38','#00BFC4','#444444'),
                     labels = c('Topo, 25% (0.54)','Topo, 50% (0.56)','Topo, 75% (0.58)','Baseline (0.48)'),
                     guide=guide_legend(override.aes = list(size = c(1.5,1.5,1.5,1.5),alpha=c(2,2,2,2)))) +
  geom_abline(slope = 1,intercept = 0,color='black',linetype='dashed',lwd=.5)+ ## replaced with baseline model
  roc_theme +
  labs(x='1 - Specificity',y='Sensitivity',color='')+
  theme(legend.position = c(.78,.2),legend.key.width = unit(25,"pt"),
        legend.background = element_rect(fill = alpha("white", .3)),
        legend.text = element_text(size=10))
topo_roc_plot
ggsave('topo_roc_plot.pdf',topo_roc_plot,width=4.5,height=4)


##### model performance PRC  #####
get_all_prc = function(mod_list) {
  if (grepl('topo',as.name(mod_list))) {dat = lioy1_topo}
  if (grepl('gapr',as.name(mod_list))) {dat = lioy1_gapr}
  if (grepl('half',as.name(mod_list))) {dat = lioy1_gapr_half}
  m_list = lapply(stringr::str_sub(mod_list,1,-7),function(x) return(get(x))) ## make a list
  res = lapply(1:length(m_list),function(test_id){
    current_fold = m_list[[test_id]]
    pred = mclapply(1:10,function(rf){
      model = lapply(1:20,function(m){
        if(is.null(current_fold[[rf]][[m]])==F) {
          return(mutate(current_fold[[rf]][[m]]$results,ntree = m))}
      })
      results = do.call(rbind.data.frame,model) %>% arrange(desc(ROC))
      best_model = current_fold[[rf]][[results[1,'ntree']]]
      test = dat[testing_idx[,test_id],]
      pred = predict(best_model,newdata=test,type='prob')
      return(pred)
    },mc.cores=5)
    prc_ls = lapply(pred,function(mod){
      x = get.pr.curve(mod,testing = dat[testing_idx[,test_id],])
      return(x)
    })  
    return(prc_ls)
  })
  return(do.call(c,res))
}

plot_all_prc = function(peak_type){
  options(warn=-1)
  if (peak_type=='gapr') {mod_list = list(gapr_models_25,gapr_models_50,gapr_models_75,rand_gapr)}
  if (peak_type=='topo') {mod_list = list(topo_models_25,topo_models_50,topo_models_75,rand_topo)}
  
  get_all_list_25 = get_all_prc(mod_list[[1]])
  df_25 = do.call(rbind,mclapply(1:length(get_all_list_25),function(i){
    tmp = data.frame(get_all_list_25[[i]]$curve[,c(1,2)])
    return(mutate(tmp,model=glue('25_{i}')))
  },mc.cores=5))
  auc_25 = sapply(get_all_list_25,function(x) x$auc.integral)
  ci_25 = ci(auc_25)[1:3] %>% round(3)
  mid_auc_25 = order(auc_25)[20]
  mid_auc_25 = filter(df_25,model==glue('25_{mid_auc_25}'))
  
  get_all_list_50 = get_all_prc(mod_list[[2]])
  df_50 = do.call(rbind,lapply(1:length(get_all_list_50),function(i){
    tmp = data.frame(get_all_list_50[[i]]$curve[,c(1,2)])
    return(mutate(tmp,model=glue('50_{i}')))
  }))
  auc_50 = sapply(get_all_list_50,function(x) x$auc.integral)
  ci_50 = ci(auc_50)[1:3] %>% round(3)
  mid_auc_50 = order(auc_50)[20]
  mid_auc_50 = filter(df_50,model==glue('50_{mid_auc_50}'))
  
  get_all_list_75 = get_all_prc(mod_list[[3]])
  df_75 = do.call(rbind,lapply(1:length(get_all_list_75),function(i){
    tmp = data.frame(get_all_list_75[[i]]$curve[,c(1,2)])
    return(mutate(tmp,model=glue('75_{i}')))
  }))
  auc_75 = sapply(get_all_list_75,function(x) x$auc.integral)
  ci_75 = ci(auc_75)[1:3] %>% round(3)
  mid_auc_75 = order(auc_75)[20]
  mid_auc_75 = filter(df_75,model==glue('75_{mid_auc_75}'))
  
  get_all_list_rand = get_all_prc(mod_list[[4]])
  df_rand = do.call(rbind,lapply(1:length(get_all_list_rand),function(i){
    tmp = data.frame(get_all_list_rand[[i]]$curve[,c(1,2)])
    return(mutate(tmp,model=glue('rand_{i}')))
  }))
  auc_rand = sapply(get_all_list_rand,function(x) x$auc.integral)
  ci_rand = ci(auc_rand)[1:3] %>% round(3)
  mid_auc_rand = order(auc_rand)[20]
  mid_auc_rand = filter(df_rand,model==glue('rand_{mid_auc_rand}'))
  df = mutate(rbind(df_25,df_50,df_75,df_rand),perc=as.factor(c(rep(paste0(peak_type,"_25"),dim(df_25)[1]),
                                                          rep(paste0(peak_type,"_50"),dim(df_50)[1]),
                                                          rep(paste0(peak_type,"_75"),dim(df_75)[1]),
                                                          rep('baseline',dim(df_rand)[1]))))
  # levels(df$perc) = c("gapr_25" ,"gapr_50", "baseline")
  gc()
  return(list(df,mid_auc_25,mid_auc_50,mid_auc_75,mid_auc_rand))
}

plot_info_gapr = plot_all_prc('gapr')
gapr_prc_plot = ggplot() + 
  geom_line(data=plot_info_gapr[[1]],aes(X1,X2,color=perc,group=model),
            alpha=ifelse(plot_info_gapr[[1]]$perc=='baseline',.1,.3),lwd=.15) +
  geom_line(data=plot_info_gapr[[2]],aes(X1,X2),alpha=.7,lwd=1.5,color='#F8766D') +
  geom_line(data=plot_info_gapr[[3]],aes(X1,X2),alpha=.7,lwd=1.5,color='#00BA38') +
  geom_line(data=plot_info_gapr[[4]],aes(X1,X2),alpha=.7,lwd=1.5,color='#00BFC4') +
  geom_line(data=plot_info_gapr[[5]],aes(X1,X2),alpha=.7,lwd=1.5,color='#444444') +
  scale_color_manual(breaks = unique(plot_info_gapr[[1]]$perc),
                     values = c('#F8766D','#00BA38','#00BFC4','#444444'),
                     labels = c('GapR, 25% (0.68)','GapR, 50% (0.68)','GapR, 75% (0.68)','Baseline (0.44)'),
                     guide=guide_legend(override.aes = list(size = c(1.5,1.5,1.5,1.5),alpha=c(2,2,2,2)))) +
  pr_theme +
  labs(x='Recall',y='Precision',color='')+
  theme(legend.position = c(.78,.2),legend.key.width = unit(25,"pt"),
        legend.background = element_rect(fill = alpha("white", .3)),
        legend.text = element_text(size=10)) 
gapr_prc_plot
setwd('~/Documents/R_scripts/Ecoli/rf_plots/')
ggsave('gapr_prc_plot.pdf',gapr_prc_plot,width=4.5,height=4)

plot_info_topo = plot_all_prc('topo')
topo_prc_plot = ggplot() + 
  geom_line(data=plot_info_topo[[1]],aes(X1,X2,color=perc,group=model),
            alpha=ifelse(plot_info_topo[[1]]$perc=='baseline',.1,.3),lwd=.15) +
  geom_line(data=plot_info_topo[[2]],aes(X1,X2),alpha=.7,lwd=1.5,color='#F8766D') +
  geom_line(data=plot_info_topo[[3]],aes(X1,X2),alpha=.7,lwd=1.5,color='#00BA38') +
  geom_line(data=plot_info_topo[[4]],aes(X1,X2),alpha=.7,lwd=1.5,color='#00BFC4') +
  geom_line(data=plot_info_topo[[5]],aes(X1,X2),alpha=.7,lwd=1.5,color='#444444') +
  scale_color_manual(breaks = unique(plot_info_topo[[1]]$perc),
                     values = c('#F8766D','#00BA38','#00BFC4','#444444'),
                     labels = c('Topo, 25% (0.41)','Topo, 50% (0.43)','Topo, 75% (0.44)','Baseline (0.44)'),
                     guide=guide_legend(override.aes = list(size = c(1.5,1.5,1.5,1.5),alpha=c(2,2,2,2)))) +
  pr_theme +
  labs(x='Recall',y='Precision',color='')+
  theme(legend.position = c(.78,.2),legend.key.width = unit(25,"pt"),
        legend.background = element_rect(fill = alpha("white", .1)),
        legend.text = element_text(size=10)) 
topo_prc_plot
ggsave('topo_prc_plot.pdf',topo_prc_plot,width=4.5,height=4)

##### model performance testing on other datasets #####
get_cross_auroc = function(mod_list,test_dat) {
  m_list = lapply(stringr::str_sub(mod_list,1,-7),function(x) return(get(x))) ## make a list
  res = lapply(1:length(m_list),function(test_id){
    current_fold = m_list[[test_id]]
    pred = mclapply(1:10,function(rf){
      model = lapply(1:20,function(m){
        if (is.null(current_fold[[rf]][[m]])==F){
          return(mutate(current_fold[[rf]][[m]]$results,ntree = m))}
      })
      results = do.call(rbind.data.frame,model) %>% arrange(desc(ROC))
      best_model = current_fold[[rf]][[results[1,'ntree']]]
      
      pred = predict(best_model,newdata=test_dat,type='prob')
      return(pred)
    },mc.cores=5)
    auc = lapply(pred,function(mod){
      x = get.roc.curve(mod,testing = test_dat[testing_idx[,test_id],])
      x$auc %>% round(4) 
    }) %>% unlist() %>% mean()
    return(auc)
  }) %>% unlist() %>% mean() %>% round(3)
  return(res)
}
get_cross_auroc(gapr_models_25,monica1_gapr)


get_cross_auprc = function(mod_list) {
  if (grepl('topo',as.name(mod_list))) {dat = lioy1_topo}
  if (grepl('gapr',as.name(mod_list))) {dat = lioy1_gapr}
  if (grepl('half',as.name(mod_list))) {dat = lioy1_gapr_half}
  m_list = lapply(stringr::str_sub(mod_list,1,-7),function(x) return(get(x))) ## make a list
  res = lapply(1:length(m_list),function(test_id){
    current_fold = m_list[[test_id]]
    pred = mclapply(1:10,function(rf){
      model = lapply(1:20,function(m){
        if(is.null(current_fold[[rf]][[m]])==F){
          return(mutate(current_fold[[rf]][[m]]$results,ntree = m))}
      })
      results = do.call(rbind.data.frame,model) %>% arrange(desc(ROC))
      best_model = current_fold[[rf]][[results[1,'ntree']]]
      test = dat[testing_idx[,test_id],]
      pred = predict(best_model,newdata=test,type='prob')
      return(pred)
    },mc.cores=5)
    auc = lapply(pred,function(mod){
      x = get.pr.curve(mod,testing = dat[testing_idx[,test_id],])
      x$auc.integral %>% round(4) 
    }) %>% unlist() %>% mean()
    return(auc)
  }) %>% unlist() %>% mean() %>% round(3)
  return(res)
}





