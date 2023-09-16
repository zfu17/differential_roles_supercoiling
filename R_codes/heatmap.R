## This script plots heatmaps for all HiC datasets
## Author: Ziqi Fu
## zfu17@jhu.edu
## Dec 13, 2022

library(ComplexHeatmap)
library(dplyr)
library(draw)
library(data.table)
library(parallel)
library(circlize)
{
  toMat = function(f){
    temp = fread(f)%>%as.data.frame()
    mat = matrix(NA,nrow = max(temp[,1:2]),ncol=max(temp[,1:2]))
    for (i in seq_along(temp[,1])){
      # mat[temp[i,2],temp[i,1]] = NA
      mat[temp[i,1],temp[i,2]] = log2(1+temp[i,3])
    }
    mat = mat*750000/sum(mat,na.rm=T) #normalization
    return(mat)
  } 
}

##### Inputting HiC matrices ICED #####
setwd('/Users/ziqi_fu/Documents/R_scripts/Ecoli/hic_matrices/')
filenames = c('lioy_LB1.matrix','lioy_LB2.matrix',
              'monica_LB1.matrix','monica_LB2.matrix',
               'CC_HpaII_rep1_5000.matrix','CC_HpaII_rep2_5000.matrix')
              # 'CC_HpaII_rep1_500.matrix','CC_HpaII_rep2_500.matrix',
              # 'CC_MluCI_500.matrix','CC_HpaII_MluCI_500.matrix')
# filenames_hires = c('CC_HpaII_rep1_5000.matrix','CC_HpaII_rep2_5000.matrix')
# hic_matDF_list =  list()
hic_matDF_list = mclapply(filenames,function(x) {toMat(x)},mc.cores=2) # toMat returns log10(1+count)
# hic_matDF_list_hires = mclapply(filenames_hires,function(x) {toMat(x)},mc.cores=2)
names(hic_matDF_list) = filenames
# names(hic_matDF_list_hires) = filenames_hires


##### generate heatmap #####
hm <- function(x,col_fun) {
  return (Heatmap(x,
                  name = ' ',
                  col = col_fun,
                  na_col = "white",
                  show_column_names = FALSE,
                  show_row_names = FALSE,
                  width = unit(8, "cm"), 
                  height = unit(8, "cm"),
                  cluster_columns=FALSE,
                  cluster_rows = FALSE
  ))
}
colorRamp2(c(.5,5), c("aquamarine","red")) ##(0,1,5) for Lioy

hm_lioy1 = hm(hic_matDF_list[[1]],colorRamp2(c(.5,5), c("aquamarine","red")))
hm_lioy2 = hm(hic_matDF_list[[2]],colorRamp2(c(2,6), c("aquamarine","red")))
hm_guo1 = hm(hic_matDF_list[[3]],colorRamp2(c(1.5,2,5), c("aquamarine",'#d17f7f',"red")))
hm_guo2 = hm(hic_matDF_list[[4]],colorRamp2(c(2,2,7), c("aquamarine",'#eccccc',"#ad1919")))


hm_list = list(hm_lioy1,hm_lioy2,hm_guo1,hm_guo2)
names(hm_list) = names(hic_matDF_list)[1:4]

setwd('~/Documents/R_scripts/Ecoli/figures/heatmaps/')

mclapply(names(hm_list)[1:4],function(x){
  pdf(file=paste0(x,'_freq.pdf'),width = 5,height=5)
  draw(hm_list[[x]])
  dev.off()
},mc.cores=4)

hm_cc1 = hm(hic_matDF_list[[5]],colorRamp2(c(.9,3.6), c("aquamarine","red")))
pdf(file='cc1_freq.pdf',width = 5,height=5)
draw(hm_cc1)
dev.off()
hm_cc2 = hm(hic_matDF_list[[6]],colorRamp2(c(1.1,3.6), c("aquamarine","red")))
pdf(file='cc2_freq.pdf',width = 5,height=5)
draw(hm_cc2)
dev.off()

##### generate heatmap high_res #####
col_fun = colorRamp2(c(0,5, 10), c("white","aquamarine","red"))
hm <- function(x) {
  return (Heatmap(x,
                  name = ' ',
                  col = col_fun,
                  na_col = "white",
                  show_column_names = FALSE,
                  show_row_names = FALSE,
                  width = unit(15, "cm"), 
                  height = unit(15, "cm"),
                  cluster_columns=FALSE,
                  cluster_rows = FALSE
  ))
}

hm_list_hires = mclapply(names(hic_matDF_list_hires), function(x) {hm(hic_matDF_list_hires[[x]])},mc.cores=4)
names(hm_list_hires) = names(hic_matDF_list_hires)

setwd('~/Documents/R_scripts/Ecoli/figures/heatmaps/')

lapply(names(hm_list_hires),function(x){
  draw(hm_list_hires[[x]])
  drawExport(paste0(x,'_freq.pdf'))
})



