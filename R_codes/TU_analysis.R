## This script is for TU analysis
## Author: Ziqi Fu
## zfu17@jhu.edu
## check convergence https://biocyc.org/ECOLI/NEW-IMAGE?type=LOCUS-POSITION&object=NIL&orgids=ECOLI&chromosome=COLI-K12&bp-range=1/50000

##### packages #####
library(ggplot2)
library(tidyquant)
library(tidyverse)
library(data.table)
library(reshape2)
library(scales)
library(glue)
library(draw)

##### functions do-not-modify #####
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
mytheme = theme_bw(base_size = 12)+
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

##### Inputting reproducible peaks (filter applied)
# see peaks processing R scripts
setwd('/Users/ziqi_fu/Documents/R_scripts/Ecoli/peaks')
gapr_df = fread('gapr_peaks.csv')  
topo_df = fread('topo_peaks.csv')  
## we filter out topo peaks that are 100% contained in a transcription unit
topo_df = filter(topo_df,tu_perc!=1) #82 removed
gapr_df = mutate(gapr_df,midpoint = ceiling((gapr_df$start+gapr_df$end)/2))
topo_df = mutate(topo_df,midpoint = ceiling((topo_df$start+topo_df$end)/2))

##### loading TU #####
tu_raw = fread('/Users/ziqi_fu/Documents/R_scripts/Ecoli/TU_annotation.csv')
tu_raw = mutate(tu_raw,signed.end=ifelse(strand==0,end,start),
                signed.start=ifelse(strand==0,start,end),width=abs(start-end)+1)

##### TU statistics #####
tu_length_plot = tu_raw %>% ggplot() +
  geom_density(aes(x=log2(width)),fill='#66A182',alpha=.7) +
  mytheme +
  labs(x='Log2 normalized TU width')
# ggsave('~/Documents/R_scripts/Ecoli/figures/tu_length_plot.pdf',tu_length_plot,width=4.5,height=3)

tu_rpk_plot=tu_raw %>% ggplot() +
  geom_density(aes(x=log2(1+newRPK)),fill='#CC79A7',alpha=.7) +
  mytheme +
  labs(x='Log2 normalized transcription activity (RPK)')
# ggsave('~/Documents/R_scripts/Ecoli/figures/tu_rpk_plot.pdf',tu_rpk_plot,width=4.5,height=3)

summary(tu_raw$width)
summary(tu_raw$newRPK)

##### TU based analysis processing for GapR #####
gapr_df = mutate(gapr_df,
                 closest_TU_end = sapply(gapr_df$midpoint,function(x) which.min(abs(x-tu_raw$signed.end))),
                 dis_to_end = sapply(gapr_df$midpoint,function(x) min(abs(x-tu_raw$signed.end))))
gapr_df = mutate(gapr_df,act_level=tu_raw[closest_TU_end,]$newRPK)    

##### TU based processing for Topo #####
topo_df = mutate(topo_df,
                 closest_TU_start = sapply(topo_df$midpoint,function(x) which.min(abs(x-tu_raw$signed.start))),
                 dis_to_start = sapply(topo_df$midpoint,function(x) min(abs(x-tu_raw$signed.start))))
topo_df = mutate(topo_df,act_level=tu_raw[closest_TU_start,]$newRPK)      

##### 2D scatter plot #####
quantile(tu_raw$newRPK,0.75) #29.34

gapr_tes_plot = ggplot(gapr_df) +
  geom_hex(aes(x=log(1+act_level),y=(dis_to_end))) +
  geom_vline(xintercept = log2(1+29.34),color='orange',lwd=1.5)+ ##29.34 is the 75% percentile of newRPK
  geom_hline(yintercept = discut,color='orange',lwd=1.5) +
  theme_bw(base_size = 10) +
  labs(x='TU Activity level in log2(1+RPK)',
       y='GapR peak: Distance to closest TES') +
  annotate('rect',xmin=-.1,xmax=log2(1+29.34),ymin=discut,ymax=5200,
           fill='purple',alpha=.2)+
  theme(legend.position = c(.9,.8),legend.key.size = unit(10,'pt'),
        legend.background=element_rect(fill = alpha("white", 0.4))) +
  scale_y_continuous(breaks=seq(0, 5200,1000)) +
  scale_x_continuous(breaks=seq(0, 10,2))
gapr_tes_plot
ggsave('~/Documents/R_scripts/Ecoli/TU_analysis_plots/gapr_tes_plot.pdf',gapr_tes_plot,height=3,width=4.5)

topo_tss_plot = ggplot(topo_df) +
  geom_hex(aes(x=log2(1+act_level),y=dis_to_start)) +
  geom_vline(xintercept = log2(1+29.34),color='orange',lwd=1.5)+
  geom_hline(yintercept = discut,color='orange',lwd=1.5) +
  theme_bw(base_size = 10) +
  labs(x='TU Activity level in log2(1+RPK)',
       y='Topo peak: Distance to closest TSS') +
  annotate('rect',xmin=-.1,xmax=log2(1+29.34),ymin=discut,ymax=8000,
           fill='purple',alpha=.2)+
  theme(legend.position = c(.9,.8),legend.key.size = unit(10,'pt'),
        legend.background=element_rect(fill = alpha("white", 0.4))) +
  scale_y_continuous(breaks=seq(0, 8000,1000)) +
  scale_x_continuous(breaks=seq(0, 18,3))
topo_tss_plot
ggsave('~/Documents/R_scripts/Ecoli/TU_analysis_plots/topo_tss_plot.pdf',topo_tss_plot,height=3,width=4.5)

##### what are the distant genes #####
## gapr
gapr_df %>% filter(dis_to_end>=discut) %>% View()


##### control: random midpoint analysis ##### 
set.seed(19980415)

random_df1 = data.frame(midpoint = sample(1:4641652,dim(gapr_df)[1],replace = F)) # 4641652 is the length of the E.coli chromosome

random_df1 = mutate(random_df1,
                 closest_TU_end = sapply(midpoint,function(x) which.min(abs(x-tu_raw$signed.end))),
                 dis_to_end = sapply(midpoint,function(x) min(abs(x-tu_raw$signed.end))),
                 closest_TU_start = sapply(midpoint,function(x) which.min(abs(x-tu_raw$signed.start))),
                 dis_to_start = sapply(midpoint,function(x) min(abs(x-tu_raw$signed.start))))
random_df1 = mutate(random_df1,act_level_end=tu_raw[closest_TU_end,]$newRPK,
                   act_level_start=tu_raw[closest_TU_start,]$newRPK)    

random_df2 = data.frame(midpoint = sample(1:4641652,dim(topo_df)[1],replace = F)) # 4641652 is the length of the E.coli chromosome

random_df2 = mutate(random_df2,
                    closest_TU_end = sapply(midpoint,function(x) which.min(abs(x-tu_raw$signed.end))),
                    dis_to_end = sapply(midpoint,function(x) min(abs(x-tu_raw$signed.end))),
                    closest_TU_start = sapply(midpoint,function(x) which.min(abs(x-tu_raw$signed.start))),
                    dis_to_start = sapply(midpoint,function(x) min(abs(x-tu_raw$signed.start))))
random_df2 = mutate(random_df2,act_level_end=tu_raw[closest_TU_end,]$newRPK,
                    act_level_start=tu_raw[closest_TU_start,]$newRPK) 


##### plot the control #####
discut =  2500
gapr_tes_plot_random = ggplot(random_df1) +
  geom_hex(aes(x=log(1+act_level_end),y=(dis_to_end))) +
  geom_vline(xintercept = log2(1+29.34),color='orange',lwd=1.5)+ ##29.34 is the 75% percentile of newRPK
  geom_hline(yintercept = discut,color='orange',lwd=1.5) +
  theme_bw(base_size = 10) +
  labs(x='TU Activity level in log2(1+RPK)',
       y='Ramdom peak: Distance to closest TES') +
  annotate('rect',xmin=-.1,xmax=log2(1+29.34),ymin=discut,ymax=7500,
           fill='purple',alpha=.2)+
  theme(legend.position = c(.9,.8),legend.key.size = unit(10,'pt'),
        legend.background=element_rect(fill = alpha("white", 0.4))) +
  scale_y_continuous(breaks=seq(0, 7000,1000)) +
  scale_x_continuous(breaks=seq(0, 10,2))
gapr_tes_plot_random
ggsave('~/Documents/R_scripts/Ecoli/TU_analysis_plots/gapr_tes_plot_random.pdf',gapr_tes_plot_random,height=3,width=4.5)

topo_tss_plot_random = ggplot(random_df2) +
  geom_hex(aes(x=log2(1+act_level_start),y=dis_to_start)) +
  geom_vline(xintercept = log2(1+29.34),color='orange',lwd=1.5)+
  geom_hline(yintercept = discut,color='orange',lwd=1.5) +
  theme_bw(base_size = 10) +
  labs(x='TU Activity level in log2(1+RPK)',
       y='Ramdom peak: Distance to closest TSS') +
  annotate('rect',xmin=-.1,xmax=log2(1+29.34),ymin=discut,ymax=7500,
           fill='purple',alpha=.2)+
  theme(legend.position = c(.9,.8),legend.key.size = unit(10,'pt'),
        legend.background=element_rect(fill = alpha("white", 0.4))) +
  # coord_cartesian(ylim=c(0,5500))+
  scale_y_continuous(breaks=seq(0, 7000,1000)) +
  scale_x_continuous(breaks=seq(0, 18,3))
topo_tss_plot_random
ggsave('~/Documents/R_scripts/Ecoli/TU_analysis_plots/topo_tss_plot_random.pdf',topo_tss_plot_random,height=3,width=4.5)


##### distance to closest TUs #####
gapr_control = data.frame(dis=c(gapr_df$dis_to_end,
                 # topo_df$dis_to_start,
                 random_df1$dis_to_end),
                 # random_df$dis_to_start),
           type=c(rep('GapR',length(gapr_df$dis_to_end)),
                  # rep('Topo',length(topo_df$dis_to_start)),
                  rep('Control',length(random_df1$dis_to_end)))) %>%
                  # rep('rand.topo',length(random_df$dis_to_start)))) %>% 
  ggplot(aes(x=log2(dis),fill=type)) +
  geom_density(alpha=.3,position = "identity")+
  mytheme+
  labs(fill='',x='Distance to closest TES, log2',y='density')+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 8))+
  theme(legend.position =c(.2,.8))+
  theme(legend.background=element_rect(fill = alpha("white", .3)))+
  annotate(geom='text',x=2,y=.05,
           label=glue('p={wilcox.test((gapr_df$dis_to_end),random_df1$dis_to_end)$p.value%>%scientific(3)}'))
gapr_control

topo_control = data.frame(dis=c(topo_df$dis_to_start,
                 random_df2$dis_to_start),
           type=c(rep('Topo',length(topo_df$dis_to_start)),
                  rep('Control',length(random_df2$dis_to_start)))) %>%
  # rep('rand.topo',length(random_df$dis_to_start)))) %>% 
  ggplot(aes(x=log2(dis),fill=type)) +
  geom_density(alpha=.3,position = "identity")+
  mytheme+
  labs(fill='',x='Distance to closest TSS, log2',y='density')+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 8))+
  theme(legend.position =c(.2,.8))+
  theme(legend.background=element_rect(fill = alpha("white", .3))) +
  scale_fill_manual(values=c(3,6))+
  annotate(geom='text',x=2,y=.05,
           label=glue('p={wilcox.test(topo_df$dis_to_start,random_df2$dis_to_start)$p.value%>%round(3)}'))
topo_control

wilcox.test((gapr_df$dis_to_end),(random_df1$dis_to_end))
wilcox.test((topo_df$dis_to_start),(random_df2$dis_to_start))
wilcox.test((random_df1$dis_to_end),(random_df2$dis_to_start))

ggsave('~/Documents/R_scripts/Ecoli/TU_analysis_plots/gapr_control.pdf',gapr_control,width=4.5,height=3)
ggsave('~/Documents/R_scripts/Ecoli/TU_analysis_plots/topo_control.pdf',topo_control,width=4.5,height=3)

summary(gapr_df$dis_to_end)
(dim(filter(gapr_df,dis_to_end>2500))[1] / dim(gapr_df)[1]*100) %>% round(3)
(dim(filter(topo_df,dis_to_start>2500))[1] / dim(topo_df)[1]*100) %>% round(3)
(length(random_df1$dis_to_end[random_df1$dis_to_end > 2500]) / dim(random_df1)[1]*100) %>% round(3)
(length(random_df2$dis_to_start[random_df2$dis_to_start > 2500]) / dim(random_df2)[1]*100) %>% round(3)

##### plot the distant sites #####
setwd('~/Documents/R_scripts/Ecoli/TU_analysis_plots')

peaks_dt = filter(random_df1,dis_to_end>discut& act_level_end<29.34)
peaks_dt = filter(random_df2,dis_to_end>discut& act_level_start<29.34)
peaks_dt = filter(gapr_df,dis_to_end>discut& act_level<29.34)
peaks_dt = tu_raw
{
  drawSettings(pageWidth = 10, pageHeight = 2, units = "cm")
  drawPage()
  drawBox(x = 0, y = 0, width = 15, height = .6,units = "cm")
  for (i in 1:length(peaks_dt[[1]])){
    # pos = -7.5+peaks_dt$midpoint[i]/4641652*15
    pos = -7.5+peaks_dt$signed.end[i]/4641652*15 #TSS TES
    drawLine(x = c(pos,pos), y = c(-.3,.3),lineColor = "green",lineWidth = .005)
  }
  #setwd("~/Desktop")
  drawExport("tes_all.pdf")
}


##### updated analysis #####
setwd('~/Documents/R_scripts/Ecoli/TU_analysis_plots')

### gapr
gapr_tu_df = 
  do.call(rbind.data.frame,
                     mclapply(1:dim(gapr_df)[1],function(g){
  peak = gapr_df[g,]
  front = filter(tu_raw,start-1000<peak$end & strand == 0)
  back = filter(tu_raw,end+1000>peak$start & strand == 1)
  all = rbind(front,back)
  all = mutate(all,dis = abs(peak$midpoint-signed.end))
  all = arrange(all,dis)
  # return(all)
  return(all[1,])
},mc.cores=4))

gapr_df = cbind.data.frame(gapr_df,gene=gapr_tu_df$`Unnamed: 0.1`,
                           tes=gapr_tu_df$signed.end,
                           activity=gapr_tu_df$newRPK) 

gapr_df = gapr_df %>% mutate(dis_to_tu = abs(midpoint-tes))

discut =  1000
gapr_plot = ggplot(gapr_df) +
  geom_hex(aes(x=log2(1+activity),y=(dis_to_tu))) +
  geom_vline(xintercept = log2(1+29.34),color='orange',lwd=1.5)+ ##29.34 is the 75% percentile of newRPK
  geom_hline(yintercept = discut,color='orange',lwd=1.5) +
  theme_bw(base_size = 16) +
  labs(x='TU Activity level in log2(1+RPK)',
       y='Distance to closest TES') +
  annotate('rect',xmin=-.1,xmax=log2(1+29.34),ymin=discut,ymax=8000,
           fill='purple',alpha=.2)+
  theme(legend.position = c(.8,.8),legend.key.size = unit(10,'pt'),
        legend.background=element_rect(fill = alpha("white", 0.4))) +
  scale_y_continuous(breaks=seq(0, 10000,1000)) +
  scale_x_continuous(breaks=seq(0, 16,2))+
  coord_cartesian(xlim=c(0,16))
gapr_plot
ggsave('gapr_tu.pdf',gapr_plot,width=5,height=3.5)

## storing the distance table 
xlsx::write.xlsx(data.frame(peak=gapr_df$midpoint,
                gapr_start = gapr_df$start,gapr_end = gapr_df$end,
                begin=gapr_df$midpoint-5000,end=gapr_df$midpoint+10000,dis=gapr_df$dis_to_tu,gene=gapr_df$gene,
                activity=log2(1+gapr_df$activity))%>%arrange(desc(dis)),'~/Documents/R_scripts/Ecoli/TU_analysis_plots/gapr_table.xlsx',
                row.names = F)

### distance to tu vs location on chromosome
gapr_position_plot = gapr_df %>% filter(dis_to_tu>discut) %>% ggplot(aes(x=dis_to_tu,y=midpoint,color=ifelse(activity<29.34,'low','high'))) + 
  geom_point(size=2)+
  labs(color='TU activity',y='Genomic coordinate',x='Distance to closest TES') +
  geom_rug(sides="l",length=unit(.05,'npc'))+
  theme_bw(base_size = 16)+
  theme(legend.position = c(.75,.8),
        legend.background=element_rect(fill = alpha("white", 0.3)))+
  coord_cartesian(xlim=c(800,6000))+
  scale_y_continuous(breaks=seq(0, 5e6,1e6),labels = glue('{0:5}Mb')) +
  scale_x_continuous(breaks=seq(1000,6000,1000),labels = glue('{1:6}Kb'))+
  geom_hline(aes(yintercept=(3925744)),color='navy',size=1) +
  geom_hline(aes(yintercept=(1609157)),color='navy',size=1)
gapr_position_plot
ggsave('gapr_position.pdf',gapr_position_plot,width=3.5,height=3.5)

gapr_position_plot_all = gapr_df %>% ggplot(aes(x=dis_to_tu,y=midpoint,color=ifelse(activity<29.34,'low','high'))) + 
  geom_point(size=1)+
  labs(color='TU activity',y='Genomic coordinate',x='Distance to closest TES') +
  geom_rug(sides="l",length=unit(.05,'npc'))+
  theme_bw(base_size = 16)+
  theme(legend.position = c(.75,.8),
        legend.background=element_rect(fill = alpha("white", 0.3)))+
  coord_cartesian(xlim=c(-100,6000))+
  scale_y_continuous(breaks=seq(0, 5e6,1e6),labels = glue('{0:5}Mb')) +
  scale_x_continuous(breaks=seq(0,6000,1000),labels = glue('{0:6}Kb')) +
  geom_hline(aes(yintercept=(3925744)),color='navy',size=1) +
  geom_hline(aes(yintercept=(1609157)),color='navy',size=1)
gapr_position_plot_all
ggsave('gapr_position_all.pdf',gapr_position_plot_all,width=5,height=3.5)

### random GapR 
rand_gapr = fread('/Users/ziqi_fu/Documents/R_scripts/Ecoli/peaks/gapr_peaks.csv')  
set.seed(19980415)
rand_gapr$start = sample(1:4641652,dim(gapr_df)[1],replace = F)
rand_gapr$end = rand_gapr$start + rand_gapr$width
rand_gapr$midpoint = ceiling(0.5*(rand_gapr$start+rand_gapr$end))

randgapr_tu_df = 
  do.call(rbind.data.frame,
          mclapply(1:dim(rand_gapr)[1],function(g){
            peak = rand_gapr[g,]
            front = filter(tu_raw,start-1000<peak$end & strand == 0)
            back = filter(tu_raw,end+1000>peak$start & strand == 1)
            all = rbind(front,back)
            all = mutate(all,dis = abs(peak$midpoint-signed.end))
            all = arrange(all,dis)
            # return(all)
            return(all[1,])
          },mc.cores=4)
  )

rand_gapr = cbind.data.frame(rand_gapr,gene=randgapr_tu_df$`Unnamed: 0.1`,
                           tes=randgapr_tu_df$signed.end,
                           activity=randgapr_tu_df$newRPK) 

rand_gapr = rand_gapr %>% mutate(dis_to_tu = abs(midpoint-tes))

discut =  1000
rand_gapr_plot = ggplot(rand_gapr) +
  geom_hex(aes(x=log2(1+activity),y=(dis_to_tu))) +
  geom_vline(xintercept = log2(1+29.34),color='orange',lwd=1.5)+ ##29.34 is the 75% percentile of newRPK
  geom_hline(yintercept = discut,color='orange',lwd=1.5) +
  theme_bw(base_size = 16) +
  labs(x='TU Activity level in log2(1+RPK)',
       y='Distance to closest TES') +
  annotate('rect',xmin=-.1,xmax=log2(1+29.34),ymin=discut,ymax=8000,
           fill='purple',alpha=.2)+
  theme(legend.position = c(.8,.8),legend.key.size = unit(10,'pt'),
        legend.background=element_rect(fill = alpha("white", 0.4))) +
  scale_y_continuous(breaks=seq(0, 10000,1000)) +
  scale_x_continuous(breaks=seq(0, 16,2))+
  coord_cartesian(xlim=c(0,16))
rand_gapr_plot
ggsave('rand_gapr_tu.pdf',rand_gapr_plot,width=5,height=3.5)


xlsx::write.xlsx(data.frame(peak=rand_gapr$midpoint,
                gapr_start = rand_gapr$start,gapr_end = rand_gapr$end,
                begin=rand_gapr$midpoint-5000,end=rand_gapr$midpoint+10000,dis=rand_gapr$dis_to_tu,gene=rand_gapr$gene,
                activity=log2(1+rand_gapr$activity))%>%arrange(desc(dis)),'rand_gapr_table.xlsx',row.names = F)

rand_gapr_position_plot = rand_gapr %>% filter(dis_to_tu > discut) %>% ggplot(aes(x=dis_to_tu,y=midpoint,color=ifelse(activity<29.34,'low','high'))) + 
  geom_point(size=2)+
  labs(color='TU activity',y='Genomic coordinate',x='Distance to closest TES') +
  geom_rug(sides="l",length=unit(.05,'npc'))+
  theme_bw(base_size = 16)+
  theme(legend.position = c(.75,.8),
        legend.background=element_rect(fill = alpha("white", 0.3)))+
  coord_cartesian(xlim=c(800,6000))+
  scale_y_continuous(breaks=seq(0, 5e6,1e6),labels = glue('{0:5}Mb')) +
  scale_x_continuous(breaks=seq(1000,6000,1000),labels = glue('{1:6}Kb'))+
  geom_hline(aes(yintercept=(3925744)),color='navy',size=1) +
  geom_hline(aes(yintercept=(1609157)),color='navy',size=1)
rand_gapr_position_plot
ggsave('rand_gapr_position.pdf',rand_gapr_position_plot,width=3.5,height=3.5)

rand_gapr_position_plot_all = rand_gapr %>% ggplot(aes(x=dis_to_tu,y=midpoint,color=ifelse(activity<29.34,'low','high'))) + 
  geom_point(size=1)+
  labs(color='TU activity',y='Genomic coordinate',x='Distance to closest TES') +
  geom_rug(sides="l",length=unit(.05,'npc'))+
  theme_bw(base_size = 16)+
  theme(legend.position = c(.75,.8),
        legend.background=element_rect(fill = alpha("white", 0.3)))+
  coord_cartesian(xlim=c(-100,6000))+
  scale_y_continuous(breaks=seq(0, 5e6,1e6),labels = glue('{0:5}Mb')) +
  scale_x_continuous(breaks=seq(0,6000,1000),labels = glue('{0:6}Kb'))+
  geom_hline(aes(yintercept=(3925744)),color='navy',size=1) +
  geom_hline(aes(yintercept=(1609157)),color='navy',size=1)
rand_gapr_position_plot_all
ggsave('rand_gapr_position_all.pdf',rand_gapr_position_plot_all,width=5,height=3.5)

### topo

topo_tu_df = 
  do.call(rbind.data.frame,
          mclapply(1:dim(topo_df)[1],function(g){
            peak = topo_df[g,]
            front = filter(tu_raw,end-1000<peak$end & strand == 1)
            back = filter(tu_raw,start+1000>peak$start & strand == 0)
            all = rbind(front,back)
            all = mutate(all,dis = abs(peak$midpoint-signed.start))
            all = arrange(all,dis)
            # return(all)
            return(all[1,])
          },mc.cores=4))

topo_df = cbind.data.frame(topo_df,gene=topo_tu_df$`Unnamed: 0.1`,
                           tss=topo_tu_df$signed.start,
                           activity=topo_tu_df$newRPK) 

topo_df = topo_df %>% mutate(dis_to_tu = abs(midpoint-tss))

discut =  1000
topo_plot = ggplot(topo_df) +
  geom_hex(aes(x=log2(1+activity),y=(dis_to_tu))) +
  geom_vline(xintercept = log2(1+29.34),color='orange',lwd=1.5)+ ##29.34 is the 75% percentile of newRPK
  geom_hline(yintercept = discut,color='orange',lwd=1.5) +
  theme_bw(base_size = 16) +
  labs(x='TU Activity level in log2(1+RPK)',
       y='Distance to closest TSS') +
  annotate('rect',xmin=-.1,xmax=log2(1+29.34),ymin=discut,ymax=10000,
           fill='purple',alpha=.2)+
  theme(legend.position = c(.8,.8),legend.key.size = unit(10,'pt'),
        legend.background=element_rect(fill = alpha("white", 0.4))) +
  scale_y_continuous(breaks=seq(0, 10000,1000)) +
  scale_x_continuous(breaks=seq(0, 18,2))+
  coord_cartesian(xlim=c(0,18))
topo_plot
ggsave('topo_tu.pdf',topo_plot,width=5,height=3.5)

## storing the distance table 
xlsx::write.xlsx(data.frame(peak=topo_df$midpoint,
                            topo_start = topo_df$start,topo_end = topo_df$end,
                            begin=topo_df$midpoint-5000,end=topo_df$midpoint+10000,dis=topo_df$dis_to_tu,gene=topo_df$gene,
                            activity=log2(1+topo_df$activity))%>%arrange(desc(dis)),'~/Documents/R_scripts/Ecoli/TU_analysis_plots/topo_table.xlsx',
                 row.names = F)

### visualize dis vs position
topo_position_plot = topo_df %>% filter(dis_to_tu>discut) %>% ggplot(aes(x=dis_to_tu,y=midpoint,color=ifelse(activity<29.34,'low','high'))) + 
  geom_point(size=2)+
  labs(color='TU activity',y='Genomic coordinate',x='Distance to closest TSS') +
  geom_rug(sides="l",length=unit(.05,'npc'))+
  theme_bw(base_size = 16)+
  theme(legend.position = c(.75,.8),
        legend.background=element_rect(fill = alpha("white", 0.3)))+
  coord_cartesian(xlim=c(800,10000))+
  scale_y_continuous(breaks=seq(0, 5e6,1e6),labels = glue('{0:5}Mb')) +
  scale_x_continuous(breaks=seq(1000,10000,2000),labels = glue('{seq(1,10,2)}Kb')) +
  geom_hline(aes(yintercept=(3925744)),color='navy',size=1) +
  geom_hline(aes(yintercept=(1609157)),color='navy',size=1)
topo_position_plot
ggsave('topo_position.pdf',topo_position_plot,width=3.5,height=3.5)

topo_position_plot_all = topo_df %>% ggplot(aes(x=dis_to_tu,y=midpoint,color=ifelse(activity<29.34,'low','high'))) + 
  geom_point(size=1)+
  labs(color='TU activity',y='Genomic coordinate',x='Distance to closest TSS') +
  geom_rug(sides="l",length=unit(.05,'npc'))+
  theme_bw(base_size = 16)+
  theme(legend.position = c(.75,.8),
        legend.background=element_rect(fill = alpha("white", 0.3)))+
  coord_cartesian(xlim=c(-200,10000))+
  scale_y_continuous(breaks=seq(0, 5e6,1e6),labels = glue('{0:5}Mb')) +
  scale_x_continuous(breaks=seq(0,10000,2000),labels = glue('{seq(0,10,2)}Kb'))+  
  geom_hline(aes(yintercept=(3925744)),color='navy',size=1) +
  geom_hline(aes(yintercept=(1609157)),color='navy',size=1)
ggsave('topo_position_all.pdf',topo_position_plot_all,width=5,height=3.5)


### random topo 
rand_topo = fread('/Users/ziqi_fu/Documents/R_scripts/Ecoli/peaks/topo_peaks.csv')  
set.seed(415)
rand_topo$start = sample(1:4641652,dim(rand_topo)[1],replace = F)
rand_topo$end = rand_topo$start + rand_topo$width
rand_topo$midpoint = ceiling(0.5*(rand_topo$start+rand_topo$end))
## we filter out topo peaks that are 100% contained in a transcription unit
joined.tu.region = do.call(c,apply(tu_raw,1,function(r) r[3]:r[4]))
rand_topo$tu_perc = mclapply(1:dim(rand_topo)[1],function(r){
  rd = rand_topo[r,] %>% unlist() %>% as.vector() 
  return(intersect(rd[3]:rd[4],joined.tu.region) %>%length() / (1+as.numeric(rd[5])))
},mc.cores=5) %>% unlist()
tu_perc_rand_topo = rand_topo$tu_perc ## saving purpose
rand_topo = filter(rand_topo,tu_perc<1)


randtopo_tu_df = 
  do.call(rbind.data.frame,
          mclapply(1:dim(rand_topo)[1],function(g){
            peak = rand_topo[g,]
            front = filter(tu_raw,start-1000<peak$end & strand == 0)
            back = filter(tu_raw,end+1000>peak$start & strand == 1)
            all = rbind(front,back)
            all = mutate(all,dis = abs(peak$midpoint-signed.end))
            all = arrange(all,dis)
            # return(all)
            return(all[1,])
          },mc.cores=4)
  )

rand_topo = cbind.data.frame(rand_topo,gene=randtopo_tu_df$`Unnamed: 0.1`,
                             tss=randtopo_tu_df$signed.start,
                             activity=randtopo_tu_df$newRPK) 

rand_topo = rand_topo %>% mutate(dis_to_tu = abs(midpoint-tss))

discut =  1000
rand_topo_plot = ggplot(rand_topo) +
  geom_hex(aes(x=log2(1+activity),y=(dis_to_tu))) +
  geom_vline(xintercept = log2(1+29.34),color='orange',lwd=1.5)+ ##29.34 is the 75% percentile of newRPK
  geom_hline(yintercept = discut,color='orange',lwd=1.5) +
  theme_bw(base_size = 16) +
  labs(x='TU Activity level in log2(1+RPK)',
       y='Distance to closest TSS') +
  annotate('rect',xmin=-.1,xmax=log2(1+29.34),ymin=discut,ymax=10000,
           fill='purple',alpha=.2)+
  theme(legend.position = c(.8,.8),legend.key.size = unit(10,'pt'),
        legend.background=element_rect(fill = alpha("white", 0.4))) +
  scale_y_continuous(breaks=seq(0, 10000,1000)) +
  scale_x_continuous(breaks=seq(0, 18,2))+
  coord_cartesian(xlim=c(0,18))
rand_topo_plot
ggsave('rand_topo_tu.pdf',rand_topo_plot,width=5,height=3.5)


xlsx::write.xlsx(data.frame(peak=rand_topo$midpoint,
                            topo_start = rand_topo$start,topo_end = rand_topo$end,
                            begin=rand_topo$midpoint-5000,end=rand_topo$midpoint+10000,dis=rand_topo$dis_to_tu,gene=rand_topo$gene,
                            activity=log2(1+rand_topo$activity))%>%arrange(desc(dis)),'rand_topo_table.xlsx',row.names = F)


### visualize dis vs position
rand_topo_position_plot = rand_topo %>% filter(dis_to_tu>discut) %>% ggplot(aes(x=dis_to_tu,y=midpoint,color=ifelse(activity<29.34,'low','high'))) + 
  geom_point(size=2)+
  labs(color='TU activity',y='Genomic coordinate',x='Distance to closest TSS') +
  geom_rug(sides="l",length=unit(.05,'npc'))+
  theme_bw(base_size = 16)+
  theme(legend.position = c(.75,.8),
        legend.background=element_rect(fill = alpha("white", 0.3)))+
  coord_cartesian(xlim=c(800,10000))+
  scale_y_continuous(breaks=seq(0, 5e6,1e6),labels = glue('{0:5}Mb')) +
  scale_x_continuous(breaks=seq(1000,10000,2000),labels = glue('{seq(1,10,2)}Kb'))+
  geom_hline(aes(yintercept=(3925744)),color='navy',size=1) +
  geom_hline(aes(yintercept=(1609157)),color='navy',size=1)
rand_topo_position_plot
ggsave('rand_topo_position.pdf',rand_topo_position_plot,width=3.5,height=3.5)

rand_topo_position_plot_all = rand_topo %>% ggplot(aes(x=dis_to_tu,y=midpoint,color=ifelse(activity<29.34,'low','high'))) + 
  geom_point(size=1)+
  labs(color='TU activity',y='Genomic coordinate',x='Distance to closest TSS') +
  geom_rug(sides="l",length=unit(.05,'npc'))+
  theme_bw(base_size = 16)+
  theme(legend.position = c(.75,.8),
        legend.background=element_rect(fill = alpha("white", 0.3)))+
  coord_cartesian(xlim=c(-200,10000))+
  scale_y_continuous(breaks=seq(0, 5e6,1e6),labels = glue('{0:5}Mb')) +
  scale_x_continuous(breaks=seq(0,10000,2000),labels = glue('{seq(0,10,2)}Kb'))+
  geom_hline(aes(yintercept=(3925744)),color='navy',size=1) +
  geom_hline(aes(yintercept=(1609157)),color='navy',size=1)
ggsave('rand_topo_position_all.pdf',rand_topo_position_plot_all,width=5,height=3.5)


##### compare densities #####
##### distance to closest TUs #####
gapr_control = data.frame(dis=c(gapr_df$dis_to_tu,
                                # topo_df$dis_to_start,
                                rand_gapr$dis_to_tu),
                          # random_df$dis_to_start),
                          type=c(rep('GapR',length(gapr_df$dis_to_tu)),
                                 # rep('Topo',length(topo_df$dis_to_start)),
                                 rep('Control',length(rand_gapr$dis_to_tu)))) %>%
  # rep('rand.topo',length(random_df$dis_to_start)))) %>% 
  ggplot(aes(x=log2(dis),fill=type)) +
  geom_density(alpha=.3,position = "identity")+
  mytheme+
  labs(fill='',x='Distance to closest TES, log2',y='density')+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 8))+
  theme(legend.position =c(.2,.8))+
  theme(legend.background=element_rect(fill = alpha("white", .3)))+
  annotate(geom='text',x=2,y=.05,
           label=glue('p={wilcox.test((gapr_df$dis_to_tu),rand_gapr$dis_to_tu)$p.value%>%scientific(3)}'))
gapr_control

topo_control = data.frame(dis=c(topo_df$dis_to_tu,
                                # topo_df$dis_to_start,
                                rand_topo$dis_to_tu),
                          # random_df$dis_to_start),
                          type=c(rep('topo',length(topo_df$dis_to_tu)),
                                 # rep('Topo',length(topo_df$dis_to_start)),
                                 rep('Control',length(rand_topo$dis_to_tu)))) %>%
  # rep('rand.topo',length(random_df$dis_to_start)))) %>% 
  ggplot(aes(x=log2(dis),fill=type)) +
  geom_density(alpha=.3,position = "identity")+
  mytheme+
  labs(fill='',x='Distance to closest TSS, log2',y='density')+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 8))+
  theme(legend.position =c(.2,.8))+
  theme(legend.background=element_rect(fill = alpha("white", .3)))+
  annotate(geom='text',x=3,y=.05,
           label=glue('p={wilcox.test((topo_df$dis_to_tu),rand_topo$dis_to_tu)$p.value%>%scientific(3)}'))
topo_control

ggsave('~/Documents/R_scripts/Ecoli/TU_analysis_plots/gapr_control.pdf',gapr_control,width=4.5,height=3)
ggsave('~/Documents/R_scripts/Ecoli/TU_analysis_plots/topo_control.pdf',topo_control,width=4.5,height=3)

summary(gapr_df$dis_to_end)
(dim(filter(gapr_df,dis_to_end>2500))[1] / dim(gapr_df)[1]*100) %>% round(3)
(dim(filter(topo_df,dis_to_start>2500))[1] / dim(topo_df)[1]*100) %>% round(3)
(length(random_df1$dis_to_end[random_df1$dis_to_end > 2500]) / dim(random_df1)[1]*100) %>% round(3)
(length(random_df2$dis_to_start[random_df2$dis_to_start > 2500]) / dim(random_df2)[1]*100) %>% round(3)


##### eCDF for distance #####
setwd('~/Documents/R_scripts/Ecoli/figures/TU_analysis_plots/')
ggplot(topo_df,aes(dis_to_tu)) + stat_ecdf(geom='step')
cdf_df = list(PS=gapr_df$dis_to_tu,NS=topo_df$dis_to_tu,rand.PS=rand_gapr$dis_to_tu,rand.NS=rand_topo$dis_to_tu)
cdf_df = cbind.data.frame(dis=unlist(cdf_df), peak=lapply(names(cdf_df),function(i){rep(i,length(cdf_df[[i]]))}) %>% unlist()) 
cdf_plot_log2 = cdf_df %>% ggplot(aes((log2(1+dis)),color=peak)) +
  stat_ecdf(geom='step',size=1.2) +
  labs(x='Distance to closest TES/TSS, log2 (bp)',y='Probability',color='') +
  mytheme+
  theme(legend.position = c(.3,.4),
        legend.background = element_rect(fill = alpha("white", .3))) +
  scale_y_continuous(breaks=seq(0,1,0.1))
cdf_plot_log2
ggsave('cdf_plot_log2.pdf',cdf_plot_log2,width=4.5,height=3)
cdf_plot = cdf_df %>% ggplot(aes(dis,color=peak)) +
  stat_ecdf(geom='step',size=1.2) +
  labs(x='Distance to closest TES/TSS (bp)',y='Probability',color='') +
  mytheme+
  theme(legend.position = c(.7,.4),
        legend.background = element_rect(fill = alpha("white", .3))) +
  scale_y_continuous(breaks=seq(0,1,0.1))
cdf_plot
ggsave('cdf_plot.pdf',cdf_plot,width=4.5,height=3)

##### percentage calculation #####
cutoff = 2500
#gapR
(dim(filter(gapr_df,dis_to_tu < cutoff))[1]/dim(gapr_df)[1]) %>% round(3) * 100
#random gapR
(dim(filter(rand_gapr,dis_to_tu < cutoff))[1]/dim(rand_gapr)[1]) %>% round(3) * 100
#topo
(dim(filter(topo_df,dis_to_tu < cutoff))[1]/dim(topo_df)[1]) %>% round(3) * 100
#rand topo
(dim(filter(rand_topo,dis_to_tu < cutoff))[1]/dim(rand_topo)[1]) %>% round(3) * 100

cutoff = 2500
dim(filter(gapr_df, dis_to_tu > cutoff & activity <29.34))[1] /  dim(filter(gapr_df,dis_to_tu > cutoff))[1]
dim(filter(gapr_df, dis_to_tu > cutoff & activity <29.34))[1] /  dim(gapr_df)[1] * 100
dim(filter(topo_df, dis_to_tu > cutoff & activity <29.34))[1] /  dim(filter(topo_df,dis_to_tu > cutoff))[1]
dim(filter(topo_df, dis_to_tu > cutoff & activity <29.34))[1] /  dim(topo_df)[1] * 100


##### intergenic #####
rand_gapr$tu_perc = mclapply(1:dim(rand_gapr)[1],function(r){
  rd = rand_gapr[r,] %>% unlist() %>% as.vector() 
  return(intersect(rd[3]:rd[4],joined.tu.region) %>%length() / (1+as.numeric(rd[5])))
},mc.cores=5) %>% unlist()

sanity_check_garp = mclapply(1:dim(gapr_df)[1],function(r){
  rd = gapr_df[r,] %>% unlist() %>% as.vector() 
  return(intersect(rd[3]:rd[4],joined.tu.region) %>%length() / (1+as.numeric(rd[5])))
},mc.cores=5) %>% unlist()

sanity_check_topo = mclapply(1:dim(topo_df)[1],function(r){
  rd = topo_df[r,] %>% unlist() %>% as.vector() 
  return(intersect(rd[3]:rd[4],joined.tu.region) %>%length() / (1+as.numeric(rd[5])))
},mc.cores=5) %>% unlist()


joined.tu.region = do.call(c,apply(tu_raw,1,function(r) r[3]:r[4]))

joined.gapr.region = do.call(c,apply(gapr_df,1,function(r) r[3]:r[4]))
joined.topo.region = do.call(c,apply(fread('topo_peaks.csv'),1,function(r) r[3]:r[4]))

length(intersect(joined.tu.region,unique(c(joined.gapr.region,joined.topo.region)))) / length(joined.tu.region)  #0.446 TU overlap with PS or NS
# length(intersect(joined.tu.region,intersect(joined.gapr.region,joined.topo.region))) / length(joined.tu.region)  #0.040 TU overlap with both PS and NS
length(intersect(joined.tu.region,setdiff(joined.gapr.region,joined.topo.region))) / length(joined.tu.region) #28.8% NS only 
.287391/0.446149


c(joined.gapr.region,joined.topo.region) %>% length()

length(joined.gapr.region)
length(intersect(joined.tu.region,c(joined.gapr.region,joined.topo.region))) / length(joined.tu.region) #32.8















