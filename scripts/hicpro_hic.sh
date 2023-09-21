## this script converts HiC-Pro outputs to .hic format (for UCSC genome browser)
## Author: Ziqi Fu
## zfu17@jhu.edu

## processing of the LB1 files
sh /Volumes/ZF_San_Disk/HiC-Pro/bin/utils/hicpro2juicebox.sh \
  -i /Volumes/ZF_San_Disk/Guo_HiC_processed/hic_results/data/LB1/Ecoli-LB1_D18-10876_Ecoli.bwt2pairs.validPairs \
  -g /Users/ziqi_fu/HiC-Pro-3.1.0/annotation/E.coli_MG1655.sizes \
  -j /Users/ziqi_fu/juicebox/juicer_tools_1.22.01.jar \
  -r /Users/ziqi_fu/HiC-Pro-3.1.0/annotation/E.coli_MG1655_HpaII.bed \
  -o /Volumes/ZF_San_Disk/Guo_HiC_processed/hic_results/data/LB1


## processing of the LB2 files
sh /Volumes/ZF_San_Disk/HiC-Pro/bin/utils/hicpro2juicebox.sh \
  -i /Volumes/ZF_San_Disk/Guo_HiC_processed/hic_results/data/LB2/Ecoli-LB2_D18-10877_Ecoli.bwt2pairs.validPairs \
  -g /Users/ziqi_fu/HiC-Pro-3.1.0/annotation/E.coli_MG1655.sizes \
  -j /Users/ziqi_fu/juicebox/juicer_tools_1.22.01.jar \
  -r /Users/ziqi_fu/HiC-Pro-3.1.0/annotation/E.coli_MG1655_HpaII.bed \
  -o /Volumes/ZF_San_Disk/Guo_HiC_processed/hic_results/data/LB2


## processing of the M9 files, didn't use in our study
sh /Volumes/ZF_San_Disk/HiC-Pro/bin/utils/hicpro2juicebox.sh \
  -i /Volumes/ZF_San_Disk/Guo_HiC_processed/hic_results/data/M9/Ecoli-M9_D18-10878_Ecoli.bwt2pairs.validPairs \
  -g /Users/ziqi_fu/HiC-Pro-3.1.0/annotation/E.coli_MG1655.sizes \
  -j /Users/ziqi_fu/juicebox/juicer_tools_1.22.01.jar \
  -r /Users/ziqi_fu/HiC-Pro-3.1.0/annotation/E.coli_MG1655_HpaII.bed \
  -o /Volumes/ZF_San_Disk/Guo_HiC_processed/hic_results/data/M9
