#from UCSC Xena TCGA PANCAN datasets
dmut = read.table2('GDC-PANCAN.mutect2_snv.tsv.gz')
dexp = read.table2('GDC-PANCAN.htseq_fpkm.annotated.txt.gz')
info = read.table2('GDC-PANCAN.project_info.tsv.gz')
#from Thorsson et al. Immunity 2018
dciber = read.table2('TCGA.Kallisto.fullIDs.cibersort.relative.tsv')
dciber$SampleID2 = substr(gsub('\\.','-',dciber$SampleID),1,16)

samps.NSCLC = subset(info,disease_code.project%in%c('LUAD','LUSC')&program_code.project=='TCGA')$sample
samps.NSCLC = intersect(intersect(colnames(dexp),samps.NSCLC), dmut$Sample_ID)

samps.NSCLC.KRAS_G12C = subset(dmut, Sample_ID%in%samps.NSCLC & gene=='KRAS' & Amino_Acid_Change=='p.G12C')$Sample_ID
samps.NSCLC.KRAS_WT = setdiff(samps.NSCLC, subset(dmut, Sample_ID%in%samps.NSCLC & gene=='KRAS')$Sample_ID)
sig.CD8T = c("CD8A", "CD8B", "CD3D", "CD3E", "CD3G", "GZMB", "PRF1", "KLRG1", "CX3CR1")

samps.NSCLC.EGFR_mut = subset(dmut, Sample_ID%in%samps.NSCLC & gene=='EGFR' & Amino_Acid_Change%in%c('p.E746_A750del','p.L858R','p.T790M'))$Sample_ID
	samps.NSCLC.EGFR_Exon19del = subset(dmut, Sample_ID%in%samps.NSCLC & gene=='EGFR' & Amino_Acid_Change%in%c('p.E746_A750del'))$Sample_ID
	samps.NSCLC.EGFR_L858R = subset(dmut, Sample_ID%in%samps.NSCLC & gene=='EGFR' & Amino_Acid_Change%in%c('p.L858R'))$Sample_ID
	samps.NSCLC.EGFR_T790M = subset(dmut, Sample_ID%in%samps.NSCLC & gene=='EGFR' & Amino_Acid_Change%in%c('p.T790M'))$Sample_ID

samps.NSCLC.EGFR_WT = setdiff(samps.NSCLC, subset(dmut, Sample_ID%in%samps.NSCLC & gene=='EGFR')$Sample_ID)

dexp.CD8T = dexp[sig.CD8T,c(samps.NSCLC.EGFR_Exon19del,samps.NSCLC.EGFR_L858R,samps.NSCLC.KRAS_G12C)]
dexp.CD8T.labels = c()
dexp.CD8T.labels[colnames(dexp.CD8T)%in%samps.NSCLC.EGFR_Exon19del] = 'EGFR Exon19del'
dexp.CD8T.labels[colnames(dexp.CD8T)%in%samps.NSCLC.EGFR_L858R] = 'EGFR L858R'
dexp.CD8T.labels[colnames(dexp.CD8T)%in%samps.NSCLC.KRAS_G12C] = 'KRAS G12C'


#MDSC gene signature
sig.MDSC = unique(xlsx::read.xlsx('data_MDSC_signature_Alshetaiwi2020SciImmunol_mouse2humanBySynGO.xlsx',sheetIndex=1)[-1,1])
sig.MDSC.ssgsea = gsva(ssgseaParam(data.matrix(dexp[,samps.NSCLC]), list(MDSC=sig.MDSC)))

df = data.frame(
	   ID = samps.NSCLC,
	   t(dexp[c('CD274','PDCD1'),samps.NSCLC]), 
	   MDSC = sig.MDSC.ssgsea['MDSC',samps.NSCLC],
	   dciber[match(samps.NSCLC,dciber$SampleID2),],
	   KRASmut = sapply(samps.NSCLC,function(x) if(x%in%samps.NSCLC.KRAS_G12C) return('KRAS G12C') else return(NA) ), 
	   EGFRmut = sapply(samps.NSCLC,function(x) if(x%in%samps.NSCLC.EGFR_Exon19del) return('EGFR Exon19del') else if(x%in%samps.NSCLC.EGFR_L858R) return('EGFR L858R') else if(x%in%samps.NSCLC.EGFR_T790M) return('EGFR T790M') else return(NA) )
)
df = subset(df,!is.na(KRASmut)|!is.na(EGFRmut))
df$mut = with(df, ifelse(!is.na(KRASmut), KRASmut, EGFRmut))
df$mut = factor(df$mut)
immune_cols <- c(grep('cells|cyte|phage|phils',colnames(dciber),value=T),'MDSC')
library(dplyr); library(tidyr)
df_long <- df %>% pivot_longer(cols = all_of(immune_cols), names_to = "immune", values_to = "fraction")
	#restrict to those of interest
	df_long = subset(df_long, immune%in%c('Macrophages.M0','Macrophages.M1','Macrophages.M2','T.cells.regulatory..Tregs.','MDSC'))
	df_long = df_long %>% mutate(immune = recode(immune, "Macrophages.M0" = "Macrophage M0", 'Macrophages.M1'='Macrophage M1', 'Macrophages.M2'='Macrophage M2', 'T.cells.regulatory..Tregs.'='Treg'))

library(ggpubr); library(ggplot2)
set.seed(123)
p = ggboxplot(df_long,x='mut',y='fraction',color='mut',size=0.7)+geom_jitter(aes(colour=mut),width=0.25)+scale_colour_manual(values=c("#00BA38","#619CFF","#F8766D"))+ stat_compare_means( comparisons=list(c("EGFR Exon19del", "EGFR L858R"), c("EGFR Exon19del", "KRAS G12C"), c("EGFR L858R", "KRAS G12C")), label = "p.signif", bracket.size=.7, size=4, family='Arial')+theme_my(base_size=10,base_family='Arial',linewidth=0.7)+facet_wrap(~immune,scales='free_y',nrow=1) + theme(axis.text.x=element_text(angle=90,vjust=.5,hjust=1))+labs(x='',y='Estimate of Relative Distribution')

cairo_pdf('fig_NSCLC_M2_MDSC_Treg.pdf',width=10,height=5)
print(p)
dev.off()



