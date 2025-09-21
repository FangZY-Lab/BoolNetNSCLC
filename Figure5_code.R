##EGFRmut, KRASmut + one of inhibitor/CTL
d.EGFRmut_EGFRi50 = read.csv(unz('EGFRmut_EGFRi50_1.0 simulation.zip','values.csv'))
d.EGFRmut_CTL25 = read.csv(unz('EGFRmut_CTL25_1.0 simulation.zip','values.csv'))
d.KRASmut_KRASi50 = read.csv(unz("KRASmut_KRASi50_1.0 simulation.zip",'values.csv'))
d.KRASmut_CTL50 = read.csv(unz("KRASmut_CTL50_1.0 simulation.zip",'values.csv'))

dtmp = list(d.EGFRmut_EGFRi50,d.EGFRmut_CTL25,d.KRASmut_KRASi50,d.KRASmut_CTL50)
names(dtmp) = c('EGFRmut_EGFRi50','EGFRmut_CTL50','KRASmut_KRASi50','KRASmut_CTL50')
dtmp = lapply(dtmp,function(x){ x$iter=1:nrow(x); x })
dtmp = lapply(dtmp,function(x){
	with(x, reshape2::melt(data.frame(iter=iter, EGFR=EGFR, RAS=RAS),id.vars='iter',variable.name='Protein',value.name='value'))
})
df = do.call(rbind,lapply(names(dtmp),function(x){ data.frame(cell_drug=x,dtmp[[x]]) }))
df$cell = gsub('_.*','',df$cell_drug)
df$drug = gsub('.*_','',df$cell_drug)
df$cell = gsub('EGFRmut','EGFR-mutant NSCLC',df$cell)
df$cell = gsub('KRASmut','KRAS-mutant NSCLC',df$cell)
df$drug = gsub('CTL.*','+ ICI',df$drug)
df$drug = gsub('EGFRi50','+ EGFR inhibitor',df$drug)
df$drug = gsub('KRASi50','+ KRAS inhibitor',df$drug)
df$drug = factor(df$drug, levels=c("+ EGFR inhibitor", "+ KRAS inhibitor", "+ ICI"))
df = subset(df,iter<=5000)
pA = ggplot(df,aes(iter,value))+geom_line(aes(colour=Protein),linewidth=0.8)+facet_wrap(vars(cell,drug),axes='all',axis.labels='all')+theme_my(base_size=10,linewidth=0.8)+labs(x='Simulation step',y='Activity level')+scale_colour_manual(values=c('#D4562E','#682487','#4485C7'))+theme(axis.text.y = element_text(margin = margin(l=1,unit="cm")),legend.position='bottom',legend.title=element_blank(),legend.key.spacing.x=unit(2,'cm'))


dtmp = list(d.EGFRmut_EGFRi50,d.EGFRmut_CTL25,d.KRASmut_KRASi50,d.KRASmut_CTL50)
names(dtmp) = c('EGFRmut_EGFRi50','EGFRmut_CTL25','KRASmut_KRASi50','KRASmut_CTL50')
dtmp = lapply(dtmp,function(x){ x$iter=1:nrow(x); x })
dtmp = lapply(dtmp,function(x){
	x = with(x, reshape2::melt(data.frame(iter=iter, Proliferation=proliferation, Growth=growth, Apoptosis=apoptosis),id.vars='iter',variable.name='Function',value.name='value'))
	x$Function = factor(x$Function,levels=c("Proliferation","Growth","Apoptosis"))
	x
})
df = do.call(rbind,lapply(names(dtmp),function(x){ data.frame(cell_drug=x,dtmp[[x]]) }))
df$cell = gsub('_.*','',df$cell_drug)
df$drug = gsub('.*_','',df$cell_drug)
df$cell = gsub('EGFRmut','EGFR-mutant NSCLC',df$cell)
df$cell = gsub('KRASmut','KRAS-mutant NSCLC',df$cell)
df$drug = gsub('CTL.*','+ ICI',df$drug)
df$drug = gsub('EGFRi50','+ EGFR inhibitor',df$drug)
df$drug = gsub('KRASi50','+ KRAS inhibitor',df$drug)
df$drug = factor(df$drug, levels=c("+ EGFR inhibitor", "+ KRAS inhibitor", "+ ICI"))
df = subset(df,iter<=5000)
pB = ggplot(df,aes(iter,value))+geom_line(aes(colour=Function),linewidth=0.8)+facet_wrap(vars(cell,drug),axes='all',axis.labels='all')+theme_my(base_size=10,linewidth=0.8)+labs(x='Simulation step',y='Activity level')+scale_colour_manual(values=c('#D4562E','#682487','#4485C7'))+theme(axis.text.y = element_text(margin = margin(l=1,unit="cm")),legend.position='bottom',legend.title=element_blank(),legend.key.spacing.x=unit(3,'cm'))

p = ggarrange(pA+theme(plot.margin=unit(c(1,1,0.1,1),"cm")), pB+theme(plot.margin=unit(c(1,1,0.1,1),"cm")), labels = c("A", "B"), nrow=2, ncol=1)

cairo_pdf('fig_treatment_single_drug.pdf',width=8,height=10)
p
dev.off()


