##EGFRmut,KRASmut two drugs
d.EGFRmut_EGFRi50CTL25 = read.csv(unz("EGFRmut_EGFRi50CTL25_1.0 simulation.zip",'values.csv'))
d.KRASmut_KRASi50CTL50 = read.csv(unz("KRASmut_KRASi50CTL50_1.0 simulation.zip",'values.csv'))

dtmp = list(d.EGFRmut_EGFRi50CTL25,d.KRASmut_KRASi50CTL50)
names(dtmp) = c('EGFRmut_EGFRi50CTL25','KRASmut_KRASi50CTL50')
dtmp = lapply(dtmp,function(x){ x$iter=1:nrow(x); x })
dtmp = lapply(dtmp,function(x){
	with(x, reshape2::melt(data.frame(iter=iter, EGFR=EGFR, RAS=RAS),id.vars='iter',variable.name='Protein',value.name='value'))
})
df = do.call(rbind,lapply(names(dtmp),function(x){ data.frame(cell_drug=x,dtmp[[x]]) }))
df$cell = gsub('_.*','',df$cell_drug)
df$drug = gsub('.*_','',df$cell_drug)
df$cell = gsub('EGFRmut','EGFR-mutant NSCLC',df$cell)
df$cell = gsub('KRASmut','KRAS-mutant NSCLC',df$cell)
df$drug = gsub('EGFRi50CTL25','+ EGFR inhibitor + ICI',df$drug)
df$drug = gsub('KRASi50CTL50','+ KRAS inhibitor + ICI',df$drug)
df = subset(df,iter<=5000)
pA = ggplot(df,aes(iter,value))+geom_line(aes(colour=Protein),linewidth=0.8)+facet_wrap(vars(cell,drug),axes='all',axis.labels='all')+theme_my(base_size=10,linewidth=0.8)+labs(x='Simulation step',y='Activity level')+scale_colour_manual(values=c('#D4562E','#682487','#4485C7'))+theme(axis.text.y = element_text(margin = margin(l=.7,unit="cm")),legend.position='bottom',legend.title=element_blank(),legend.key.spacing.x=unit(2,'cm'))


dtmp = list(d.EGFRmut_EGFRi50CTL25,d.KRASmut_KRASi50CTL50)
names(dtmp) = c('EGFRmut_EGFRi50CTL25','KRASmut_KRASi50CTL50')
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
df$drug = gsub('EGFRi50CTL25','+ EGFR inhibitor + ICI',df$drug)
df$drug = gsub('KRASi50CTL50','+ KRAS inhibitor + ICI',df$drug)
df = subset(df,iter<=5000)
pB = ggplot(df,aes(iter,value))+geom_line(aes(colour=Function),linewidth=0.8)+facet_wrap(vars(cell,drug),axes='all',axis.labels='all')+theme_my(base_size=10,linewidth=0.8)+labs(x='Simulation step',y='Activity level')+scale_colour_manual(values=c('#D4562E','#682487','#4485C7'))+theme(axis.text.y = element_text(margin = margin(l=.7,unit="cm")),legend.position='bottom',legend.title=element_blank(),legend.key.spacing.x=unit(3,'cm'))

	#load ODE result
res.ode = readRDS(r'(result_ode_simulation.rds)')

dtmp = res.ode[c("EGFRmut","KRASmut","EGFRmut_EGFRi50CTL25","KRASmut_KRASi50CTL50")]
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
df$drug = gsub('EGFRi50CTL25','+ EGFR inhibitor + ICI',df$drug)
df$drug = gsub('KRASi50CTL50','+ KRAS inhibitor + ICI',df$drug)
df$drug = gsub('^EGFRmut$','+ None',df$drug)
df$drug = gsub('^KRASmut$','+ None',df$drug)
df = subset(df,iter<=5000)
df = dplyr::rename(df,Drug=drug)
pC1 = ggplot(subset(df,grepl('EGFR',cell)),aes(iter,value))+geom_line(aes(linetype=Drug,colour=Function),linewidth=0.8)+facet_wrap(vars(cell),axes='all',axis.labels='all')+theme_my(base_size=10,linewidth=0.8)+labs(x='Simulation time',y='Activity level')+scale_colour_manual(values=c('#D4562E','#682487','#4485C7'))+theme(axis.text.y = element_text(margin = margin(l=.7,unit="cm")),legend.position='bottom',legend.title=element_blank(),legend.key.spacing.x=unit(4,'cm'))+guides(linetype=guide_legend(nrow=2),colour='none')+lims(y=c(0,1.5))
pC2 = ggplot(subset(df,grepl('KRAS',cell)),aes(iter,value))+geom_line(aes(linetype=Drug,colour=Function),linewidth=0.8)+facet_wrap(vars(cell),axes='all',axis.labels='all')+theme_my(base_size=10,linewidth=0.8)+labs(x='Simulation time',y='Activity level')+scale_colour_manual(values=c('#D4562E','#682487','#4485C7'))+theme(axis.text.y = element_text(margin = margin(l=.7,unit="cm")),legend.position='bottom',legend.title=element_blank(),legend.key.spacing.x=unit(4,'cm'),axis.title.y=element_blank())+guides(linetype=guide_legend(nrow=2),colour='none')+lims(y=c(0,1.5))
pC = ggarrange(pC1+theme(plot.margin=unit(c(.1,0,.1,1),"cm")), pC2+theme(plot.margin=unit(c(.1,1,.1,0.3),"cm")))

p = ggarrange(pA+theme(plot.margin=unit(c(.1,1,.1,1),"cm")), pB+theme(plot.margin=unit(c(.1,1,.1, 1),"cm")), pC, labels = c("A", "B", "C"), nrow=3, ncol=1, vjust=0, hjust=0)+theme(plot.margin=unit(c(0.5, 2, 0.1, 0.1),"cm"))


cairo_pdf('fig_treatment_two_drug.pdf',width=8,height=8)
p
dev.off()



