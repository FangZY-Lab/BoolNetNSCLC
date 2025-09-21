##drug resistance
d.KRASmut_KRASi50to10 = read.csv(unz('KRASmut_KRASi50_iter17666KRASi10_1.0 simulation.zip','values.csv'))
s.KRASmut_KRASi50to10 = read.csv(unz('KRASmut_KRASi50_iter17666KRASi10_1.0 simulation.zip','settings.csv'))
t.resist = which(s.KRASmut_KRASi50to10$KRAS.inhibtor==10)[1]
d.KRASmut_KRASi50to10 = d.KRASmut_KRASi50to10[-(1:(t.resist-5000)),]
s.KRASmut_KRASi50to10 = s.KRASmut_KRASi50to10[-(1:(t.resist-5000)),]

d.KRASmut_CTL50to10 = read.csv(unz('KRASmut_CTL50_iter13915CTL10_1.0 simulation.zip','values.csv'))
s.KRASmut_CTL50to10 = read.csv(unz('KRASmut_CTL50_iter13915CTL10_1.0 simulation.zip','settings.csv'))
t.resist = which(s.KRASmut_CTL50to10$CTL==10)[1]
d.KRASmut_CTL50to10 = d.KRASmut_CTL50to10[-(1:(t.resist-10000)),]
s.KRASmut_CTL50to10 = s.KRASmut_CTL50to10[-(1:(t.resist-10000)),]

d.KRASmut_KRASi50CTL50to10 = read.csv(unz('KRASmut_KRASi50CTL50_iter6012KRASi10_iter11021CTL10_1.0 simulation.zip','values.csv'))
s.KRASmut_KRASi50CTL50to10 = read.csv(unz('KRASmut_KRASi50CTL50_iter6012KRASi10_iter11021CTL10_1.0 simulation.zip','settings.csv'))
t.resist = which(s.KRASmut_KRASi50CTL50to10$KRAS.inhibtor==10)[1]
d.KRASmut_KRASi50CTL50to10 = d.KRASmut_KRASi50CTL50to10[-(1:(t.resist-5000)),]
s.KRASmut_KRASi50CTL50to10 = s.KRASmut_KRASi50CTL50to10[-(1:(t.resist-5000)),]

t.resist2 = which(s.KRASmut_KRASi50CTL50to10$CTL==10)[1] #only for this one
t.resist = 5000 #common for all


dtmp = list(d.KRASmut_KRASi50to10, d.KRASmut_CTL50to10, d.KRASmut_KRASi50CTL50to10)
names(dtmp) = c('KRASmut_KRASi50to10', 'KRASmut_CTL50to10', 'KRASmut_KRASi50CTL50to10')
dtmp = lapply(dtmp,function(x){ x$iter=1:nrow(x); x })
dtmp = lapply(dtmp,function(x){
	with(x, reshape2::melt(data.frame(iter=iter, EGFR=EGFR, RAS=RAS),id.vars='iter',variable.name='Protein',value.name='value'))
})
df = do.call(rbind,lapply(names(dtmp),function(x){ data.frame(cell_drug=x,dtmp[[x]]) }))
df$cell = gsub('_.*','',df$cell_drug)
df$drug = gsub('to10','',gsub('.*_','',df$cell_drug))
df$cell = gsub('KRASmut','KRAS-mutant NSCLC',df$cell)
df$drug = gsub('KRASi50CTL50','+ KRAS inhibitor + ICI',df$drug)
df$drug = gsub('KRASi50','+ KRAS inhibitor',df$drug)
df$drug = gsub('CTL50','+ ICI',df$drug)
df$drug = factor(df$drug,levels=c('+ KRAS inhibitor','+ ICI','+ KRAS inhibitor + ICI'))
df = subset(df,iter<=15000)
pA = ggplot(df,aes(iter,value))+geom_line(aes(colour=Protein),linewidth=0.8)+facet_wrap(vars(cell,drug),axes='all',axis.labels='all')+theme_my(base_size=10,linewidth=0.8)+labs(x='Simulation step',y='Activity level')+scale_colour_manual(values=c('#D4562E','#682487','#4485C7'))+theme(axis.text.y = element_text(margin = margin(l=.7,unit="cm")),legend.position='bottom',legend.title=element_blank(),legend.key.spacing.x=unit(2,'cm'))
pA = pA+geom_vline(xintercept=c(5000,10000),color='cyan',linewidth=1,linetype=2)


dtmp = list(d.KRASmut_KRASi50to10, d.KRASmut_CTL50to10, d.KRASmut_KRASi50CTL50to10)
names(dtmp) = c('KRASmut_KRASi50to10', 'KRASmut_CTL50to10', 'KRASmut_KRASi50CTL50to10')
dtmp = lapply(dtmp,function(x){ x$iter=1:nrow(x); x })
dtmp = lapply(dtmp,function(x){
	x = with(x, reshape2::melt(data.frame(iter=iter, Proliferation=proliferation, Growth=growth, Apoptosis=apoptosis),id.vars='iter',variable.name='Function',value.name='value'))
	x$Function = factor(x$Function,levels=c("Proliferation","Growth","Apoptosis"))
	x
})
df = do.call(rbind,lapply(names(dtmp),function(x){ data.frame(cell_drug=x,dtmp[[x]]) }))
df$cell = gsub('_.*','',df$cell_drug)
df$drug = gsub('to10','',gsub('.*_','',df$cell_drug))
df$cell = gsub('KRASmut','KRAS-mutant NSCLC',df$cell)
df$drug = gsub('KRASi50CTL50','+ KRAS inhibitor + ICI',df$drug)
df$drug = gsub('KRASi50','+ KRAS inhibitor',df$drug)
df$drug = gsub('CTL50','+ ICI',df$drug)
df$drug = factor(df$drug,levels=c('+ KRAS inhibitor','+ ICI','+ KRAS inhibitor + ICI'))
df = subset(df,iter<=15000)
pB = ggplot(df,aes(iter,value))+geom_line(aes(colour=Function),linewidth=0.8)+facet_wrap(vars(cell,drug),axes='all',axis.labels='all')+theme_my(base_size=10,linewidth=0.8)+labs(x='Simulation step',y='Activity level')+scale_colour_manual(values=c('#D4562E','#682487','#4485C7'))+theme(axis.text.y = element_text(margin = margin(l=.7,unit="cm")),legend.position='bottom',legend.title=element_blank(),legend.key.spacing.x=unit(3,'cm'))
pB = pB+geom_vline(xintercept=c(5000,10000),color='cyan',linewidth=1,linetype=2)

p = ggarrange(pA+theme(plot.margin=unit(c(1,1,0.1,1),"cm")), pB+theme(plot.margin=unit(c(1,1,0.1,1),"cm")), labels = c("A", "B"), nrow=2, ncol=1)

cairo_pdf('fig_treatment_resistance.pdf',width=8,height=7)
p
dev.off()




