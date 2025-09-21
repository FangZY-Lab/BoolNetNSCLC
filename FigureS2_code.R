library(ggplot2)
library(gridExtra)
library(ggpubr)

##Robustness: drug efficacy
plist = list()
for(x in c(15,20,25,30)){
	d.EGFRmut_EGFRi50CTL25 = read.csv(unz(sprintf("result__Robust_on_CTL__EGFRmutNSCLC__CTL%d.zip",x),'values.csv'))
	d.KRASmut_KRASi50CTL50 = read.csv(unz(sprintf("result__Robust_on_CTL__KRASmutNSCLC__CTL%d.zip",2*x),'values.csv'))
	
	dtmp = list(d.EGFRmut_EGFRi50CTL25,d.KRASmut_KRASi50CTL50)
	names(dtmp) = c('EGFRmut_EGFRi50CTL25','KRASmut_KRASi50CTL50')
	dtmp = lapply(dtmp,function(x){ x$iter=1:nrow(x); x })
	dtmp = lapply(dtmp,function(x){
		with(x, reshape2::melt(data.frame(iter=iter, EGFR=EGFR, RAS=RAS, CTL=CTL),id.vars='iter',variable.name='Protein',value.name='value'))
	})
	df = do.call(rbind,lapply(names(dtmp),function(x){ data.frame(cell_drug=x,dtmp[[x]]) }))
	df$cell = gsub('_.*','',df$cell_drug)
	df$drug = gsub('.*_','',df$cell_drug)
	df$cell = gsub('EGFRmut','EGFR-mutant NSCLC',df$cell)
	df$cell = gsub('KRASmut','KRAS-mutant NSCLC',df$cell)
	df$drug = gsub('EGFRi50CTL25','+ EGFR inhibitor + ICI',df$drug)
	df$drug = gsub('KRASi50CTL50','+ KRAS inhibitor + ICI',df$drug)
	df = subset(df,iter<=5000)
	pA = ggplot(df,aes(iter,value))+geom_line(aes(colour=Protein),linewidth=0.8)+facet_wrap(vars(cell,drug),axes='all',axis.labels='all')+theme_my(base_size=10,linewidth=0.8)+labs(x='Simulation step',y='Activity level')+scale_colour_manual(values=c('#D4562E','#682487','#4485C7'))+theme(axis.text.y = element_text(margin = margin(l=.7,unit="cm")),legend.position='bottom',legend.title=element_blank(),legend.key.spacing.x=unit(1,'cm'))+lims(y=c(0,100))


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
	pB = ggplot(df,aes(iter,value))+geom_line(aes(colour=Function),linewidth=0.8)+facet_wrap(vars(cell,drug),axes='all',axis.labels='all')+theme_my(base_size=10,linewidth=0.8)+labs(x='Simulation step',y='Activity level')+scale_colour_manual(values=c('#D4562E','#682487','#4485C7'))+theme(axis.text.y = element_text(margin = margin(l=.7,unit="cm")),legend.position='bottom',legend.title=element_blank(),legend.key.spacing.x=unit(2,'cm'))+lims(y=c(0,100))

	p = ggarrange(pA+theme(plot.margin=unit(c(.1,1,.1,1),"cm")), pB+theme(plot.margin=unit(c(.1,1,.1, 1),"cm")), ncol=1)
	plist[[as.character(x)]] = p
}

p = ggarrange(plist[[1]],plist[[2]],plist[[3]],plist[[4]],ncol=2,nrow=2,labels=LETTERS[1:4])
cairo_pdf('fig_robust_CTL_efficacy.pdf',width=10,height=10)
print(p)
dev.off()


