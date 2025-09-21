library(ggplot2)
library(ggpubr)
library(gridExtra)

d.Norm = read.csv(unz('normal_lung_1.0 simulation.zip','values.csv'))
d.Norm$iter = 1:nrow(d.Norm)
d.EGFRmut = read.csv(unz('EGFRmut_1.0 simulation.zip','values.csv'))
d.EGFRmut$iter = 1:nrow(d.EGFRmut)
d.KRASmut = read.csv(unz('KRASmut_1.0 simulation.zip','values.csv'))
d.KRASmut$iter = 1:nrow(d.KRASmut)

df.Norm = with(d.Norm, rbind(data.frame(iter,Protein='EGFR',value=EGFR),data.frame(iter,Protein='RAS',value=RAS)))
df.EGFRmut = with(d.EGFRmut, rbind(data.frame(iter,Protein='EGFR',value=EGFR),data.frame(iter,Protein='RAS',value=RAS)))
df.KRASmut = with(d.KRASmut, rbind(data.frame(iter,Protein='EGFR',value=EGFR),data.frame(iter,Protein='RAS',value=RAS)))
df = rbind(data.frame(cell='Normal Lung',df.Norm), data.frame(cell='EGFR-mutant NSCLC',df.EGFRmut), data.frame(cell='KRAS-mutant NSCLC',df.KRASmut))
df = subset(df,iter<=5000)
df$cell = factor(df$cell,levels=c('Normal Lung','EGFR-mutant NSCLC','KRAS-mutant NSCLC'))
pA = ggplot(df,aes(iter,value))+geom_line(aes(colour=Protein),linewidth=0.8)+facet_wrap(~cell,axes='all',axis.labels='all')+theme_my(base_size=10,linewidth=0.8)+labs(x='Simulation step',y='Activity level')+scale_colour_manual(values=c('#D4562E','#682487'))+theme(axis.text.y = element_text(margin = margin(l=1,unit="cm")),legend.position='bottom',legend.title=element_blank(),legend.key.spacing.x=unit(2,'cm'))


df.Norm = with(d.Norm, rbind(data.frame(iter,Function='Proliferation',value=proliferation),data.frame(iter,Function='Growth',value=growth),data.frame(iter,Function='Apoptosis',value=apoptosis)))
df.EGFRmut = with(d.EGFRmut, rbind(data.frame(iter,Function='Proliferation',value=proliferation),data.frame(iter,Function='Growth',value=growth),data.frame(iter,Function='Apoptosis',value=apoptosis)))
df.KRASmut = with(d.KRASmut, rbind(data.frame(iter,Function='Proliferation',value=proliferation),data.frame(iter,Function='Growth',value=growth),data.frame(iter,Function='Apoptosis',value=apoptosis)))
df = rbind(data.frame(cell='Normal Lung',df.Norm), data.frame(cell='EGFR-mutant NSCLC',df.EGFRmut), data.frame(cell='KRAS-mutant NSCLC',df.KRASmut))
df = subset(df,iter<=5000)
df$cell = factor(df$cell,levels=c('Normal Lung','EGFR-mutant NSCLC','KRAS-mutant NSCLC'))
df$Function = factor(df$Function,levels=c("Proliferation","Growth","Apoptosis"))
pB = ggplot(df,aes(iter,value))+geom_line(aes(colour=Function),linewidth=0.8)+facet_wrap(~cell,axes='all',axis.labels='all')+theme_my(base_size=10,linewidth=0.8)+labs(x='Simulation step',y='Activity level')+scale_colour_manual(values=c('#D4562E','#682487','#4485C7'))+theme(axis.text.y = element_text(margin = margin(l=1,unit="cm")),legend.position='bottom',legend.title=element_blank(),legend.key.spacing.x=unit(3,'cm'))

p = ggarrange(pA+theme(plot.margin=unit(c(1,1,0.1,1),"cm")), pB+theme(plot.margin=unit(c(1,1,0.1,1),"cm")), labels = c("A", "B"), nrow=2, ncol=1)

cairo_pdf('fig_Normal_vs_NSCLC.pdf',width=8,height=7)
p
dev.off()


