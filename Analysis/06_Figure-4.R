# -------------------------------------------------------------------------
#'This scripts generates Figure 4 in the manuscript
#
# -------------------------------------------------------------------------


# Import Data -------------------------------------------------------------

#Import Decay summary Data
summary_decay<-read_csv('Output/Figure4/summary_decay.csv')%>%
  mutate_at(c('Number of MVA Doses','alpha_val'),as.factor)

#Import Efficacy Deca data
combined_df<-read_csv('Output/Figure4/Efficacy_decay.csv')

#Import Immunogenicity experimental Data
datatab <- read_csv("Output/Data/ELISA_decay.csv") %>%
  mutate(Study=ifelse(Study==23,5,Study),Time=Week/52) %>%
  unite(Group,Dose2,`Dose 3`,Vaccinia_EXP,sep='_',remove=FALSE)%>%
  mutate_at(c('Group','Formulation','Study'),as.factor)%>%
  mutate(Group=as.integer(Group))%>%
  filter(Group==1|Group==5|Group==7|Group==9)%>%
  mutate(Group=factor(Group))

#Import Decay Immunogenicity data
q_df<-read_csv('Output/Figure4/Decay_quantiles.csv')%>%
  mutate(Group=factor(Group))

# Fig4A -------------------------------------------------------------------
# Generate Plot
plot_decay<-ggplot(data=datatab,
                   aes(x=Time,y=GMT,ymin=LCI,ymax=UCI,color=Group,fill=Group))+
  geom_point(aes(shape=label),size=2,alpha=0.8)+
  geom_line(aes(shape=label),alpha=0.3)+
  geom_errorbar(aes(shape=label),alpha=0.8)+
  labs(y='GMT',x='Time (years)')+
  scale_y_continuous(trans='log10')+
  theme_classic()+
  geom_line(data=q_df,size=1,linetype=1) +
  geom_ribbon(data=q_df,alpha=0.1,linetype=0) +
  scale_shape_manual(name='Study',
                     breaks=levels(as.factor(datatab$label)),
                     values=1:nlevels(as.factor(datatab$label)),
                     labels=levels(as.factor(datatab$label)),
                     guide=guide_legend(override.aes = list(color='black',alpha=1)))+
  theme(axis.title = element_text(size=12),
        axis.text=element_text(size=10),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12))+
  scale_color_manual('Number of \nMVA-BN Doses',
                     breaks=c(1,5,9,7),
                     values=c(colD2),
                     labels=c('1-Dose','2-Dose (Day-28)','2-Dose (Day-730)','3-Dose'))+
  scale_fill_manual('Number of \nMVA-BN Doses',
                    breaks=c(1,5,9,7),
                    values=c(colD2),
                    labels=c('1-Dose','2-Dose (Day-28)','2-Dose (Day-730)','3-Dose'))+
  guides(fill='none',color='none')

#Save Figure
save_file='Output\\Figures\\longterm_antibody_decay_years.jpg'
ggsave(save_file,plot_decay,width=10,height=5,dpi=600)
save_file='Output\\Figures\\longterm_antibody_decay_years.eps'
ggsave(save_file,plot_decay,width=10,height=5,dpi=600)


# Fig4b -------------------------------------------------------------------

#Generate Plot
plot_GMT<-ggplot(data=summary_decay,aes(color=`Number of MVA Doses`,x=`Final Dose Day`)) +
  geom_point(aes(y=I_GMT),position=position_dodge2(width=0.1),size=2.5) +
  geom_errorbar(aes(ymin=I_LCI,ymax=I_UCI),position=position_dodge2(width=10,preserve='single'),size=0.8) +
  scale_y_continuous(trans = 'log10',limits=c(3.5,4500)) +
  scale_x_continuous(trans = 'log10',
                     breaks=c(1,7,28,730),
                     labels=c('0','7','28','730'))+
  labs(x='Interval between first and last dose',y='Peak ELISA endpoint GMT')+
  scale_color_manual('Number of\nMVA-BN Doses',
                     breaks=c('0','1','2'),
                     values=colD,
                     labels=c('1-Dose','2-Doses','3-Doses')) +
  theme_classic()+
  theme(axis.text = element_text(size=10),
        axis.title = element_text(size=12),
        legend.text=element_text(size=10),
        legend.title = element_text(size=12))

#Save Figure
ggsave('Output//Figure4//Fig4b.pdf',
       plot=plot_GMT,
       width=5,height=4,dpi=1000)


# Fig 4C ------------------------------------------------------------------

plot_longGMT<-ggplot(data=summary_decay,aes(color=`Number of MVA Doses`,x=`Final Dose Day`,alpha=alpha_val)) +
  geom_point(aes(y=GMT),position=position_dodge2(width=0.1),size=2.5) +
  geom_errorbar(aes(ymin=LCI,ymax=UCI),position=position_dodge2(width=10,preserve='single'),size=0.8) +
  ylab('Predicted ELISA endpoint GMT\n at 1 year') +
  xlab('Interval between first and last dose')+
  scale_y_continuous(trans = 'log10',limits=c(3.5,4500)) +
  scale_x_continuous(trans = 'log10',
                     breaks=c(1,7,28,730),
                     labels=c('0','7','28','730'))+
  scale_color_manual('Number of\nMVA-BN Doses',
                     breaks=c('0','1','2'),
                     values=colD,
                     labels=c('1-Dose','2-Doses','3-Doses')) +
  scale_alpha_manual(breaks=c('0','1'),
                     values=c(0.5,1),
                     guide='none') +
  theme_classic()+
  theme(axis.text = element_text(size=10),
        axis.title = element_text(size=12),
        legend.text=element_text(size=10),
        legend.title = element_text(size=12))

#Save Figure
ggsave('Output//Figure4//Fig4b.pdf',
       plot=plot_longGMT,
       width=5,height=5,dpi=300)

# Fig4D -------------------------------------------------------------------

plot_longterm<-ggplot(data=combined_df,
                      aes(x=time,y=Efficacy,ymin=LCI,ymax=UCI,color=Dose,fill=Dose))+
  geom_line(size=1)+
  geom_ribbon(size=1,linetype=0,alpha=0.3)+
  xlab('Time (years)') +
  scale_y_continuous(limits=c(0,100),expand = c(0,0))+
  theme_classic()+
  scale_x_continuous(n.breaks=6,limits=c(0,10),expand = c(0,0.2))+
  annotate('rect', xmin=2, xmax=10, ymin=0, ymax=100, alpha=.1, fill='black')+
  facet_grid(cols=vars(Dose))+
  scale_color_manual(name='Number of Doses',
                     breaks=levels(as.factor(combined_df$Dose)),
                     values=c(colD2),
                     labels=c('1-Dose','2-Dose (Day-28)','2-Dose(Day-730)','3-Dose'))+
  scale_fill_manual(name='Number of Doses',
                    breaks=levels(as.factor(combined_df$Dose)),
                    values=c(colD2),
                    labels=c('1-Dose','2-Dose (Day-28)','2-Dose(Day-730)','3-Dose'))+
  theme(strip.background = element_blank(),
        strip.text=element_text(size=12),
        axis.text=element_text(size=11),
        axis.title = element_text(size=12),
        legend.text = element_text(size=11),
        legend.title = element_text(size=12))

ggsave('Output/Figures/Manuscript_efficacy_decay.jpg',
       plot=plot_longterm,
       width=10,height=5,dpi=600)
ggsave('Output/Figures/Manuscript_efficacy_decay.pdf',
       plot=plot_longterm,
       width=10,height=5,dpi=600)


# Fig4 --------------------------------------------------------------------
plot_decay_base<-ggarrange(plot_GMT,plot_longGMT,
                           labels=c("B","C"),
                           ncol=2,common.legend = TRUE,legend='right')
decay_model_full<-ggarrange(plot_decay,plot_decay_base,plot_longterm,
                            labels=c("A","","D"),
                            nrow=3)
show(decay_model_full)
ggsave('Output/Figure4/Manuscript_decay_full.jpg',
       plot=decay_model_full,
       height=11,width=9,dpi=600)
ggsave('Output/Figure4/Manuscript_decay_full.eps',
       plot=decay_model_full,
       height=11,width=9,dpi=600)

