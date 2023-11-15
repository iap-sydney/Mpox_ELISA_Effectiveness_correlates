#Load efficacy study data
df_studies<-read_csv('Data/Efficacy_studies.csv')

# Create Figure 1A --------------------------------------------------------

#Load estimated individual study effectiveness
df_output_indiv<-read_csv('Output/Figure1/Indiv_eff.csv')

#Load combined effectiveness estimates
df_output_combined<-read_csv('Output/Figure1/combined_eff.csv')

#Create Figure Labels
Vacc_labs<-c('1st Gen','MVA-BN\n1-Dose','MVA-BN\n2-Dose')

#Create Plot
eff_plot<-ggplot(df_output_indiv,aes(x=Vaccine,y=X50.,ymin=X2.5.,ymax=X97.5.,shape=Study,group=Study,color=Study))+
  geom_point(position=position_dodge(width=0.4),size=2,alpha=1)+
  geom_errorbar( position=position_dodge(width=0.4),alpha=1)+
  geom_line(position=position_dodge2(width=0.4),alpha=0.7)+
  geom_point(data=df_output_combined,shape=16,size=2)+
  geom_errorbar(data=df_output_combined)+
  scale_y_continuous(limits=c(0,100))+
  labs(y='Vaccine Effectiveness (%)') +
  scale_x_discrete(labels=Vacc_labs) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text = element_text(size=10),
        axis.title = element_text(size=12),
        legend.text = element_text(size=10))+
  scale_shape_manual(name='Study',
                     breaks=c(df_studies$Author,'Combined'),
                     values=c(1:nlevels(as.factor(df_output_indiv$Study)),16),
                     labels=c(df_studies$Label,'Combined Estimate'))+
  scale_color_manual(name='Study',
                     breaks=c(df_studies$Author,'Combined'),
                     values=c(col_scale,'black'),
                     labels=c(df_studies$Label,'Combined Estimate'),
                     guide=guide_legend(override.aes = list(shape=c(1:nlevels(as.factor(df_output_indiv$Study)),16))))

#Save Figure
ggsave('Output//Figure1//Fig1A.pdf',
       plot=eff_plot,
       height=5,width = 8,dpi=600)


# Create Figure 1B --------------------------------------------------------
#load Study relative efficacy data
REtib<-read_csv('Output/Figure1/REtib.csv')%>%
  mutate(Vaccine='\n')

#load Combined RE data
REquant<-read_csv('Output/Figure1/RE_comb.csv')%>%
  mutate(Vaccine='\n')

#Create plot
plot_RE<-ggplot(REtib,aes(x=Vaccine,y=X2,ymin=X1,ymax=X3,shape=Study,color=Study))+
  geom_point(position=position_dodge(width=0.5),size=2,alpha=1)+
  geom_errorbar(position=position_dodge(width=0.5),alpha=1)+
  geom_point(data=REquant,size=2)+
  geom_errorbar(data=REquant)+
  scale_y_continuous(limits=c(0,1.5))+
  labs(y='Risk Ratio\n(MVA-BN 2-Dose vs 1-Dose)',x='')+
  scale_shape_manual(name='Study',
                     breaks=c(df_studies$Author,'Combined'),
                     values=c(1:nlevels(as.factor(df_output_indiv$Study)),16),
                     labels=c(df_studies$Label,'Combined Estimate'))+
  scale_color_manual(name='Study',
                     breaks=c(df_studies$Author,'Combined'),
                     values=c(col_scale,'black'),
                     labels=c(df_studies$Label,'Combined Estimate'))+
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title = element_text(size=12),
        axis.text = element_text(size=10),
        legend.text = element_text(size=10),
        plot.margin = unit(c(5.5,10,5.5,10),'points'))

#Save Figure
ggsave('Output/Figure1/Fig1B.pdf',
       plot=plot_RE,
       height=5,width = 8,dpi=600)

#Combine Figures
plot_efficacy<-ggarrange(eff_plot,plot_RE,labels=c('A','B'),ncol=2,common.legend = TRUE,widths=c(5,3),legend='right')
ggsave('Output/Figure1/Fig1.jpg',
       plot=plot_efficacy,
       height=4,width = 6,dpi=1000)
ggsave('Output/Figure1/Fig1.pdf',
       plot=plot_efficacy,
       height=4,width = 6,dpi=1000)
