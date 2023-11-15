
# Supplementary Figures ---------------------------------------------------



# Fig S3A -----------------------------------------------------------------
df_formulation <- read_csv('Output/Samples/ELISA_effect.csv')
plot_ELISA <- ggplot(df_formulation,aes(x=`samples_ELISA.mu_f`))+
  geom_density()+
  theme_classic()+
  ylab('Posterior Density')+
  xlab(TeX('ELISA fixed effect, $y_E$'))

ggsave('Output//Figures//Sup_S3a.png',width=10,height=5)


# Fig S3B -----------------------------------------------------------------
df_formulation <- read_csv('Output/Samples/antibody_means.csv')
plot_formulation <- ggplot(df_formulation,aes(x=Formulation))+
  geom_density()+
  theme_classic()+
  ylab('Posterior Density')+
  xlab(TeX('Formulation fixed effect, $y_f$'))

ggsave('Output//Figures//Sup_S3b.png',width=10,height=5)


# Fig S3 ------------------------------------------------------------------

plot_FE<-ggarrange(plot_ELISA,plot_formulation,labels=c("A",'B'))
ggsave('Output/Figures/Supp_S3.jpg',
       plot = plot_FE,
       width=6,height=3,dpi=1000)
ggsave('Output/Figures/Supp_S3.pdf',
       plot = plot_FE,
       width=6,height=3,dpi=1000)


# FigS4 -------------------------------------------------------------------

summary_decay<-read_csv('Output/Figure4/summary_decay.csv')%>%
  mutate_at(c('Number of MVA Doses','alpha_val'),as.factor)
plot_prop<-ggplot(data=summary_decay,aes(color=`Number of MVA Doses`,x=`Final Dose Day`,alpha=alpha_val)) +
  geom_point(aes(y=percent),position=position_dodge2(width=0.1),size=3) +
  geom_errorbar(aes(ymin=p_LCI,ymax=p_UCI),position=position_dodge2(width=10,preserve='single'),size=0.8) +
  ylab('Percentage of slow decay antibodies (%)') +
  scale_y_continuous(limits = c(0,100)) +
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

ggsave('Output//Figures//supp_4.jpg',
       plot=plot_prop,
       width=5,height=4,dpi=300)
ggsave('Output//Figures//supp_4.pdf',
       plot=plot_prop,
       width=5,height=4,dpi=300)


# Generate supp table 4 ---------------------------------------------------

#Import Decay Immunogenicity data
q_df<-read_csv('Output/Figure4/Decay_quantiles.csv')%>%
  mutate(Group=factor(Group),Time=Time*12)%>%
  mutate_at(c('GMT','UCI','LCI'),round)%>%
  mutate(CI=paste('(',LCI,'.',UCI,')',sep=''))%>%
  mutate(vals=paste(GMT,CI))%>%
  select(Time,vals,Group)%>%
  filter(Time %in% c(1,2,3,6,12,18,24))%>%
  pivot_wider(values_from = vals,names_from = Time)%>%
  mutate(Group=c('1-Dose','2-Dose (Day-28)','3-Dose','2-Dose (Day-730)'))%>%
  write_csv('Output/TableS4.csv')

c(1,5,9,7)
c('1-Dose','2-Dose (Day-28)','2-Dose (Day-730)','3-Dose')


# Generate linear model fit -----------------------------------------------

# Load Data ---------------------------------------------------------------

#import sampled means
mean_df<-read_csv('Output/Samples/antibody_means.csv')

#import sampled efficacies
Efficacies_df<-read_csv('Output/Samples/Combined-Effectiveness-sample.csv')


#Import Effectiveness quantiles with GMT
eff_df<-read_csv('Output/Figures/Efficacy_quantiles_lin.csv')
eff_df$LCI<-log(1-eff_df$LCI)
eff_df$median<-log(1-eff_df$median)
eff_df$UCI<-log(1-eff_df$UCI)


# pair samples from effectiveness and immunogencity -----------------------

#MVA-BN 1 dose set
dose1df<-data.frame(mean_df$`MVA-BN 1-dose`,log(1-Efficacies_df$`MVA-BN 1-dose`/100))
colnames(dose1df)<-c('GMT','Efficacy')
dose1df$Vaccine='MVA-BN 1-dose'

#MVA-BN 2 dose set
dose2df<-data.frame(mean_df$`MVA-BN 2-dose`,log(1-Efficacies_df$`MVA-BN 2-dose`/100))
colnames(dose2df)<-c('GMT','Efficacy')
dose2df$Vaccine='MVA-BN 2-dose'

#Historic Vaccination set
doseddf<-data.frame(mean_df$`Historic Vaccination`,log(1-Efficacies_df$`First Generation`/100))
colnames(doseddf)<-c('GMT','Efficacy')
doseddf$Vaccine='First Generation'

#combine sampled data
big_df=rbind(dose1df,dose2df,doseddf)


# Figure 3 ----------------------------------------------------------------

#generate figure
lin_plot<-ggplot(eff_df,aes(x=GMT,y=median)) +
  geom_line(color='black')+
  geom_ribbon(aes(ymin=LCI,ymax=UCI),color='black',fill='black',alpha=0.3,linetype=2)+
  geom_density_2d(data=big_df,aes(y=Efficacy,color=Vaccine),
                  contour_var = 'ndensity',
                  bins=5,size=0.5)+
  scale_color_manual(breaks = levels(as.factor(big_df$Vaccine)),
                     values= colV,
                     label=levels(as.factor(big_df$Vaccine)))+
  theme_classic()+
  scale_y_continuous('log(OR)') +
  scale_x_continuous(breaks=1:4,labels=10^(1:4))+
  theme(legend.text = element_text(size=10),
        legend.title = element_text(size=12),
        axis.title = element_text(size=12),
        axis.text = element_text(size=10)
  )

#Save Figure
ggsave("Output//Figures//lin_curve_fit.jpg",
       plot=logistic_plot,
       width=5,height=4,dpi=1000)
ggsave("Output//Figures//lin_curve_fit.pdf",
       plot=logistic_plot,
       width=5,height=4,dpi=1000)
