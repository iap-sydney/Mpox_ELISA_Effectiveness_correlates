
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
       plot=summary_decay,
       width=5,height=4,dpi=1000)
ggsave('Output//Figures//supp_4.pdf',
       plot=summary_decay,
       width=5,height=4,dpi=1000)