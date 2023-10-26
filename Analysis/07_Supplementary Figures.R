
# Supplementary Figures ---------------------------------------------------



# Fig S3A -----------------------------------------------------------------
df_formulation <- read_csv('Output/Samples/ELISA_effect.csv')
plot_ELISA <- ggplot(df_formulation,aes(x=`samples_ELISA.mu_f`))+
  geom_density()+
  theme_classic()+
  ylab('Posterior Density')+
  xlab(TeX('ELISA fixed effect, $y_E$'))

ggsave('Output//Figures//Elisa_effect.png',width=10,height=5)


# Fig S3B -----------------------------------------------------------------
df_formulation <- read_csv('Output/Samples/antibody_means.csv')
plot_formulation <- ggplot(df_formulation,aes(x=Formulation))+
  geom_density()+
  theme_classic()+
  ylab('Posterior Density')+
  xlab(TeX('Formulation fixed effect, $y_f$'))

ggsave('Output//Figures//formulation_effect.png',width=10,height=5)


# Fig S3 ------------------------------------------------------------------

plot_FE<-ggarrange(plot_ELISA,plot_formulation,labels=c("A",'B'))
ggsave('Output/Figures/Supp_fixed_effects.jpg',
       plot = plot_FE,
       width=6,height=3,dpi=1000)
ggsave('Output/Figures/Supp_fixed_effects.eps',
       plot = plot_FE,
       width=6,height=3,dpi=1000)
