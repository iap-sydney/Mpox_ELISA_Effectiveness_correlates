
# Neutralising antibody model ---------------------------------------------


#Read the data from the csv file.
datatab <- read_csv("Output/Data/ELISA_full_set_update.csv",show_col_types = FALSE)%>%
  mutate(ELISA=replace_na(ELISA,0))%>%
  mutate(Dose=if_else(Dose=='Dryvax','1st Gen',Dose))
datatab$col_group<-paste(datatab$Dose,datatab$Formulation)

#Convert strings to factors for indexing
datatab$Dose <- as.factor(datatab$Dose)
datatab$Study <- as.factor(datatab$Study)
datatab$ELISA <- as.factor(datatab$ELISA)
datatab$Formulation <- as.factor(datatab$Formulation)
datatab$col_group <- as.factor(datatab$col_group)
meandata<-datatab%>%
  group_by(Dose)%>%
  summarise(meandose=mean(loggmt1))

Tx<-sum(datatab$loggmt1^2)
#Define the model data
data <- list(
  N = nrow(datatab),
  y_n = log10(datatab$GMT),
  s_n= datatab$std,
  n_n=datatab$Number,
  num_doses=nlevels(datatab$Dose),
  num_studies=nlevels(datatab$Study),
  num_form=nlevels(datatab$Formulation),
  study_ind=as.integer(datatab$Study),
  dose_ind=as.integer(datatab$Dose),
  form_ind=as.integer(datatab$Formulation)
)

#fit the model data using stan
fit_unpooled <- stan(
  file = "Stan/neut_model_regress.stan",  # Stan program
  data = data,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 2000,          # number of warmup iterations per chain
  iter = 7000,            # total number of iterations per chain
  seed=10,
  cores=4,
  control = list(max_treedepth=12,stepsize=0.01)
)


# Extract samples ---------------------------------------------------------

samples<-rstan::extract(fit_unpooled)
means<-samples$mu_n
quants<-t(apply(means,2,quantile,probs=c(0.025,0.5,0.975)))
df<-data.frame(quants)
df$Study<-datatab$label
df$Dose<-datatab$Dose
df$Formulation <- gsub('Dryvax','Historic Vaccination',datatab$Formulation)
phi <- samples$sigma_s
mu_1 <- rnorm(length(phi),0,phi)+samples$mu_d[,1]
mean_quant <- quantile(samples$mu_d[,1],probs=c(0.025,0.5,0.975))
mean_quant1 <- 10**quantile(samples$mu_d[,1],probs=c(0.025,0.5,0.975))
mean_quant2 <- 10**quantile(samples$mu_d[,2],probs=c(0.025,0.5,0.975))
mean_quant3 <- 10**quantile(samples$mu_d[,3],probs=c(0.025,0.5,0.975))
mean_quant4 <- 10**quantile(samples$mu_d[,2]+samples$mu_f,probs=c(0.025,0.5,0.975))
mean_quant5 <- 10**quantile(samples$mu_d[,3]+samples$mu_f,probs=c(0.025,0.5,0.975))
mean_quants <- as_tibble(t(data.frame(mean_quant1,mean_quant4,mean_quant2,mean_quant5,mean_quant3)))%>%
  mutate(Dose=gsub('Historic Vaccination First Generation','1st Gen',gsub(' [A-Z]*$','',levels(datatab$col_group))))%>%
  mutate(col_group=levels(datatab$col_group))%>%
  mutate(logmin=log10(`2.5%`),logmed=log10(`50%`),logmax=log10(`97.5%`))

colForm<-c('darkorange2','dodgerblue','dodgerblue4','mediumpurple','mediumpurple4')
datatab$Dose <- gsub('Historic Vaccination','1st Gen',datatab$Dose)

# Plot Figures ------------------------------------------------------------

col<-c("dodgerblue2","mediumpurple3","darkorange2")

dose_names<-as_labeller(list('Historic Vaccination'='1st Gen',
                              'MVA-BN 1-dose'='1dose',
                              'MVA-BN 1-dose'='2dose'))

plot1<-ggplot(datatab,aes(x=label,y=GMT,shape=ELISA,color=col_group,group=Arm)) +
  geom_point(position=position_dodge(width=0.5),size=2.5) +
  facet_grid(cols=vars(Dose),scales='free',space='free',labeller = dose_names) + 
  geom_errorbar(aes(ymin=LCI,ymax=UCI), position=position_dodge(width=0.5)) + 
  geom_hline(data=mean_quants,aes(yintercept=`50%`,color=col_group),alpha=1) + 
  geom_hline(data=mean_quants,aes(yintercept=`2.5%`,color=col_group),alpha=1,linetype='dashed') + 
  geom_hline(data=mean_quants,aes(yintercept=`97.5%`,color=col_group),alpha=1,linetype='dashed') +
  theme_classic() +
  labs(y='ELISA endpoint GMT',x='Immunogenicity Study')+
  scale_shape_manual(values=c(15,16,17),
                     labels = c('Unspecified', "450nm/0.3",'400nm/0.35')) +
  scale_color_manual(values=colForm) +
  scale_y_continuous(trans='log10') +
  scale_x_discrete(guide = guide_axis(angle = 65))+
  theme(strip.background = element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        axis.title.x = element_blank(),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14),
        strip.text = element_text(size=14))+
  guides(color='none')
show(plot1)

ggsave('Output//Figures//Mean_plots.jpg',width=10,height=6,dpi=2000)
ggsave('Output//Figures//Mean_plots.eps',width=10,height=6,dpi=2000)
ggsave('Output//Figures//Mean_plots.pdf',width=10,height=6,dpi=2000)


# Output Study Means ------------------------------------------------------

mean_vaccine <- samples$mu_d
mead_df <- data.frame(mean_vaccine)
colnames(mead_df)<-c('Historic Vaccination','MVA-BN 1-dose','MVA-BN 2-dose')
write.csv(mead_df,'Output/Data/antibody_means.csv')

# Model without study effect ----------------------------------------------

data <- list(
  N = nrow(datatab),
  y_n = log10(datatab$GMT),
  s_n= datatab$std,
  n_n=datatab$Number,
  num_doses=nlevels(datatab$Dose),
  num_form=nlevels(datatab$Formulation),
  dose_ind=as.integer(datatab$Dose),
  form_ind=as.integer(datatab$Formulation)
)

fit_pooled <- stan(
  file = "Stan/neut_model_pooled.stan",  # Stan program
  data = data,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 2000,          # number of warmup iterations per chain
  iter = 7000,            # total number of iterations per chain
  seed=10,
  cores=4,
  control = list(max_treedepth=12,stepsize=0.01)
)

loo1<-loo(fit_unpooled)
#loo2<-loo(fit_pooled)
#loo_compare(loo1,loo2)


# Mean Plots --------------------------------------------------------------

df1dose <- datatab[as.integer(datatab$Dose)==1,]
df2dose <- datatab[as.integer(datatab$Dose)==2,]
dfdrydose <- datatab[as.integer(datatab$Dose)==3,]

shapes<-c(1,2,3)
col<-c('red2','springgreen4','dodgerblue2')

plot1<-ggplot(df1dose,aes(x=label,y=GMT,shape=ELISA,color=Formulation)) +
  geom_point(position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin=LCI,ymax=UCI), position=position_dodge(width=0.5)) +
  geom_hline(aes(yintercept=mean_quant1[2]),alpha=0.5) +
  geom_hline(aes(yintercept=mean_quant1[1]),alpha=0.3,linetype='dashed') +
  geom_hline(aes(yintercept=mean_quant1[3]),alpha=0.3,linetype='dashed') +
  scale_y_continuous(limits=c(1,1500),trans='log10') +
  scale_shape_manual(name='ELISA',
                     breaks=levels(df1dose$ELISA),
                     values=c(1,2,3),
                     labels = c('Unspecified', "450nm/0.3",'400nm/0.35')) +
  scale_color_manual(name='Formulation',
                     breaks=levels(df1dose$Formulation),
                     values=col,
                     labels = c('Freeze Dried', "First Generation",'Liquid Frozen')) +
  theme_classic() +
  scale_x_discrete(guide = guide_axis(angle = 45))
show(plot1)

plot2<-ggplot(df2dose,aes(x=label,y=GMT,shape=ELISA,color=Formulation)) +
  geom_point(position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin=LCI,ymax=UCI), position=position_dodge(width=0.5)) +
  geom_hline(aes(yintercept=mean_quant2[2]),alpha=0.5) +
  geom_hline(aes(yintercept=mean_quant2[1]),alpha=0.3,linetype='dashed') +
  geom_hline(aes(yintercept=mean_quant2[3]),alpha=0.3,linetype='dashed') +
  scale_y_continuous(limits=c(1,1500),trans='log10') +
  scale_shape_manual(name='ELISA',
                     breaks=levels(df1dose$ELISA),
                     values=c(1,2,3),
                     labels = c('Unspecified', "450nm/0.3",'400nm/0.35')) +
  scale_color_manual(name='Formulation',
                     breaks=levels(df1dose$Formulation),
                     values=col,
                     labels = c('Freeze Dried', "First Generation",'Liquid Frozen')) +
  theme_classic() +
  scale_x_discrete(guide = guide_axis(angle = 45))
show(plot2)

plot3<-ggplot(dfdrydose,aes(x=label,y=GMT,shape=ELISA,color=Formulation)) +
  geom_point(position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin=LCI,ymax=UCI), position=position_dodge(width=0.5)) +
  geom_hline(aes(yintercept=mean_quant3[2]),alpha=0.5) +
  geom_hline(aes(yintercept=mean_quant3[1]),alpha=0.3,linetype='dashed') +
  geom_hline(aes(yintercept=mean_quant3[3]),alpha=0.3,linetype='dashed') +
  scale_y_continuous(limits=c(1,1500),trans='log10') +
  scale_shape_manual(name='ELISA',
                     breaks=levels(df1dose$ELISA),
                     values=c(2,2,3),
                     labels = c('Unspecified', "450nm/0.3",'400nm/0.35')) +
  scale_color_manual(name='Formulation',
                     breaks=levels(df1dose$Formulation),
                     values=col,
                     labels = c('Freeze Dried', "First Generation",'Liquid Frozen')) +
  theme_classic()+
  scale_x_discrete(guide = guide_axis(angle = 45))
show(plot3)

figure<-ggarrange(plot1,plot2,plot3,ncol=3,common.legend = TRUE, label.y='GMT') + 
  scale_shape_manual(name='ELISA',
                     breaks=levels(df1dose$ELISA),
                     values=c(2,3,4),
                     labels = c('Unspecified', "450nm/0.3",'400nm/0.35')) +
  scale_color_manual(name='Formulation',
                     breaks=levels(df1dose$Formulation),
                     values=col,
                     labels=levels(df$Formulation))
show(figure)

ggsave('Output//Figures//Mean_plots2.png',scale=1.5)


# Formulation Effect Plot -------------------------------------------------
samples<-rstan::extract(fit_unpooled)
df_formulation <- data.frame(samples$mu_f)
plot_formulation <- ggplot(df_formulation,aes(x=samples$mu_f))+
  geom_density()+
  theme_classic()+
  ylab('Posterior Density')+
  xlab(TeX('Formulation fixed effect, $y_f$'))
show(plot_formulation)

ggsave('Output//Figures//formulation_effect.png',width=10,height=5)



# ELISA effect ------------------------------------------------------------

df_ELISA <- filter(datatab,as.numeric(ELISA)> 1)
data <- list(
  N = nrow(df_ELISA),
  y_n = log10(df_ELISA$GMT),
  s_n= df_ELISA$std,
  n_n=df_ELISA$Number,
  num_doses=nlevels(as.factor(df_ELISA$Dose)),
  num_studies=nlevels(df_ELISA$Study),
  num_form=nlevels(df_ELISA$ELISA)-1,
  study_ind=as.integer(df_ELISA$Study),
  dose_ind=as.integer(as.factor(df_ELISA$Dose)),
  form_ind=as.integer(df_ELISA$ELISA)-1
)
fit_ELISA<- stan(
  file = "Stan/neut_model_regress.stan",  # Stan program
  data = data,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 2000,          # number of warmup iterations per chain
  iter = 7000,            # total number of iterations per chain
  seed=10,
  cores=4,
  control = list(max_treedepth=12,stepsize=0.01)
)


samples_ELISA<-rstan::extract(fit_ELISA)
df_formulation <- data.frame(samples_ELISA$mu_f)
plot_ELISA <- ggplot(df_formulation,aes(x=samples_ELISA$mu_f))+
  geom_density()+
  theme_classic()+
  ylab('Posterior Density')+
  xlab(TeX('ELISA fixed effect, $y_E$'))
show(plot_ELISA)

ggsave('Output//Figures//Elisa_effect.png',width=10,height=5)

plot_FE<-ggarrange(plot_ELISA,plot_formulation,labels=c("A",'B'))
show(plot_FE)
ggsave('Output/Figures/Supp_fixed_effects.jpg',plot = plot_FE,width=6,height=3,dpi=1000)
ggsave('Output/Figures/Supp_fixed_effects.eps',plot = plot_FE,width=6,height=3,dpi=1000)
