df_data<-read_csv('Output/Data/Eff.csv') %>%
  mutate(Study=gsub('\\d\\d-','',Study))%>%
  replace(is.na(.),0)

df_studies<-read_csv('Data/Efficacy_studies.csv')

df_data=mutate(df_data,Indiv=if_else(Dose==0,'Unvaccinated',paste(Vaccine,Study,Dose)))
df_data$Indiv<-gsub('Unvaccinated \\w*','Unvaccinated',df_data$Indiv)
df_data$Vaccinated=df_data$Vaccine
df_data$Vaccinated[df_data$Vaccine!='Unvaccinated']<-'Vaccinated'

df_data$Vaccinated=as.factor(df_data$Vaccinated)
df_data$Indiv=as.factor(df_data$Indiv)
df_data$index=as.factor(df_data$Study)
df_data$group=as.factor(df_data$group)
df_data$Vacc_dose=as.factor(df_data$Vacc_dose)

df_2dose=filter(df_data,Dose==2)%>%
  mutate(Indiv=fct_drop(Indiv))

df_unvacc=df_data[df_data$Vaccine=='Unvaccinated',]
df_vacc=df_data[df_data$Vaccine!='Unvaccinated',]

data <- list(
    J = nrow(df_vacc),
    Num_pop = df_data$Pop,
    num_cases=df_data$Case,
    num_groups=nlevels(df_data$group),
    group_index=as.integer(df_data$group),
    num_vaccines=nlevels(df_data$Vacc_dose)-1,
    Vaccine_ind=as.integer(df_data$Vacc_dose),
    num_vacc_studies=nlevels(df_data$index),
    vacc_study_ind=as.integer(df_data$index)
)
  
eff_fit <- stan(
    file = "Stan/Efficacy_regress.stan",  # Stan program
    data = data,    # named list of data
    chains = 4,             # number of Markov chains
    warmup = 5000,          # number of warmup iterations per chain
    iter = 35000,            # total number of iterations per chain
    seed=10,
    cores=4,
    control = list(max_treedepth=12,stepsize=0.5,adapt_delta=0.9)
)
eff_fit_pooled <- stan(
  file = "Stan/Efficacy_pooled.stan",  # Stan program
  data = data,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 10000,          # number of warmup iterations per chain
  iter = 40000,            # total number of iterations per chain
  seed=10,
  cores=4,
  control = list(max_treedepth=12,stepsize=0.5,adapt_delta=0.9)
)
  
samples<-rstan::extract(eff_fit)
Effs<-(samples$p_s)
quants<-t(apply(Effs,2,quantile,probs=c(0.025,0.5,0.975)))
df_output<-data.frame(quants)
df_output$Study=levels(df_data$index)

plot1<-ggplot(df_output,aes(x=Study,y=X50.))
plot1<-plot1+geom_point(position=position_dodge(width=0.5))
plot1<-plot1+geom_errorbar(aes(ymin=X2.5.,ymax=X97.5.), position=position_dodge(width=0.5))
show(plot1)

samples<-rstan::extract(eff_fit)
Effs<-100*(1-exp(samples$V))
quants<-t(apply(Effs,2,quantile,probs=c(0.025,0.5,0.975)))
df_output_combined<-data.frame(quants)
mean_quant=quantile(samples$E,probs=c(0.025,0.5,0.975))
df_output_combined$Vaccine=c('First Generation','MVA-BN 1-dose','MVA-BN 2-dose')
df_output_combined$Study='Combined'

RE_comb<-quantile(100*(1-exp(samples$V[2,]-samples$V[3,])),probs=c(0.025,0.5,0.975))

Eff_df<-data.frame(Effs)
colnames(Eff_df)<-c('First Generation','MVA-BN 1-dose','MVA-BN 2-dose')
write_csv(Eff_df,"Output//Data//Efficacy.csv")


data2 <- list(
  J = nrow(df_vacc),
  Num_pop = df_data$Pop,
  num_cases=df_data$Case,
  num_groups=nlevels(df_data$group),
  group_ind=as.integer(df_data$group),
  num_vaccines=nlevels(df_data$Indiv)-1,
  Vaccine_ind=as.integer(df_data$Indiv)
)


eff_fit_individual <- stan(
  file = "Stan/Efficacy.stan",  # Stan program
  data = data2,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 10000,          # number of warmup iterations per chain
  iter = 40000,            # total number of iterations per chain
  seed=10,
  cores=4,
  control = list(max_treedepth=12,stepsize=0.5,adapt_delta=0.9)
)

samples2<-rstan::extract(eff_fit_individual)
Effs<-100*(1-exp(samples2$V))
quants<-t(apply(Effs,2,quantile,probs=c(0.025,0.5,0.975)))
df_output_indiv<-data.frame(quants)
df_output_indiv$Study<-gsub(' \\d','',c(levels(df_data$Indiv)[1:15]))
df_output_indiv$Study<-gsub('\\w*\\W','',df_output_indiv$Study)
df_output_indiv$Vaccine<-c(gsub('MVA-BN1','MVA-BN 1-dose',gsub(' [A-Z,a-z]*','',levels(df_data$Indiv)[1:15])))
df_output_indiv$Vaccine<-gsub('Dryvax1','First Generation',df_output_indiv$Vaccine)
df_output_indiv$Vaccine<-as.numeric(as.factor(gsub('MVA-BN2','MVA-BN 2-dose',df_output_indiv$Vaccine)))
df_output=rbind(df_output_indiv,df_output_combined)

dose2_studies<-levels(df_2dose$Indiv)
studies<-levels(df_data$Indiv)
N2<-length(dose2_studies)
REtab<-data.frame(matrix(nrow=N2,ncol=3))
for (i in 1:N2){
  study2<-dose2_studies[i]
  j<-match(study2,studies)
  study1<-gsub('2','1',study2)
  k<-match(study1,studies)
  RE<-exp(samples2$V[,j]-samples2$V[,k])
  REtab[i,]=quantile(RE,probs=c(0.025,0.5,0.975))
}
RE<-exp(samples$V[,3]-samples$V[,2])
#REtab[N2+1,]
REquant=c(quantile(RE,probs=c(0.025,0.5,0.975)))
REtib<-tibble(REtab)%>%
  mutate(Vaccine='MVA-BN\n',Study=gsub(' 2','',c(gsub('MVA-BN ','',dose2_studies))))
Com_RE<-tibble(REquant)%>%
  mutate(Vaccine='MVA-BN\n',Study='Combined')

#col_scale<- qualitative_hcl(11)

Vacc_labs<-c('First Generation','MVA-BN 1-dose','MVA-BN 2-doses')
Vacc_labs<-c('1st Gen','MVA-BN\n1-Dose','MVA-BN\n2-Dose')
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
show(eff_plot)

ggsave('Output//Figures//Efficacy_comparison.png',plot=eff_plot,height=5,width = 8,dpi=600)
ggsave('Output//Figures//Efficacy_comparison.pdf',plot=eff_plot,height=5,width = 8,dpi=600)


plot_RE<-ggplot(REtib,aes(x=Vaccine,y=X2,ymin=X1,ymax=X3,shape=Study,color=Study))+
  geom_point(position=position_dodge(width=0.5),size=2,alpha=1)+
  geom_errorbar(position=position_dodge(width=0.5),alpha=1)+
  geom_point(aes(x='MVA-BN\n',y=REquant[2],shape='Combined',color='Combined'),size=2)+
  geom_errorbar(aes(x='MVA-BN\n',ymin=REquant[1],ymax=REquant[3],color='Combined'))+
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
        legend.text = element_text(size=10))
show(plot_RE)
ggsave('Output/Figures/Relative Efficacy.png',height=5,width=4)

plot_efficacy<-ggarrange(eff_plot,plot_RE,labels=c('A','B'),ncol=2,common.legend = TRUE,widths=c(5,3),legend='right')
show(plot_efficacy)
ggsave('Output/Figures/Manuscript_efficacy.jpg',
       plot=plot_efficacy,
       height=4,width = 6,dpi=1000)
ggsave('Output/Figures/Manuscript_efficacy.eps',
       plot=plot_efficacy,
       height=4,width = 6,dpi=1000)

Efficacies_df=data.frame(Effs)
colnames(Efficacies_df)<-levels(df_data$Indiv)[1:(nlevels(df_data$Indiv)-1)]
write.csv(Efficacies_df,'Output/data/efficacies_vaccines.csv')



# Fit universal efficacy model --------------------------------------------
data_one <- list(
  J = nrow(df_vacc),
  Num_pop = df_data$Pop,
  num_cases=df_data$Case,
  num_groups=nlevels(df_data$group),
  group_index=as.integer(df_data$group),
  num_vaccines=1,
  Vaccine_ind=rep(1,nrow(df_vacc)),
  num_vacc_studies=nlevels(df_data$index),
  vacc_study_ind=as.integer(df_data$index)
)

eff_fit_one <- stan(
  file = "Stan/Efficacy_regress.stan",  # Stan program
  data = data_one,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 10000,          # number of warmup iterations per chain
  iter = 40000,            # total number of iterations per chain
  seed=10,
  cores=4,
  control = list(max_treedepth=12,stepsize=0.5,adapt_delta=0.9)
)

eff_fit_one_pooled <- stan(
  file = "Stan/Efficacy_pooled.stan",  # Stan program
  data = data_one,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 10000,          # number of warmup iterations per chain
  iter = 40000,            # total number of iterations per chain
  seed=10,
  cores=4,
  control = list(max_treedepth=12,stepsize=0.5,adapt_delta=0.9)
)

samples<-rstan::extract(eff_fit_one)
Effs<-(samples$p_s)
quants<-t(apply(Effs,2,quantile,probs=c(0.025,0.5,0.975)))
df_output<-data.frame(quants)
df_output$Study=levels(df_data$index)

plot1<-ggplot(df_output,aes(x=Study,y=X50.))
plot1<-plot1+geom_point(position=position_dodge(width=0.5))
plot1<-plot1+geom_errorbar(aes(ymin=X2.5.,ymax=X97.5.), position=position_dodge(width=0.5))
show(plot1)

loo1<-waic(extract_log_lik(eff_fit))
loo2<-waic(extract_log_lik(eff_fit_one))
loo3<-waic(extract_log_lik(eff_fit_individual))
loo4<-waic(extract_log_lik(eff_fit_pooled))
loo5<-waic(extract_log_lik(eff_fit_one_pooled))
print(loo_compare(loo1,loo2,loo3,loo4,loo5))

elppd1<-elpd(extract_log_lik(eff_fit))
elppd2<-elpd(extract_log_lik(eff_fit_one))
elppd3<-elpd(extract_log_lik(eff_fit_individual))
elppd4<-elpd(extract_log_lik(eff_fit_pooled))
elppd5<-elpd(extract_log_lik(eff_fit_one_pooled))
print(loo_compare(elppd1,elppd2,elppd3,elppd4,elppd5))
