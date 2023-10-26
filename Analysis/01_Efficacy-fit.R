# -------------------------------------------------------------------------
#'Code used to estimate the effectiveness of the different vaccines including
#'comparison of one and two doses of MVA-BN
#
# -------------------------------------------------------------------------



# Import Data -------------------------------------------------------------

#Load the Effectiveness data
df_data<-read_csv('Output/Data/Eff.csv') %>%
  mutate(Indiv=if_else(Dose==0,'Unvaccinated',paste(Vaccine,Study,Dose)),
         Vaccinated=if_else(Dose==0,'Unvaccinated','Vaccinated'),
         index=Study)%>%
  replace(is.na(.),0)%>%
  mutate_at(c('Vaccinated','Indiv','index','group','Vacc_dose'),as.factor)

#Load efficacy study data
df_studies<-read_csv('Data/Efficacy_studies.csv')

#Create Data Frame for 2 dose dataset
df_2dose=filter(df_data,Dose==2)%>%
  mutate(Indiv=fct_drop(Indiv))


# Analyse combined estimates of effectiveness -----------------------------

#Create Data list for Stan
data <- list(
    J = nrow(df_data),
    Num_pop = df_data$Pop,
    num_cases=df_data$Case,
    num_groups=nlevels(df_data$group),
    group_index=as.integer(df_data$group),
    num_vaccines=nlevels(df_data$Vacc_dose)-1,
    Vaccine_ind=as.integer(df_data$Vacc_dose),
    num_vacc_studies=nlevels(df_data$index),
    vacc_study_ind=as.integer(df_data$index)
)
 
#Run Stan algorithm for sampling 
eff_fit <- stan(
    file = "Stan/Effectiveness_combined.stan",  # Stan program
    data = data,    # named list of data
    chains = 4,             # number of Markov chains
    warmup = warmup,          # number of warmup iterations per chain
    iter = sampling,            # total number of iterations per chain
    seed=10,
    cores=4,
    control = list(max_treedepth=12,stepsize=0.5,adapt_delta=0.9)
)

#Extract samples from RStan Object
samples<-rstan::extract(eff_fit)
#Determine posterior samples of effectiveness and summary statistics
Effs<-100*(1-exp(samples$V))
quants<-t(apply(Effs,2,quantile,probs=c(0.025,0.5,0.975)))
#Create Data Frame for combined estimates
df_output_combined<-data.frame(quants)
mean_quant=quantile(samples$E,probs=c(0.025,0.5,0.975))
df_output_combined$Vaccine=c('First Generation','MVA-BN 1-dose','MVA-BN 2-dose')
df_output_combined$Study='Combined'
write_csv(df_output_combined,'Output/Figure1/combined_eff.csv')

#Calculate posterior samples and summaries for RE_combined
RE<-exp(samples$V[,3]-samples$V[,2])
REquant=c(quantile(RE,probs=c(0.025,0.5,0.975)))
Com_RE<-tibble(REquant)%>%
  mutate(Vaccine='MVA-BN\n',
         Study='Combined',
         y=c('X1','X2','X3'))%>%
  pivot_wider(names_from = 'y',values_from = REquant)%>%
  write_csv('Output/Figure1/RE_comb.csv')

#Save Posterior samples of effectiveness.
Eff_df<-data.frame(Effs)
colnames(Eff_df)<-c('First Generation','MVA-BN 1-dose','MVA-BN 2-dose')
Eff_df$`Risk Ratio`<-RE
write_csv(Eff_df,"Output//Samples//Combined-Effectiveness-sample.csv")


# Analyse Data for individual Effectiveness -------------------------------

#Create Data list for Stan
data2 <- list(
  J = nrow(df_data),
  Num_pop = df_data$Pop,
  num_cases=df_data$Case,
  num_groups=nlevels(df_data$group),
  group_ind=as.integer(df_data$group),
  num_vaccines=nlevels(df_data$Indiv)-1,
  Vaccine_ind=as.integer(df_data$Indiv)
)

#Run RStan Model
eff_fit_individual <- stan(
  file = "Stan/Effectiveness_individual.stan",  # Stan program
  data = data2,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = warmup,          # number of warmup iterations per chain
  iter = sampling,            # total number of iterations per chain
  seed=10,
  cores=4,
  control = list(max_treedepth=12,stepsize=0.5,adapt_delta=0.9)
)

#Extract Samples
samples2<-rstan::extract(eff_fit_individual)
#Calculate per study effectiveness
Effs<-100*(1-exp(samples2$V))
quants<-t(apply(Effs,2,quantile,probs=c(0.025,0.5,0.975)))
df_output_indiv<-data.frame(quants)
df_output_indiv$Study<-gsub(' \\d','',c(levels(df_data$Indiv)[1:15]))
df_output_indiv$Study<-gsub('\\w*\\W','',df_output_indiv$Study)
df_output_indiv$Vaccine<-c(gsub('MVA-BN1','MVA-BN 1-dose',gsub(' [A-Z,a-z]*','',levels(df_data$Indiv)[1:15])))
df_output_indiv$Vaccine<-gsub('Dryvax1','First Generation',df_output_indiv$Vaccine)
df_output_indiv$Vaccine<-as.numeric(as.factor(gsub('MVA-BN2','MVA-BN 2-dose',df_output_indiv$Vaccine)))
write_csv(df_output_indiv,'Output/Figure1/Indiv_eff.csv')

#Calculate per study Risk Ratio for 1 to 2 doses MVA-BN
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
REtib<-tibble(REtab)%>%
  mutate(Vaccine='MVA-BN\n',
         Study=gsub(' 2','',c(gsub('MVA-BN ','',dose2_studies))))%>%
  write_csv('Output/Figure1/REtib.csv')


#Save samples for individual effectiveness
Efficacies_df=data.frame(Effs)
colnames(Efficacies_df)<-levels(df_data$Indiv)[1:(nlevels(df_data$Indiv)-1)]
write.csv(Efficacies_df,'Output/samples/individual-effectiveness.csv')
