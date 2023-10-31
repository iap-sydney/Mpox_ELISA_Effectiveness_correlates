# -------------------------------------------------------------------------
#'Script for analysing the immunogenicity data achieved from vaccination
#
# -------------------------------------------------------------------------


# Import Data -------------------------------------------------------------


#Read the immunogenicity data from the csv file.
datatab <- read_csv("Output/Data/ELISA_full_set_update.csv",show_col_types = FALSE)%>%
  mutate(ELISA=replace_na(ELISA,0),
         col_group=paste(Dose,Formulation))%>%
  mutate_at(c('col_group','Formulation','Study','ELISA','Dose'),as.factor)


#Obtain mean estimates
meandata<-datatab%>%
  group_by(Dose)%>%
  summarise(meandose=mean(loggmt1))


# Fit Immunogenicity Model ------------------------------------------------

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
  file = "Stan/Immunogenicity_model.stan",  # Stan program
  data = data,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = warmup,          # number of warmup iterations per chain
  iter = sampling,            # total number of iterations per chain
  seed=10,
  cores=4,
  control = list(max_treedepth=12,stepsize=0.01)
)


# Extract samples ---------------------------------------------------------

samples<-rstan::extract(fit_unpooled)

#Extract GMTs by vaccine
mean_vaccine <- samples$mu_d
mean_df <- data.frame(mean_vaccine)
colnames(mean_df)<-c('Historic Vaccination','MVA-BN 1-dose','MVA-BN 2-dose')
#Add Formulation sample
mean_df$Formulation<-samples$mu_f
#save posterior sample
write.csv(mean_df,'Output/samples/antibody_means.csv')

#Generate quantiles for dose plus formulation effect
mean_quant <- quantile(samples$mu_d[,1],probs=c(0.025,0.5,0.975))
mean_quant1 <- 10**quantile(samples$mu_d[,1],probs=c(0.025,0.5,0.975))
mean_quant2 <- 10**quantile(samples$mu_d[,2],probs=c(0.025,0.5,0.975))
mean_quant3 <- 10**quantile(samples$mu_d[,3],probs=c(0.025,0.5,0.975))
mean_quant4 <- 10**quantile(samples$mu_d[,2]+samples$mu_f,probs=c(0.025,0.5,0.975))
mean_quant5 <- 10**quantile(samples$mu_d[,3]+samples$mu_f,probs=c(0.025,0.5,0.975))
#Save quantiles as csv for Figure 2
mean_quants <- as_tibble(t(data.frame(mean_quant1,mean_quant4,mean_quant2,mean_quant5,mean_quant3)))%>%
  mutate(Dose=gsub('1st Gen First Generation','1st Gen',gsub(' [A-Z]*$','',levels(datatab$col_group))))%>%
  mutate(col_group=levels(datatab$col_group))%>%
  mutate(logmin=log10(`2.5%`),logmed=log10(`50%`),logmax=log10(`97.5%`))%>%
  write_csv('Output/Figure2/Indiv_eff.csv')


# ELISA effect ------------------------------------------------------------

#Restrict data to those with known ELISA
df_ELISA <- filter(datatab,as.numeric(ELISA)> 1)

#State Data for Stan
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

#Fit Stan Model accounting for ELISA
fit_ELISA<- stan(
  file = "Stan/Immunogenicity_model.stan",  # Stan program
  data = data,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = warmup,          # number of warmup iterations per chain
  iter = sampling,            # total number of iterations per chain
  seed=10,
  cores=4,
  control = list(max_treedepth=12,stepsize=0.01)
)

#Save ELISA samples
samples_ELISA<-rstan::extract(fit_ELISA)
df_formulation <- data.frame(samples_ELISA$mu_f)
write_csv(df_formulation,'output/samples/ELISA_effect.csv')

