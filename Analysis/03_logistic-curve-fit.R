# -------------------------------------------------------------------------
#'This script combines the immunogenicity data and the effectiveness data 
#'using a logistic model
#
# -------------------------------------------------------------------------

# Import Efficacy Data ----------------------------------------------------
df_eff<-read_csv('Output/Data/Eff.csv',show_col_types = FALSE) %>%
    mutate(index=Study)%>%
    replace(is.na(.),0)%>%
    mutate_at(c('index','group','Vacc_dose'),as.factor)

# Import Immunogenicity data ----------------------------------------------

datatab <- read_csv("Output/Data/ELISA_full_set_update.csv",show_col_types = FALSE)%>%
  mutate(ELISA=replace_na(ELISA,0))%>%
  mutate_at(c('Vaccine','Study','Formulation','Dose'),as.factor)


# Fit the logistic model to both datasets ---------------------------------

#Generate random sample used for evaluating expectations quickly
sample_size=100
sample=rnorm(sample_size)

# Define the model data
data <- list(
  N = nrow(datatab),
  y_n = log10(datatab$GMT),
  s_n= datatab$std,
  n_n=datatab$Number,
  num_studies=nlevels(datatab$Study),
  num_form=nlevels(datatab$Formulation),
  study_ind=as.integer(datatab$Study),
  dose_ind=as.integer(datatab$Dose),
  form_ind=as.integer(datatab$Formulation),
  J = nrow(df_eff),
  Num_pop = df_eff$Pop,
  num_cases=df_eff$Case,
  num_groups=nlevels(df_eff$group),
  group_index=as.integer(df_eff$group),
  num_vaccines=nlevels(df_eff$Vacc_dose)-1,
  Vaccine_ind=as.integer(df_eff$Vacc_dose),
  num_vacc_studies=nlevels(df_eff$index),
  vacc_study_ind=as.integer(df_eff$index),
  sample_size=sample_size,
  sample=sample
)

#Fit the logistic model
full_fit <- stan(
  file = "Stan/lofgistic_model.stan",  # Stan program
  data = data,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = warmup,          # number of warmup iterations per chain
  iter = sampling,            # total number of iterations per chain
  seed=10,
  cores=4,
  control = list(max_treedepth=12,stepsize=0.5,adapt_delta=0.9)
)


# Analyse model outputs ---------------------------------------------------

#Extract samples
samples<-rstan::extract(full_fit)
k <- samples$k
A <- samples$A
N <- 61
full_fit_df <- data.frame(k,A)
write.csv(full_fit_df,'Output/samples/logistic_curve_parameters.csv')

#Estimate Effectiveness for varying means with fixed standard deviation
eff_df <- data.frame(row.names = c('LCI','median','UCI'))
E <- rep(0,1000)
GMTs <- seq(1,4,length.out=N)
for (i in 1:N){
  for (j in 1:length(E)){
    E[j] <- mean(1/(1+exp(-k[j]*(0.5*sample+GMTs[i]-2.5)-A[j])))
  }
  eff_df[,i] <- quantile(E,probs = c(0.025,0.5,0.975))
}
eff_df <- as.data.frame(t(eff_df))
eff_df$GMT <- GMTs
write.csv(eff_df,'Output/Figure3/Efficacy_quantiles.csv')



