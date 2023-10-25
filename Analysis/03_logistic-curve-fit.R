
# Import Efficacy Data ----------------------------------------------------

df_eff<-read_csv('Output/Data/Eff.csv')%>%
  mutate(Study=gsub('\\d\\d-','',Study))%>%
  replace(is.na(.),0)


# Remove Problem Datasets -------------------------------------------------

remove_study <- FALSE

if (remove_study){
  df_eff <- filter(df_eff,Study!='Deputy')
}


# Convert categorical variables to factors --------------------------------

df_eff$index <- as.factor(df_eff$Study)
df_eff$group <- as.factor(df_eff$group)
df_eff$Vacc_dose <- as.factor(df_eff$Vacc_dose)

df_2dose<-filter(df_eff,Dose==2)


# Create subsets for vaccination status -----------------------------------

df_unvacc=df_eff[df_eff$Vacc_dose=='Unvaccinated',]
df_vacc=df_eff[df_eff$Vacc_dose!='Unvaccinated',]



# Import Immunogenicity data ----------------------------------------------

datatab <- read_csv("Data/ELISA_full_set_update.csv",show_col_types = FALSE)%>%
  mutate(ELISA=replace_na(ELISA,0))

#Convert strings to factors for indexing
datatab$Vaccine <- as.factor(datatab$Vaccine)
datatab$Study <- as.factor(datatab$Study)
datatab$Formulation <- as.factor(datatab$Formulation)

sample=rnorm(5000)

# Define the model data
data <- list(
  N = nrow(datatab),
  y_n = log10(datatab$GMT),
  s_n= datatab$Std,
  n_n=datatab$Number,
  num_studies=nlevels(datatab$Study),
  num_form=nlevels(datatab$Formulation),
  study_ind=as.integer(datatab$Study),
  dose_ind=as.integer(datatab$Vaccine),
  form_ind=as.integer(datatab$Formulation),
  J = nrow(df_eff),
  Num_pop = df_eff$Pop,
  num_cases=df_eff$Case,
  num_groups=nlevels(df_eff$group),
  group_index=as.integer(df_eff$group),
  num_vaccines=nlevels(df_eff$Vacc_dose)-1,
  num_2dose=nrow(df_2dose),
  Vaccine_ind=as.integer(df_eff$Vacc_dose),
  num_vacc_studies=nlevels(df_eff$index),
  vacc_study_ind=as.integer(df_eff$index),
  sample=sample
)



# Fit the logistic model --------------------------------------------------

full_fit <- stan(
  file = "Stan/efficacy_model.stan",  # Stan program
  data = data,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 2000,          # number of warmup iterations per chain
  iter = 7000,            # total number of iterations per chain
  seed=10,
  cores=4,
  control = list(max_treedepth=12,stepsize=0.5,adapt_delta=0.9)
)


# Analyse model outputs ---------------------------------------------------

samples<-rstan::extract(full_fit)
k <- samples$k
A <- samples$A
N <- 61
full_fit_df <- data.frame(k,A)
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

write.csv(full_fit_df,'Output/data/full_fit_parameters.csv')
write.csv(eff_df,'Output/data/Efficacy_quantiles.csv')


# Plot Efficacy Quantiles -------------------------------------------------

plot<-ggplot(data=eff_df,aes(x=GMT,y=median)) +
  geom_line() +
  geom_line(aes(y=LCI)) +
  geom_line(aes(y=UCI))
show(plot)


# Fit Saturated Model -----------------------------------------------------

sat_fit <- stan(
  file = "Stan/Saturated_model.stan",  # Stan program
  data = data,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 2000,          # number of warmup iterations per chain
  iter = 7000,            # total number of iterations per chain
  seed=10,
  cores=4,
  control = list(max_treedepth=12,stepsize=0.2,adapt_delta=0.95)
)


# Fit Flat Model ----------------------------------------------------------

flat_fit <- stan(
  file = "Stan/Flat_model.stan",  # Stan program
  data = data,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 2000,          # number of warmup iterations per chain
  iter = 7000,            # total number of iterations per chain
  seed=10,
  cores=4,
  control = list(max_treedepth=12,stepsize=0.5,adapt_delta=0.9)
)


# Compare Models ----------------------------------------------------------

loo1 <- waic(extract_log_lik(sat_fit))
loo2 <- waic(extract_log_lik(full_fit))
loo3 <- waic(extract_log_lik(flat_fit))
print(loo_compare(loo1,loo2,loo3))

elppd1<-elpd(extract_log_lik(sat_fit))
elppd2<-elpd(extract_log_lik(full_fit))
elppd3<-elpd(extract_log_lik(flat_fit))
print(loo_compare(elppd1,elppd2,elppd3))

# Analyse Saturated model output ------------------------------------------

samples<-rstan::extract(sat_fit)
k <- samples$k
A <- samples$A
M <- samples$M
write_csv(data.frame(k,A,M),"Output//Data//sat_fit.csv")
N <- 61
eff_df <- data.frame(row.names = c('LCI','median','UCI'))
GMTs <- seq(1,4,length.out=N)
for (i in 1:N){
  for (j in 1:length(E)){
    E[j] <- mean(1/(M[j]+exp(-k[j]*(0.5*sample+GMTs[i]-2.5)-A[j])))
  }
  eff_df[,i] <- quantile(E,probs = c(0.025,0.5,0.975))
}
eff_df <- as.data.frame(t(eff_df))
eff_df$GMT <- GMTs

write.csv(eff_df,'Output/data/Efficacy_sat_quantiles.csv')

plot<-ggplot(data=eff_df,aes(x=GMT,y=median)) +
  geom_line() +
  geom_line(aes(y=LCI)) +
  geom_line(aes(y=UCI))
show(plot)


# Model Fit ---------------------------------------------------------------
indx<-0
tab<-tribble(~LOR,~Group,~Vacc_dose,~Study)
for (i in 1:3){
  for (j in 1:length(levels(df_eff$group))){
    group_s<-levels(df_eff$group)[j]
    vacc_dose=levels(df_eff$Vacc_dose)[i]
    selection<-filter(df_eff,group==group_s , Vacc_dose==vacc_dose)
    if (nrow(selection)==1){
      cases<-selection$Case
      non_cases<-selection$Pop-selection$Case
      odds<-(cases+1)/(non_cases+1)
      study<-as.numeric(selection$index)
      unvacc<-filter(df_eff,group==group_s,Vacc_dose=='Unvaccinated')
      cases<-unvacc$Case
      non_cases<-unvacc$Pop-unvacc$Case
      un_odds<-(cases+1)/(non_cases+1)
      LOR=log(odds/un_odds)
      indx<-indx+1
      tab<-add_row(tab,LOR=LOR,
                   Group=group_s,
                   Vacc_dose=i,
                   Study=study)
    }
  }
}

indx<-0
tab2<-tribble(~LOR,~Vacc_dose,~Study)
for (i in 1:3){
  for (j in 1:length(levels(df_eff$index))){
    study<-levels(df_eff$index)[j]
    vacc_dose=levels(df_eff$Vacc_dose)[i]
    selection<-filter(df_eff,index==study , Vacc_dose==vacc_dose)
    if (nrow(selection)>0){
      cases<-selection$Case
      non_cases<-selection$Pop-selection$Case
      odds<-(cases+1)/(non_cases+1)
      unvacc<-filter(df_eff,index==study,Vacc_dose=='Unvaccinated')
      cases<-unvacc$Case
      non_cases<-unvacc$Pop-unvacc$Case
      un_odds<-(cases+1)/(non_cases+1)
      LOR=mean(log(odds/un_odds))
      indx<-indx+1
      tab2<-add_row(tab2,LOR=LOR,
                   Vacc_dose=i,
                   Study=j)
    }
  }
}

samples<-rstan::extract(full_fit)

V<-samples$V
ps<-samples$p_s
N<-length(V[,1])
S<-nrow(tab2)
R<-0

for (i in 1:N){
  EOR<-rep(0,S)
  for (j in 1:S){
    EOR[j]<-V[i,tab2[j,]$Vacc_dose]+ps[i,tab2[j,]$Study]
  }
  res<-EOR-tab2[,1]
  R[i]=var(EOR)/(var(EOR)+var(res))
}

samples2<-rstan::extract(flat_fit)
V<-samples2$V
ps<-samples2$p_s
N<-length(V)
S<-nrow(tab2)
R2<-0

for (i in 1:N){
  EOR<-rep(0,S)
  for (j in 1:S){
    EOR[j]<-V[i]+ps[i,tab2[j,]$Study]
  }
  res<-EOR-tab2[,1]
  R2[i]=var(EOR)/(var(EOR)+var(res))
}


# Test Statistics ---------------------------------------------------------

#Calculating a test statistic for the full model.
tab<-tribble(~LOR,~vacc_ind,~unvacc_ind, ~unvacc_pop, ~vacc_pop)
df_eff$position<-seq(1,nrow(df_eff),by=1)
for (j in 1:length(levels(df_eff$group))){
    group_s<-levels(df_eff$group)[j]
    vacc_dose=levels(df_eff$Vacc_dose)[1]
    selection<-filter(df_eff,group==group_s , Vacc_dose==vacc_dose)
    if (nrow(selection)==1){
      cases<-selection$Case
      vacc_ind<-selection$position
      case_pop<-selection$Pop
      non_cases<-case_pop-selection$Case
      case_ind<-selection
      odds<-(cases+1)/(non_cases+1)
      study<-as.numeric(selection$index)
      unvacc<-filter(df_eff,group==group_s,Vacc_dose=='Unvaccinated')
      cases<-unvacc$Case
      unvacc_pop<-unvacc$Pop
      non_cases<-unvacc_pop-selection$Case
      unvacc_ind<-unvacc$position
      un_odds<-(cases+1)/(non_cases+1)
      LOR=log(odds/un_odds)
      tab<-add_row(tab,LOR=LOR,
                   vacc_ind=vacc_ind,
                   unvacc_ind=unvacc_ind,
                   vacc_pop=case_pop,
                   unvacc_pop=unvacc_pop)
    }
}

samples<-rstan::extract(full_fit)
sims<-samples$sim_cases
Tdata<-mean(tab$LOR)
Tdatamin<-min(tab$LOR)
Trep<-0
Trepmin<-0
Trepmax<-0
for (i in 1:N){
  LBOR<-0
  for (j in 1:nrow(tab)){
    vacc_case<-sims[i,tab[j,]$vacc_ind]
    odds_vacc<-(vacc_case+1)/(tab[j,]$vacc_pop-vacc_case+1)
    unvacc_case<-sims[i,tab[j,]$unvacc_ind]
    odds_unvacc<-(unvacc_case+1)/(tab[j,]$unvacc_pop-unvacc_case+1)
    LBOR[j]<-log(odds_vacc/odds_unvacc)
  }
  Trep[i]<-mean(LBOR)
  Trepmin[i]<-min(LBOR)
  Trepmax[i]<-max(LBOR)
}


#Calculating the test statistic in the flat model case
samples2<-rstan::extract(flat_fit)

tab<-tribble(~LOR,~vacc2_ind,~unvacc2_ind, ~unvacc2_pop, ~vacc2_pop,~vacc1_ind,~unvacc1_ind, ~unvacc1_pop, ~vacc1_pop,~weight)
df_eff$position<-seq(1,nrow(df_eff),by=1)

vacc_dose=levels(df_eff$Vacc_dose)[3]
selection<-filter(df_eff, Vacc_dose==vacc_dose)

for (j in 1:length(levels(df_eff$group))){
  group_s<-levels(df_eff$group)[j]
  selection2<-filter(selection,group==group_s)
  if (nrow(selection2)==1){
    study<-selection2$Study
    cases<-selection2$Case
    vacc2_ind<-selection2$position
    case2_pop<-selection2$Pop
    non_cases<-case2_pop-selection2$Case
    odds2<-(cases+1)/(non_cases+1)
    unvacc2<-filter(df_eff,group==group_s,Vacc_dose=='Unvaccinated')
    cases<-unvacc2$Case
    unvacc2_pop<-unvacc2$Pop
    non_cases<-unvacc2_pop-unvacc2$Case
    unvacc2_ind<-unvacc2$position
    un_odds2<-(cases+1)/(non_cases+1)
    unvacc<-filter(df_eff,group==group_s,Vacc_dose=='MVA-BN_1')
    if (nrow(unvacc)==1){
      cases<-unvacc$Case
      case1_pop<-unvacc$Pop
      non_cases<-case1_pop-unvacc$Case
      vacc1_ind<-unvacc$position
      unvacc1_ind<-unvacc2_ind
      unvacc1_pop<-unvacc2_pop
      odds1<-(cases+1)/(non_cases+1)
      LOR=log(odds2/odds1)
    }else{
      unvacc<-filter(df_eff,Study==study,Vacc_dose=='MVA-BN_1')
      cases<-unvacc$Case
      case1_pop<-unvacc$Pop
      non_cases<-unvacc_pop-unvacc$Case
      vacc1_ind<-unvacc$position
      odds1<-(cases+1)/(non_cases+1)
      group_n<-unvacc$group
      unvacc1<-filter(df_eff,group==group_n,Vacc_dose=='Unvaccinated')
      cases<-unvacc1$Case
      unvacc1_pop<-unvacc1$Pop
      non_cases<-unvacc1_pop-unvacc1$Case
      unvacc1_ind<-unvacc1$position
      un_odds1<-(cases+1)/(non_cases+1)
      
      LOR=log(odds2/un_odds2/(odds1/un_odds1))
    }
    weighting<-filter(df_eff,Study==study,Vacc_dose=='MVA-BN_1')
    weight<-1/nrow(weighting)
    tab<-add_row(tab,LOR=LOR,
                 vacc1_ind=vacc1_ind,
                 unvacc1_ind=unvacc1_ind,
                 vacc1_pop=case1_pop,
                 unvacc1_pop=unvacc1_pop,
                 vacc2_ind=vacc2_ind,
                 unvacc2_ind=unvacc2_ind,
                 vacc2_pop=case2_pop,
                 unvacc2_pop=unvacc2_pop,
                 weight=weight)
  }
}

sims<-samples2$sim_cases
Tdata<-weighted.mean(tab$LOR,tab$weight)
Tmean<-mean(tab$LOR)
Tdatamin<-min(tab$LOR)
Trep<-0
Trepmin<-0
Trepmax<-0
for (i in 1:N){
  LBOR<-0
  for (j in 1:nrow(tab)){
    vacc2_case<-sims[i,tab[j,]$vacc2_ind]
    odds2_vacc<-(vacc2_case+1)/(tab[j,]$vacc2_pop-vacc2_case+1)
    unvacc2_case<-sims[i,tab[j,]$unvacc2_ind]
    odds2_unvacc<-(unvacc2_case+1)/(tab[j,]$unvacc2_pop-unvacc2_case+1)
    vacc1_case<-sims[i,tab[j,]$vacc1_ind]
    odds1_vacc<-(vacc1_case+1)/(tab[j,]$vacc1_pop-vacc1_case+1)
    unvacc1_case<-sims[i,tab[j,]$unvacc1_ind]
    odds1_unvacc<-(unvacc1_case+1)/(tab[j,]$unvacc1_pop-unvacc1_case+1)
    LBOR[j]<-log(odds2_vacc/odds2_unvacc/(odds1_vacc/odds1_unvacc))
  }
  Trep[i]<-weighted.mean(LBOR,tab$weight)
  Trepmin[i]<-min(LBOR)
  Trepmax[i]<-max(LBOR)
}
