# -------------------------------------------------------------------------
#' This script models the decay of antibody titers overtime and compares
#' the results of different decay models.
#
# -------------------------------------------------------------------------

# Load antibody data ------------------------------------------------------

datatab <- read_csv("Output/Data/ELISA_decay.csv") %>%
  mutate(Study=ifelse(Study==23,5,Study)) %>%
  unite(Group,Dose2,`Dose 3`,Vaccinia_EXP,sep='_',remove=FALSE)%>%
  mutate_at(c('Group','Formulation','Study'),as.factor)

# Model the decay using Stan ----------------------------------------------

#Define Data
data <- list(
  N = nrow(datatab),
  y_n = log10(datatab$GMT),
  time=datatab$Week,
  s_n=datatab$std,
  n_n=datatab$Number,
  num_doses=nlevels(datatab$Group),
  dose_ind=as.integer(datatab$Group),
  num_form = 2, 
  num_studies= nlevels(datatab$Study),
  study_ind = as.integer(datatab$Study),
  form_ind=as.integer(datatab$Formulation)
)

#fit1 Best fit model with decay rates fixed across regimens 
fit1 <- stan(
  file = "Stan/Decay_Model.stan",  # Stan program
  data = data,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = warmup,          # number of warmup iterations per chain
  iter = sampling,            # total number of iterations per chain
  seed=10,
  cores=4,
  control = list(max_treedepth=14,stepsize=1),
  algorithm = 'NUTS'
)

#fit2 Two phase model with varying decay rates
fit2 <- stan(
  file = "Stan/Decay_Model_large.stan",  # Stan program
  data = data,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = warmup,          # number of warmup iterations per chain
  iter = sampling,            # total number of iterations per chain
  seed=10,
  cores=4,
  control = list(max_treedepth=14,stepsize=1),
  algorithm = 'NUTS'
)

#fit3 One phase model with consistent decay rate
fit3 <- stan(
  file = "Stan/Decay_Model_linear.stan",  # Stan program
  data = data,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = warmup,          # number of warmup iterations per chain
  iter = sampling,            # total number of iterations per chain
  seed=10,
  cores=4,
  control = list(max_treedepth=14,stepsize=1),
  algorithm = 'NUTS'
)

#fit4 Single Phase decay model with different decay rates
fit4 <- stan(
  file = "Stan/Decay_Model_linear_large.stan",  # Stan program
  data = data,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = warmup,          # number of warmup iterations per chain
  iter = sampling,            # total number of iterations per chain
  seed=10,
  cores=4,
  control = list(max_treedepth=14,stepsize=1),
  algorithm = 'NUTS'
)

#fit5 No decay model
fit5 <- stan(
  file = "Stan/Decay_Model_Flat.stan",  # Stan program
  data = data,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = warmup,          # number of warmup iterations per chain
  iter = sampling,            # total number of iterations per chain
  seed=10,
  cores=4,
  control = list(max_treedepth=14,stepsize=1),
  algorithm = 'NUTS'
)


# Extract Data from main fit ----------------------------------------------

#Extract data
output<-rstan::extract(fit1)
#extract decay rates
D1<-output$decay1
D2<-output$decay2
#Extract proportion of short-term antibodies
tc<-output$tc
#extract initial mean and sd
mu<-output$mu_d
sd<-output$sigma_d

#Save Samples
write_csv(data.frame(mu),"Output\\Samples\\Initial_GMT_decay.csv")
write_csv(data.frame(tc),"Output\\Samples\\Decay_proportions.csv")
write_csv(data.frame(sd),"Output\\Samples\\Decay_sd.csv")
decay_data<-tibble(data.frame(D1,D2))%>%
  write_csv('Output/Samples/Decay_rates.csv')

# Generate decay titer quantiles ------------------------------------------
times=seq(0,104)
pred_decay<-data.frame()
q=array(dim=c(105,12))
for (j in 1:105){
  titer1<-10**mu[,1]*(tc[,1]*exp(-D1*(times[j]))+(1-tc[,1])*exp(-D2*(times[j])))
  titer2<-10**mu[,5]*(tc[,5]*exp(-D1*(times[j]))+(1-tc[,5])*exp(-D2*(times[j])))
  titer3<-10**mu[,7]*(tc[,7]*exp(-D1*(times[j]))+(1-tc[,7])*exp(-D2*(times[j])))
  titer4<-10**mu[,9]*(tc[,9]*exp(-D1*(times[j]))+(1-tc[,9])*exp(-D2*(times[j])))
  q[j,1:3]<-quantile(titer1,probs=c(0.025,0.5,0.975))
  q[j,4:6]<-quantile(titer2,probs=c(0.025,0.5,0.975))
  q[j,7:9]<-quantile(titer3,probs=c(0.025,0.5,0.975))
  q[j,10:12]<-quantile(titer4,probs=c(0.025,0.5,0.975))
}
q_df<-data.frame(q)
q_df$Time<-times/52
write_csv(q_df,'Output/Samples/GMT_decay_quants.csv')

#Reorganise data
q1<-select(q_df,Time,X1,X2,X3)%>%
  rename(LCI=X1,GMT=X2,UCI=X3)%>%
  mutate(Group=1)
q2<-select(q_df,Time,X4,X5,X6)%>%
  rename(LCI=X4,GMT=X5,UCI=X6)%>%
  mutate(Group=5)
q3<-select(q_df,Time,X7,X8,X9)%>%
  rename(LCI=X7,GMT=X8,UCI=X9)%>%
  mutate(Group=7)
q4<-select(q_df,Time,X10,X11,X12)%>%
  rename(LCI=X10,GMT=X11,UCI=X12)%>%
  mutate(Group=9)
q_df<-rbind(q1,q2,q3,q4)%>%
  write_csv('Output/Samples/Decay_quantiles.csv')

# Generate Summaries for Figure 4 -----------------------------------------

#extract percentile summaries
percent_quant<-as_tibble(t(apply((1-tc)*100,2,quantile,
                                 probs=c(0.025,0.5,0.975))))
names(percent_quant)<-c('p_LCI','percent','p_UCI')
percent_quant$Group<-levels(datatab$Group)

#Extract titer summaries
titer_quant<-data.frame(t(apply(10**mu,2,quantile,
                                probs=c(0.025,0.5,0.975))))
names(titer_quant)<-c('I_LCI','I_GMT','I_UCI')
titer_quant$Group<-levels(datatab$Group)

#Calculate titers at 1 year post vaccination
Long_titer<- array(dim=c(20000,9))
for (i in 1:nlevels(datatab$Group)){
  Long_titer[,i]=(1-tc[,i])*(10**mu[,i])*exp(-D2*52)
}
write_csv(data.frame(Long_titer),"Output\\Samples\\GMT_1year.csv")

#Extract long-titer summaries
long_titer_quant<-as_tibble(t(apply(Long_titer,2,quantile,
                                    probs=c(0.025,0.5,0.975))))
names(long_titer_quant)<-c('LCI','GMT','UCI')
long_titer_quant$Group<-levels(datatab$Group)

#Group summaries together
summary_decay<-full_join(percent_quant,titer_quant,by=join_by(Group))%>%
  full_join(long_titer_quant,by=join_by(Group))%>%
  separate_wider_delim(cols=Group,delim='_',
                       names=c('Dose2','Dose3','Vaccinia Experience'))%>%
  mutate_at(c('Dose2','Dose3'),as.numeric)%>%
  mutate(`Number of MVA Doses`=factor((Dose2>0)+(Dose3>0)))%>%
  mutate(`Final Dose Day`=if_else(Dose2==0,1,
                                  (abs(Dose2-Dose3)+(`Dose3`+Dose2))/2))%>%
  mutate(Dose2=factor(Dose2))%>%
  filter(`Vaccinia Experience`==0)%>%
  mutate(alpha_val=c('1','0','0','1','1','1','1'))%>%
  write_csv('output/Figure4/summary_decay.csv')



# Compare Models using WAIC -----------------------------------------------

#Compare WAIC
loo1<-waic(extract_log_lik(fit1))
loo2<-waic(extract_log_lik(fit2))
loo3<-waic(extract_log_lik(fit3))
loo4<-waic(extract_log_lik(fit4))
loo5<-waic(extract_log_lik(fit5))
print(loo_compare(loo1,loo2,loo3,loo4,loo5))

#Compare elppd
elppd1<-elpd(extract_log_lik(fit1))
elppd2<-elpd(extract_log_lik(fit2))
elppd3<-elpd(extract_log_lik(fit3))
elppd4<-elpd(extract_log_lik(fit4))
elppd5<-elpd(extract_log_lik(fit5))
print(loo_compare(elppd1,elppd2,elppd3,elppd4,elppd5))


# Calculate time till 1 dose ----------------------------------------------
mean_df<-read_csv('Output/Samples/antibody_means.csv')
mu_h<-mean_df$`Historic Vaccination`

samples<-rstan::extract(fit1)
k1<-samples$decay1
k2<-samples$decay2
f3<-samples$tc[,7]
f2<-samples$tc[,5]
f1<-samples$tc[,1]
x3<-samples$mu_d[,7]
x2<-samples$mu_d[,5]
x1<-samples$mu_d[,1]
gmt0<-x2-x1
gmt2<-x3-x1
gmt1<-x2-mu_h
t<-0
t2<-0
t3<-0
for (i in 1:length(f1)){
  t[i]<-nleqslv(10,function (x) decay_titer(x,gmt0[i],k1[i],k2[i],f2[i]))$x
  t3[i]<-nleqslv(10,function (x) decay_titer(x,gmt2[i],k1[i],k2[i],f3[i]))$x
  t2[i]<-nleqslv(10,function (x) decay_titer(x,gmt1[i],k1[i],k2[i],f2[i]))$x
}

times<-tibble(t,t2,t3)%>%
  write_csv('Output/Samples/Decay_time_to_1dose.csv')
  


# Fold Differences --------------------------------------------------------

x<-samples$mu_d
f<-samples$tc
k1<-samples$decay1
k2<-samples$decay2
FD7_28<-10^(x[,5]-x[,8])
FD700_28<-10^(x[,9]-x[,5])
FD700_700<-10^(x[,7]-x[,9])
FD700_1<-10^(x[,7]-x[,1])

gmt11<-decay_titer(52,x[,1],k1,k2,f[,1])
gmt12<-decay_titer(52,x[,5],k1,k2,f[,5])
gmt13<-decay_titer(52,x[,7],k1,k2,f[,7])
gmt1730<-decay_titer(52,x[,9],k1,k2,f[,9])
FD1_28_0<-10^(gmt11-gmt12)
FD1_700_28<-10^(gmt13-gmt12)
FD1_700_700<-10^(gmt13-gmt1730)

fds<-tibble(data.frame(FD1_28_0,FD1_700_28,FD1_700_700))%>%
  write_csv('output/samples/fold_drops_decay1.csv')

