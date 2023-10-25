
# Load antibody data ------------------------------------------------------

datatab <- read_csv("Output/Data/ELISA_decay.csv") %>%
  mutate(Study=ifelse(Study==23,5,Study)) %>%
  unite(Group,Dose2,`Dose 3`,Vaccinia_EXP,sep='_',remove=FALSE)%>%
  mutate_at(c('Group','Formulation','Study'),as.factor)
  

# plot the raw antibody data ----------------------------------------------

plot=ggplot(datatab, aes(Week, GMT, colour = Group)) + 
  geom_point() + 
  scale_y_continuous('GMT', trans = 'log10') +
  theme_classic()
show(plot)
ggsave('Output/Figures/GMT_decay.png',width=10,height=5)


#plot2=ggplot(datatab, aes(Week, fold_drop, colour = Group)) + 
#  geom_point() +
#  scale_y_continuous(trans = 'log10')
#show(plot2)

#ggsave('Output/Figures/fold_drop_decay.png',width=10,height=5)


# Define data for Stan Model ----------------------------------------------

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


# Fit Stan Model ----------------------------------------------------------

#fit1 Best fit model with decay rates fixed across regimens 
fit1 <- stan(
  file = "Stan/Stan_Model.stan",  # Stan program
  data = data,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 15000,          # number of warmup iterations per chain
  iter = 20000,            # total number of iterations per chain
  seed=10,
  cores=4,
  control = list(max_treedepth=14,stepsize=1),
  algorithm = 'NUTS'
)

#fit2 Two phase model with varying decay rates
fit2 <- stan(
  file = "Stan/Stan_Model_large.stan",  # Stan program
  data = data,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 15000,          # number of warmup iterations per chain
  iter = 20000,            # total number of iterations per chain
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
  warmup = 15000,          # number of warmup iterations per chain
  iter = 20000,            # total number of iterations per chain
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
  warmup = 15000,          # number of warmup iterations per chain
  iter = 20000,            # total number of iterations per chain
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
  warmup = 15000,          # number of warmup iterations per chain
  iter = 20000,            # total number of iterations per chain
  seed=10,
  cores=4,
  control = list(max_treedepth=14,stepsize=1),
  algorithm = 'NUTS'
)
# Print Diagnostic plots --------------------------------------------------

pars=c('tc[1]','tc[2]','tc[3]','tc[4]','tc[5]','tc[6]','tc[7]','tc[8]','tc[9]','decay2','decay1','lp__')

trace=mcmc_trace(fit1,pars=pars)
show(trace)
post=mcmc_dens(fit1,pars=pars)
show(post)

# Display decay plots -----------------------------------------------------

output=as.data.frame(fit1)
output<-output[c(1:11,length(output))]
write_csv(output, "Output\\Data\\decay.csv")

output<-rstan::extract(fit1)
D1<-output$decay1
D2<-output$decay2
tc<-output$tc
mu<-output$mu_d
sd<-output$sigma_d
write_csv(data.frame(mu),"Output\\Data\\mu_d.csv")
write_csv(data.frame(sd),"Output\\Data\\sigma_d.csv")

Long_titer<- array(dim=c(20000,9))
for (i in 1:nlevels(datatab$Group)){
  Long_titer[,i]=(1-tc[,i])*(10**mu[,i])*exp(-D2*52)
}
write_csv(data.frame(Long_titer),"Output\\Data\\long_titer.csv")

times=seq(0,104)
pred_decay<-data.frame()
for (i in 1:nlevels(datatab$Group)){
  q=array(dim=c(105,3))
  for (j in 1:105){
    #titer<-10**(-D1[,i]*(times[j]-2)+(D1[,i]-D2)*(abs(times[j]-2-tc)+(times[j]-2-tc))/2)
    titer<-10**mu[,i]*(tc[,i]*exp(-D1*(times[j]))+(1-tc[,i])*exp(-D2*(times[j])))
    q[j,]<-quantile(titer,probs=c(0.025,0.5,0.975))
  }
  q_df<-data.frame(q)
  q_df$Time<-times
  group<-toString(levels(datatab$Group)[i])
  plot_decay<-ggplot(data=q_df,aes(x=Time)) +
    geom_line(color='red',aes(y=X2)) +
    geom_ribbon(aes(ymin=X1,ymax=X3),color='red',fill='blue',alpha=0.2) +
    geom_point(data=datatab[as.integer(datatab$Group)==i,],aes(x=Week,y=GMT,color=Study)) +
    geom_errorbar(data=datatab[as.integer(datatab$Group)==i,],aes(x=Week,ymin=LCI,ymax=UCI,color=Study)) +
    labs(y='GMT',x='Weeks',title=group)+
    scale_y_continuous(trans='log10')+
    theme_classic()
  show(plot_decay)
  save_file=glue('Output\\Figures\\{group}_decay.png')
  ggsave(save_file,plot_decay,width=10,height=5)
}


# Plot longterm decay -----------------------------------------------------

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
q_df$Time<-times
write_csv(q_df,'Output/Data/GMT_decay.csv')
plot_decay<-ggplot(data=q_df,aes(x=Time,color='black')) +
  geom_point(data=datatab[as.integer(datatab$Group)==1,],aes(x=Week,y=GMT,shape=label),color='red',alpha=0.5) +
  geom_line(data=datatab[as.integer(datatab$Group)==1,],aes(x=Week,y=GMT,shape=label),color='red',alpha=0.5) +
  geom_errorbar(data=datatab[as.integer(datatab$Group)==1,],aes(x=Week,ymin=LCI,ymax=UCI),color='red',alpha=0.5)+
  geom_point(data=datatab[as.integer(datatab$Group)==5,],aes(x=Week,y=GMT,shape=label),color='blue',alpha=0.5) +
  geom_line(data=datatab[as.integer(datatab$Group)==5,],aes(x=Week,y=GMT,shape=label),color='blue',alpha=0.5) +
  geom_errorbar(data=datatab[as.integer(datatab$Group)==5,],aes(x=Week,ymin=LCI,ymax=UCI),color='blue',alpha=0.5) +
  labs(y='GMT',x='Weeks')+
  scale_y_continuous(trans='log10')+
  theme_classic()+
  geom_line(color='black',aes(y=X2),linetype=1) +
  geom_ribbon(aes(ymin=X1,ymax=X3),color='black',alpha=0,linetype=2) +
  geom_line(color='black',aes(y=X5),linetype=1) +
  geom_ribbon(aes(ymin=X4,ymax=X6),color='black',alpha=0,linetype=2) +
  scale_shape_manual(breaks=levels(as.factor(datatab$label)),
                     values=1:nlevels(as.factor(datatab$label)),
                     labels=levels(as.factor(datatab$label)),
                     guide=guide_legend(override.aes = list(color='black',
                                                            alpha=1)))
show(plot_decay)
save_file='Output\\Figures\\longterm_antibody_decay.jpg'
ggsave(save_file,plot_decay,width=10,height=5,dpi=600)


# Plot percentages --------------------------------------------------------


percent_quant<-as_tibble(t(apply((1-tc)*100,2,quantile,probs=c(0.025,0.5,0.975))))
names(percent_quant)<-c('LCI','GMT','UCI')
percent_quant$Group<-levels(datatab$Group)
percent_quant<-separate_wider_delim(percent_quant,cols=Group,delim='_',names=c('Dose2','Dose3','Vaccinia Experience'))%>%
  mutate_at(c('Dose2','Dose3'),as.numeric)%>%
  mutate(`Number of MVA Doses`=factor((Dose2>0)+(Dose3>0)))%>%
  mutate(`Final Dose Day`=if_else(Dose2==0,1,(abs(Dose2-Dose3)+(`Dose3`+Dose2))/2))%>%
  mutate(Dose2=factor(Dose2))%>%
  filter(`Vaccinia Experience`==0)
plot_prop<-ggplot(data=percent_quant,aes(color=`Number of MVA Doses`,x=`Final Dose Day`)) +
  geom_point(aes(y=GMT),position=position_dodge2(width=0.1),size=3) +
  geom_errorbar(aes(ymin=LCI,ymax=UCI),position=position_dodge2(width=10,preserve='single'),size=0.8) +
  ylab('Percentage of slow decay antibodies (%)') +
  scale_y_continuous(limits = c(0,100)) +
  scale_x_continuous(trans = 'log10')+
  scale_color_manual('Number of\nMVA-BN Doses',
                     breaks=c('0','1','2'),
                     values=colD,
                     labels=c('1-Dose','2-Doses','3-Doses')) +
  theme_classic()+
  theme(axis.text = element_text(size=10),
        axis.title = element_text(size=12),
        legend.text=element_text(size=10),
        legend.title = element_text(size=12))
show(plot_prop)

ggsave('Output//Figures//Percentage_longterm.jpg',width=5,height=4,dpi=1000)

# Plot initial antibodies -------------------------------------------------

titer_quant<-data.frame(t(apply(10**mu,2,quantile,probs=c(0.025,0.5,0.975))))
names(titer_quant)<-c('LCI','GMT','UCI')
titer_quant$Group<-levels(datatab$Group)
titer_quant<-separate_wider_delim(titer_quant,cols=Group,delim='_',names=c('Dose2','Dose3','Vaccinia Experience'))%>%
  mutate_at(c('Dose2','Dose3'),as.numeric)%>%
  mutate(`Number of MVA Doses`=factor((Dose2>0)+(Dose3>0)))%>%
  mutate(`Final Dose Day`=if_else(Dose2==0,1,(abs(Dose2-Dose3)+(`Dose3`+Dose2))/2))%>%
  mutate(Dose2=factor(Dose2))%>%
  filter(`Vaccinia Experience`==0)

write_csv(titer_quant,'Output/Data/Initial_Antibody_titers.csv')

plot_GMT<-ggplot(data=titer_quant,aes(color=`Number of MVA Doses`,x=`Final Dose Day`)) +
  geom_point(aes(y=GMT),position=position_dodge2(width=0.1),size=2.5) +
  geom_errorbar(aes(ymin=LCI,ymax=UCI),position=position_dodge2(width=10,preserve='single'),size=0.8) +
  ylab('GMT') +
  #facet_grid(cols=vars(`Number of MVA Doses`),scales='free',space='free',shrink=TRUE)+
  scale_y_continuous(trans = 'log10',limits=c(3.5,4500)) +
  scale_x_continuous(trans = 'log10',
                     breaks=c(1,7,28,730),
                     labels=c('0','7','28','730'))+
  labs(x='Interval between first and last dose',y='Peak ELISA endpoint GMT')+
  scale_color_manual('Number of\nMVA-BN Doses',
                       breaks=c('0','1','2'),
                       values=colD,
                       labels=c('1-Dose','2-Doses','3-Doses')) +
  theme_classic()+
  theme(axis.text = element_text(size=10),
        axis.title = element_text(size=12),
        legend.text=element_text(size=10),
        legend.title = element_text(size=12))
show(plot_GMT)

ggsave('Output//Figures//Initial_GMT.jpg',width=5,height=4,dpi=1000)



# Plot long-term antiboy titers -------------------------------------------


long_titer_quant<-as_tibble(t(apply(Long_titer,2,quantile,probs=c(0.025,0.5,0.975))))
names(long_titer_quant)<-c('LCI','GMT','UCI')
long_titer_quant$Group<-levels(datatab$Group)
long_titer_quant<-separate_wider_delim(long_titer_quant,cols=Group,delim='_',names=c('Dose2','Dose3','Vaccinia Experience'))%>%
  mutate_at(c('Dose2','Dose3'),as.numeric)%>%
  mutate(`Number of MVA Doses`=factor((Dose2>0)+(Dose3>0)))%>%
  mutate(`Final Dose Day`=if_else(Dose2==0,1,(abs(Dose2-Dose3)+(`Dose3`+Dose2))/2))%>%
  mutate(Dose2=factor(Dose2))%>%
  filter(`Vaccinia Experience`==0)%>%
  mutate(alpha_val=c('1','0','0','1','1','1','1'))


plot_longGMT<-ggplot(data=long_titer_quant,aes(color=`Number of MVA Doses`,x=`Final Dose Day`,alpha=alpha_val)) +
  geom_point(aes(y=GMT),position=position_dodge2(width=0.1),size=2.5) +
  geom_errorbar(aes(ymin=LCI,ymax=UCI),position=position_dodge2(width=10,preserve='single'),size=0.8) +
  ylab('Predicted ELISA endpoint GMT\n at 1 year') +
  xlab('Interval between first and last dose')+
  scale_y_continuous(trans = 'log10',limits=c(3.5,4500)) +
  scale_x_continuous(trans = 'log10',
                     breaks=c(1,7,28,730),
                     labels=c('0','7','28','730'))+
  scale_color_manual('Number of\nMVA-BN Doses',
                     breaks=c('0','1','2'),
                     values=colD,
                     labels=c('1-Dose','2-Doses','3-Doses')) +
  scale_alpha_manual(breaks=c('0','1'),
                     values=c(0.5,1),
                     guide='none')+
  theme_classic()+
  theme(axis.text = element_text(size=10),
        axis.title = element_text(size=12),
        legend.text=element_text(size=10),
        legend.title = element_text(size=12))
show(plot_longGMT)

ggsave('Output//Figures//initial_longterm_titer.jpg',width=10,height=5)

plot_s_l<-ggarrange(plot_GMT,plot_longGMT,labels=c("A","B"),ncol=2,common.legend = TRUE,legend='right')
show(plot_s_l)

ggsave('Output//Figures//GMT_fit.jpg',width=8,height=4,dpi=1000)
ggsave('Output//Figures//GMT_fit.eps',width=8,height=4,dpi=1000)

# plot Study Bias ---------------------------------------------------------

mu_s=rstan::extract(fit1)$mu_j
mu_df=data.frame(t(apply(mu_s,2,quantile,probs=c(0.025,0.5,0.975))))
mu_df$Study=levels(datatab$Study)

plot_bias<-ggplot(mu_df, aes(x=Study)) +
  geom_point(aes(y=X50.)) +
  geom_errorbar(aes(ymin=X2.5.,ymax=X97.5.))
show(plot_bias)



#  Fit decay model for 1 dose data only -----------------------------------

dose1_decay=filter(datatab,Group=='0_0_0')%>%
  mutate_at(c('Group','Formulation','Study'),fct_drop)

data1 <- list(
  N = nrow(dose1_decay),
  y_n = log10(dose1_decay$GMT),
  time=dose1_decay$Week,
  s_n=dose1_decay$std,
  n_n=dose1_decay$Number,
  num_doses=1,
  dose_ind=as.integer(dose1_decay$Group),
  num_form = 2, 
  num_studies= nlevels(dose1_decay$Study),
  study_ind = as.integer(dose1_decay$Study),
  form_ind=as.integer(dose1_decay$Formulation)
)

fit1dose <- stan(
  file = "Stan/Stan_Model.stan",  # Stan program
  data = data1,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 15000,          # number of warmup iterations per chain
  iter = 20000,            # total number of iterations per chain
  seed=10,
  cores=4,
  control = list(max_treedepth=14,stepsize=1),
  algorithm = 'NUTS'
)




# Fit decay model for 2 dose only -----------------------------------------


dose2_decay=filter(datatab,Group=='28_0_0')%>%
  mutate_at(c('Group','Formulation','Study'),fct_drop)

data2 <- list(
  N = nrow(dose2_decay),
  y_n = log10(dose2_decay$GMT),
  time=dose2_decay$Week,
  s_n=dose2_decay$std,
  n_n=dose2_decay$Number,
  num_doses=1,
  dose_ind=as.integer(dose2_decay$Group),
  num_form = 2, 
  num_studies= nlevels(dose2_decay$Study),
  study_ind = as.integer(dose2_decay$Study),
  form_ind=as.integer(dose2_decay$Formulation)
)

fit2dose <- stan(
  file = "Stan/Stan_Model.stan",  # Stan program
  data = data2,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 15000,          # number of warmup iterations per chain
  iter = 20000,            # total number of iterations per chain
  seed=10,
  cores=4,
  control = list(max_treedepth=14,stepsize=1),
  algorithm = 'NUTS'
)


output<-rstan::extract(fit1dose)
D11<-output$decay1
D21<-output$decay2
tc<-output$tc[,1]
mu<-output$mu_d[,1]
sd<-output$sigma_d

Long_titer<- array(dim=20000)
Long_titer<-(1-tc)*(10**mu)*exp(-D21*52)+(tc)*(10**mu)*exp(-D11*52)
write_csv(data.frame(Long_titer),"Output\\Data\\long_titer_1dose.csv")

times=seq(0,106)
pred_decay<-data.frame()
q=array(dim=c(107,3))
for (j in 1:105){
  titer<-10**mu*(tc*exp(-D11*(times[j]))+(1-tc)*exp(-D21*(times[j])))
  q[j,]<-quantile(titer,probs=c(0.025,0.5,0.975))
}
q_df<-data.frame(q)
q_df$Time<-times

output<-rstan::extract(fit2dose)
D12<-output$decay1
D22<-output$decay2
tc<-output$tc[,1]
mu<-output$mu_d[,1]
sd<-output$sigma_d

Long_titer<- array(dim=20000)
Long_titer<-(1-tc)*(10**mu)*exp(-D22*52)+(tc)*(10**mu)*exp(-D12*52)
write_csv(data.frame(Long_titer),"Output\\Data\\long_titer_2dose.csv")


for (j in 1:105){
  titer<-10**mu*(tc*exp(-D12*(times[j]))+(1-tc)*exp(-D22*(times[j])))
  q[j,]<-quantile(titer,probs=c(0.025,0.5,0.975))
}
q2_df<-data.frame(q)
q2_df$Time<-times

plot_decay<-ggplot(data=q_df,aes(x=Time)) +
  geom_line(color='red',aes(y=X2)) +
  geom_ribbon(aes(ymin=X1,ymax=X3),color='red',fill='blue',alpha=0.2) +
  geom_point(data=dose1_decay,aes(x=Week,y=GMT)) +
  geom_line(data=dose1_decay,aes(x=Week,y=GMT),alpha=0.2) +
  geom_errorbar(data=dose1_decay,aes(x=Week,ymin=LCI,ymax=UCI)) +
  geom_line(data=q2_df,color='red',aes(y=X2)) +
  geom_ribbon(data=q2_df, aes(ymin=X1,ymax=X3),color='red',fill='blue',alpha=0.2) +
  geom_point(data=dose2_decay,aes(x=Week,y=GMT)) +
  geom_line(data=dose2_decay,aes(x=Week,y=GMT),alpha=0.2) +
  geom_errorbar(data=dose2_decay,aes(x=Week,ymin=LCI,ymax=UCI)) +
  labs(y='GMT',x='Weeks')+
  scale_y_continuous(trans='log10')+
  theme_classic()
show(plot_decay)
ggsave('Output\\Figures\\dose1_decay.png',plot_decay,width=10,height=5)


loo1<-waic(extract_log_lik(fit1))
loo2<-waic(extract_log_lik(fit2))
loo3<-waic(extract_log_lik(fit3))
loo4<-waic(extract_log_lik(fit4))
loo5<-waic(extract_log_lik(fit5))
print(loo_compare(loo1,loo2,loo3,loo4,loo5))

elppd1<-elpd(extract_log_lik(fit1))
elppd2<-elpd(extract_log_lik(fit2))
elppd3<-elpd(extract_log_lik(fit3))
elppd4<-elpd(extract_log_lik(fit4))
elppd5<-elpd(extract_log_lik(fit5))
print(loo_compare(elppd1,elppd2,elppd3,elppd4,elppd5))


# Calculate time till 1 dose ----------------------------------------------
mean_df<-read_csv('Output/Data/antibody_means.csv')
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


# Fold Differences --------------------------------------------------------

x<-samples$mu_d
f<-samples$tc
k1<-samples$decay1
k2<-samples$decay2
FD7_28<-10^(x[,5]-x[,8])
FD700_28<-10^(x[,9]-x[,5])
FD700_700<-10^(x[,7]-x[,9])
FD700_1<-10^(x[,7]-x[,1])

gmt11<-decay_titer(4,x[,1],k1,k2,f[,1])
gmt12<-decay_titer(52,x[,5],k1,k2,f[,5])
gmt13<-decay_titer(52,x[,7],k1,k2,f[,7])
gmt1730<-decay_titer(52,x[,9],k1,k2,f[,9])
FD1_28_0<-10^(gmt11-gmt12)
FD1_700_28<-10^(gmt13-gmt12)
FD1_700_700<-10^(gmt13-gmt1730)

