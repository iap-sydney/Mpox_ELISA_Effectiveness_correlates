##R-script to determine the long-term efficacy estimates.

col<-c('red2','dodgerblue2','springgreen4')
# Import decay data -------------------------------------------------------


decay<-read.csv("Output\\Data\\decay.csv",sep=',')
d2<-decay$decay2*52
d1<-decay$decay1*52

mu<-read_csv("Output\\Data\\mu_d.csv")
sd<-read_csv("Output\\Data\\sigma_d.csv")

f1<-1-decay$tc.1.
m1<-mu$X1
s1<-sd$X1

f2<-1-decay$tc.5.
m2<-mu$X5
s2<-sd$X5

f3<-1-decay$tc.7.
m3<-mu$X7
s3<-sd$X7

f4<-1-decay$tc.9.
m4<-mu$X9
s4<-sd$X9


# Import logistic efficacy model data -------------------------------------

full_fit_df<-read_csv("Output//Data//full_fit_parameters.csv")
k<-full_fit_df$k
A<-full_fit_df$A

sat_fit_df<-read_csv("Output//Data//sat_fit.csv")
k_sat<-sat_fit_df$k
A_sat<-sat_fit_df$A
M_sat<-sat_fit_df$M

# Generate miscellaneous  -------------------------------------------------

sample<-rnorm(1000)
N<-length(m1)
e<-exp(1)
times<-seq(0,10,by=0.1)
quant<-array(dim=c(3,length(times)))
quant_sat<-array(dim=c(3,length(times)))

# Simulate effectiveness of the vaccine -----------------------------------

for (l in 1:length(times)){
  E<-array(dim=N)
  E_sat<-array(dim=N)
  sss<-1/(1+exp(-k*(matrix(s2)%*%matrix(sample,nrow=1)+m2+log10(f2*exp(-d2*times[l])+(1-f2)*exp(-d1*times[l]))-2.5)-A))
  E<-rowMeans(sss)
  #for (j in 1:N){
  #  E[j]<-mean(1/(1+exp(-k[j]*(sample*s2[j]+m2[j]+log10(f2*exp(-d2*times[l])+(1-f2)*exp(-d1*times[l]))-2.5)-A[j])))
    #E_sat[j]<-mean(1/(M_sat[j]+exp(-k_sat[j]*(sample*s2[j]+m2[j]+log10(f2*exp(-d2*times[l])+(1-f2)*exp(-d1*times[l]))-2.5)-A_sat[j])))
  #}
  #E<-1/(1+exp(-k*(m2+log10(f2*exp(-d2*times[l])+(1-f2)*exp(-d1*times[l]))-2.5)-A))
  E_sat<-1/(M_sat+exp(-k_sat*(m2+log10(f2*exp(-d2*times[l])+(1-f2)*exp(-d1*times[l]))-2.5)-A_sat))
  quant[,l]<-quantile(E,probs=c(0.025,0.5,0.975))
  quant_sat[,l]<-quantile(E_sat,probs=c(0.025,0.5,0.975))
}

# Amalgamate data into dataframe ------------------------------------------

longterm_df=data.frame(t(quant*100))
colnames(longterm_df)<-c('LCI','Efficacy','UCI')
longterm_df$time=times

longterm_sat_df=data.frame(t(quant_sat*100))
colnames(longterm_sat_df)<-c('LCI','Efficacy','UCI')
longterm_sat_df$time=times

# Simulate effectiveness of the vaccine -----------------------------------

for (l in 1:length(times)){
  E<-array(dim=N)
  E_sat<-array(dim=N)
  sss<-1/(1+exp(-k*(matrix(s1)%*%matrix(sample,nrow=1)+m1+log10(f1*exp(-d2*times[l])+(1-f1)*exp(-d1*times[l]))-2.5)-A))
  E<-rowMeans(sss)
  #E<-1/(1+exp(-k*(m1+log10(f1*exp(-d2*times[l])+(1-f1)*exp(-d1*times[l]))-2.5)-A))
  E_sat<-1/(M_sat+exp(-k_sat*(m1+log10(f1*exp(-d2*times[l])+(1-f1)*exp(-d1*times[l]))-2.5)-A_sat))
  quant[,l]<-quantile(E,probs=c(0.025,0.5,0.975))
  quant_sat[,l] <- quantile(E_sat,probs=c(0.025,0.5,0.975))
}

# Amalgamate data into dataframe ------------------------------------------

longterm1_df=data.frame(t(quant*100))
colnames(longterm1_df)<-c('LCI','Efficacy','UCI')
longterm1_df$time=times

longterm1_sat_df=data.frame(t(quant_sat*100))
colnames(longterm1_sat_df)<-c('LCI','Efficacy','UCI')
longterm1_sat_df$time=times

# Plot long-term Efficacy 2 dose -------------------------------------------------

plot_longterm2<-ggplot(longterm_df,aes(x=time,y=Efficacy,ymin=LCI,ymax=UCI)) +
  geom_line()+
  geom_ribbon(color='black',linetype=2,alpha=0.2)+
  xlab('Time (years)') +
  scale_y_continuous(limits=c(0,100),expand = c(0,0))+
  theme_classic()+
  scale_x_continuous(limits=c(0,10),expand = c(0,0))
show(plot_longterm2)
ggsave("Output//Figures//longterm_efficacy2.png",width=10,height=5)

plot_longterm<-ggplot(longterm_sat_df,aes(x=time,y=Efficacy,ymin=LCI,ymax=UCI)) +
  geom_line()+
  geom_ribbon(color='black',linetype=2,alpha=0)+
  xlab('Time (years)') +
  scale_y_continuous(limits=c(0,100))
show(plot_longterm)
ggsave("Output//Figures//longterm_sat_efficacy2.png",width=10,height=5)

# Plot long-term Efficacy -------------------------------------------------

plot_longterm1<-ggplot(longterm1_df,aes(x=time,y=Efficacy,ymin=LCI,ymax=UCI)) +
  geom_line(data=longterm1_df)+
  geom_ribbon(data=longterm1_df,color='black',linetype=2,alpha=0.2)+
  xlab('Time (years)') +
  scale_y_continuous(limits=c(0,100),expand = c(0,0))+
  theme_classic()+
  scale_x_continuous(limits=c(0,10),expand = c(0,0))
show(plot_longterm1)
ggsave("Output//Figures//longterm_efficacy1.png",width=10,height=5)


plot_longterm_sat<-ggplot(longterm1_sat_df,aes(x=time,y=Efficacy,ymin=LCI,ymax=UCI)) +
  geom_line()+
  geom_ribbon(color='red',fill='blue',alpha=0.3)+
  xlab('Time (years)') +
  scale_y_continuous(limits=c(0,100))
show(plot_longterm_sat)
ggsave("Output//Figures//longterm_sat_efficacy1.png",width=10,height=5)


# 3 dose long-term --------------------------------------------------------


for (l in 1:length(times)){
  E<-array(dim=N)
  E_sat<-array(dim=N)
  sss<-1/(1+exp(-k*(matrix(s3)%*%matrix(sample,nrow=1)+m3+log10(f3*exp(-d2*times[l])+(1-f3)*exp(-d1*times[l]))-2.5)-A))
  E<-rowMeans(sss)
    #E<-1/(1+exp(-k*(m3+log10(f3*exp(-d2*times[l])+(1-f3)*exp(-d1*times[l]))-2.5)-A))
  E_sat<-1/(M_sat+exp(-k_sat*(m3+log10(f3*exp(-d2*times[l])+(1-f3)*exp(-d1*times[l]))-2.5)-A_sat))
  quant[,l]<-quantile(E,probs=c(0.025,0.5,0.975))
  quant_sat[,l]<-quantile(E_sat,probs=c(0.025,0.5,0.975))
}

# Amalgamate data into dataframe ------------------------------------------

longterm_df3=data.frame(t(quant*100))
colnames(longterm_df3)<-c('LCI','Efficacy','UCI')
longterm_df3$time=times

longterm_sat_df3=data.frame(t(quant_sat*100))
colnames(longterm_sat_df3)<-c('LCI','Efficacy','UCI')
longterm_sat_df3$time=times

# 2 dose boost long-term --------------------------------------------------------


for (l in 1:length(times)){
  E<-array(dim=N)
  E_sat<-array(dim=N)
  sss<-1/(1+exp(-k*(matrix(s4)%*%matrix(sample,nrow=1)+m4+log10(f4*exp(-d2*times[l])+(1-f4)*exp(-d1*times[l]))-2.5)-A))
  E<-rowMeans(sss)
  #E<-1/(1+exp(-k*(m3+log10(f3*exp(-d2*times[l])+(1-f3)*exp(-d1*times[l]))-2.5)-A))
  E_sat<-1/(M_sat+exp(-k_sat*(m4+log10(f4*exp(-d2*times[l])+(1-f4)*exp(-d1*times[l]))-2.5)-A_sat))
  quant[,l]<-quantile(E,probs=c(0.025,0.5,0.975))
  quant_sat[,l]<-quantile(E_sat,probs=c(0.025,0.5,0.975))
}

# Amalgamate data into dataframe ------------------------------------------

longterm_df4=data.frame(t(quant*100))
colnames(longterm_df4)<-c('LCI','Efficacy','UCI')
longterm_df4$time=times

longterm_sat_df4=data.frame(t(quant_sat*100))
colnames(longterm_sat_df4)<-c('LCI','Efficacy','UCI')
longterm_sat_df4$time=times

# Plot long-term Efficacy -------------------------------------------------

plot_longterm3<-ggplot(longterm_df3,aes(x=time,y=Efficacy,ymin=LCI,ymax=UCI)) +
  geom_line()+
  geom_ribbon(color='black',linetype=2,alpha=0.3)+
  xlab('Time (years)') +
  scale_y_continuous(limits=c(0,100),expand = c(0,0))+
  theme_classic()+
  scale_x_continuous(limits=c(0,10),expand = c(0,0))
show(plot_longterm3)
ggsave("Output//Figures//longterm_efficacy3.png",width=10,height=5)

plot_longterm<-ggplot(longterm_sat_df3,aes(x=time,y=Efficacy,ymin=LCI,ymax=UCI)) +
  geom_line()+
  geom_ribbon(color='red',fill='blue',alpha=0.3)+
  xlab('Time (years)') +
  scale_y_continuous(limits=c(0,100))+
  scale_x_continuous(limits=c(0,10))
show(plot_longterm)
ggsave("Output//Figures//longterm_sat_efficacy3.png",width=10,height=5)


# Aggregate longterm efficacy plots ---------------------------------------

plot_long_eff<-ggarrange(plot_longterm1,plot_longterm2,plot_longterm3,
                        labels=c('A','B','C'),
                        ncol=1, nrow=3)
show(plot_long_eff)

longterm1_df$Dose='1-Dose'
longterm_df$Dose='2-Dose (Day-28)'
longterm_df3$Dose='3-Dose'
longterm_df4$Dose='2-Dose (Day-730)'

combined_df<-rbind(longterm1_df,longterm_df,longterm_df3,longterm_df4)
write_csv(combined_df,'Output/data/Longterm_Efficacy.csv')

plot_longterm<-ggplot(data=combined_df,
                      aes(x=time,y=Efficacy,ymin=LCI,ymax=UCI,color=Dose,fill=Dose))+
  geom_line(size=1)+
  geom_ribbon(size=1,linetype=0,alpha=0.3)+
  xlab('Time (years)') +
  scale_y_continuous(limits=c(0,100),expand = c(0,0))+
  theme_classic()+
  scale_x_continuous(n.breaks=6,limits=c(0,10),expand = c(0,0.2))+
  annotate('rect', xmin=2, xmax=10, ymin=0, ymax=100, alpha=.1, fill='black')+
  facet_grid(cols=vars(Dose))+
  scale_color_manual(name='Number of Doses',
                     breaks=levels(as.factor(combined_df$Dose)),
                     values=c(colD2),
                     labels=c('1-Dose','2-Dose (Day-28)','2-Dose(Day-730)','3-Dose'))+
  scale_fill_manual(name='Number of Doses',
                     breaks=levels(as.factor(combined_df$Dose)),
                     values=c(colD2),
                     labels=c('1-Dose','2-Dose (Day-28)','2-Dose(Day-730)','3-Dose'))+
  theme(strip.background = element_blank(),
        strip.text=element_text(size=12),
        axis.text=element_text(size=11),
        axis.title = element_text(size=12),
        legend.text = element_text(size=11),
        legend.title = element_text(size=12))
show(plot_longterm)

ggsave('Output/Figures/Manuscript_efficacy_decay.jpg',width=10,height=5,dpi=600)
ggsave('Output/Figures/Manuscript_efficacy_decay.eps',width=10,height=5,dpi=600)


# Plot Antibody Decay -----------------------------------------------------

q_df<-read_csv('Output/Data/GMT_decay.csv')%>%
  mutate(Time=Time/52)
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
  mutate(Group=factor(Group))


datatab <- read_csv("Output/Data/ELISA_decay.csv") %>%
  mutate(Study=ifelse(Study==23,5,Study),Time=Week/52) %>%
  unite(Group,Dose2,`Dose 3`,Vaccinia_EXP,sep='_',remove=FALSE)%>%
  mutate_at(c('Group','Formulation','Study'),as.factor)%>%
  mutate(Group=as.integer(Group))%>%
  filter(Group==1|Group==5|Group==7|Group==9)%>%
  mutate(Group=factor(Group))

plot_decay<-ggplot(data=datatab,
                   aes(x=Time,y=GMT,ymin=LCI,ymax=UCI,color=Group,fill=Group))+
  geom_point(aes(shape=label),size=2,alpha=0.8)+
  geom_line(aes(shape=label),alpha=0.3)+
  geom_errorbar(aes(shape=label),alpha=0.8)+
  labs(y='GMT',x='Time (years)')+
  scale_y_continuous(trans='log10')+
  theme_classic()+
  geom_line(data=q_df,size=1,linetype=1) +
  geom_ribbon(data=q_df,alpha=0.1,linetype=0) +
  scale_shape_manual(name='Study',
                     breaks=levels(as.factor(datatab$label)),
                     values=1:nlevels(as.factor(datatab$label)),
                     labels=levels(as.factor(datatab$label)),
                     guide=guide_legend(override.aes = list(color='black',alpha=1)))+
  theme(axis.title = element_text(size=12),
        axis.text=element_text(size=10),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12))+
  scale_color_manual('Number of \nMVA-BN Doses',
                     breaks=c(1,5,9,7),
                     values=c(colD2),
                     labels=c('1-Dose','2-Dose (Day-28)','2-Dose (Day-730)','3-Dose'))+
  scale_fill_manual('Number of \nMVA-BN Doses',
                    breaks=c(1,5,9,7),
                    values=c(colD2),
                    labels=c('1-Dose','2-Dose (Day-28)','2-Dose (Day-730)','3-Dose'))+
  guides(fill='none',color='none')
show(plot_decay)
  

save_file='Output\\Figures\\longterm_antibody_decay_years.jpg'
ggsave(save_file,plot_decay,width=10,height=5,dpi=600)
save_file='Output\\Figures\\longterm_antibody_decay_years.eps'
ggsave(save_file,plot_decay,width=10,height=5,dpi=600)

Decay_comb<-ggarrange(plot_decay,plot_longterm,
                      labels=c('A','B'),
                      nrow=2)
show(Decay_comb)
ggsave('Output/Figures/Manuscript_decay_comb.jpg',plot=Decay_comb,height=7,width=9,dpi=600)
ggsave('Output/Figures/Manuscript_decay_comb.eps',plot=Decay_comb,height=7,width=9,dpi=600)


plot_decay_base<-ggarrange(plot_GMT,plot_longGMT,labels=c("B","C"),ncol=2,common.legend = TRUE,legend='right')
decay_model_full<-ggarrange(plot_decay,plot_decay_base,plot_longterm,
                            labels=c("A","","D"),
                            nrow=3)
show(decay_model_full)
ggsave('Output/Figures/Manuscript_decay_full.jpg',plot=decay_model_full,height=11,width=9,dpi=600)
ggsave('Output/Figures/Manuscript_decay_full.eps',plot=decay_model_full,height=11,width=9,dpi=600)



# 2 year efficacy ---------------------------------------------------------

time<-2
E1<-1/(1+exp(-k*(m1+log10(f1*exp(-d2*time)+(1-f1)*exp(-d1*time))-2.5)-A))
E2<-1/(1+exp(-k*(m2+log10(f2*exp(-d2*time)+(1-f2)*exp(-d1*time))-2.5)-A))
diff<-1-(1-E2)/(1-E1)
quantile(diff,probs=c(0.025,0.5,0.975))
