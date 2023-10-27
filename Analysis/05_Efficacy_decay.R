# -------------------------------------------------------------------------
#' R-script that calculates the long-term efficacy estimates 
#
# -------------------------------------------------------------------------


# Import decay data -------------------------------------------------------

#Load Decay Data
decay<-read.csv("Output\\Samples\\Decay_rates.csv",sep=',')
d2<-decay$D2*52
d1<-decay$D1*52

# Load means and standard deviations
mu<-read_csv("Output\\samples\\Initial_GMT_decay.csv")
sd<-read_csv("Output\\Samples\\Decay_sd.csv")

#load proportion data
tc<-read_csv("Output\\Samples\\Decay_proportions.csv")

#select 1 dose
f1<-1-tc$X1
m1<-mu$X1
s1<-sd$X1

#select 2 dose normal
f2<-1-tc$X5
m2<-mu$X5
s2<-sd$X5

#select 2 dose delayed
f3<-1-tc$X7
m3<-mu$X7
s3<-sd$X7

#select 3 dose
f4<-1-tc$X9
m4<-mu$X9
s4<-sd$X9


# Import logistic efficacy model data -------------------------------------

full_fit_df<-read_csv("Output//Samples//logistic_curve_parameters.csv")
k<-full_fit_df$k
A<-full_fit_df$A

# Generate miscellaneous  -------------------------------------------------

sample<-rnorm(1000)
N<-length(m1)
e<-exp(1)
times<-seq(0,10,by=0.1)
quant<-array(dim=c(3,length(times)))


# Simulate 2 dose effectiveness -------------------------------------------

#Simulate effectiveness over time
for (l in 1:length(times)){
  E<-array(dim=N)
  sss<-1/(1+exp(-k*(matrix(s2)%*%matrix(sample,nrow=1)+m2+log10(f2*exp(-d2*times[l])+(1-f2)*exp(-d1*times[l]))-2.5)-A))
  E<-rowMeans(sss)
  quant[,l]<-quantile(E,probs=c(0.025,0.5,0.975))
}

#Create dataframe with quantiles
longterm_2dose=data.frame(t(quant*100))
colnames(longterm_2dose)<-c('LCI','Efficacy','UCI')
longterm_2dose$time=times



# Simulate the effectiveness of 1 dose of vaccine -------------------------

#Simulate effectiveness
for (l in 1:length(times)){
  E<-array(dim=N)
  sss<-1/(1+exp(-k*(matrix(s1)%*%matrix(sample,nrow=1)+m1+log10(f1*exp(-d2*times[l])+(1-f1)*exp(-d1*times[l]))-2.5)-A))
  E<-rowMeans(sss)
  quant[,l]<-quantile(E,probs=c(0.025,0.5,0.975))
}

#Create dataframe with quantiles
longterm_1dose=data.frame(t(quant*100))
colnames(longterm_1dose)<-c('LCI','Efficacy','UCI')
longterm_1dose$time=times



# 3 dose long-term --------------------------------------------------------

#Simulate the effectiveness
for (l in 1:length(times)){
  E<-array(dim=N)
  sss<-1/(1+exp(-k*(matrix(s3)%*%matrix(sample,nrow=1)+m3+log10(f3*exp(-d2*times[l])+(1-f3)*exp(-d1*times[l]))-2.5)-A))
  E<-rowMeans(sss)
  quant[,l]<-quantile(E,probs=c(0.025,0.5,0.975))
}

# Amalgamate data into dataframe ------------------------------------------

longterm_3dose=data.frame(t(quant*100))
colnames(longterm_3dose)<-c('LCI','Efficacy','UCI')
longterm_3dose$time=times


# 2 dose boost long-term --------------------------------------------------------

#Simulate effectiveness
for (l in 1:length(times)){
  E<-array(dim=N)
  sss<-1/(1+exp(-k*(matrix(s4)%*%matrix(sample,nrow=1)+m4+log10(f4*exp(-d2*times[l])+(1-f4)*exp(-d1*times[l]))-2.5)-A))
  E<-rowMeans(sss)
  quant[,l]<-quantile(E,probs=c(0.025,0.5,0.975))
}

# Amalgamate data into dataframe ------------------------------------------

longterm_boost=data.frame(t(quant*100))
colnames(longterm_boost)<-c('LCI','Efficacy','UCI')
longterm_boost$time=times


#Combine Data frames
longterm_1dose$Dose='1-Dose'
longterm_2dose$Dose='2-Dose (Day-28)'
longterm_3dose$Dose='3-Dose'
longterm_boost$Dose='2-Dose (Day-730)'

combined_df<-rbind(longterm_1dose,longterm_2dose,longterm_3dose,longterm_boost)
write_csv(combined_df,'Output/Figure4/Efficacy_decay.csv')



# Calculate efficacy after 2 years ----------------------------------------
time<-2
E1<-1/(1+exp(-k*(m1+log10(f1*exp(-d2*time)+(1-f1)*exp(-d1*time))-2.5)-A))
E2<-1/(1+exp(-k*(m2+log10(f2*exp(-d2*time)+(1-f2)*exp(-d1*time))-2.5)-A))
diff<-1-(1-E2)/(1-E1)
quantile(diff,probs=c(0.025,0.5,0.975))
