# This scripts imports the original data file and creates additional columns

# Import Efficacy Data ----------------------------------------------------

eff<-read_csv('Data/Efficacy_data_regress_single.csv')%>%
  mutate(Vacc_dose=if_else(Dose==0,'Unvaccinated',paste(Vaccine,Dose,sep='_')),
         group=paste(Study,Group,sep='_'))%>%
  write_csv(file='Output/Data/Eff.csv')
