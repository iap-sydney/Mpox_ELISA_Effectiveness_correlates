
# -------------------------------------------------------------------------
# This script loads the immunogenicity and formats the full dataset into
# the categories to be used in each component of the analysis
#
# -------------------------------------------------------------------------


# Load Clinical Trial Data ------------------------------------------------

ELISA_study<-read_csv("Data/ELISA_studies.csv")

# Load Immunogenicity Data -------------------------------------------------

#Clean ELISA data and adjust estimates for consistent treatment of values 
# below LOD
ELISA<-read_csv("Data/ELISA_data_full_set.csv")%>%
  mutate(loggmt=log10(GMT),
         std = (log10(UCI)-log10(`LCI`))/2*sqrt(`Number`)/qt(0.975,df=`Number`))%>%
  mutate(SS=(`Number`-1)*std+`Number`*loggmt^2)%>%
  mutate(loggmt1=loggmt-Seronegative*log10(LOD)/`Number`,
         SS1=SS-Seronegative*(log10(LOD)^2))%>%
  mutate(GMT=10^loggmt1,
         std=(SS1-`Number`*(loggmt1^2))/(`Number`-1))%>%
  mutate(LCI=10^(loggmt1+qt(0.025,df=`Number`)*std/sqrt(`Number`)),
         UCI=10^(loggmt1+qt(0.975,df=`Number`)*std/sqrt(`Number`)))%>%
  filter(!is.na(std))%>%
  mutate(`POX-MVA`=Study)%>%
  left_join(ELISA_study,by=join_by(`POX-MVA`))%>%
  mutate(label=paste(Author,' (',Year,')',sep=''))%>%
  mutate(Formulation=if_else(Formulation=='Dryvax','Historic Vaccination',Formulation),
         label=if_else(label=='NA (NA)',Identifier,label))

#Filter two-dose data from full dataset
ELISA_2dose <- filter(ELISA,Day==42, Dose2==28,
                      (`Dose 3`>42 | is.na(`Dose 3`)),
                      Vaccinia_EXP==0)%>%
  mutate(Dose='MVA-BN 2-dose')

#Filter one-dose data from full dataset
ELISA_1dose <- filter(ELISA,Day==28,
                      (Dose2>=28 | Dose2==0),
                      Vaccinia_EXP==0)%>%
  mutate(Dose='MVA-BN 1-dose')

#Filter the vaccinia experienced data from full dataset
ELISA_exp <- filter(ELISA,Day==0,
                    Vaccinia_EXP==1)%>%
  mutate(Dose='Historic Vaccination')%>%
  mutate(Formulation='First Generation')

#Combine immunogenicity data into dataset for analysis
ELISA_immunogenecity <- rbind(ELISA_1dose,ELISA_2dose,ELISA_exp)%>% 
  write_csv('Output/Data/ELISA_full_set_update.csv')


# Process data for decay curves -------------------------------------------

Decay<-filter(ELISA,!is.na(std))%>%
  mutate(`Dose 3`=replace_na(`Dose 3`,0))%>%
  mutate(Final_day=(abs(Dose2-`Dose 3`)+(`Dose 3`+Dose2))/2)%>%
  mutate(Final_day=ifelse(Day==Final_day,0,Final_day))%>%
  mutate(Dose2=Final_day)%>%
  mutate(Dose2=if_else(`Dose 3`>0,28,Dose2))%>%
  mutate(peak=(abs(Final_day-14)+(Final_day+14))/2+14)%>%
  mutate(Week=(Day-peak)/7)%>%
  filter(Week>=0)%>%
  write_csv('Output/Data/ELISA_decay.csv')

