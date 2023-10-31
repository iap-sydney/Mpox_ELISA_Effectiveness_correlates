# -------------------------------------------------------------------------
#'Script contains the information needed to generate Figure 3.
#
# -------------------------------------------------------------------------


# Load Data ---------------------------------------------------------------

#import sampled means
mean_df<-read_csv('Output/Samples/antibody_means.csv')

#import sampled efficacies
Efficacies_df<-read_csv('Output/Samples/Combined-Effectiveness-sample.csv')

#Import Effectiveness quantiles with GMT
eff_df<-read_csv('Output/Figure3/Efficacy_quantiles.csv')
eff_df$LCI<-eff_df$LCI*100
eff_df$median<-eff_df$median*100
eff_df$UCI<-eff_df$UCI*100


# pair samples from effectiveness and immunogencity -----------------------

#MVA-BN 1 dose set
dose1df<-data.frame(mean_df$`MVA-BN 1-dose`,Efficacies_df$`MVA-BN 1-dose`)
colnames(dose1df)<-c('GMT','Efficacy')
dose1df$Vaccine='MVA-BN 1-dose'

#MVA-BN 2 dose set
dose2df<-data.frame(mean_df$`MVA-BN 2-dose`,Efficacies_df$`MVA-BN 2-dose`)
colnames(dose2df)<-c('GMT','Efficacy')
dose2df$Vaccine='MVA-BN 2-dose'

#Historic Vaccination set
doseddf<-data.frame(mean_df$`Historic Vaccination`,Efficacies_df$`First Generation`)
colnames(doseddf)<-c('GMT','Efficacy')
doseddf$Vaccine='First Generation'

#combine sampled data
big_df=rbind(dose1df,dose2df,doseddf)


# Figure 3 ----------------------------------------------------------------

#generate figure
logistic_plot<-ggplot(eff_df,aes(x=GMT,y=median,color=Vaccine)) +
  geom_line(color='black')+
  geom_ribbon(aes(ymin=LCI,ymax=UCI),color='black',fill='black',alpha=0.3,linetype=2)+
  geom_density_2d(data=big_df,aes(y=Efficacy,color=Vaccine),
                  contour_var = 'ndensity',
                  bins=5,size=0.5)+
  scale_color_manual(breaks = levels(as.factor(big_df$Vaccine)),
                     values= colV,
                     label=levels(as.factor(big_df$Vaccine)))+
  theme_classic()+
  scale_y_continuous('Effectiveness (%)',limits=c(0,100)) +
  scale_x_continuous(breaks=1:4,labels=10^(1:4))+
  theme(legend.text = element_text(size=10),
        legend.title = element_text(size=12),
        axis.title = element_text(size=12),
        axis.text = element_text(size=10)
        )

#Save Figure
ggsave("Output//Figure3//Logistic_curve_fit.jpg",
       plot=logistic_plot,
       width=5,height=4,dpi=1000)
ggsave("Output//Figure3//Logistic_curve_fit.pdf",
       plot=logistic_plot,
       width=5,height=4,dpi=1000)
