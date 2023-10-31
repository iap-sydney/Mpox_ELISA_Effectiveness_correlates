# -------------------------------------------------------------------------
#'Script for generating Figure 2 from the manuscript
#
# -------------------------------------------------------------------------

#Load immunogenicity Data 
datatab <- read_csv("Output/Data/ELISA_full_set_update.csv",show_col_types = FALSE)%>%
  mutate(ELISA=replace_na(ELISA,0),
         col_group=paste(Dose,Formulation),
         Dose=if_else(Dose=='Historic Vaccination','1st Gen',Dose))%>%
  mutate_at(c('col_group','Formulation','Study','ELISA','Dose'),as.factor)

#Load Estimated Means
mean_quants<-read_csv('Output/Figure2/Indiv_eff.csv')%>%
  mutate(Dose=if_else(Dose=='Historic Vaccination First Generation',
                      '1st Gen',Dose))

# Plot Figures ------------------------------------------------------------
#Define Colours
colForm<-c('darkorange2','dodgerblue','dodgerblue4','mediumpurple','mediumpurple4','darkorange2')

#Create Plot
plot_immunogenicity<-ggplot(datatab,aes(x=label,y=GMT,shape=ELISA,color=col_group,group=Arm)) +
  geom_point(position=position_dodge(width=0.5),size=2.5) +
  facet_grid(cols=vars(Dose),scales='free',space='free') + 
  geom_errorbar(aes(ymin=LCI,ymax=UCI), position=position_dodge(width=0.5)) + 
  geom_hline(data=mean_quants,aes(yintercept=`50%`,color=col_group),alpha=1) + 
  geom_hline(data=mean_quants,aes(yintercept=`2.5%`,color=col_group),alpha=1,linetype='dashed') + 
  geom_hline(data=mean_quants,aes(yintercept=`97.5%`,color=col_group),alpha=1,linetype='dashed') +
  theme_classic() +
  labs(y='ELISA endpoint GMT',x='Immunogenicity Study')+
  scale_shape_manual(values=c(15,16,17),
                     labels = c('Unspecified', "450nm/0.3",'400nm/0.35')) +
  scale_color_manual(values=colForm) +
  scale_y_continuous(trans='log10') +
  scale_x_discrete(guide = guide_axis(angle = 65))+
  theme(strip.background = element_blank(),
        axis.text=element_text(size=10),
        axis.title=element_text(size=12),
        axis.title.x = element_blank(),
        legend.position = c(0.88,0.25),
        legend.text=element_text(size=10),
        legend.title=element_text(size=12),
        strip.text = element_text(size=12))+
  guides(color='none')

#Save Figure
ggsave('Output//Figure2//Fig2.jpg',
       plot=plot_immunogenicity,
       width=8,height=4.5,dpi=2000)
ggsave('Output//Figure2//Fig2.eps',
       plot=plot_immunogenicity,
       width=8,height=4.5,dpi=2000)
ggsave('Output//Figure2//Fig2.pdf',
       plot=plot_immunogenicity,
       width=8,height=4.5,dpi=2000)
