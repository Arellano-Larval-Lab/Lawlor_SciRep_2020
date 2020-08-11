# # Lawlor and Arellano 2020
# Scientific Reports

# Script 4. chemistry treatments over time figure



# load data
#----------------------
load(here("data_processed","growthchem.Rdata"))
load(here("data_processed","growthchem_means.Rdata"))
#----

# plot temperatures over time
#-------------------------------------
temps <- growthchem %>%
  
  # add target temps
  mutate(Target_Temp = case_when(Cup %in% c("A1","A2","A3","A4","A5",
                                            "B1","B2","B3","B4","B5") ~ "~14",
                                 Cup %in% c("C1","C2","C3","C4","C5",
                                            "D1","D2","D3","D4","D5") ~ "~18",
                                 Cup %in% c("E1","E2","E3","E4","E5",
                                            "F1","F2","F3","F4","F5") ~ "~22",
                                 Cup %in% c("G1","G2","G3","G4","G5",
                                            "H1","H2","H3","H4","H5") ~ "~26",
                                 Cup %in% c("I1","I2","I3","I4","I5",
                                            "J1","J2","J3","J4","J5") ~ "~29")) %>%
  
  # plot
  ggplot(aes(x=day,y=Temp,group=Cup,color=Target_Temp))+
  geom_line()+
  ylim(11,32)+
  
  scale_color_discrete(name="Target\nTemperature (°C)",
                       limits=c("~29","~26","~22","~18","~14"))+
  ylab("Temperature (°C)")+
  xlab("Day")+
  coord_cartesian(xlim=c(0,18),expand = F)+
  theme_classic()+
  theme(legend.position = "right")
#----



# plot salinities over time
#------------------------------------
sals <- growthchem %>% 
  
  # add targets
  mutate(Target_Sal = case_when(Cup %in% c("I4","G1","F2","C1","B2") ~ "12",
                                Cup %in% c("J5","H1","F5","C2","B4") ~ "15",
                                Cup %in% c("I1","G3","F3","D3","B1") ~ "18",
                                Cup %in% c("J2","G2","E5","C4","B5") ~ "21",
                                Cup %in% c("I5","H3","E4","D4","A3") ~ "24",
                                Cup %in% c("J1","H2","E3","C5","A1") ~ "27",
                                Cup %in% c("I2","G5","E1","D5","A5") ~ "30",
                                Cup %in% c("I3","G4","E2","D2","A4") ~ "33",
                                Cup %in% c("J3","H5","F4","C3","B3") ~ "36",
                                Cup %in% c("J4","H4","F1","D1","A2") ~ "39")) %>%
  ggplot(aes(x=day,y=Sal,group=Cup,color=Target_Sal))+
  geom_line()+
  ylim(11,41)+
  scale_color_discrete(name="Target Salinity\n(PSU)",
                       limits=c("39","36","33","30","27","24","21","18","15","12"))+
  ylab("Salinity (PSU)")+
  xlab("Day")+
  coord_cartesian(xlim=c(0,18),expand = F)+
  theme_classic()+
  theme(legend.position = "right")
#----



# plot pH over time by tube color (bubbled air treatment)
#-------------------------------------
pHs <- growthchem %>%
  mutate(Color = recode(Color,
                        "yellow" = "400",
                        "blue" = "800",
                        "green" = "1200",
                        "red" = "1600")) %>%
  ggplot(aes(x=day,y=pH,group=Cup,color=Color))+
  geom_line()+
  scale_color_discrete(name = expression(atop(pCO[2]*" of                ", "bubbled air (ppm)")),
                       limits = c("400","800","1200","1600")) +
  theme(axis.title = element_text(hjust=0))+
  #scale_color_manual(values = c("yellow" = "gold","blue" = "cornflowerblue","green" = "forestgreen","red"="tomato2"))+
  ylab("pH")+
  xlab("Day") +
  coord_cartesian(xlim=c(0,18),expand = F)+
  theme_classic()+
  theme(legend.position = "right")
#----


# merge and save treatments over time (fig S1)
#------------------------------------
FigS1 <- egg::ggarrange(temps,sals,pHs,
                        ncol = 1,
                        labels=c("a","b","c"))

ggsave(FigS1,filename = here::here("figures","FigS1.eps"),
       device = cairo_ps,
       scale=1.6,
       height = 7,
       width = 5)
#----


# plot trt averages with sterr (fig 2)
#--------------------------------------
pal <- pnw_palette("Bay",n=50,type="continuous")
Trts <- ggplot(growthchem_means,aes(x=Sal,y=Temp,group=Cup,color=pH))+
  geom_point(size=3)+
  #scale_color_viridis()+
  ylab("Temperature (°C)")  +  xlab("Salinity (PSU)")+
  scale_color_gradientn(colours=rev(pal),
                        guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))+
  geom_errorbar(aes(ymin = Temp-Tempsd,ymax = Temp+Tempsd),size=.75,width=.25)  +
  geom_errorbarh(aes(xmin = Sal-Salsd,xmax = Sal+Salsd),size=.75,height=.25)+
  geom_text_repel(aes(label=Cup),hjust=0.5, vjust=0.5,color="black",size=2.5,fontface="bold")+
  ylim(12,31)+xlim(11,41)+
  theme_classic() +
  theme(
   # legend.key.width  = unit(.4,'cm'),
  #  legend.key.height  = unit(1,'cm'),
    panel.border = element_rect(color="black",fill="transparent"),
    #legend.margin = margin(l=-.05,unit="in"),
    plot.margin = margin(0,0,0,0,"in"),
    legend.position = "right",
    legend.margin = margin(l=-.05,unit="in")
  ) +
  guides(color=guide_colorbar(barwidth = unit(.13,"in"),
                             barheight = unit(4,"in"),
                             ticks.colour = "black",
                             frame.colour = "black",
                             title.position = "top"))



# save
ggsave(filename = here::here("figures","fig2.2.eps"),
       device = cairo_ps,
       height = 5.5,
       width = 8)
#----



# plot trt averages with sterr, other orientations (fig S2)
#--------------------------------------
# plot with temp by pH, color = sal
pal <- pnw_palette("Bay",n=50,type="continuous")
Trts2 <- ggplot(growthchem_means,aes(x=Temp,y=pH,group=Cup,color=Sal))+
  geom_point(size=3)+
  #scale_color_viridis()+
  xlab("Temperature (°C)")  +  ylab("pH")+
  scale_color_gradientn(colours=rev(pal),name= "Salinity\n(PSU)",
                        guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))+
  geom_errorbarh(aes(xmin = Temp-Tempsd,xmax = Temp+Tempsd),size=.75,height=.01)  +
  geom_errorbar(aes(ymin = pH-pHsd,ymax = pH+pHsd),size=.75,width=.25)+
  geom_text_repel(aes(label=Cup),hjust=0.5, vjust=0.5,color="black",size=2.5,fontface="bold")+
  xlim(12,31)+ylim(7.55,8.15)+
  theme_classic()+
  theme(
    legend.key.width  = unit(.4,'cm'),
    legend.key.height  = unit(1,'cm'),
    panel.border = element_rect(color="black",fill="transparent")
    
  )



# plot with Sal by pH, color = temp
pal <- pnw_palette("Bay",n=50,type="continuous")
Trts3 <- ggplot(growthchem_means,aes(x=pH,y=Sal,group=Cup,color=Temp))+
  geom_point(size=3)+
  ylab("Salinity (PSU)")  +  xlab("pH")+
  scale_color_gradientn(colours=rev(pal),name = "Temperature\n(°C)",
                        guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))+
  geom_errorbar(aes(ymin = Sal-Salsd,ymax = Sal+Salsd),size=.75,width=.01)  +
  geom_errorbarh(aes(xmin = pH-pHsd,xmax = pH+pHsd),size=.75,height=.5)+
  geom_text_repel(aes(label=Cup),hjust=0.5, vjust=0.5,color="black",size=2.5,fontface="bold")+
  ylim(11,41)+xlim(7.55,8.15)+
  theme_classic() +
  theme(
    legend.key.width  = unit(.4,'cm'),
    legend.key.height  = unit(1,'cm'),
    panel.border = element_rect(color="black",fill="transparent")
    
  )


FigS2 <- egg::ggarrange(Trts2,Trts3,
                        labels=c("a","b"))


ggsave(FigS2, filename = here::here("figures","figS2.eps"),
       device = cairo_ps,
       width = 8,
       height = 11,
       scale = .8)
#----





