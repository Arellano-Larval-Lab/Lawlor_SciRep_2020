# # Lawlor and Arellano 2020
# Scientific Reports

# Script 5. individual effects figures


# individual effects - growth

# load data
#--------------------------------
load(here("data_processed","regressions.Rdata"))
load(here("data_processed","growthchem_means.Rdata"))

# merge together
growthchem_means <- regressions %>%
  select(Cup,Gr) %>%
  right_join(growthchem_means,by="Cup")
#----

# Spline of Salinity
#-----------------------------------------
p1 <- growthchem_means %>%
  ggplot( aes(x=Sal,y=Gr))+
  geom_point()+
  geom_smooth(method="loess",span=1,color=c("cornflowerblue"))+
  xlab("Salinity (PSU)") + 
  #ylab("Growth Rate (μm/day)")+
  ylab(element_blank())+
  theme_classic()+
  theme(panel.background = element_rect(fill="transparent",color=NA),
        panel.border = element_rect(color="black",fill="transparent")
  )
#----


# Spline of Temp
#-----------------------------------------
p2 <- growthchem_means %>%
  ggplot(aes(x=Temp,y=Gr))+
  geom_point()+
  geom_smooth(method="loess", span=1,color=c("cornflowerblue"))+
  xlab("Temperature (°C)") + 
  #ylab("Growth Rate (μm/day)")+
  ylab(element_blank())+
  theme_classic()+
  theme(panel.background = element_rect(fill="transparent",color=NA),
        panel.border = element_rect(color="black",fill="transparent")
  )
#----

#Spline of pH
#-----------------------------------------
p3 <- growthchem_means %>%
  ggplot( aes(x=pH,y=Gr))+
  geom_point()+
  geom_smooth(method="loess",span=1,color=c("cornflowerblue"))+
  xlab("pH") + 
  #ylab("Growth Rate (μm/day)")+
  ylab(element_blank())+
  theme_classic()+
  theme(panel.background = element_rect(fill="transparent",color=NA),
        panel.border = element_rect(color="black",fill="transparent")
  )
#----

# Merge and Annotate axis
#-----------------------------------------
FigS3 <- egg::ggarrange(p1,p2,p3,ncol=1) %>%
  ggpubr::annotate_figure(left = ggpubr::text_grob( expression(paste("Growth Rate (",mu,"m/day)")),rot=90))
#----

### Fig S3 save ####
#-----------------------------------------
ggsave(FigS3, file= here::here("figures","figS3.eps"),
       device=cairo_ps,
       height=8,
       width = 4.25,
       scale=1.2)
