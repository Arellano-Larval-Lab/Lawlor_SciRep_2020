# # Lawlor and Arellano 2020
# Scientific Reports

# Script 5. heatmap figure growth

# libraries 
#----------------------------------------
library(here)
library(tidyverse)
library(akima)
library(PNWColors)
#----



# load data
#---------------------------------------
load(here("data_processed","regressions.Rdata"))
load(here("data_processed","growthchem_means.Rdata"))
#----


# merge data objects
#-------------------------------------------
growthchem_means <- regressions %>%
  select(Cup,Gr) %>%
  right_join(growthchem_means,by="Cup")
#----



# make interpolation df (Akima package)
#-----------------------------------------
interpdf <-interp2xyz(interp(x=growthchem_means$Sal, 
                             y=growthchem_means$Temp, 
                             z=growthchem_means$Gr, 
                             duplicate="mean"), data.frame=TRUE)
#----


# make base interp heatmap
#-----------------------------------------
pal <- pnw_palette("Bay",100,type="continuous")

heatmap <- interpdf %>%
  filter(!is.na(z)) %>%
  as_tibble() %>%
  ggplot(aes(x = x, y = y, z = z, fill = z)) + 
  #geom_tile() + 
  geom_raster()+
  ylab("Temperature (°C)")+
  xlab("Salinity (PSU)")+
  coord_cartesian(xlim = c(12, 40), ylim = c(10, 30))+
  geom_contour(color = "white", alpha = 0.05) + 
  scale_fill_gradientn(colours=pal,name="Growth Rate (µm/day)",
                       limits=c(0,max(regressions$Gr)),
                       guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) + 
  theme_classic()# +
#  theme(legend.position = c(.1,.9),
#        legend.justification = c(0,1),
#        legend.background = element_rect(color="black",fill=alpha("white",.7)))
# this makes a great heatmap of temp vs sal vs growth
#----


# add treatment points and make pretty
#---------------------------------------
Fig5a <- heatmap +
  geom_errorbar(data=growthchem_means,inherit.aes=F,aes(x=Sal,ymin = Temp-Tempsd,ymax = Temp+Tempsd),size=.5,width=.1)+
  geom_errorbarh(data=growthchem_means,inherit.aes=F,aes(y=Temp,xmin = Sal-Salsd,xmax = Sal+Salsd),size=.5,height=.1)+
  geom_point(data=growthchem_means,inherit.aes=F,aes(x=Sal,y=Temp,group=Cup,fill=Gr),shape=21,color="black",size=4)+

  theme(legend.background = element_rect(fill="transparent",color=NA),
        legend.key = element_rect(colour = "transparent", fill = alpha("red", 0)),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.key.width  = unit(2,'cm'),
        legend.key.height  = unit(.5,'cm'),
        panel.border = element_rect(color="black",fill="transparent")
        
  )+
  guides(fill = guide_colorbar(title.position = "top")) +
  labs(title="Interpolated growth rates")


