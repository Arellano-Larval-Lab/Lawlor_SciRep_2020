# # Lawlor and Arellano 2020
# Scientific Reports

# Script 7. larval growth model selection (GAM)

# GAM library
#---------------
library(mgcv)
#----


# load files and make df
#---------------
load(here("data_processed","growthchem_means.Rdata"))
load(here("data_processed","regressions.Rdata"))

growthgam <- growthchem_means %>% 
  left_join(regressions,by=c("Cup")) %>%
  select(Cup,Gr,Temp,Sal,pH,pCO2,Ar,DIC)
#----


# Make list of GAM formulas to test
#-----------------------------------------
# pH, pCO2, and Ar never used together
# Ar and Sal never used together

gam_trials <- list( 
  # list all possible models here to test each. 
  # pH, pCO2, and Ar never used together
  # Ar and Sal never used together bc correlations
  
  # Null model
  Gr ~ NULL,
  # ind. main effects,
  Gr ~ s(Temp),
  Gr ~ s(Sal),
  Gr ~ s(pH),
  Gr ~ s(pCO2),
  Gr ~ s(Ar),
  Gr ~ s(DIC),
  
  # two variables
  Gr ~ s(Temp) + s(Sal),
  Gr ~ s(Temp) + s(pH),
  Gr ~ s(Temp) + s(pCO2),
  Gr ~ s(Temp) + s(Ar),
  Gr ~ s(Sal) + s(pH),
  Gr ~ s(Sal) + s(pCO2),
  
  # three variables
  Gr ~ s(Temp) + s(Sal) + s(pH),
  Gr ~ s(Temp) + s(Sal) + s(pCO2),
  #Gr ~ s(Temp) + s(Sal) + s(Ar),
  
  # add interactions
  Gr ~ te(Temp, Sal),
  Gr ~ te(Temp, Sal) + s(pH),
  Gr ~ te(Temp, Sal) + s(pCO2),
  # Gr ~ te(Temp, Sal) + s(Ar),
  Gr ~ te(Temp, pH),
  Gr ~ te(Temp, pH) + s(Sal),
  Gr ~ te(Temp, pCO2),
  Gr ~ te(Temp, pCO2) + s(Sal),                  
  Gr ~ te(Temp, Ar),
  #  Gr ~ te(Temp, Ar) + s(Sal), 
  Gr ~ te(Sal, pH),
  Gr ~ te(Sal, pH) + s(Temp),
  Gr ~ te(Sal, pCO2),
  Gr ~ te(Sal, pCO2) + s(Temp), 
  
  # interactions while one is linear
  Gr ~ te(Sal,by=Temp),
  Gr ~ te(Sal, by=Temp) + s(pH),
  Gr ~ te(Sal, by=Temp) + s(pCO2),
  #Gr ~ te(Sal, by=Temp) + s(Ar),
  
  Gr ~ te(Temp, by=Sal),
  Gr ~ te(Temp, by=Sal) + s(pH),
  Gr ~ te(Temp, by=Sal) + s(pCO2),
  #  Gr ~ te(Temp, by=Sal) + s(Ar),
  
  Gr ~ te(Temp,by=pH),
  Gr ~ te(Temp,by=pH) + s(Sal),
  Gr ~ te(Temp,by=pCO2),
  Gr ~ te(Temp,by=pCO2) + s(Sal),
  Gr ~ te(Temp,by=Ar),
  #  Gr ~ te(Temp, by=DIC),
  #  Gr ~ te(Temp, DIC),
  # Gr ~ te(Temp,by=Ar) + s(Sal),
  
  Gr ~ te(Sal,by=pH),
  Gr ~ te(Sal,by=pH) + s(Temp),
  Gr ~ te(Sal,by=pCO2),
  Gr ~ te(Sal,by=pCO2) + s(Temp),
  # Gr ~ te(Sal,by=Ar),
  # Gr ~ te(Sal,by=Ar) + s(Temp),
  
  # two interactions? 
  Gr ~ te(Temp, pH) + te(Temp, Sal),
  Gr ~ te(Temp, pCO2) + te(Temp, Sal),
  #Gr ~ te(Temp, Ar) + te(Temp, Sal),
  Gr ~ te(Sal, pH) + te(Temp, Sal),
  Gr ~ te(Sal, pCO2) + te(Temp, Sal)
)
#----


# run the whole list and rank models
#-------------------------------------
growth_gam_select <- purrr::map_df(gam_trials, ~{
  gam <- gam(.x, 
             data= growthgam,  # always use same df
             method="REML")  # always use REML smoothing method
  data.frame(cbind(formula = as.character(list(formula(.x))),
                   round(broom::glance(gam),2),
                   R2 = round(summary(gam)$r.sq,2),
                   AICc = round(MuMIn::AICc(gam),2)))
})

# add rank columns
growth_gam_select_table <- growth_gam_select %>% mutate('AICc Rank' = dense_rank(AICc),
                                                        'BIC Rank' = dense_rank(BIC)) %>%
  select(formula, df, logLik, AICc, 'AICc Rank', BIC, 'BIC Rank', R2)

# save for Supp table S2. Add annotations in inkscape
write.csv(growth_gam_select_table,file = here::here("data_processed","Growth_GAM_selection_table.csv"))
#----


# choose best model
#----------------------------------
# choose final model (AKA lowest AIC/highest R2)
growthgam_final <- gam(Gr ~ te(Temp,Sal), 
                       data= growthgam,  # always use same df
                       method="REML") 


summary(growthgam_final)$p.table
summary(growthgam_final)$s.table
summary(growthgam_final)

save(growthgam_final,file = here("data_processed","growthgam_final.Rdata"))
saveRDS(growthgam_final, file = here("data_processed","growthgam_finalRDS.rds"))


# view it
vis.gam(growthgam_final, 
        plot.type = "contour", 
        color = "terrain", 
        main = "tensor product")

#----

mod <- readRDS(here("data_processed",'growthgam_finalRDS.rds'))


predict(mod,list(Temp=27,Sal=30))


# make interpolated heatmap (fig 5a)
#----------------------------------
growthchem_means <- regressions %>%
  select(Cup,Gr) %>%
  right_join(growthchem_means,by="Cup")

interpdf <-interp2xyz(interp(x=growthchem_means$Sal, 
                             y=growthchem_means$Temp, 
                             z=growthchem_means$Gr, 
                             duplicate="mean"), data.frame=TRUE)

# plot interpolated heatmap
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
  scale_fill_gradientn(colours=pal,name="Growth Rate\n(µm/day)",
                       limits=c(0,max(regressions$Gr)),
                       guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) + 
  theme_classic()# +
#  theme(legend.position = c(.1,.9),
#        legend.justification = c(0,1),
#        legend.background = element_rect(color="black",fill=alpha("white",.7)))
# this makes a nonlinearly interpolated heatmap of temp vs sal vs growth


# add points and format
fig5a <- heatmap +
  geom_errorbar(data=growthchem_means,inherit.aes=F,aes(x=Sal,ymin = Temp-Tempsd,ymax = Temp+Tempsd),size=.5,width=.1)+
  geom_errorbarh(data=growthchem_means,inherit.aes=F,aes(y=Temp,xmin = Sal-Salsd,xmax = Sal+Salsd),size=.5,height=.1)+
  geom_point(data=growthchem_means,inherit.aes=F,aes(x=Sal,y=Temp,group=Cup,fill=Gr),shape=21,color="black",size=4)+
  
  theme(legend.background = element_rect(fill="transparent",color=NA),
        legend.key = element_rect(colour = "transparent", fill = alpha("red", 0)),
        legend.position = "right",
        legend.direction = "vertical",
        #legend.key.height  = unit(2,'cm'),
        #legend.key.width  = unit(.35,'cm'),
        panel.border = element_rect(color="black",fill="transparent"),
        legend.margin=margin(l = -.05,b=0, unit='in'),
        plot.margin = margin(0,0,0,0,unit="pt")
  )+
  guides(fill = guide_colorbar(title.position = "top",
                               ticks.colour = "black",
                               frame.colour = "black",
                               #title.hjust = .5,
                               barheight  = unit(2.5,'in'),
                               barwidth   = unit(.11,'in')
  )) 

 # guides(fill = FALSE)+

?guide_colorbar
#----

?unit
# predict growth rates with GAM 
#------------------------------------
# Simulate and predict gam
# make function to predict from GAM model
larvagrow_gam <- function(GrTemp,GrSal){
  
  LarvaGr <- predict(growthgam_final,list(Temp=GrTemp,Sal=GrSal))
  return(LarvaGr)
}

# save fuction
save(larvagrow_gam, file=here("data_processed","growfunction.Rdata"))


# simulate df to be predicted by GAM
gamsim <- data.frame(matrix(ncol=3,nrow=1120))
colnames(gamsim) <- c("temp","sal","Gr")
gamsim$sal <- rep(rep(seq(from=12, to=39.5,by=.5)),each=20)
gamsim$temp <- rep(11:30,times=56)
# populate with predicted growthrate
gamsim$Gr <- larvagrow_gam(gamsim$temp,gamsim$sal)
range(gamsim$Gr) # view range of growth rates
gamsim$Gr[gamsim$Gr<0] <- 0 # change non-zero values to 0
#----


# Plot predicted values (fig 5b)
#--------------------------------------
pal <- pnw_palette("Bay",100)
heatmap_gam <- gamsim %>%
  ggplot(aes(x = sal, y = temp, z = Gr, fill = Gr)) + 
  #  geom_tile() + 
  geom_raster() +
  ylab("Temperature (°C)")+
  xlab("Salinity (PSU)")+
  coord_cartesian(xlim = c(12, 40), ylim = c(10, 30))+
  geom_contour(color = "white", alpha =.15) + 
  scale_fill_gradientn(colours=pal,name="Model Predicted\nGrowth Rate\n(µm/day)",
                       limits=c(0,max(growthchem_means$Gr)),
                       guide = guide_colorbar(order=1,ticks.colour = "black",frame.colour = "black")) + 
  
  theme_classic() 


save(heatmap_gam,file=here("data_processed","heatmap_gam.Rdata"))

# Add real data and errors and format
# add predicted value and diff to
growthchem_means <- growthchem_means %>%
  mutate(gampredict = larvagrow_gam(growthchem_means$Temp,growthchem_means$Sal),
         gamdiff = abs(gampredict - Gr))


fig5b <- heatmap_gam +
  geom_errorbar(data=growthchem_means,inherit.aes=F,aes(x=Sal,ymin = Temp-Tempsd,ymax = Temp+Tempsd,color=gamdiff),size=.5,width=.1)+
  geom_errorbarh(data=growthchem_means,inherit.aes=F,aes(y=Temp,xmin = Sal-Salsd,xmax = Sal+Salsd,color=gamdiff),size=.5,height=.1)+
  geom_point(data=growthchem_means,inherit.aes=F,aes(x=Sal,y=Temp,group=Cup,fill=Gr,color=gamdiff),shape=21,size=4)+
  scale_color_gradient(low=c("black"),high = c("red"),name="Absolute difference\n|Predicted-Real|\n(µm/Day)",guide = guide_legend(order=0))+
  annotate("text",label="Adj. R-sq:  0.92",lineheight=.9,x=12,y=11,hjust=0,vjust=0,color="white") +
  
  # add border
  theme(panel.background = element_rect(color="black",fill="transparent",size=.5)) +
  
  theme(legend.position = "right",
        legend.background  = element_blank(),
      #  legend.box.background = element_rect(color="black"),
      legend.margin=margin(l = -.05,b=0, unit='in',
                           r = 0) ,
      plot.margin = margin(0,0,0,0,unit="pt")
        )+
  guides(fill = guide_colorbar(order=1,
                               title.position = "top",
                               ticks.colour = "black",
                               frame.colour = "black",
                               #title.hjust = .5,
                               barheight  = unit(2.5,'in'),
                               barwidth   = unit(.11,'in')),
         color = guide_legend(order = 2,
                              #title.hjust = .5,
                              title.position = "top",
                              #keywidth  = unit(1,"in")
                              ))
#----



quartz()
# merge and save
#-----------------------------------------
fig5 <- fig5a / fig5b + 
  plot_annotation(tag_levels = 'a')  & 
  theme(legend.justification = "left",
        plot.margin = margin(0,0,0,0,unit="pt"))

ggsave(fig5,
       filename = here::here("figures","Fig5.3.pdf"),
       dpi=300,
       device = cairo_pdf,
       width = 8.5,
       height = 10.5)
#----


f <- formula(survmod)
paste(deparse(f, width.cutoff = 500), collapse="")
