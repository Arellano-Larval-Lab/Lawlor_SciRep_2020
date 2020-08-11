# # Lawlor and Arellano 2020
# Scientific Reports

# Script 8. present and future PLDs plot

# Load data
#-----------------------------
load(here("data_processed","fieldsitemeans.Rdata")) # measured field valies
load(here("data_processed","growfunction.Rdata"))
load(here("data_processed","site_order.Rdata")) #north-south site order
load(here("data_processed","heatmap_gam.Rdata"))
load(here("data_processed","survival_sim.Rdata"))
load(here("data_processed","growthgam_final.Rdata"))
load(here("data_processed","sitevals.Rdata"))


#----


SSM_hotvals <- sitevals %>%
  filter(temp14 > 13)

# Make a df of PLD and GR in 2014 and 295
#-----------------------------------------------
startsize <- 156
endsize <- 260


# calculate PLD and Gr for SSM 2014 and 2095 points
SSM_points <- SSM_hotvals %>%
  
  # calculate growth rates
  mutate(Gr14 = as.numeric(larvagrow_gam(temp14,sal14)),
         Gr95 = as.numeric(larvagrow_gam(temp95,sal95))) %>%
  
  # calculate PLDs
  mutate(PLD14 = as.numeric((endsize - startsize) / Gr14),
         PLD95 = as.numeric((endsize - startsize) / Gr95)) %>%
  
  ungroup() %>%
  arrange(site) %>%
  mutate(NSorder  = c(1:9))





field_points <- fieldsitemeans %>% 
  mutate(fieldGr = as.numeric(larvagrow_gam(meantemp,meansal)),
         fieldPLD = as.numeric((endsize - startsize) / fieldGr))
#----



# make change in temperature plot (fig 7a)
#-------------------------------------------
fig7a <-  heatmap_gam +
  
  # add contours for survival
  metR::geom_text_contour(inherit.aes=F, data=survival_sim,aes(x=sal,y=temp,z=surv),color="grey30")+
  stat_contour(inherit.aes=F,data=survival_sim,aes(x=sal,y=temp,z=surv),color="grey30",alpha=.5,size=.25)+
  
  # add points for 2014 SSM
  geom_point(data=SSM_points,inherit.aes=F,aes(x=sal14,y=temp14),shape=17,color="darkblue",fill="cornflowerblue")+
  
  # add points for 2095 SSM
  geom_point(data=SSM_points,inherit.aes=F,aes(x=sal95,y=temp95),shape=19,color="tomato2",fill="tomato2") +
  
  # connect
  geom_segment(data=SSM_points,inherit.aes=F,aes(x=sal14,y=temp14,xend=sal95,yend=temp95,group=site),
               color="black",arrow = arrow(length=unit(0.01,"npc")),size=.25)+
  
  # label SSM sites
  geom_text_repel(data=SSM_points,inherit.aes = F,aes(group=site,sal14,temp14,label=NSorder)) + 
  
  # add field points
  geom_point(inherit.aes=F,data=field_points,aes(x=meansal,y=meantemp),shape=5,size=2,stroke=1.5)+
  
  # label
  geom_text_repel(inherit.aes=F,data=field_points,aes(x=meansal,y=meantemp,label=paste(site,year)),
                  xlim=c(32,NA),ylim = c(NA,19),force=5,vjust=0,hjust=0)  +
  
  # theme
  theme(panel.background = element_rect(color="black",
                                       fill="transparent",
                                       size=.5),
        legend.position = "bottom",
        legend.direction = "horizontal")+
   
   # add fake points to make a color scale
   geom_point(inherit.aes=F,
              data=data.frame(x=c(23,23,21),y=c(24,24,21),color=c("Modeled\n2014 Site","Modeled\n2095 Site","Measured\nField Site")),
              aes(x=x,y=y,color=color),
              alpha=0) + 
  scale_color_discrete(limits =c("Modeled\n2014 Site","Modeled\n2095 Site","Measured\nField Site"))+
  # guides
  guides(fill = guide_colorbar(order=1,
                               title = "Model Predicted Growth Rate (µm/day)",
                               title.position = "top",
                               barwidth = unit(3.5,"in"),
                               barheight = unit(.12,"in"),
                               title.hjust = 0,
                               frame.colour = "black",
                               ticks.colour = "black"),
         color = guide_legend(override.aes = 
                                list(pch = c(17,19,5), 
                                     color = c("darkblue", "tomato2","black"),
                                     alpha=c(1,1,1),
                                     size=c(1.5,1.5,2),
                                     stroke = c(1,1,2)),
                              title = NULL,
                              title.position = "top",
                              keyheight = unit(.12,"in"))) +
  
  # change margins
  theme(plot.margin = margin(0,0,0,0,"in"),
        legend.margin = margin(t=-0,0,0,0,"in"),
        legend.background =   element_blank(),
        #legend.box.background = element_rect(color="black")
       )
 #----
 fig7a


# prepare df for fig7b
#----------------------------------------------
fieldsite_order <- c("Fidalgo Bay 2014", "Fidalgo Bay 2017", "Liberty Bay 2017","Liberty Bay 2018","Liberty Bay 2019")

fieldsitemeans_reformat <- field_points %>%
  select(site,year,meantemp,meansal,fieldGr,fieldPLD) %>%
  rename(
    sal14 = meansal,
    temp14 = meantemp,
    Gr14 = fieldGr,
    PLD14 = fieldPLD) %>%
  mutate(site = paste(site,year)) %>%
  select(-year) 


full_sites <- plyr::rbind.fill(SSM_points,fieldsitemeans_reformat) %>% 
  mutate(NSorder = c(1:nrow(.))) %>%
  mutate(site = forcats::fct_inorder(site))

full_sites %>% glimpse()
#----



# make figure7b
#---------------------------------------------
fig7b <- full_sites %>%
  ggplot(aes(x=site,y=PLD14)) +
  # add 2014 sites and field sites
  geom_point(color=c(rep("darkblue",time=9),rep("black",time=5)),
             shape=c(rep(17,time=9),rep(5,time=5)),
             size = c(rep(1.5,time=9),rep(2,time=5)),
             stroke = c(rep(1,time=9),rep(1.5,time=5)))+
  
  # add 2014 mean
  stat_summary(inherit.aes=F,data = (full_sites %>% filter(site %in% site_order) ),geom="segment",
               aes(y=mean(PLD14),yend=mean(PLD14),x=0,xend=9.5),
               fun=mean,color="darkblue",size=1)+
  # 2014 mean label
 # stat_summary(data = (full_sites %>% filter(site %in% site_order)),
 #              geom="text", color="darkblue",
 #              aes(x=.1,y=mean(PLD14)+.25,label=round(mean(PLD14),2)),
 #              fun=mean,vjust=0,hjust=0)+
  

  # add 2095 points
  geom_point(aes(x=site,y=PLD95),
             color="tomato2",
             shape=19,
             size = 1.5,
             stroke =1) + 
  
  # add 2095 mean
  stat_summary(inherit.aes=F,data = (full_sites %>% filter(site %in% site_order) ),geom="segment",
               aes(y=mean(PLD95),yend=mean(PLD95),x=0,xend=9.5),
               fun=mean,color="tomato2",size=1) + 
  # 2095 mean label
#  stat_summary(data = (full_sites %>% filter(site %in% site_order)),
#               geom="text", color="tomato2",
#               aes(x=.1,y=mean(PLD95)+.25,label=round(mean(PLD95),2)),
#               fun=mean,vjust=0,hjust=0)+
  
  # add segments 
  geom_segment(inherit.aes=F,aes(x=site,y=(PLD14-1),xend=site,yend=(PLD95+1),group=site),
               color="black",arrow = arrow(length=unit(0.01,"npc")),size=.15)+
  
  # theme stuff
  theme_classic()+
  theme(axis.text.x = element_text(angle=45,hjust=1))+
  geom_vline(xintercept = 9.5,linetype="dotted")+
  
  # axis labels
  scale_y_continuous(name="Predicted PLD (days)",
                     limits =c(0,60))+
  
  # add labels
  annotate(geom="text", hjust=.5,
           x=5,y=60,label="Salish Sea Model Values")+
  annotate(geom="text", hjust=.5,
           x=12,y=60,label="Field Values") +
  
  # add panel box
  theme(panel.background = element_rect(color="black",
                                        fill="transparent",
                                        size=.5),
        axis.title.x=element_blank()) +
  
  # margins
  theme(plot.margin = margin(0,0,0,0,"in"),
        plot.background = element_blank())
#----

quartz()

# merge and save
#-----------------------------
fig7 <- 
fig7a / fig7b +
  plot_annotation(tag_levels = 'a')+
  plot_layout(heights = c(.55,.45)) &
  theme(plot.margin = margin(0,0,0,0))

ggsave(fig7,
       filename=here::here("figures","fig7.2.pdf"),
       dpi=300, 
       device = cairo_pdf,
       height = 10.5,
       width = 7.25)

width = 7.25
height = 9.5