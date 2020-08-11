# # Lawlor and Arellano 2020
# Scientific Reports

# Script 2. make figure of size and eyespots through time


# Load data
load(here("data_processed","growth_trimmed.Rdata"))
load(here("data_processed","mort_trimmed.Rdata"))
load(here("data_processed","growthsummary.Rdata"))

# Jitter plot of sizes over time (Fig 3a)
#------------------------------------
fig3a <- growth_trimmed %>%
  
  filter(!(day==0 & Cup != "A1")) %>% # first sample is repeated for every cup. remove all but one
  
  mutate(Eyes = case_when(Eyes == "" ~ "Larva without\neyespot",
                          Eyes == "y" ~ "Larva with\neyespot")) %>%
  
  ggplot(aes(x=day,y=um,fill=Eyes,size=Eyes)) +
  geom_jitter(shape=21,alpha =.7) +
  scale_fill_manual(values= rev(c("black","red"))) +
  scale_size_manual(values=rev(c(1.5, 2))) +
  xlab("Day") + 
  ylab(expression(paste("Size (", mu,"m)")))+
  coord_cartesian(ylim = c(125,360),xlim=c(-1.5,18.5),expand = F)+
  geom_hline(yintercept = 260, linetype="dashed",size=.5)+
  theme_classic() +
  theme(legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.background = element_rect(color="transparent",fill="transparent"),
        legend.spacing.y = unit(2,"pt"),
        legend.title = element_blank(),
        legend.text = element_text(lineheight = .8),
        legend.key.height=unit(1.5,"line"),
        legend.key.width=unit(.25,"line")) +
  theme(panel.background = element_rect(color="black",fill="transparent",size=.5))
#----



# Average size over time (fig 3b)
#----------------------------------
fig3b <- 
  growthsummary %>%  
  mutate(MostEyed = case_when(Eyes < 25 ~ "<25% eyed",
                              Eyes >= 25 ~ ">25% eyed")) %>%
  
  # start plot
  ggplot( aes(x=day, y=um, group=Cup)) + 
  
  geom_line() +
  geom_errorbar(aes(ymin=um-se, ymax=um+se), width=.1) +
  geom_point(shape=21,aes(fill=MostEyed, size = MostEyed)) +
  
  # control size and fill
  scale_fill_manual(values= c("red","black"),
                    breaks = c(">25% eyed","<25% eyed")) +
  scale_size_manual(values=c(2, 1.5),
                    breaks = c(">25% eyed","<25% eyed")) +
  
  # draw hline
  geom_hline(yintercept = 260, linetype="dashed",size=.5)+
  
  # labs
  labs(x="Day",y = element_blank())+
  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())+
  
  # theme stuff
  theme_classic() +
  theme(legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.background = element_rect(color="transparent",fill="transparent"),
        legend.spacing.y = unit(2,"pt"),
        legend.title = element_blank(),
        legend.text = element_text(lineheight = .8),
        legend.key.height=unit(1.5,"line"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.key.width=unit(.25,"line"))+
  
  # coords
  coord_cartesian(ylim = c(125,360),xlim=c(-1.5,18.5),expand = F) +
  theme(panel.background = element_rect(color="black",fill="transparent",size=.5))

#----




# combine and export
#----------------------------
fig3 <- fig3a + fig3b + 
  plot_annotation(tag_levels = 'a') &
  theme(plot.margin = margin(0,0,0,0,unit = "pt"))

ggsave(filename=here::here("figures","Fig3.pdf"),
       device=cairo_pdf,
       width=8.5,
       height=5.25)

rm(fig3a,fig3b)
