# # Lawlor and Arellano 2020
# Scientific Reports

# Script 3. make figure of every cup growth rate 



# load data
#------------------------------------------
load(here("data_processed","growthsummary.Rdata"))
load(here("data_processed","regressions.Rdata"))
load(here("data_processed","sampledays.Rdata"))
load(here("data_processed","successdf.Rdata"))

#----



# make figure of growth regressions
#---------------------------------------------
pal <- pnw_palette("Bay",50,type="continuous") # create palette


fig4 <- ggplot(growthsummary,aes(x=day,y=um)) +
  geom_point() +
  #color backgrounds
  geom_rect(aes(fill = Gr, 
                xmin = -Inf, xmax=Inf, 
                ymin = -Inf, ymax = Inf),
            regressions,
            inherit.aes = FALSE) +
  scale_fill_gradientn(colors = pal, name = "Growth Rate\n(µm/day)",
                       guide = guide_colorbar(frame.colour = "black", 
                                              ticks.colour = "black",
                                              order=1))+
  
  # draw dotted line
  geom_hline(yintercept = 260,linetype="dashed",size=.3)+
  
  # plot size measurements with SE
  geom_errorbar(aes(ymin=um-se,ymax=um+se),width=.1)+
  geom_line() +
  geom_point(aes(),size=1.5,shape=18) +
  
  # draw regression line
  geom_abline(data=regressions, 
              aes(slope=Gr,
                  intercept=intercept),
              linetype="solid",color="white",alpha=.75)+
  
  # draw death line
  geom_vline(data=   sampledays[sampledays$days.before.mort<=sampledays$days.before.comp
                                & is.na(successdf$successtime),],
             aes( xintercept = sampleday), color="black",size=.4,linetype="dotted")+
  
  # add competence points
  geom_point(data=growthsummary[growthsummary$Eyes>=25,],
             inherit.aes=F,
             aes(x=day,y=um,color=Eyes),
             size=2.5,shape=16)+
  scale_color_continuous(low="grey90",high = "black",name="% Eyed",guide= guide_legend(),limits=c(5,100),breaks=c(25,50,75,100))+
  #add black background on competence points
  geom_point(data=growthsummary[growthsummary$Eyes>=25,],inherit.aes=F,aes(x=day,y=um),size=2.5,shape=21,color="black",fill="transparent")+
  # add competence tim
  geom_text(data=successdf,inherit.aes = F, x =1, y = 300,hjust="center",size=3.5,
            aes(label= ifelse(is.na(successtime),"",paste(successtime))))+ 
  
  # facet
  facet_wrap(~Cup_f,nrow=5) +
  xlab("Day") + 
  ylab(expression(paste("Size (",mu,"m)"))) +
  
  # theme stuff
  theme_bw()+
  theme(strip.background = element_rect(fill="white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.box.margin = margin(0,0,0,30),
        plot.margin = margin(35,0,5,5))+
  
  # legend positions
  guides(color = guide_legend(order = 2))
#----


# save figure
#----------------------------------------------
ggsave(fig4, filename = here::here("figures","fig4.pdf"),
       device=cairo_pdf,
       width = 11,
       height = 6.7,
       units=c("in"),
       dpi=300)
# still need to add extra annotations to this figure in inkscape





