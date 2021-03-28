# # Lawlor and Arellano 2020
# Scientific Reports

# Script 6. habitat suitability model selection (logistic regression)



# install libraries
#-------------------------------------
library(MuMIn) # for AICc
library(PNWColors)
#----

# load data
#-------------------------------------
load(here("data_processed","growthsummary.Rdata"))
load(here("data_processed","growthchem_means.Rdata"))
#----

# find treatments that "survived" (made it to >25% competence without dying)
#-------------------------------------
survtrts <- growthchem_means %>%
  select(Cup_f, Temp, Sal,pH,pCO2,Ar,DIC) %>% 
  mutate(survival = 0) %>% # add row for susvival
  mutate(survival = case_when(Cup_f %in% as.character(unique(growthsummary$Cup_f[growthsummary$Eyes >25])) ~ 1,
                              TRUE ~0)) # change survival to 1 for treatments that made it > 25% eyes
#----


# make list of all possible model combinations to test GLMs
#----------------------------------------
formulas <- list(# list all models here
  # again, always pH OR pCO2 OR Ar
  # never Ar and Salinity
  
  # NULL model
  survival ~ NULL,
  
  # Single predictor
  survival ~ Temp,
  survival ~ Sal,
  survival ~ pH,
  survival ~ pCO2,
  survival ~ Ar,
  survival ~ DIC,
  
  # two predictors
  survival ~ Temp + Sal,
  survival ~ Temp + pH,
  survival ~ Temp + pCO2,
  survival ~ Temp + Ar,
  survival ~ Sal + pH,
  survival ~ Sal + pCO2,
  
  # three additive
  survival ~ Temp + Sal + pH,
  survival ~ Temp + Sal + pCO2,
  #survival ~ Temp + Sal + Ar,
  
  
  #start quadratic models
  # add in arguments for quadratic temp and sal (because they span huge range)
  # single quadratic 
  survival ~ Temp + I(Temp^2),
  survival ~  Sal + I(Sal^2),
  
  # one linear one quadratic 
  survival ~ Temp + I(Temp^2) + Sal,
  survival ~ Temp + I(Temp^2) + pH,
  survival ~ Temp + I(Temp^2) + pCO2,
  survival ~ Temp + I(Temp^2) + Ar,
  survival ~  Sal + I(Sal^2) + pH,
  survival ~  Sal + I(Sal^2) + pCO2,
  survival ~ Sal + I(Sal^2) + Temp,
  
  survival ~ Temp + I(DIC^2) + DIC,
  
  
  
  # two linear one quadratic 
  survival ~  Temp + I(Temp^2) + Sal + pH,
  survival ~  Temp + I(Temp^2) + Sal + pCO2,
  survival ~  Sal + I(Sal^2) + Temp + pH,
  survival ~  Sal + I(Sal^2) + Temp + pCO2,
  survival ~  DIC + I(DIC^2) + Temp + pH,
  survival ~  DIC + I(DIC^2) + Temp + pCO2,
  
  # two quadratic 
  survival ~  Sal + I(Sal^2) + Temp + I(Temp^2),
  
  # two quadratic one linear
  survival ~  Sal + I(Sal^2) + Temp + I(Temp^2) + pH,
  survival ~  Sal + I(Sal^2) + Temp + I(Temp^2) + pCO2,
  survival ~  DIC + I(DIC^2) + Temp + I(Temp^2) + pCO2,
  
  # add interactions between temp/sal
  survival ~ Temp + Sal + Temp*Sal,
  survival ~ Temp + Sal + Temp*Sal + pH,
  survival ~ Temp + Sal + Temp*Sal + pCO2,
  survival ~ Temp + pH + Temp*pH,
  survival ~ Temp + pH + Temp*pH + Sal,
  survival ~ Temp + pCO2 + Temp*pCO2,
  survival ~ Temp + pCO2 + Temp*pCO2 + Sal,
  survival ~ Temp + Ar + Temp*Ar,
  survival ~ Sal + pH + Sal*pH,
  survival ~ Sal + pH + Sal*pH + Temp,
  survival ~ Sal + pCO2 + Sal*pCO2,
  survival ~ Sal + pCO2 + Sal*pCO2 + Temp,
  
  
  
  # quadratic with interaction
  survival ~ Temp + Sal + Temp*Sal + I(Temp^2),
  survival ~ Temp + Sal + Temp*Sal + I(Sal^2),
  survival ~ Temp + Sal + Temp*Sal + pH + I(Temp^2),
  survival ~ Temp + Sal + Temp*Sal + pH + I(Sal^2),
  survival ~ Temp + Sal + Temp*Sal + pCO2 + I(Temp^2),
  survival ~ Temp + Sal + Temp*Sal + pCO2 + I(Sal^2),
  survival ~ Temp + DIC + Temp*DIC + pCO2 + I(DIC^2),
  
  survival ~ Temp + I(Temp^2) + Sal + I(Sal^2) + Temp*Sal
  
)
#----




# test all the models in the list and compare
#----------------------------------
surv_selection <- purrr::map_df(formulas, ~{
  mod <- glm(.x, data= survtrts, family="binomial")
  data.frame(formula = format(.x), 
             AICc = round(AICc(mod),2), # AICc from MuMIn package adjusts AIC for small sample size
             BIC = round(BIC(mod),2),
             PseudoR = round(DescTools::PseudoR2(mod,which=c("McFaddenAdj")),2) # mcfaddens R2
  )
})

warnings() 
# this gives 10 warnings that probabilities of 1 or 0 occurred (for the top 10 models, presumably)
# this means our models are predicting somewhere that treatments may have 0% or 100% probability of surviving
# this warning message is explained in Venables & Ripley (2002, pp. 197–8)
# see excerpt here: https://stackoverflow.com/questions/8596160/why-am-i-getting-algorithm-did-not-converge-and-fitted-prob-numerically-0-or
# and Ripley's expansion on this here: http://math.yorku.ca/Who/Faculty/Monette/S-news/0027.html
# basically, because there is a small sample size of successful treatments and they are all
# lumped together in multidimensional space, our model predicts that there will be a 100% chance
# of survival in these conditions. This probably isn't true, biologically, but is a limitation
# of the sample size and number of model predictors. This suggests that if we predict probability for 
# temperatures further outside of the tested range will be inflated. It also suggests that
# our model is predicting perfect or near-perfect separation, which could be a red flag, but just
# looking at our own survivability results in Fig 4, we can see that near-perfect separation actually
# does exist in our dataset (you could draw a simple parabola around survivable treatments), so 
# we probably shouldn't change this: https://stats.stackexchange.com/questions/11109/how-to-deal-with-perfect-separation-in-logistic-regression
# for our puroses, this is ok, since we are not predicting above these treatment levels,
# and in fact, perfect separation in our experimental dataset actually did occur, but 
# should be kept in mind in the future if applying this model to anything else. 
# Still, the coefficients and deviations look okay, and this model does actually give the 
# best prediction of our experimental data, so we will keep these models despite the warning.  
#----


# rank models (supp table 1)
#---------------------------------------
# add rank columns ordered by formula order
surv_selection_table <- 
  surv_selection %>%  
  mutate('AICc Rank' = dense_rank(AICc),
         'BIC Rank' = dense_rank(BIC)) %>% 
  select(formula, AICc, 'AICc Rank', BIC, 'BIC Rank', PseudoR)

# save file - this will be supp table. added notations in inkscape. 
write.csv(surv_selection_table,file = here::here("data_processed","Surv_selection_table.csv"))
#----


# set best model
#----------------------------------------

survmod <- glm( survival ~ Sal + I(Sal^2) + Temp,
                family="binomial",
                data=survtrts)
# again, this gives the probabilities 0 or 1 warning.. see above for why we use anyway.
summary(survmod)
AICc(survmod)
glance(survmod)
#----



# make surv model prediction figure (fig 6a)
#------------------------------------
predicted.data <- data.frame(
  probability.survive = survmod$fitted.values,
  survive = survtrts$survival)
predicted.data <- predicted.data[order(predicted.data$probability.survive,decreasing=F),]
predicted.data$rank <- 1:nrow(predicted.data) 

# create classification table
theProbs <- fitted(survmod)
table(theProbs > .5,survtrts$survival) # 33 of 34 dead treatments were predicted correctly, and 15 out of 16 living predicted correctly

survtrts$survival <- as.logical(survtrts$survival)
quartz()
pal <- pnw_palette("Winter",n=2,type="discrete")
fig6a <- ggplot(data=predicted.data,aes(x=rank, y=probability.survive))+
  geom_point(aes(color=factor(survive)),alpha=1,shape=4,stroke=2, show.legend = F)+
  scale_color_manual(values=pal,labels=c("Unsuitable","Suitable"))+
  scale_y_continuous(breaks=c(0,1),labels = c("0%","100%")) +
  xlab("Index")+
  ylab("Predicted probability of habitat suitability")+
  geom_hline(yintercept = .50,linetype="dashed",size=.25)+
  labs(color="Experimental \nOutcome")+
  theme_classic()+
  theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    #  legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.position = c(.05,.95),
    legend.justification = c(0,1),
    legend.background = element_rect(color="black",size=.25),
    legend.key.height=unit(1.5,"line"),
  )+
  annotate("curve",
           x=30,y=.96,
           xend=38.25,yend=.915,
           curvature=-.18,
           size=.25)+
  annotate("label",
           x=30,y=.96,
           hjust=1,
           fill=pal[2],
           color="black",
           label="Suitable\nTreatments",
           lineheight=.8,
           size=3.5)+
  annotate("curve",
           x=40,y=.23,
           xend=33.5,yend=.19,
           curvature=.18,
           size=.25)+
  annotate("label",
           x=40,y=.23,
           hjust=0,
           fill=pal[1],
           color="grey90",
           label="Unsuitable\nTreatments",
           lineheight=.8,
           size=3.5)+
  annotate("text",
           x=1.7,y=.52,
           hjust=0,
           vjust=0,
           label="Predicted\nSuitable",
           size=3.5,
           lineheight=.8)+
  annotate("text",
           x=1.7,y=.48,
           hjust=0,
           vjust=1,
           label="Predicted\nUnsuitable",
           size=3.5,
           lineheight=.8) +
  annotate("segment",
           x=1,y=.48,
           xend=1,yend=.41,
           arrow = arrow(length = unit(0.015, "npc"),
                         type="closed",
                         angle=25))+
  annotate("segment",
           x=1,y=.52,
           xend=1,yend=.59,
           arrow = arrow(length = unit(0.015, "npc"),
                         type="closed",
                         angle=25)) +
  theme(panel.background = element_rect(color="black",fill="transparent"),
        plot.margin = margin(0,0,0,0,"in"))
#----

fig6a


# make survival prediction heatmap (fig 6b)
#-------------------------------------
# create dataset to simulate survival across conditions
survival_sim <- data.frame(matrix(ncol=3,nrow=560))
colnames(survival_sim) <- c("temp","sal","surv")
survival_sim$sal <- rep(12:39,each=20)
survival_sim$temp <- rep(11:30,times=28)

# create function that predicts probability of survivalable treatment (>25% of larvae will be succeed)
larvasurvive <- function(survtemp,survsal){
  survprob <- predict(survmod,list(Temp=survtemp,Sal=survsal),type="response")
  return(round(survprob*100,digits=2))
}

#predict survival probability using "larva survive" function
survival_sim$surv <- larvasurvive(survival_sim$temp,survival_sim$sal)
survival_sim$surv_threshold <- survival_sim$surv>.5

save(survival_sim,file=here("data_processed","survival_sim.Rdata"))
# Make figure 6b

pal <- pnw_palette("Winter",n=50,type="continuous")
fig6b<- survival_sim %>%
  ggplot(aes(x = sal, y = temp, z = surv, fill = surv)) + 
  # geom_tile() + 
  geom_raster()+
  #scale_fill_gradientn(colors=alpha(c("green","red"),.05)) +
  scale_fill_gradient(low ="#2d2926",high="#81a9ad", 
                      name="Likelihood of \nhabitat suitability",
                      guide = guide_colorbar(frame.colour = "black",ticks.colour = "black" ),
  ) +
  #scale_fill_gradientn(colors=pal) +
  ylab("Temperature (°C)")+
  xlab("Salinity (PSU)")+
  coord_cartesian(xlim = c(12, 40), ylim = c(10, 30))+
  theme_classic()+  
  stat_contour(inherit.aes=F,data=survival_sim,aes(x=sal,y=temp,z=surv),color="#ececec",alpha=.5,size=.25)+
  metR::geom_text_contour(inherit.aes=F, data=survival_sim,aes(x=sal,y=temp,z=surv),alpha=.75,color="#ececec")+
  theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    legend.background = element_rect(fill = "transparent") # get rid of legend bg
  )+
  annotate("text",label="McFadden's Adj. Pseudo R-sq:  0.67",x=12,y=11,hjust=0,vjust=0,color="white") +
  theme(legend.position = c(.05,.95),
        #legend.position = "top",
        legend.justification = c(0,1),
        legend.background = element_rect(fill=alpha("white",.95),color="black",size=.5),
        legend.title = element_text(size=8),
        legend.key.width  = unit(.5,'cm'),
        #legend.key.height = unit(.5,"cm"),
        legend.text = element_text(size=7),
        panel.background = element_rect(color="black",fill="transparent"),
        plot.margin = margin(0,0,0,0,"in")) +
  guides(fill = guide_colorbar(barwidth = unit(.15,"in"),
                               barheight = unit(2.75,"in"),
                               title.hjust = 0,
                               frame.colour = "black",
                               ticks.colour = "black"))
#----
 
 
 
 # merge and save
 #-----------------------
 Fig6_inset <- fig6a / fig6b + 
  plot_annotation(tag_levels = 'a') +
  plot_layout(heights = c(.45,.55)) &
  theme(plot.margin = margin(0,0,0,0))
quartz() 



ggsave(Fig6_inset,filename = here::here("figures","fig6.pdf"),
       device = cairo_pdf(),
       dpi=300,
       width = 7.25,
       height = 9.5)




