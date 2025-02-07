library(tidyverse)
library(lubridate)
library(reshape2)
library(viridis)
library(ggridges)
library(igraph)
library(patchwork)
library(brms)
library(ggbeeswarm)
library(modelr)
library(tidybayes)
library(ggalluvial)
library(ggplot2)
theme_set(theme_classic(base_size = 16))

setwd("~/Desktop/SPRF/Laikipia Fission/Drafts/Code/Laikipia_GH")

## LOAD DATA FRAMES
## PLOT 1 DATA
group.size<-read.csv("Population_History_Laikipia.csv")
group.dispersals<-read.csv("Laikipia_Dispersals.csv")

## NETWORK DATA
network.measures<-read.csv("Network_Measures_Laikipia.csv")

## DYADIC GROOMING DATA
ff.dyads<-read.csv("FF_Dyads_Laikipia.csv")
fm.dyads<-read.csv("FM_Dyads_Laikipia.csv")

## AGGRESSION DATA
agg.by.id<-read.csv("Agg_Data_ID.csv")
agg.by.dyad<-read.csv("Agg_Data_Dyadic.csv")
agg.by.event<-read.csv("Agg_Data_Events.csv")
agg.by.event.intra<-read.csv("Agg_Data_Intra.csv")
agg.by.event.inter<-read.csv("Agg_Data_Inter.csv")

## INTER-GROUP GROOMING DATA
intergroup.grooming<-read.csv("Intergroup_Grooming.csv")

## MAKE LARGER POPULATION HISTORY FIGURE - PLOT 1a
colors<-c("gray","#440154FF","#21908CFF","#FDE725FF")
group.size$year.date<-as.Date(group.size$year.date, format="%m/%d/%y")
group.size$group<-as.factor(group.size$group)
group.size$group<-factor(group.size$group, levels=c("OGS","PHG", "ENK", "YNT"))
plot1a<-ggplot(group.size,
              aes(x = year.date, alluvium = group, y=groupsize,
                  fill = group)) +
  scale_fill_manual(values=colors) +
  scale_color_manual(values=colors) +
  geom_flow(aes(fill=group, color=group),alpha=0.6, curve_type="sigmoid", curve_range=0.1, width=0) +
  theme(legend.position = "none") + ylim(0,80) + ylab("# of adults and subadults\n(both sexes)") + theme(plot.title=element_text(hjust=0.5, vjust=0.1, face='bold')) +
  annotate("rect", xmin=as.Date("2013-11-01"), xmax=as.Date("2017-12-05"), ymin=0, ymax=80, alpha=0.6, fill="darkgray", color="black") +
  annotate("text", x=as.Date("2015-11-01"), y=65, label= "Study period",size=7) + xlab("Year")
plot1a

## MAKE POPULATION HISTORY FIGURES FOR STUDY PERIOD - PLOT 1b, 1c
colors<-c("#21908CFF","#FDE725FF","#440154FF")
group.dispersals$group<-factor(group.dispersals$group,levels=c("ENK","YNT","PHG"))
group.dispersals$period<-factor(group.dispersals$period,levels=c("Stable","Dispersal 1","Dispersal 2", "Fission"))

plot1b<-ggplot(group.dispersals[group.dispersals$sex=="f",],
              aes(x = period, stratum = group, alluvium = as.character(id),
                  fill = group, label = group)) +
  scale_fill_manual(values=colors) +
  geom_stratum(alpha=0.6) +
  geom_flow(alpha=0.3, curve_type = "quintic") +
  geom_text(stat = "stratum", size = 3) + xlab("Study period") +
  theme(legend.position = "none") + ylim(0,40) + ylab("# of adults and subadults\n(sex-specific)") +
  ggtitle("females") + theme(plot.title=element_text(hjust=0.5, vjust=0.1, face='bold'))

plot1c<-ggplot(group.dispersals[group.dispersals$sex=="m",],
              aes(x = period, stratum = group, alluvium = as.character(id),
                  fill = group, label = group)) +
  scale_fill_manual(values=colors) +
  geom_stratum(alpha=0.6, na.rm=F) +
  geom_flow(alpha=0.3, curve_type = "quintic") +
  geom_text(stat = "stratum", size = 3) + xlab("Study period") +
  theme(legend.position = "none") + ylim(0,20) + ylab("") +
  ggtitle("males") + theme(plot.title=element_text(hjust=0.5, vjust=0.1, face='bold'))

plot1a/(plot1b + plot1c) + plot_annotation(tag_levels = "a")

## MAKE NETWORK MEASURE PLOTS - PLOT 2a,2b
network.measures$Period<-factor(network.measures$Period, levels=c("Stable","Dispersal 1", "Dispersal 2", "Fission"))
plot2a<-ggplot(data=network.measures, aes(x=Period, y=Density, group=Group, color=Group)) + 
  geom_line(lwd=2.2) + geom_point(size=4) + ylim(0,1) + xlab("") +
  scale_color_viridis(discrete = T) + theme(legend.position = c(0.65,0.9)) 

plot2b<-ggplot(data=network.measures, aes(x=Period, y=Modularity, group=Group, color=Group)) + 
  geom_line(lwd=2.2) + geom_point(size=4) + ylim(0,0.5) + geom_hline(yintercept = 0.3, cex=1, lty=2) +
  scale_color_viridis(discrete = T) + theme(legend.position = "none")
plot2a/plot2b + plot_annotation(tag_levels = "a")

## MODEL FF GROOMING DYADS
ff.dyads$outcome<-factor(ff.dyads$outcome, levels=c("stayed together", "split apart", "briefly together", "never together"))
brm_ff_groom<-brm(data=ff.dyads, 
             Grooming ~  outcome + period + group.combo + scale(rval) + scale(RankDiff) + offset(log(Obs.effort))+ 
               (1|mm(id1,id2)), chains=4, cores=2, iter=2000, family="hurdle_gamma",
             control = list(max_treedepth=12,adapt_delta=0.99))
summary(brm_ff_groom, prob=0.89)

brm_ff_intx<-brm(data=ff.dyads, 
              Grooming ~  outcome*period + group.combo + scale(rval) + scale(RankDiff) + offset(log(Obs.effort))+ 
                (1|mm(id1,id2)), chains=4, cores=2, iter=2000, family="hurdle_gamma",
              control = list(max_treedepth=12,adapt_delta=0.99))
summary(brm_ff_intx, prob=0.89)
conditional_effects(brm_ff_intx, "outcome:period")

## VISUALIZE FF EFFECTS
palette1<-c("#D1B79EFF","#7EBAC2FF") 
palette2<-c("#403369FF", "#5C5992FF", "#AE93BEFF", "#B4DAE5FF")
## MAKE GROUP COMBO PLOT
grid = ff.dyads %>%
  data_grid(rval=mean(ff.dyads$rval), RankDiff=mean(ff.dyads$RankDiff), group.combo, period, outcome, id1, id2, Obs.effort=60)

means = grid %>%
  add_epred_draws(brm_ff_groom, ndraws=250, allow_new_levels=T)

check<-means %>%
  group_by(group.combo, .draw) %>%
  dplyr::summarise(Mean.Pred=mean(.epred))

overall<-check %>%
  group_by(group.combo) %>%
  dplyr::summarise(Grand.Mean=mean(Mean.Pred), Min=min(Mean.Pred), Max=max(Mean.Pred), 
                   Pred.25=quantile(Mean.Pred, probs=0.25), Pred.75=quantile(Mean.Pred, probs=0.75),
                   LL.Pred=quantile(Mean.Pred, probs=0.055), UL.Pred=quantile(Mean.Pred, probs=0.945))

plot3a<-ggplot() + 
  geom_violin(aes(y = Mean.Pred, x= group.combo, fill=group.combo), data = check, alpha = 0.9, trim=F) +
  scale_fill_manual(values=palette1) +
  geom_point(aes(y = Grand.Mean, x= group.combo), data = overall, alpha = 1, size=2, color = "black") +
  geom_errorbar(aes(ymin = Pred.25, ymax=Pred.75, x= group.combo), data = overall, width=0, linewidth=1.6, color = "black") +
  geom_errorbar(aes(ymin = LL.Pred, ymax=UL.Pred, x= group.combo), data = overall, width=0, linewidth=0.8, color = "black") +
  ylab("female-female grooming\n(minutes/hour)") + xlab("current group") + theme(legend.position ="none") + ylim(0,0.15)
plot3a

## MAKE OUTCOME PLOT
check<-means %>%
  group_by(outcome, .draw) %>%
  dplyr::summarise(Mean.Pred=mean(.epred))

overall<-check %>%
  group_by(outcome) %>%
  dplyr::summarise(Grand.Mean=mean(Mean.Pred), Min=min(Mean.Pred), Max=max(Mean.Pred), 
                   Pred.25=quantile(Mean.Pred, probs=0.25), Pred.75=quantile(Mean.Pred, probs=0.75),
                   LL.Pred=quantile(Mean.Pred, probs=0.055), UL.Pred=quantile(Mean.Pred, probs=0.945))

plot3b<-ggplot() + 
  geom_violin(aes(y = Mean.Pred, x= outcome, fill=outcome), data = check, alpha = 0.9, trim=F) +
  scale_fill_manual(values=palette2) +
  geom_point(aes(y = Grand.Mean, x= outcome), data = overall, alpha = 1, size=2, color = "black") +
  geom_errorbar(aes(ymin = Pred.25, ymax=Pred.75, x= outcome), data = overall, width=0, linewidth=1.6, color = "black") +
  geom_errorbar(aes(ymin = LL.Pred, ymax=UL.Pred, x= outcome), data = overall, width=0, linewidth=0.8, color = "black") +
  ylab("") + xlab("eventual outcome") + theme(legend.position ="none") + ylim(0,0.15)
plot3b

## MODEL FM GROOMING DYADS
fm.dyads$outcome<-factor(fm.dyads$outcome, levels=c("stayed together", "split apart", "briefly together", "never together"))
brm_fm_groom<-brm(data=fm.dyads, 
               Grooming~  outcome + period + group.combo + scale(Male.Rank)*scale(Fem.Rank) + offset(log(Obs.effort)) +     
                 (1|mm(id1,id2)), chains=4, cores=2, iter=4000, family="hurdle_gamma",
               control = list(adapt_delta=0.99))
summary(brm_fm_groom, prob=0.89)

## VISUALIZE FM EFFECTS
## MAKE GROUP COMBO PLOT
grid = fm.dyads %>%
  data_grid(group.combo, outcome, id1, id2, period, Male.Rank=mean(fm.dyads$Male.Rank), Fem.Rank=mean(fm.dyads$Fem.Rank), Obs.effort=60)

means = grid %>%
  add_epred_draws(brm_fm_groom, ndraws=250, allow_new_levels=T)

check<-means %>%
  group_by(group.combo, .draw) %>%
  dplyr::summarise(Mean.Pred=mean(.epred))

overall<-check %>%
  group_by(group.combo) %>%
  dplyr::summarise(Grand.Mean=mean(Mean.Pred), Min=min(Mean.Pred), Max=max(Mean.Pred), 
                   Pred.25=quantile(Mean.Pred, probs=0.25), Pred.75=quantile(Mean.Pred, probs=0.75),
                   LL.Pred=quantile(Mean.Pred, probs=0.055), UL.Pred=quantile(Mean.Pred, probs=0.945))

plot3c<-ggplot() + 
  geom_violin(aes(y = Mean.Pred, x= group.combo, fill=group.combo), data = check, alpha = 0.9, trim=F) +
  scale_fill_manual(values=palette1) +
  geom_point(aes(y = Grand.Mean, x= group.combo), data = overall, alpha = 1, size=2, color = "black") +
  geom_errorbar(aes(ymin = Pred.25, ymax=Pred.75, x= group.combo), data = overall, width=0, linewidth=1.6, color = "black") +
  geom_errorbar(aes(ymin = LL.Pred, ymax=UL.Pred, x= group.combo), data = overall, width=0, linewidth=0.8, color = "black") +
  ylab("female-male grooming\n(minutes/hour)") + xlab("current group") + theme(legend.position ="none") + ylim(0,0.2)
plot3c

## MAKE OUTCOME PLOT
check<-means %>%
  group_by(outcome, .draw) %>%
  dplyr::summarise(Mean.Pred=mean(.epred))

overall<-check %>%
  group_by(outcome) %>%
  dplyr::summarise(Grand.Mean=mean(Mean.Pred), Min=min(Mean.Pred), Max=max(Mean.Pred), 
                   Pred.25=quantile(Mean.Pred, probs=0.25), Pred.75=quantile(Mean.Pred, probs=0.75),
                   LL.Pred=quantile(Mean.Pred, probs=0.055), UL.Pred=quantile(Mean.Pred, probs=0.945))

plot3d<-ggplot() + 
  geom_violin(aes(y = Mean.Pred, x= outcome, fill=outcome), data = check, alpha = 0.8, trim=F) +
  scale_fill_manual(values=palette2) +
  geom_point(aes(y = Grand.Mean, x= outcome), data = overall, alpha = 1, size=2, color = "black") +
  geom_errorbar(aes(ymin = Pred.25, ymax=Pred.75, x= outcome), data = overall, width=0, linewidth=1.6, color = "black") +
  geom_errorbar(aes(ymin = LL.Pred, ymax=UL.Pred, x= outcome), data = overall, width=0, linewidth=0.8, color = "black") +
  ylab("") + xlab("eventual outcome") + theme(legend.position ="none") + ylim(0,0.2)
plot3d

## COMBINE GROOMING FIGURES
layout<-"AABBBB
         CCDDDD"
plot3a+plot3b+
  plot3c+plot3d + plot_layout(design=layout) + plot_annotation(tag_levels = "a")

## AGGRESSION ANALYSES
## ANALYSIS 1: AGGRESSION BY INDIVIDUAL
agg.by.id$period<-factor(agg.by.id$period, levels=c("Stable","Dispersal 1", "Dispersal 2","Fission"))
brm_agg_id<-brm(data=agg.by.id, 
                N.Agg ~ period + dispersals + groupsize + offset(log(Obs.Effort)) + 
                  (1|id), chains=4, cores=2, iter=4000, family="zero_inflated_poisson",
                control = list(max_treedepth=12,adapt_delta=0.99))
summary(brm_agg_id, prob=0.89)

## MAKE DISPERSALS PLOT
palette1<-c("#D1B79EFF","#7EBAC2FF") 
palette2<-c("#403369FF", "#5C5992FF", "#AE93BEFF", "#B4DAE5FF")
grid = agg.by.id %>%
  data_grid(dispersals, groupsize, period="Stable", id, Obs.Effort=1)

means = grid %>%
  add_predicted_draws(brm_agg_id, ndraws=250, allow_new_levels=T)

check<-means %>%
  group_by(dispersals, .draw) %>%
  dplyr::summarise(Mean.Pred=mean(.prediction))

overall<-check %>%
  group_by(dispersals) %>%
  dplyr::summarise(Grand.Mean=mean(Mean.Pred), Min=min(Mean.Pred), Max=max(Mean.Pred), 
                   Pred.25=quantile(Mean.Pred, probs=0.25), Pred.75=quantile(Mean.Pred, probs=0.75),
                   LL.Pred=quantile(Mean.Pred, probs=0.055), UL.Pred=quantile(Mean.Pred, probs=0.945))

plot1a<-ggplot() + 
  geom_violin(aes(y = Mean.Pred, x= dispersals, fill=dispersals), data = check, alpha = 0.9, trim=F) +
  scale_fill_brewer(palette="Blues") +
  geom_point(aes(y = Grand.Mean, x= dispersals), data = overall, alpha = 1, size=2, color = "black") +
  geom_errorbar(aes(ymin = Pred.25, ymax=Pred.75, x= dispersals), data = overall, width=0, linewidth=1.6, color = "black") +
  geom_errorbar(aes(ymin = LL.Pred, ymax=UL.Pred, x= dispersals), data = overall, width=0, linewidth=0.8, color = "black") +
  ylab("") + xlab("Imminent dispersals?") + theme(legend.position ="none") + ylim(0,4)
plot1a

## MAKE PERIOD PLOT
grid = agg.by.id %>%
  data_grid(dispersals="no", groupsize, period, id, Obs.Effort=1)

means = grid %>%
  add_predicted_draws(brm_agg_id, ndraws=250, allow_new_levels=T)

check<-means %>%
  group_by(period, .draw) %>%
  dplyr::summarise(Mean.Pred=mean(.prediction))

overall<-check %>%
  group_by(period) %>%
  dplyr::summarise(Grand.Mean=mean(Mean.Pred), Min=min(Mean.Pred), Max=max(Mean.Pred), 
                   Pred.25=quantile(Mean.Pred, probs=0.25), Pred.75=quantile(Mean.Pred, probs=0.75),
                   LL.Pred=quantile(Mean.Pred, probs=0.055), UL.Pred=quantile(Mean.Pred, probs=0.945))

plot1b<-ggplot() + 
  geom_violin(aes(y = Mean.Pred, x= period, fill=period), data = check, alpha = 0.9, trim=F) +
  scale_fill_brewer(palette="Purples") +
  geom_point(aes(y = Grand.Mean, x= period), data = overall, alpha = 1, size=2, color = "black") +
  geom_errorbar(aes(ymin = Pred.25, ymax=Pred.75, x= period), data = overall, width=0, linewidth=1.6, color = "black") +
  geom_errorbar(aes(ymin = LL.Pred, ymax=UL.Pred, x= period), data = overall, width=0, linewidth=0.8, color = "black") +
  ylab("aggression received\n(events/hour)") + xlab("Study period") + theme(legend.position ="none") + ylim(0,4)
plot1b

layout<-"AAABB"
plot1b+plot1a+plot_annotation(tag_levels = "a") + plot_layout(design = layout)

## ANALYSIS 2: DYADIC-LEVEL
agg.by.dyad$period<-factor(agg.by.dyad$period, levels=c("Stable","Dispersal 1", "Dispersal 2","Fission"))
brm_agg<-brm(data=agg.by.dyad, 
             N.Agg ~  period + outcome + group.combo + scale(rval) + scale(RankDiff) + offset(log(Obs.effort))+ 
               (1|mm(id1,id2)), chains=4, cores=2, iter=2000, family="zero_inflated_poisson",
             control = list(max_treedepth=12,adapt_delta=0.99))
summary(brm_agg, prob=0.89)

## ANALYSIS 3: EVENT-LEVEL
agg.by.event$period<-factor(agg.by.event$period, levels=c("Stable","Dispersal 1", "Dispersal 2","Fission"))
brm_intensity<-brm(data=agg.by.event, 
                   Intensity ~  period + outcome + rval + RankDiff +
                     (1|mm(id1,id2)), chains=4, cores=2, iter=2000, family="bernoulli",
                   control = list(max_treedepth=12,adapt_delta=0.99))
summary(brm_intensity, prob=0.89)

## ANALYSIS 4: EVENT-LEVEL - INTRAGROUP
agg.by.event.intra$period<-factor(agg.by.event.intra$period, levels=c("Stable","Dispersal 1", "Dispersal 2","Fission"))
brm_winner_intra<-brm(data=agg.by.event.intra, 
                      win ~  relative.size + period + (1|mm(focal,partner)), chains=4, cores=2, iter=2000, family="bernoulli",
                      control = list(max_treedepth=12,adapt_delta=0.99))
summary(brm_winner_intra,prob=0.89)

## ANALYSIS 5: EVENT-LEVEL - INTERGROUP
agg.by.event.inter$period<-factor(agg.by.event.inter$period, levels=c("Stable","Dispersal 1", "Dispersal 2","Fission"))
brm_winner_inter<-brm(data=agg.by.event.inter,
                win ~  relative.size*period + (1|mm(focal,partner)), chains=4, cores=2, iter=2000, family="bernoulli",
                control = list(max_treedepth=12,adapt_delta=0.99))
summary(brm_winner_inter,prob=0.89)

## INTERGROUP GROOMING
intergroup.grooming$Period<-factor(intergroup.grooming$Period, levels=c("Stable","Dispersal 1", "Dispersal 2","Fission"))
intergroup.grooming$outcome<-factor(intergroup.grooming$outcome, levels=c("non-disperser", "disperser"))

brm_inter_groom<-brm(data=intergroup.grooming, 
                     Between.Dur|trials(Groom.Duration) ~ outcome*Period + outcome*sex.combo + (1|focal), chains=4, cores=2, 
                     iter=2000, family="zero_inflated_binomial",
                     control = list(max_treedepth=12,adapt_delta=0.99))
print(summary(brm_inter_groom), digits=2, prob=0.89)
print(summary(brm_group_inter), digits=2, prob=0.89)

conditional_effects(brm_inter_groom, "outcome:Period")
conditional_effects(brm_group_inter, "outcome:Period")

## PLOT INTER-GROUP EFFECT
library(modelr)
library(tidybayes)
grid = intergroup.grooming %>%
  data_grid(focal, outcome, sex.combo, Period, Groom.Duration=1)

means = grid %>%
  add_epred_draws(brm_inter_groom, ndraws=100, allow_new_levels=T)

check<-means %>%
  group_by(outcome, Period, .draw) %>%
  dplyr::summarise(Mean.Pred=mean(.epred))

range1<-means %>%
  group_by(outcome, Period, .draw) %>%
  dplyr::summarise(Mean.Pred=mean(.epred))

overall<-check %>%
  group_by(outcome, Period) %>%
  dplyr::summarise(Mean.Pred=mean(Mean.Pred))

summary.range1<-range1 %>%
  group_by(outcome, Period) %>%
  dplyr::summarise(LL=quantile(Mean.Pred,0.055), Med=quantile(Mean.Pred,0.5), UL=quantile(Mean.Pred,0.945))

palette3<-c("#BCBDDC", "#9E9AC8", "#756BB1", "#54278F")

summary.range1$Period<-as.factor(summary.range1$Period)
check$Period<-as.factor(check$Period)
overall$Period<-as.factor(overall$Period)
intergroup.grooming$Period<-as.factor(intergroup.grooming$Period)

levels(summary.range1$Period)<-c("Stable","Dispersal 1","Dispersal 2","Fission")
levels(check$Period)<-c("Stable","Dispersal 1","Dispersal 2","Fission")
levels(overall$Period)<-c("Stable","Dispersal 1","Dispersal 2","Fission")
levels(intergroup.grooming$Period)<-c("Stable","Dispersal 1","Dispersal 2","Fission")

library(ggbeeswarm)
plot1a<-ggplot() + 
  scale_color_manual(values = palette3) + labs(size="Observation hours\nper individidual") +
  geom_quasirandom(data=intergroup.grooming, aes(y=Between.Prop,x=Period,size=log(Groom.Duration)/20, color=Period),alpha=0.3, width=0.4) +
  geom_errorbar(data=summary.range1,aes(x=Period,ymin=LL,ymax=UL), width=0.1, lwd=0.8) +
  geom_point(data=summary.range1,aes(x=Period,y=Med), cex=2) + scale_size_continuous(range = c(0,3)) +
  facet_grid(cols=vars(outcome)) +
  ylab("proportion of grooming effort\ndirected outside group") + xlab("study period") + theme(legend.position ="none")
plot1a

##
check<-means %>%
  group_by(sex.combo, outcome, .draw) %>%
  dplyr::summarise(Mean.Pred=mean(.epred))

range1<-means %>%
  group_by(sex.combo, outcome, .draw) %>%
  dplyr::summarise(Mean.Pred=mean(.epred))

overall<-check %>%
  group_by(sex.combo, outcome) %>%
  dplyr::summarise(Mean.Pred=mean(Mean.Pred))

summary.range1<-range1 %>%
  group_by(sex.combo, outcome,) %>%
  dplyr::summarise(LL=quantile(Mean.Pred,0.05), Med=quantile(Mean.Pred,0.5), UL=quantile(Mean.Pred,0.95))

library(RColorBrewer)

plot1b<-ggplot() + 
  scale_color_brewer(palette =  "Accent") + labs(size="Observation hours\nper individidual") +
  geom_quasirandom(data=intergroup.grooming, aes(y=Between.Prop,x=sex.combo,size=log(Groom.Duration)/20, color=sex.combo),alpha=0.7, width=0.4) +
  geom_errorbar(data=summary.range1,aes(x=sex.combo,ymin=LL,ymax=UL), width=0.1, lwd=0.8) +
  geom_point(data=summary.range1,aes(x=sex.combo,y=Med), cex=2) + scale_size_continuous(range = c(0,3)) +
  facet_grid(cols=vars(outcome)) +
  ylab("proportion of grooming effort\ndirected outside group") + xlab("sex pairing") + theme(legend.position ="none")
plot1b

plot1a/plot1b + plot_annotation(tag_levels = "a")

