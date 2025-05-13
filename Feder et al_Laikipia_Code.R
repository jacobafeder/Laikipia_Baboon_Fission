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

## Q1 + Q2
## DYADIC GROOMING DATA
ff.dyads<-read.csv("FF_Dyads_Laikipia.csv")
fm.dyads<-read.csv("FM_Dyads_Laikipia.csv")

## SAMPLE SIZES
length(unique(c(ff.dyads$id1,ff.dyads$id2)))
length(unique(c(ff.dyads$Dyad)))
length(unique(c(fm.dyads$id1,fm.dyads$id2)))
length(unique(c(fm.dyads$Dyad)))
length(unique(group.dispersals$id[group.dispersals$sex=="f"]))
length(unique(group.dispersals$id[group.dispersals$sex=="m"]))

## Q3
## AGGRESSION DATA
agg.by.id<-read.csv("Agg_Data_ID.csv")
agg.by.dyad<-read.csv("Agg_Data_Dyadic.csv")
agg.by.event<-read.csv("Agg_Data_Events.csv")
agg.by.event.inter<-read.csv("Agg_Data_Inter.csv")

## Q4
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

## MODEL FF GROOMING DYADS
ff.dyads$outcome<-factor(ff.dyads$outcome, levels=c("stayed together", "split apart", "briefly together", "never together"))
ff.dyads$period<-factor(ff.dyads$period, levels=c("Stable", "Dispersal 1", "Dispersal 2", "Fission"))

## MODEL MAIN EFFECTS -- HURDLE GAMMA
## MODEL 1A CO-RESIDENCE
brm_ff_groom1<-brm(data=ff.dyads, bf(Grooming ~  period + group.combo  + scale(rval) + scale(RankDiff) + offset(log(Obs.effort)) + (1|mm(id1,id2)), 
                                    hu ~   period + group.combo + scale(rval) + scale(RankDiff) + log(Obs.effort) + (1|mm(id1,id2))),
              chains=4, cores=2, iter=4000, family="hurdle_gamma",
             control = list(max_treedepth=12,adapt_delta=0.99))
summary(brm_ff_groom1, prob=0.89)
conditional_effects(brm_ff_groom1, "group.combo", prob=0.89)

## MODEL 1B OUTCOME
brm_ff_groom2<-brm(data=ff.dyads, bf(Grooming ~  period + outcome  + scale(rval) + scale(RankDiff) + offset(log(Obs.effort)) + (1|mm(id1,id2)), 
                                     hu ~   period + outcome + scale(rval) + scale(RankDiff) + log(Obs.effort) + (1|mm(id1,id2))),
                   chains=4, cores=2, iter=4000, family="hurdle_gamma",
                   control = list(max_treedepth=12,adapt_delta=0.99))
summary(brm_ff_groom2, prob=0.89)
conditional_effects(brm_ff_groom2, "outcome", prob=0.89)

## MODEL 1C CO-RESIDENCE + OUTCOME 
brm_ff_groom3<-brm(data=ff.dyads, bf(Grooming ~  period + group.combo + outcome  + scale(rval) + scale(RankDiff) + offset(log(Obs.effort)) + (1|mm(id1,id2)), 
                                     hu ~   period + group.combo + outcome + scale(rval) + scale(RankDiff) + log(Obs.effort) + (1|mm(id1,id2))),
                   chains=4, cores=2, iter=4000, family="hurdle_gamma",
                   control = list(max_treedepth=12,adapt_delta=0.99))
summary(brm_ff_groom3, prob=0.89)
conditional_effects(brm_ff_groom3, "outcome", prob=0.89)
conditional_effects(brm_ff_groom3, "group.combo", prob=0.89)

loo1<-loo(brm_ff_groom1)
loo2<-loo(brm_ff_groom2)
loo3<-loo(brm_ff_groom3)
loo_compare(loo1,loo2,loo3)

model_weights(brm_ff_groom1, brm_ff_groom2, brm_ff_groom3)

## MODEL INTX EFFECTS FOR MODEL 2 -- HURDLE GAMMA
brm_ff_intx2<-brm(data=ff.dyads, bf(Grooming ~  outcome*period + scale(rval) + scale(RankDiff) + offset(log(Obs.effort)) + (1|mm(id1,id2)), 
                                   hu ~  outcome*period + scale(rval) + scale(RankDiff) + log(Obs.effort) + (1|mm(id1,id2))),
                 chains=4, cores=2, iter=4000, family="hurdle_gamma",
                 control = list(max_treedepth=12,adapt_delta=0.99))
summary(brm_ff_intx2, prob=0.89)
conditional_effects(brm_ff_intx2, "outcome", prob=0.89)
conditional_effects(brm_ff_intx2, "period")
conditional_effects(brm_ff_intx2, "rval")

loo1<-loo(brm_ff_groom2, pointwise = TRUE)
loo2<-loo(brm_ff_intx2, pointwise = TRUE)
loo_compare(loo1,loo2)

## MODEL INTX EFFECTS FOR MODEL 3 -- HURDLE GAMMA
brm_ff_intx3<-brm(data=ff.dyads, bf(Grooming ~  outcome*period + group.combo + scale(rval) + scale(RankDiff) + offset(log(Obs.effort)) + (1|mm(id1,id2)), 
                                    hu ~  outcome*period + group.combo + scale(rval) + scale(RankDiff) + log(Obs.effort) + (1|mm(id1,id2))),
              chains=4, cores=2, iter=4000, family="hurdle_gamma", save_pars = save_pars(all=T),
             control = list(max_treedepth=12,adapt_delta=0.99))
summary(brm_ff_intx3, prob=0.89)
conditional_effects(brm_ff_intx3, "outcome", prob=0.89)
conditional_effects(brm_ff_intx3, "period")
conditional_effects(brm_ff_intx3, "group.combo")
conditional_effects(brm_ff_intx3, "rval")

loo1<-loo(brm_ff_groom3, pointwise = TRUE)
loo2<-loo(brm_ff_intx3, pointwise = TRUE)
loo_compare(loo1,loo2)

## TRACK RELATEDNESS ACROSS SPLIT VS. STAYED
related.data<-ff.dyads[ff.dyads$period=="Stable" & ff.dyads$group1=="PHG" & ff.dyads$group2=="PHG",]
related.data$related<-0
related.data$related[related.data$rval>0]<-1
related.test<-brm(data=related.data, related ~ outcome + (1|mm(id1,id2)),
                  chains=4, cores=2, iter=4000, family="bernoulli",
                  control = list(max_treedepth=12,adapt_delta=0.99))
summary(related.test, prob=0.89)
conditional_effects(related.test, "outcome")
pp_check(related.test)

## VISUALIZE FF EFFECTS
palette1<-c("#D1B79EFF","#7EBAC2FF") 
palette2<-c("#403369FF", "#5C5992FF", "#AE93BEFF", "#B4DAE5FF")

## MAKE GROUP COMBO PLOT
grid = ff.dyads %>%
  data_grid(rval=mean(ff.dyads$rval), RankDiff=mean(ff.dyads$RankDiff), group.combo, period, outcome, id1, id2, Obs.effort=2400)

means = grid %>%
  add_epred_draws(brm_ff_groom1, ndraws=250, allow_new_levels=T)

check<-means %>%
  group_by(group.combo, .draw) %>%
  dplyr::summarise(Mean.Pred=mean(.epred))

overall<-check %>%
  group_by(group.combo) %>%
  dplyr::summarise(Grand.Mean=mean(Mean.Pred), Min=min(Mean.Pred), Max=max(Mean.Pred), 
                   Pred.25=quantile(Mean.Pred, probs=0.25), Pred.75=quantile(Mean.Pred, probs=0.75),
                   LL.Pred=quantile(Mean.Pred, probs=0.055), UL.Pred=quantile(Mean.Pred, probs=0.945))

plot3a<-ggplot() + 
  geom_violin(aes(y = Mean.Pred/40, x= group.combo, fill=group.combo), data = check, alpha = 0.9, trim=F) +
  scale_fill_manual(values=palette1) +
  geom_point(aes(y = Grand.Mean/40, x= group.combo), data = overall, alpha = 1, size=2, color = "black") +
  geom_errorbar(aes(ymin = Pred.25/40, ymax=Pred.75/40, x= group.combo), data = overall, width=0, linewidth=1.6, color = "black") +
  geom_errorbar(aes(ymin = LL.Pred/40, ymax=UL.Pred/40, x= group.combo), data = overall, width=0, linewidth=0.8, color = "black") +
  ylab("female-female grooming\n(minutes/hour)") + xlab("current group") + theme(legend.position ="none") + ylim(0,0.25)
plot3a

## MAKE OUTCOME PLOT
means = grid %>%
  add_epred_draws(brm_ff_groom2, ndraws=250, allow_new_levels=T)

check<-means %>%
  group_by(outcome, .draw) %>%
  dplyr::summarise(Mean.Pred=mean(.epred))

overall<-check %>%
  group_by(outcome) %>%
  dplyr::summarise(Grand.Mean=mean(Mean.Pred), Min=min(Mean.Pred), Max=max(Mean.Pred), 
                   Pred.25=quantile(Mean.Pred, probs=0.25), Pred.75=quantile(Mean.Pred, probs=0.75),
                   LL.Pred=quantile(Mean.Pred, probs=0.055), UL.Pred=quantile(Mean.Pred, probs=0.945))
?geom_violin()

plot3b<-ggplot() + 
  geom_violin(aes(y = Mean.Pred/40, x= outcome, fill=outcome), data = check, alpha = 0.9, trim=F, width=0.5, scale="width") +
  scale_fill_manual(values=palette2) +
  geom_point(aes(y = Grand.Mean/40, x= outcome), data = overall, alpha = 1, size=2, color = "black") +
  geom_errorbar(aes(ymin = Pred.25/40, ymax=Pred.75/40, x= outcome), data = overall, width=0, linewidth=1.6, color = "black") +
  geom_errorbar(aes(ymin = LL.Pred/40, ymax=UL.Pred/40, x= outcome), data = overall, width=0, linewidth=0.8, color = "black") +
  ylab("") + xlab("eventual outcome") + theme(legend.position ="none") + ylim(0,0.25)
plot3b

## MODEL FM GROOMING DYADS
fm.dyads$outcome<-factor(fm.dyads$outcome, levels=c("stayed together", "split apart", "briefly together", "never together"))
fm.dyads$period<-factor(fm.dyads$period, levels=c("Stable", "Dispersal 1", "Dispersal 2", "Fission"))

## MODEL 1D GROUP
brm_fm_groom1<-brm(data=fm.dyads, bf(Grooming~   + group.combo + period + scale(Male.Rank) + siring + offset(log(Obs.effort)) + (1|mm(id1,id2)),
                                    hu~    group.combo + period + scale(Male.Rank) + siring + log(Obs.effort) + (1|mm(id1,id2))),
                  chains=4, cores=2, iter=4000, family="hurdle_gamma",
                  control = list(max_treedepth=12,adapt_delta=0.99))
summary(brm_fm_groom1, prob=0.89)
conditional_effects(brm_fm_groom1, "group.combo")

## MODEL 1E OUTCOME
brm_fm_groom2<-brm(data=fm.dyads, bf(Grooming~ outcome + period + scale(Male.Rank) + siring + offset(log(Obs.effort)) + (1|mm(id1,id2)),
                                     hu~  outcome + period + scale(Male.Rank) + siring + log(Obs.effort) + (1|mm(id1,id2))),
                   chains=4, cores=2, iter=4000, family="hurdle_gamma",
                   control = list(max_treedepth=12,adapt_delta=0.99))
summary(brm_fm_groom2, prob=0.89)
conditional_effects(brm_fm_groom2, "outcome")

## MODEL 1F OUTCOME
brm_fm_groom3<-brm(data=fm.dyads, bf(Grooming~ group.combo + outcome + period + scale(Male.Rank) + siring + offset(log(Obs.effort)) + (1|mm(id1,id2)),
                                     hu~  group.combo + outcome + period + scale(Male.Rank) + siring + log(Obs.effort) + (1|mm(id1,id2))),
                   chains=4, cores=2, iter=4000, family="hurdle_gamma",
                   control = list(max_treedepth=12,adapt_delta=0.99))
summary(brm_fm_groom3, prob=0.89)

loo1<-loo(brm_fm_groom1)
loo2<-loo(brm_fm_groom2)
loo3<-loo(brm_fm_groom3)
loo_compare(loo1,loo2,loo3)
model_weights(brm_fm_groom1, brm_fm_groom2, brm_fm_groom3)

## VISUALIZE FM EFFECTS
## MAKE GROUP COMBO PLOT
grid = fm.dyads %>%
  data_grid(group.combo, outcome, id1, id2, period, siring="No", Male.Rank=mean(fm.dyads$Male.Rank), Obs.effort=2400)

means = grid %>%
  add_epred_draws(brm_fm_groom1, ndraws=250, allow_new_levels=T)

check<-means %>%
  group_by(group.combo, .draw) %>%
  dplyr::summarise(Mean.Pred=mean(.epred))

overall<-check %>%
  group_by(group.combo) %>%
  dplyr::summarise(Grand.Mean=mean(Mean.Pred), Min=min(Mean.Pred), Max=max(Mean.Pred), 
                   Pred.25=quantile(Mean.Pred, probs=0.25), Pred.75=quantile(Mean.Pred, probs=0.75),
                   LL.Pred=quantile(Mean.Pred, probs=0.055), UL.Pred=quantile(Mean.Pred, probs=0.945))

plot3c<-ggplot() + 
  geom_violin(aes(y = Mean.Pred/40, x= group.combo, fill=group.combo), data = check, alpha = 0.9, trim=F) +
  scale_fill_manual(values=palette1) +
  geom_point(aes(y = Grand.Mean/40, x= group.combo), data = overall, alpha = 1, size=2, color = "black") +
  geom_errorbar(aes(ymin = Pred.25/40, ymax=Pred.75/40, x= group.combo), data = overall, width=0, linewidth=1.6, color = "black") +
  geom_errorbar(aes(ymin = LL.Pred/40, ymax=UL.Pred/40, x= group.combo), data = overall, width=0, linewidth=0.8, color = "black") +
  ylab("female-male grooming\n(minutes/hour)") + xlab("current group") + theme(legend.position ="none") + ylim(0,0.2)
plot3c

## MAKE OUTCOME PLOT
means = grid %>%
  add_epred_draws(brm_fm_groom2, ndraws=250, allow_new_levels=T)

check<-means %>%
  group_by(outcome, .draw) %>%
  dplyr::summarise(Mean.Pred=mean(.epred))

overall<-check %>%
  group_by(outcome) %>%
  dplyr::summarise(Grand.Mean=mean(Mean.Pred), Min=min(Mean.Pred), Max=max(Mean.Pred), 
                   Pred.25=quantile(Mean.Pred, probs=0.25), Pred.75=quantile(Mean.Pred, probs=0.75),
                   LL.Pred=quantile(Mean.Pred, probs=0.055), UL.Pred=quantile(Mean.Pred, probs=0.945))
?scale_y_log10
plot3d<-ggplot() + 
  geom_violin(aes(y = Mean.Pred/40, x= outcome, fill=outcome), data = check, alpha = 0.8, trim=F) +
  scale_fill_manual(values=palette2) +
  geom_point(aes(y = Grand.Mean/40, x= outcome), data = overall, alpha = 1, size=2, color = "black") +
  geom_errorbar(aes(ymin = Pred.25/40, ymax=Pred.75/40, x= outcome), data = overall, width=0, linewidth=1.6, color = "black") +
  geom_errorbar(aes(ymin = LL.Pred/40, ymax=UL.Pred/40, x= outcome), data = overall, width=0, linewidth=0.8, color = "black") +
  ylab("") + xlab("eventual outcome") + theme(legend.position ="none") + ylim(0,0.2)
plot3d

## COMBINE GROOMING FIGURES
layout<-"AABBBB
         CCDDDD"
plot3a+plot3b+
  plot3c+plot3d + plot_layout(design=layout) + plot_annotation(tag_levels = "a")
ggsave(file="Laikipia_Fission_Males_Females.jpg", units="cm", width=28, height=16, dpi=300)

## Q2 CAPTURE NETWORK-WIDE STRUCTURAL CHANGES
## USE MODEL PREDICTIONS TO MAKE POSTERIOR NETWORKS
## AND EXTRACT NETWORK METRICS
newdata <- ff.dyads
newdata$Obs.effort<-2400

newdata2 <- fm.dyads
newdata2$Obs.effort<-2400

pp1<-predict(brm_ff_groom2, newdata = newdata, ndraws=1000, summary=F)
pp2<-predict(brm_fm_groom3, newdata = newdata2, ndraws=1000, summary=F)
newdata<-newdata[colnames(newdata) %in% c("id1","id2", "period", "group1", "group2")]
newdata2<-newdata2[colnames(newdata2) %in% c("id1","id2", "period", "group1", "group2")]

## GET PHG METRICS
density1.PHG<-c()
density2.PHG<-c()
density3.PHG<-c()
density4.PHG<-c()

modularity1.PHG<-c()
modularity2.PHG<-c()
modularity3.PHG<-c()
modularity4.PHG<-c()

cv1.PHG<-c()
cv2.PHG<-c()
cv3.PHG<-c()
cv4.PHG<-c()

## PREP ENK METRICS
density1.ENK<-c()
density2.ENK<-c()
density3.ENK<-c()
density4.ENK<-c()

modularity1.ENK<-c()
modularity2.ENK<-c()
modularity3.ENK<-c()
modularity4.ENK<-c()

cv1.ENK<-c()
cv2.ENK<-c()
cv3.ENK<-c()
cv4.ENK<-c()

## PREP YNT METRICS
density4.YNT<-c()
modularity4.YNT<-c()
cv4.YNT<-c()

i<-3
for (i in 1:1000) {
  newdata$weight<-pp1[i,]
  newdata2$weight<-pp2[i,]
  edge.list<-rbind(newdata, newdata2)
  PHG.list<-edge.list[edge.list$group1=="PHG" & edge.list$group2=="PHG",]
  ENK.list<-edge.list[edge.list$group1=="ENK" & edge.list$group2=="ENK",]
  YNT.list<-edge.list[edge.list$group1=="YNT" & edge.list$group2=="YNT",]
  
  ## STABLE PHG
  stable.list<-PHG.list[PHG.list$period=="Stable",]
  graph1<-graph_from_data_frame(stable.list[,-c(1,4,5)], directed=F)
  isolates<-which(igraph::strength(graph1)==0)
  graph1<-delete_vertices(graph1,isolates)
  density1.PHG[i]<-sum(E(graph1)$weight>0)/length(E(graph1))
  modularity1.PHG[i]<-modularity(cluster_louvain(graph1))
  cv1.PHG[i]<-sd(E(graph1)$weight)/mean(E(graph1)$weight)
  ## DISPERSAL 1 PHG
  dispersal1.list<-PHG.list[PHG.list$period=="Dispersal 1",]
  graph2<-graph_from_data_frame(d = dispersal1.list[-c(1,4,5)], directed = F)
  isolates<-which(igraph::strength(graph2)==0)
  graph2<-delete_vertices(graph2,isolates)
  sort(E(graph2)$weight)
  density2.PHG[i]<-sum(E(graph2)$weight>0)/length(E(graph2))
  modularity2.PHG[i]<-modularity(cluster_louvain(graph2))
  cv2.PHG[i]<-sd(E(graph2)$weight)/mean(E(graph2)$weight)
  ## DISPERSAL 2 PHG
  dispersal2.list<-PHG.list[PHG.list$period=="Dispersal 2",]
  graph3<-graph_from_data_frame(d = dispersal2.list[,-c(1,4,5)], directed = F)
  isolates<-which(igraph::strength(graph3)==0)
  graph3<-delete_vertices(graph3,isolates)
  density3.PHG[i]<-sum(E(graph3)$weight>0)/length(E(graph3))
  modularity3.PHG[i]<-modularity(cluster_louvain(graph3))
  cv3.PHG[i]<-sd(E(graph3)$weight)/mean(E(graph3)$weight)
  ## FISSION PHG
  fission.list<-PHG.list[PHG.list$period=="Fission",]
  graph4<-graph_from_data_frame(d = fission.list[,-c(1,4,5)], directed = F)
  isolates<-which(igraph::strength(graph4)==0)
  graph4<-delete_vertices(graph4,isolates)
  density4.PHG[i]<-sum(E(graph4)$weight>0)/length(E(graph4))
  modularity4.PHG[i]<-modularity(cluster_louvain(graph4))
  cv4.PHG[i]<-sd(E(graph4)$weight)/mean(E(graph4)$weight)
  
  ## STABLE ENK
  stable.list<-ENK.list[ENK.list$period=="Stable",]
  graph1<-graph_from_data_frame(stable.list[,-c(1,4,5)], directed=F)
  isolates<-which(igraph::strength(graph1)==0)
  graph1<-delete_vertices(graph1,isolates)
  density1.ENK[i]<-sum(E(graph1)$weight>0)/length(E(graph1))
  modularity1.ENK[i]<-modularity(cluster_louvain(graph1))
  cv1.ENK[i]<-sd(E(graph1)$weight)/mean(E(graph1)$weight)
  ## DISPERSAL 1 ENK
  dispersal1.list<-ENK.list[ENK.list$period=="Dispersal 1",]
  graph2<-graph_from_data_frame(d = dispersal1.list[,-c(1,4,5)], directed = F)
  isolates<-which(igraph::strength(graph2)==0)
  graph2<-delete_vertices(graph2,isolates)
  density2.ENK[i]<-sum(E(graph2)$weight>0)/length(E(graph2))
  modularity2.ENK[i]<-modularity(cluster_louvain(graph2))
  cv2.ENK[i]<-sd(E(graph2)$weight)/mean(E(graph2)$weight)
  ## DISPERSAL 2 ENK
  dispersal2.list<-ENK.list[ENK.list$period=="Dispersal 2",]
  graph3<-graph_from_data_frame(d = dispersal2.list[,-c(1,4,5)], directed = F)
  isolates<-which(igraph::strength(graph3)==0)
  graph3<-delete_vertices(graph3,isolates)
  density3.ENK[i]<-sum(E(graph3)$weight>0)/length(E(graph3))
  modularity3.ENK[i]<-modularity(cluster_louvain(graph3))
  cv3.ENK[i]<-sd(E(graph3)$weight)/mean(E(graph3)$weight)
  ## FISSION ENK
  fission.list<-ENK.list[ENK.list$period=="Fission",]
  graph4<-graph_from_data_frame(d = fission.list[,-c(1,4,5)], directed = F)
  isolates<-which(igraph::strength(graph4)==0)
  graph4<-delete_vertices(graph4,isolates)
  density4.ENK[i]<-sum(E(graph4)$weight>0)/length(E(graph4))
  modularity4.ENK[i]<-modularity(cluster_louvain(graph4))
  cv4.ENK[i]<-sd(E(graph4)$weight)/mean(E(graph4)$weight)
  
  ## FISSION YNT
  fission.list<-YNT.list[YNT.list$period=="Fission",]
  graph4<-graph_from_data_frame(d = fission.list[,-c(1,4,5)], directed = F)
  isolates<-which(igraph::strength(graph4)==0)
  graph4<-delete_vertices(graph4,isolates)
  density4.YNT[i]<-sum(E(graph4)$weight>0)/length(E(graph4))
  modularity4.YNT[i]<-modularity(cluster_louvain(graph4))
  cv4.YNT[i]<-sd(E(graph4)$weight)/mean(E(graph4)$weight)
}

## METRIC BOUNDS
metric.data<-data.frame(Period=c("Stable","Dispersal 1", "Dispersal 2", "Fission",
                    "Stable","Dispersal 1", "Dispersal 2", "Fission", "Fission"),
           Group=c("PHG","PHG","PHG","PHG","ENK","ENK","ENK","ENK","YNT"),
           LL.Mod=NA, Mid.Mod=NA, UL.Mod=NA,
           LL.Dens=NA, Mid.Dens=NA, UL.Dens=NA,
           LL.CV=NA, Mid.CV=NA, UL.CV=NA)

metric.data$Group<-factor(metric.data$Group,levels=c("PHG","ENK","YNT"))

metric.data[1,3:5]<-quantile(modularity1.PHG, probs = c(0.055,0.5,0.945))
metric.data[2,3:5]<-quantile(modularity2.PHG, probs = c(0.055,0.5,0.945))
metric.data[3,3:5]<-quantile(modularity3.PHG, probs = c(0.055,0.5,0.945))
metric.data[4,3:5]<-quantile(modularity4.PHG, probs = c(0.055,0.5,0.945))
metric.data[5,3:5]<-quantile(modularity1.ENK, probs = c(0.055,0.5,0.945))
metric.data[6,3:5]<-quantile(modularity2.ENK, probs = c(0.055,0.5,0.945))
metric.data[7,3:5]<-quantile(modularity3.ENK, probs = c(0.055,0.5,0.945))
metric.data[8,3:5]<-quantile(modularity4.ENK, probs = c(0.055,0.5,0.945))
metric.data[9,3:5]<-quantile(modularity4.YNT, probs = c(0.055,0.5,0.945))
metric.data[1,6:8]<-quantile(density1.PHG, probs = c(0.055,0.5,0.945))
metric.data[2,6:8]<-quantile(density2.PHG, probs = c(0.055,0.5,0.945))
metric.data[3,6:8]<-quantile(density3.PHG, probs = c(0.055,0.5,0.945))
metric.data[4,6:8]<-quantile(density4.PHG, probs = c(0.055,0.5,0.945))
metric.data[5,6:8]<-quantile(density1.ENK, probs = c(0.055,0.5,0.945))
metric.data[6,6:8]<-quantile(density2.ENK, probs = c(0.055,0.5,0.945))
metric.data[7,6:8]<-quantile(density3.ENK, probs = c(0.055,0.5,0.945))
metric.data[8,6:8]<-quantile(density4.ENK, probs = c(0.055,0.5,0.945))
metric.data[9,6:8]<-quantile(density4.YNT, probs = c(0.055,0.5,0.945))

metric.data[1,9:11]<-quantile(cv1.PHG, probs = c(0.055,0.5,0.945))
metric.data[2,9:11]<-quantile(cv2.PHG, probs = c(0.055,0.5,0.945))
metric.data[3,9:11]<-quantile(cv3.PHG, probs = c(0.055,0.5,0.945))
metric.data[4,9:11]<-quantile(cv4.PHG, probs = c(0.055,0.5,0.945))
metric.data[5,9:11]<-quantile(cv1.ENK, probs = c(0.055,0.5,0.945))
metric.data[6,9:11]<-quantile(cv2.ENK, probs = c(0.055,0.5,0.945))
metric.data[7,9:11]<-quantile(cv3.ENK, probs = c(0.055,0.5,0.945))
metric.data[8,9:11]<-quantile(cv4.ENK, probs = c(0.055,0.5,0.945))
metric.data[9,9:11]<-quantile(cv4.YNT, probs = c(0.055,0.5,0.945))

## MAKE NETWORK MEASURE PLOTS - PLOT 2a,2b
metric.data$Period<-factor(metric.data$Period, levels=c("Stable","Dispersal 1", "Dispersal 2", "Fission"))
plot2a<-ggplot(data=metric.data, aes(x=Period, y=Mid.Dens, group=Group, color=Group)) +
  geom_line(lwd=1.5, position=position_dodge(width=0.5)) +
  geom_pointrange(aes(ymin=LL.Dens, ymax=UL.Dens), cex=1, lwd=1.5, position = position_dodge(width=0.5)) + ylim(0,1) + xlab("") +
  scale_color_viridis(discrete = T) + theme(legend.position = c(0.65,0.9))  + ylab("Density")

plot2b<-ggplot(data=metric.data, aes(x=Period, y=Mid.CV, group=Group, color=Group)) +
  geom_line(lwd=1.5, position=position_dodge(width=0.5)) +
  geom_pointrange(aes(ymin=LL.CV, ymax=UL.CV), cex=1, lwd=1.5, position = position_dodge(width=0.5)) + ylim(0,4) + xlab("") +
  scale_color_viridis(discrete = T) + theme(legend.position = "none")  + ylab("CV")

plot2c<-ggplot(data=metric.data, aes(x=Period, y=Mid.Mod, group=Group, color=Group)) +
  geom_line(lwd=1.5, position=position_dodge(width=0.5)) + geom_hline(aes(yintercept=0.3),lty="dashed") +
  geom_pointrange(aes(ymin=LL.Mod, ymax=UL.Mod), cex=1, lwd=1.5, position = position_dodge(width=0.5)) + ylim(0,0.6) + xlab("") +
  scale_color_viridis(discrete = T) + theme(legend.position = "none")  + ylab("Modularity")

plot2a/plot2c + plot_annotation(tag_levels = "a")

## Q3 AGGRESSION ANALYSES
## MODEL 3a: AGGRESSION BY INDIVIDUAL
dispersers<-unique(agg.by.event.intra$focal[agg.by.event.intra$relative.size=="leaver"])
agg.by.id$dispersers<-"No"
agg.by.id$dispersers[agg.by.id$id %in% dispersers]<-"Yes"

agg.by.id$period<-factor(agg.by.id$period, levels=c("Stable","Dispersal 1", "Dispersal 2","Fission"))
brm_agg_id<-brm(data=agg.by.id, 
                N.Agg ~ period + dispersals + dispersers + groupsize + offset(log(Obs.Effort)) + 
                  (1|id), chains=4, cores=2, iter=4000, family="zero_inflated_poisson",
                control = list(max_treedepth=12,adapt_delta=0.99))
summary(brm_agg_id, prob=0.89)

conditional_effects(brm_agg_id, "dispersals")
conditional_effects(brm_agg_id, "period")
output<-summary(brm_agg_id, prob=0.89, digits=3)
print(output$fixed,digits=3)

## MAKE DISPERSALS PLOT
palette1<-c("#D1B79EFF","#7EBAC2FF") 
palette2<-c("#403369FF", "#5C5992FF", "#AE93BEFF", "#B4DAE5FF")
grid = agg.by.id %>%
  data_grid(dispersals, dispersers, groupsize, period="Stable", id, Obs.Effort=1)

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
  ylab("") + xlab("Group unstable?") + theme(legend.position ="none") + ylim(0,4)
plot1a

## MAKE PERIOD PLOT
grid = agg.by.id %>%
  data_grid(dispersals="no", dispersers="No", groupsize, period, id, Obs.Effort=1)

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

## MODEL S3 DYADIC-LEVEL
agg.by.dyad$outcome<-factor(agg.by.dyad$outcome,  levels=c("stayed together", "split apart", "briefly together", "never together"))
agg.by.dyad$period<-factor(agg.by.dyad$period, levels=c("Stable","Dispersal 1", "Dispersal 2","Fission"))
brm_agg<-brm(data=agg.by.dyad, 
             N.Agg ~  period + outcome + group.combo + scale(rval) + scale(RankDiff) + offset(log(Obs.effort))+ 
               (1|mm(id1,id2)), chains=4, cores=2, iter=4000, family="zero_inflated_poisson",
             control = list(max_treedepth=12,adapt_delta=0.99))
summary(brm_agg, prob=0.89)

output<-summary(brm_agg, prob=0.89, digits=3)
print(output$fixed,digits=3)

## MODEL 3B: EVENT-LEVEL
agg.by.event$outcome<-factor(agg.by.event$outcome,  levels=c("stayed together", "split apart", "briefly together", "never together"))
agg.by.event$period<-factor(agg.by.event$period, levels=c("Stable","Dispersal 1", "Dispersal 2","Fission"))
brm_intensity<-brm(data=agg.by.event, 
                   Intensity ~  period + outcome + rval + RankDiff +
                     (1|mm(id1,id2)), chains=4, cores=2, iter=4000, family="bernoulli",
                   control = list(max_treedepth=12,adapt_delta=0.99))
summary(brm_intensity, prob=0.89)
pp_check(brm_intensity)

 ## MODEL 3C: EVENT-LEVEL - INTERGROUP
agg.by.event.inter$outcome<-factor(agg.by.event.inter$outcome,  levels=c("stayed together", "split apart", "briefly together", "never together"))
agg.by.event.inter$period<-factor(agg.by.event.inter$period, levels=c("Stable","Dispersal 1", "Dispersal 2","Fission"))
brm_winner_inter<-brm(data=agg.by.event.inter,
                win ~  relative.size + period + (1|mm(focal,partner)), chains=4, cores=2, iter=4000, family="bernoulli",
                control = list(max_treedepth=12,adapt_delta=0.99))
summary(brm_winner_inter,prob=0.89)

## Q4
## INTERGROUP GROOMING
intergroup.grooming$Period<-factor(intergroup.grooming$Period, levels=c("Stable","Dispersal 1", "Dispersal 2","Fission"))
intergroup.grooming$outcome<-factor(intergroup.grooming$outcome, levels=c("non-disperser", "disperser"))

brm_inter_groom<-brm(data=intergroup.grooming, 
                     Between.Dur|trials(Groom.Duration) ~ outcome*Period + outcome*sex.combo + (1|focal), chains=4, cores=2, 
                     iter=4000, family="zero_inflated_binomial",
                     control = list(max_treedepth=12,adapt_delta=0.99))
print(summary(brm_inter_groom), digits=2, prob=0.89)

pp_check(brm_inter_groom)
conditional_effects(brm_inter_groom, "Period:outcome")
conditional_effects(brm_inter_groom, "outcome:sex.combo")

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

