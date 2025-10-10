##HISTORICALEXPOSURE STUDY SCRIPT CH2###
##SURVIVAL AND DEV COUNTS###

setwd("~/Desktop/purp dev data")

library(ggplot2)
library(tidyr)
library(dplyr)
library(tidyverse)
library(car)
library(lme4)
library(emmeans)
library(RColorBrewer)

#figure theme to keep things consistent
theme_box <- function(base_size = 11, base_family = '') {
  theme_classic() %+replace% 
    theme(text = element_text(size = 20)) 
}
#FIRST merge datasets across all 4 experiments in excel!! :) 
#Don't forget three tabs of survival data for the three stages

#Fertilization/Development## 

#EarlyDevCheck#----
fullfert<-read.csv("2024_WBML/HistoricalEx_Devcheck.csv")
fullfert <- fullfert %>% 
  mutate(total= rowSums(across(7:8), na.rm=T)) %>% 
  mutate (prop.fert=(Developed/total)) 

#Organize Populations with better labels
fullfert$Population<- factor(x=fullfert$Population, labels=c(
  "Cab" ="Urban, Cabrillo Beach",
  "Trp" = "Nonurban, Treasure Park",
  "Twp" = "Nonurban, Twin Pts",
  "Wh" = "Urban, White Pt"
))

fullfert.summary<- fullfert %>% 
  group_by(MatePair,Beaker, Treatment, Experiment, Population) %>% 
  summarize(mean.prop=mean(prop.fert), se.prop=(sd(prop.fert)/sqrt(length(prop.fert))))
#View (fullfert.summary)

#Reorder so urban sites are first and nonurban are second
fullfert.summary$Population=factor(fullfert.summary$Population,levels=c("Urban, Cabrillo Beach", "Urban, White Pt", "Nonurban, Treasure Park", "Nonurban, Twin Pts"))

#HEdev<-
ggplot(fullfert.summary, aes(x=Population, y=mean.prop, fill=Treatment)) +
  geom_boxplot()+ theme_classic()+
  #geom_line(aes(color= MatePair, group=MatePair))+
  labs(title = 'Development Success Across Populations', x = 'Population', y = 'Proportion of Developed Embryos')+
  scale_fill_manual(values = c("NP" = "#A13D63", "C" = "#C8E9A0"))


#lines per Mate Pair
#HEMPdev<-
ggplot(fullfert.summary, aes(x=Treatment, y=mean.prop, color=MatePair)) +
  geom_point(size=2)+ ylim(0,1)+
  geom_line(aes(color= MatePair, group=MatePair))+
  facet_wrap(~Population)+theme_classic()+ 
  labs(title = 'Development Success Across Populations', x = 'Population', y = 'Proportion of Developed Embryos')+
  scale_color_viridis_d(option="plasma")

#Stats: are these groups statically different?
#ANOVA https://statsandr.com/blog/anova-in-r/

#Stats: are these groups statically different?
#GLM
TP=glm(mean.prop~Treatment*Population, data=fullfert.summary)
summary(TP) # Diff in development success in white point (p=0.00147) 
#and interaction of treatment and white point (p=0.01735)

TPE=glm(mean.prop~Treatment*Population+Experiment, data=fullfert.summary) 
summary(TPE) #difference in dev success in white point (p=0.00199) 
#and interaction of treatment and white point (p=0.01791) 
#with slight effect of treatment (p=0.08149)

AIC(TP, TPE) #models fit the data similarly...so choose least complicated aka without experiment

#subset without low fert values
fullfert.summary.sub<- fullfert.summary %>% 
  group_by(MatePair,Beaker, Treatment, Experiment, Population, se.prop) %>% 
  filter(mean.prop >= 0.70) %>% 
  select(mean.prop)

#mate pair variation across pops in development success
ggplot(fullfert.summary.sub, aes(x=Treatment, y=mean.prop, color=MatePair)) +
  geom_point(size=2)+ ylim(0,1)+
  geom_line(aes(color= MatePair, group=MatePair))+
  facet_wrap(~Population)+theme_classic()+ 
  labs(title = 'Development Success Across Populations (subsetted)', x = 'Population', y = 'Proportion of Developed Embryos')+
  scale_color_viridis_d(option="plasma")

HE_devsub<-
  ggplot(fullfert.summary.sub, aes(x=Population, y=mean.prop, fill=Treatment)) +
  geom_boxplot()+ theme_classic()+
  #geom_line(aes(color= MatePair, group=MatePair))+
  labs(title = 'Development Success Across Populations(subsetted)', x = 'Population', y = 'Proportion of Developed Embryos')+
    scale_fill_manual(values = c("NP" = "#C77DFF", "C" = "#C8E9A0"))

ggsave("2024_WBML/figsHE/HEdev_subset.png",HE_devsub, width=18, height=10, units = "cm")

TPsub=glm(mean.prop~Treatment*Population+Experiment, data=fullfert.summary.sub)
summary(TPsub)

##what about urban vs nonurban?
#add urban column
urbfullfert.summary.sub <- fullfert.summary.sub %>% 
  mutate(urban=(Population))

#rename so Cab &Wh= urban
urbfullfert.summary.sub$urban<- factor(x=urbfullfert.summary.sub$urban, labels=c(
  "Urban, Cabrillo Beach" ="urban",
  "Urban, White Pt" = "urban",
  "Nonurban, Treasure Park" = "nonurban",
  "Nonurban, Twin Pts" = "nonurban"
))

ggplot(urbfullfert.summary.sub, aes(x=urban, y=mean.prop, fill=Treatment)) +
  geom_boxplot()+ theme_classic()+
  #geom_line(aes(color= MatePair, group=MatePair))+
  labs(title = 'Development Success Across Populations(subsetted)', x = 'Population', y = 'Proportion of Developed Embryos')+
  scale_fill_brewer(palette = "RdGy")

sub1=glm(mean.prop~Treatment*urban+Experiment, data=urbfullfert.summary.sub)
summary(sub1)

#Levene's test, more flexible if data isn't normally distributed
leveneTest(mean.prop ~urban, data = urbfullfert.summary.sub) #not sig more var in urban vs nonurban

##subset even further by removing experiment 1 (chaos experiment)
fullfert.summary.sub2<- fullfert.summary.sub %>% 
  group_by(MatePair,Beaker, Treatment, Experiment, Population, se.prop) %>% 
  filter(Experiment !="1") %>% 
  select(mean.prop)
View (fullfert.summary.sub2)

#HE_devsub2<-
ggplot(fullfert.summary.sub2, aes(x=Treatment, y=mean.prop, fill=Population)) +
  #facet_wrap(~Experiment)+
  geom_boxplot()+ theme_classic()+ 
  #geom_jitter(alpha=0.9, aes(color=Treatment))+
  labs(title = 'Development Success Across Populations (subsetted)', x = 'Population', y = 'Proportion of Developed Embryos')+
  scale_fill_brewer(palette = "PRGn")

ggsave("HEdev_subset2.png",HE_devsub2, width=20, height=17, units = "cm")

subsub=glm(mean.prop~Treatment*Population+Experiment, data=fullfert.summary.sub2)
summary(subsub)

#Survival ----
#compare survival shifts across stages

#first looking at overall shift: blastula to 7dpf Pluteus

#read in blastula file
HEblast<-read.csv("2024_WBML/HE_Survival-Blast.csv")

###remove experiment 1 samples from here
#HEblast<- HEblast %>% 
#  group_by(MatePair, Treatment, Experiment, Population) %>% 
#  filter(Experiment !="1")

#combine experiment & beaker columns to have a unique identifier
HEblast<-HEblast%>%
  unite("Beaker_Exp", Beaker:Experiment, remove=FALSE) %>% 
  mutate(
  environment = case_when(
    Population %in% c("Wh", "Cab") ~ "Urban",
    Population %in% c("Twp", "Trp") ~ "Nonurban"))


#get mean count data per mL and make summary datasheet
HEblast.summary<- HEblast %>% 
  group_by(Stage, Beaker_Exp, Treatment, Population, MatePair, environment) %>% 
  summarize(mean.per.mL =mean((Count)*2), se.prop=(sd(Count)/sqrt(length(Count)))) %>% 
  rename(mean.blast.per.mL= mean.per.mL, blast.se.prop=se.prop)

#read in pluteus file
HEplut<-read.csv("2024_WBML/HE_survival-plut.csv")

#remove experiment 1 samples from here
#HEplut<- HEplut %>% 
#  group_by(MatePair, Treatment, Experiment, Population) %>% 
#  filter(Experiment !="1")

#combine experiment & beaker columns to have a unique identifier
HEplut<-HEplut%>%
  unite("Beaker_Exp", Beaker:Experiment, remove=FALSE) %>% 
  mutate(
    environment = case_when(
      Population %in% c("Wh", "Cab") ~ "Urban",
      Population %in% c("Twp", "Trp") ~ "Nonurban"))

#get mean count data per mL and make summary datasheet
HEplut.summary<- HEplut %>% 
  group_by(Stage, Beaker_Exp, Treatment, Population, MatePair, Experiment, environment) %>% 
  summarize(mean.per.mL =mean((Count)*2), se.prop=(sd(Count)/sqrt(length(Count))))%>% 
  rename(mean.plut.per.mL= mean.per.mL, plut.se.prop=se.prop)

#merge two datasets together using unique identifier, Beaker_Exp
HE_BP<-merge(HEblast.summary,HEplut.summary, by="Beaker_Exp")
#View(HE_BP) #looks like it merged successfully since it should have 96 entries

#rename Treatment column 
HE_BP<-HE_BP %>% dplyr::select (1:16) %>% 
  rename(Treatment=Treatment.x, Population=Population.x, MatePair=MatePair.x)

#remove unwanted columns, like duplicate columns
HE_BP<-HE_BP %>% select(-c(Treatment.y, MatePair.y, Population.y, environment.y))

#add new column for survival differences between blast and plut
HE_BP <- HE_BP %>% 
  mutate (prop_survival=(mean.plut.per.mL/mean.blast.per.mL)) %>% 
  mutate (prop_se=(plut.se.prop/blast.se.prop))


#Organize Populations with better labels
HE_BP$Population<- factor(x=HE_BP$Population, labels=c(
  "Cab" ="Urban, Cabrillo Beach",
  "Trp" = "Nonurban, Treasure Park",
  "Twp" = "Nonurban, Twin Pts",
  "Wh" = "Urban, White Pt"
))

#Reorder so urban sites are first and nonurban are second
HE_BP$Population=factor(HE_BP$Population,levels=c("Urban, Cabrillo Beach", "Urban, White Pt", "Nonurban, Treasure Park", "Nonurban, Twin Pts"))

##NOW we are ready to plot after organizing the data!

HEBP_exp<-
  ggplot(HE_BP, aes(x=environment.x, y=prop_survival, fill=Treatment)) +
  geom_boxplot()+ theme_box()+
 # facet_wrap(~Population)+
 # geom_line(aes(color= MatePair, group=MatePair))+
  labs(title = 'Survival', x = 'Population', 
       y = 'Proportion Survival')+
    scale_fill_manual(values = c("NP" = "#C77DFF", "C" = "#C8E9A0"))
  
ggsave("HE_BP_propsurvival.png",HEBP_exp, width=30, height=15, units = "cm") 

#HEBP_urbsurvival<-
ggplot(HE_BP, aes(x=environment.x, y=prop_survival, fill=Treatment)) +
  geom_boxplot()+ theme_classic()+
  #geom_line(aes(color= MatePair, group=MatePair))+
  labs(title = 'Survival Between Blastula and 7dpf Pluteus Across Populations', x = 'Population', 
       y = 'Proportion Survival')+
  scale_fill_brewer(palette = "PRGn")
ggsave("HE_BP_urbanpropsurvival.png",HEBP_urbsurvival, width=30, height=15, units = "cm") 

#test for abn variation with treatment, blocking by population
survstats<-lmer(prop_survival~Treatment*Population+ (1|Experiment), data=HE_BP) 
anova(survstats) 

#lets test for the amount of variance within vs between populations for survival relative to control
#Levene's test, more flexible if data isn't normally distributed
leveneTest(prop_survival ~ Treatment*environment.x, data = HE_BP) #no difference in homogeneity

#what about survival per beaker relative to control?
# get the comparison of survival from the control to treatments for each dose in excel
#in datasheet= look at difference in treatment - control for each MP for each pop == prop_survival_2control

HESurvival2Csplit<-HE_BP %>% 
  spread("Treatment", "prop_survival") %>% 
  rename(prop_survivalC= C, prop_survivalNP=NP) 
  

##separate dataframe for control and dataframe for treatment
#Order by population and MP
HESurvivalC<-HESurvival2Csplit %>% 
  select(Beaker_Exp,Population,MatePair,prop_survivalC)%>% 
  unite("MP_Population", MatePair:Population, remove=FALSE)%>% 
  drop_na()

HESurvivalNP<-HESurvival2Csplit %>% 
  select(Beaker_Exp,Population,MatePair,prop_survivalNP) %>% 
  unite("MP_Population", MatePair:Population, remove=FALSE)%>% 
  drop_na()

#now merge back together so we can look at the difference
HESurvival2C<-left_join(HESurvivalC, HESurvivalNP, by="MP_Population")
  
#remove duplicated columns and rename, and add new column for difference between control and treatment
HESurvival2C<-HESurvival2C %>% 
  select(-c(Beaker_Exp.y,MatePair.y, Population.y)) %>% 
  rename(Beaker_Exp=Beaker_Exp.x, Population=Population.x, MatePair=MatePair.x) %>% 
  mutate (prop_survival_2C=(prop_survivalNP-prop_survivalC))

#relsurvival<-
ggplot(HESurvival2C, aes(x=Population, y=prop_survival_2C, fill=Population)) +
  geom_point(aes(color=Population))+
  geom_boxplot()+ theme_classic()+
  labs(title = 'Relative Survival to Control Across Populations', x = 'Population', 
       y = 'Relative Proportion Survival to Control Across Families')+ scale_fill_brewer(palette = "PRGn") +scale_color_brewer(palette = "PRGn")

model<-glm(prop_survival_2C~Population, data=HESurvival2C) 
anova(model)

#test for abn variation with treatment, blocking by population
survstats<-lmer(prop_survival~Treatment*Population+ (1|Experiment), data=HE_BP) 
anova(survstats) 

#lets test for the amount of variance within vs between populations for survival relative to control
#Levene's test, more flexible if data isn't normally distributed
leveneTest(prop_survival_2C ~ Population, data = HESurvival2C) #no difference in homogeneity

shapiro.test(HESurvival2C$prop_survival_2C) #greater than 0.05 so the data is not normally distributed (which I expected)

#is the proportion survival significantly less than 0 for all populations? Do test for each population separately
# Split data by population
populations <- split(HESurvival2C$prop_survival_2C, HESurvival2C$Population)
print(populations)

wilcox_tests <- lapply(populations, function(prop_survival_2C) wilcox.test(prop_survival_2C, mu = 0, alternative = "less"))
print(wilcox_tests) #Cab p=0.533, Wh p=0.186, Trp p=0.156, Twp p=0.282
#If the p-value is less than your chosen significance level (commonly 0.05), 
#you can reject the null hypothesis and conclude that the values in that population are significantly lower than 0.

##what about urban vs nonurban?
#add urban column
HE_Survival2Curb <- HESurvival2C %>% 
  mutate(urban=(Population))

#rename so Cab &Wh= urban
HE_Survival2Curb$urban<- factor(x=HE_Survival2Curb$urban, labels=c(
  "Urban, Cabrillo Beach" ="urban",
  "Nonurban, Treasure Park" = "nonurban",
  "Nonurban, Twin Pts" = "nonurban",
  "Urban, White Pt" = "urban"
))

#HE_BP_urb<-
ggplot(HE_Survival2Curb, aes(x=urban, y=prop_survival_2C, fill=urban)) +
  geom_boxplot()+ theme_box()+
  geom_hline(yintercept = 0, color="grey")+ #0=#everyone survived
  labs(title = 'Relative Survival to Control Across Urban vs Nonurban Sites', x = 'Site Classification', 
       y = 'Proportion Survival')+
  scale_color_manual(values = c("urban" = "#5C5649", "nonurban" = "#99D1EC")) +  scale_fill_manual(values = c("urban" = "#5C5649", "nonurban" = "#99D1EC"))

ggsave("HE_BPurb.png",HE_BP_urb, width=20, height=10, units = "cm") 

ggplot(HE_Survival2Curb, aes(x=Population, y=prop_survival_2C, fill=Population)) +
  geom_boxplot()+ theme_box()+
  geom_hline(yintercept = 0, color="grey")+ #0=#everyone survived
  labs(title = 'Relative Survival to Control', x = 'Population', 
       y = 'Proportion Survival')+
  scale_color_manual(values = c("Urban, Cabrillo Beach"="#5C5649", "Urban, White Pt"="#A69C8B", "Nonurban, Treasure Park"="#0F6DD2", "Nonurban, Twin Pts"="#99D1EC")) + 
  scale_fill_manual(values = c("Urban, Cabrillo Beach"="#5C5649", "Urban, White Pt"="#A69C8B", "Nonurban, Treasure Park"="#0F6DD2", "Nonurban, Twin Pts"="#99D1EC"))

urb=glm(prop_survival_2C~urban,data=HE_Survival2Curb)
summary(urb) #no effect

#lets test for the amount of variance within vs between urban vs nonurban for survival relative to control
#Levene's test, more flexible if data isn't normally distributed
leveneTest(prop_survival_2C ~ urban, data = HE_Survival2Curb)

shapiro.test(HE_Survival2Curb$prop_survival_2C) #greater than 0.05 so the data is not normally distributed (which I expected)

#is the proportion survival significantly less than 0 for all populations? Do test for each population separately
# Split data by population
city <- split(HE_Survival2Curb$prop_survival_2C, HE_Survival2Curb$urban)
print(city)

wilcox_tests <- lapply(city, function(prop_survival_2C) wilcox.test(prop_survival_2C, mu = 0, alternative = "less"))
print(wilcox_tests) #urban p=0.2496, nonurban p=0.1078


#subset experiment1 out and look at urban vs nonurban####

subHE_Survival2Curb <- subsurvival2C %>% 
  mutate(urban=(Population))

#rename so Cab &Wh= urban
subHE_Survival2Curb$urban<- factor(x=subHE_Survival2Curb$urban, labels=c(
  "Urban, Cabrillo Beach" ="urban",
  "Urban, White Pt" = "urban",
  "Nonurban, Treasure Park" = "nonurban",
  "Nonurban, Twin Pts" = "nonurban"
))

#HE_BP_urb<-
ggplot(subHE_Survival2Curb, aes(x=urban, y=prop_survival_2C, fill=urban)) +
  geom_boxplot()+ theme_classic()+
  #facet_wrap(~Experiment)+
  #geom_point(aes(color=Population))+
  #geom_line(aes(color=Population, group=Population))+
  labs(title = 'Relative Survival to Control Across Urban vs Nonurban Sites (subsetted)', x = 'Site Classification', 
       y = 'Proportion Survival')+
  scale_fill_brewer(palette = "PRGn")

ggsave("HE_BPurb.png",HE_BP_urb, width=20, height=10, units = "cm") 

urb=aov(prop_survival_2C~urban,data=subHE_Survival2Curb)
summary(urb) #no effect

#lets test for the amount of variance within vs between urban vs nonurban for survival relative to control
#Levene's test, more flexible if data isn't normally distributed
leveneTest(prop_survival_2C ~ urban, data = subHE_Survival2Curb)
#####
