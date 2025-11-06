##HISTORICALEXPOSURE STUDY SCRIPT CH2###
##ABNORMALITY DATA###

setwd("~/Desktop/purp dev data")

library(ggplot2)
library(tidyr)
library(dplyr)
library(tidyverse)
library(car)
library(lme4)
library(emmeans)
library(RColorBrewer)
library(wesanderson)

#figure theme to keep things consistent
theme_box <- function(base_size = 11, base_family = '') {
  theme_classic() %+replace% 
    theme(text = element_text(size = 20)) 
}


#Blastula####
HE_abntypeB<-read.csv("2024_WBML/HE_Blastula_abn-type.csv") #categorical data
HE_abncountB<-read.csv("2024_WBML/HE_Blastula_abn-counts.csv") #count data

###count data first####
HE_abncountB<- HE_Blastula_abn_counts %>% 
  group_by(BeakerID, Population, MatePair, Treatment, Urban) %>% 
  mutate (prop.abn=(abn_count/total))

#remove NA sections
HE_abncountB<- HE_abncountB %>%
  select_if(~ !all(is.na(.))) %>% 
  filter_all(all_vars(!is.na(.)))

#Organize Populations with better labels
HE_abncountB$Population<- factor(x=HE_abncountB$Population, labels=c(
  "Cab" ="Urban, Cabrillo Beach",
  "Trp" = "Nonurban, Treasure Park",
  "Twp" = "Nonurban, Twin Pts",
  "Wh" = "Urban, White Pt"
))

#population differences
ggplot(HE_abncountB, aes(x=Population, y=prop.abn, fill=Treatment)) +
  geom_boxplot()+ theme_classic()+
  labs(title = 'Proportion Abnormal Blastula Across Treatments', x = 'Nonylphenol Concentration',
       y = 'Proportion of Abnormal Blastula')+
  scale_fill_manual(values = c("NP" = "#C77DFF", "C" = "#C8E9A0"))

#urban vs not
b_urb<-
  ggplot(HE_abncountB, aes(x=Urban, y=prop.abn, fill=Treatment)) +
  geom_boxplot()+ theme_box()+
  #geom_line(aes(color= MatePair, group=MatePair))+
  labs(title = 'A. Blastula', x = '', y = 'Proportion Abnormalities')+
  scale_fill_manual(values = c("NP" = "#C77DFF", "C" = "#C8E9A0"))

ggsave("blastula_urban.png",b_urb, width=20, height=15, units = "cm") 

#test for abn variation with treatment, blocking by population
bstats<-lmer(prop.abn~Treatment*Urban+ (1|Experiment), data=HE_abncountB) 
anova(bstats) 

#lets test for the amount of variance within vs between populations for survival relative to control
#Levene's test, more flexible if data isn't normally distributed
leveneTest(ratio~Treatment*Urban, data=HE_abncountB) 

#since it isn't we need to do a nonparametric test for abn in this stage
library(ARTool)
HE_abncountB$Treatment <- factor(HE_abncountB$Treatment) #need to all be factors
HE_abncountB$Urban     <- factor(HE_abncountB$Urban)
HE_abncountB$Experiment <- factor(HE_abncountB$Experiment)

bstats.fix <- art(ratio ~ Treatment * Urban + (1|Experiment), data = HE_abncountB)
anova(bstats.fix)

##relative to control calcs####
##split for rel2C
HE_abncountB_split<-HE_abncountB %>% 
  spread("Treatment", "prop.abn") %>% 
  rename(prop_abnC= C, prop_abnNP=NP) 

##separate dataframe for control and dataframe for treatment
#Order by population and MP
HE_abncountB_C<-HE_abncountB_split %>% 
  select(BeakerID,Population,MatePair,Experiment,Urban,prop_abnC)%>% 
  unite("MP_Population", MatePair:Population, remove=FALSE)%>% 
  drop_na()

abncountB_NP<-HE_abncountB_split %>% 
  select(BeakerID,Population,MatePair,Experiment,Urban,prop_abnNP)%>% 
  unite("MP_Population", MatePair:Population, remove=FALSE)%>% 
  drop_na()

#now merge back together so we can look at the difference
HE_abncountB<-left_join(HE_abncountB_C, abncountB_NP, by="MP_Population")

#remove duplicated columns and rename, and add new column for difference between control and treatment
HE_abncountB<-HE_abncountB %>% 
  select(-c(BeakerID.y,MatePair.y, Population.y, Urban.y, Experiment.y)) %>% 
  rename(BeakerID=BeakerID.x, Population=Population.x, MatePair=MatePair.x, Urban=Urban.x, Experiment= Experiment.x) %>% 
  mutate (prop_abn2C=(prop_abnNP-prop_abnC))

#Reorder so urban sites are first and nonurban are second
HE_abncountB$Urban=factor(HE_abncountB$Urban,levels=c("Urban", "Nonurban"))

#rel2CAbn<-
ggplot(HE_abncountB, aes(x=Urban, y=prop_abn2C, fill=Urban)) +
  geom_boxplot()+ theme_box()+
#  facet_wrap(~Experiment)+
  geom_hline(yintercept = 0, color="gray")+
  labs(title = 'Blastula', x = 'Site Classification', 
       y = 'Relative Proportion Abnormalities')+ 
  scale_color_manual(values = c("Urban" = "#5C5649", "Nonurban" = "#99D1EC")) +  scale_fill_manual(values = c("Urban" = "#5C5649", "Nonurban" = "#99D1EC"))


ggsave("rel2CAbn.png",rel2CAbn, width=20, height=11, units = "cm") 

ggplot(HE_abncountB, aes(x=Population, y=prop_abn2C, fill=Population)) +
  geom_boxplot()+ theme_box()+
  #  facet_wrap(~Experiment)+
  geom_hline(yintercept = 0, color="gray")+
  labs(title = 'Blastula', x = 'Site Classification', 
       y = 'Relative Proportion Abnormalities')+ 
  scale_color_manual(values = c("Urban, Cabrillo Beach"="#5C5649", "Urban, White Pt"="#A69C8B", "Nonurban, Treasure Park"="#0F6DD2", "Nonurban, Twin Pts"="#99D1EC")) + 
  scale_fill_manual(values = c("Urban, Cabrillo Beach"="#5C5649", "Urban, White Pt"="#A69C8B", "Nonurban, Treasure Park"="#0F6DD2", "Nonurban, Twin Pts"="#99D1EC"))

leveneTest(prop_abn2C ~ Urban, data = HE_abncountB)

#test for abn variation with treatment, blocking by population
modelA<-glm(prop_abn2C~ Population, data=HE_abncountB) 
summary(modelA) #diff across populations & urban vs nonurban

shapiro.test(HE_abncountB$prop_abn2C) #greater than 0.05 so the data is not normally distributed (which I expected)

#is the proportion survival significantly less than 0 for all populations? Do test for each population separately
# Split data by population
city <- split(HE_abncountB$prop_abn2C, HE_abncountB$Urban)

wilcox_tests <- lapply(city, function(prop_abn2C) wilcox.test(prop_abn2C, mu = 0, alternative = "less"))
print(wilcox_tests) #urban p-value = 0.991; nonurban p-value = 0.3574

###type####
#now let's do type

#Organize Populations with better labels
HE_abntypeB$Population<- factor(x=HE_abntypeB$Population, labels=c(
  "Cab" ="Urban, Cabrillo Beach",
  "Trp" = "Nonurban, Treasure Park",
  "Twp" = "Nonurban, Twin Pts",
  "Wh" = "Urban, White Pt"
))

#add in total count column from other dataset
HE_abntypeB_full<- HE_abntypeB %>% 
  group_by(BeakerID, Population, Treatment, Urban, abnormality_type) %>%
  tally() %>%
  rename(abn_count=n)

#group with data from count information
HE_abntypeB_full<-HE_abntypeB_full %>%
  group_by(BeakerID, Population, Treatment, Urban, abnormality_type) %>% 
  left_join(., HE_abncountB %>% select(BeakerID, total))

#get ratio of abn types
HE_abntypeB_full<-HE_abntypeB_full %>%
  mutate(prop_abn=(abn_count/total)) %>% 
  group_by(BeakerID, Population, Treatment, Urban, abnormality_type)%>%
  summarize(totalprop_abn=sum(prop_abn), total=mean(total))


typesgraph_b<-
  ggplot(HE_abntypeB_full, aes(fill=abnormality_type, y=totalprop_abn, x=Treatment)) +
  geom_bar(position="stack", stat="identity")+ 
  facet_wrap("Urban")+
  theme_box()+
  labs(title = 'A. Blastula', x = 'Treatment', y = 'Proportion of Abnormality Types')+  
  scale_fill_brewer(palette="Spectral")

ggsave("blast_abntypes.png",typesgraph_b, width=20, height=15, units = "cm")

#gastrula ####
HE_abncountG<-read.csv("2024_WBML/HE_Gastrula_abn-counts.csv") #count data
HE_abntypeG<-read.csv("2024_WBML/HE_Gastrula_abn-type.csv") #categorical data

###count data first####
#remove NA sections
HE_abncountG<- HE_abncountG %>%
  select_if(~ !all(is.na(.))) %>% 
  filter_all(all_vars(!is.na(.)))

#Organize Populations with better labels
HE_abncountG$Population<- factor(x=HE_abncountG$Population, labels=c(
  "Cab" ="Urban, Cabrillo Beach",
  "Trp" = "Nonurban, Treasure Park",
  "Twp" = "Nonurban, Twin Pts",
  "Wh" = "Urban, White Pt"
))

#population differences
abnGg<-
  ggplot(HE_abncountG, aes(x=Urban, y=ratio, fill=Treatment)) +
  geom_boxplot()+ theme_box()+
  labs(title = 'B. Gastrula ', x = '',
       y = 'Proportion of Abnormalities')+
    scale_fill_manual(values = c("NP" = "#C77DFF", "C" = "#C8E9A0"))

ggsave("gastrula_urban.png",abnGg, width=20, height=15, units = "cm")   
  
#save plots as pictures
#ggsave("abn.png",urbabnGg, width=32, height=15, units = "cm")

#test for abn variation with treatment, blocking by population
gstats<-lmer(ratio~Treatment*Urban+ (1|Experiment), data=HE_abncountG) 
anova(gstats) 

#Levene's test, more flexible if data isn't normally distributed
leveneTest(ratio~Treatment*Urban, data=HE_abncountG) 

#since it isn't we need to do a nonparametric test for abn in this stage
library(ARTool)
HE_abncountG$Treatment <- factor(HE_abncountG$Treatment) #need to all be factors
HE_abncountG$Urban     <- factor(HE_abncountG$Urban)
HE_abncountG$Experiment <- factor(HE_abncountG$Experiment)

gstats.fix <- art(ratio ~ Treatment * Urban + (1|Experiment), data = HE_abncountG)
anova(gstats.fix)


##separate dataframe for coUrban##separate dataframe for control and dataframe for treatment
HE_abncountG_split<-HE_abncountG %>% 
  spread("Treatment", "ratio") %>% 
  rename(prop_abnC= C, prop_abnNP=NP)

#Order by population and MP
HE_abncountG_C<-HE_abncountG_split %>% 
  select(BeakerID,Population,MatePair,Experiment,Urban,prop_abnC)%>% 
  unite("MP_Population", MatePair:Population, remove=FALSE)%>% 
  drop_na()

abncountG_NP<-HE_abncountG_split %>% 
  select(BeakerID,Population,MatePair,Experiment,Urban,prop_abnNP)%>% 
  unite("MP_Population", MatePair:Population, remove=FALSE)%>% 
  drop_na()

#now merge back together so we can look at the difference
HE_abncountG_full<-left_join(HE_abncountG_C, abncountG_NP, by="MP_Population")

#remove duplicated columns and rename, and add new column for difference between control and treatment
HE_abncountG_full<-HE_abncountG_full %>% 
  select(-c(BeakerID.y,MatePair.y, Population.y, Urban.y, Experiment.y)) %>% 
  rename(BeakerID=BeakerID.x, Population=Population.x, MatePair=MatePair.x, Urban=Urban.x, Experiment= Experiment.x) %>% 
  mutate (prop_abn2C=(prop_abnNP-prop_abnC))

#Reorder so urban sites are first and nonurban are second
HE_abncountG_full$Urban=factor(HE_abncountG_full$Urban,levels=c("Urban", "Nonurban"))

#rel2CAbnG<-
  ggplot(HE_abncountG_full, aes(x=Urban, y=prop_abn2C, fill=Urban)) +
  geom_boxplot()+ theme_box()+
  #  facet_wrap(~Experiment)+
  geom_hline(yintercept = 0, color="gray")+
  labs(title = 'Gastrula', x = 'Site Classification', 
       y = 'Relative Proportion Abnormalities')+
  scale_color_manual(values = c("Urban" = "#5C5649", "Nonurban" = "#99D1EC")) +  scale_fill_manual(values = c("Urban" = "#5C5649", "Nonurban" = "#99D1EC"))

ggplot(HE_abncountG_full, aes(x=Population, y=prop_abn2C, fill=Population)) +
    geom_boxplot()+ theme_box()+
    geom_hline(yintercept = 0, color="gray")+
    labs(title = 'Gastrula', x = 'Site Classification', 
         y = 'Relative Proportion Abnormalities')+ 
    scale_color_manual(values = c("Urban, Cabrillo Beach"="#5C5649", "Urban, White Pt"="#A69C8B", "Nonurban, Treasure Park"="#0F6DD2", "Nonurban, Twin Pts"="#99D1EC")) + 
    scale_fill_manual(values = c("Urban, Cabrillo Beach"="#5C5649", "Urban, White Pt"="#A69C8B", "Nonurban, Treasure Park"="#0F6DD2", "Nonurban, Twin Pts"="#99D1EC"))
  

ggsave("rel2CAbnG.png",rel2CAbnG, width=20, height=11, units = "cm") 
leveneTest(prop_abn2C ~ Urban, data = HE_abncountG_full)

#test for abn variation with treatment, blocking by population
modelG<-glm(prop_abn2C~ Population, data=HE_abncountG_full) 
summary(modelG) #diff across populations & urban vs nonurban

###type####
#now let's do type

#Organize Populations with better labels
HE_abntypeG$Population<- factor(x=HE_abntypeG$Population, labels=c(
  "Cab" ="Urban, Cabrillo Beach",
  "Trp" = "Nonurban, Treasure Park",
  "Twp" = "Nonurban, Twin Pts",
  "Wh" = "Urban, White Pt"
))

#add in total count column from other dataset
HE_abntypeG_full<- HE_abntypeG %>% 
  group_by(BeakerID, Population, Treatment, Urban, abnormality_type) %>%
  tally() %>%
  rename(abn_count=n)

#group with data from count information
HE_abntypeG_full<-HE_abntypeG_full %>%
  group_by(BeakerID, Population, Treatment, Urban, abnormality_type) %>% 
  left_join(., HE_abncountB %>% select(BeakerID, total))

#get ratio of abn types
HE_abntypeG_full<-HE_abntypeG_full %>%
  mutate(prop_abn=(abn_count/total)) %>% 
  group_by(BeakerID, Population, Treatment, Urban, abnormality_type)%>%
  summarize(totalprop_abn=sum(prop_abn), total=mean(total))

typesgraph_g<-
  ggplot(HE_abntypeG_full, aes(fill=abnormality_type, y=totalprop_abn, x=Treatment)) +
  geom_bar(position="stack", stat="identity")+ 
  facet_wrap("Urban")+
  theme_box()+
  labs(title = 'B. Gastrula', x = 'Treatment', y = 'Proportion of Abnormality Types')+  
  scale_fill_brewer(palette="Spectral")

ggsave("gast_abntypes.png",typesgraph_g, width=20, height=15, units = "cm")


ggplot(HE_abntypeG_full, aes(fill=abnormality_type, y=totalprop_abn, x=Treatment)) +
  geom_bar(position="stack", stat="identity")+ 
  facet_wrap("Population")+
  theme_classic()+
  labs(title = 'Gastrula Abnormalities Observed', x = 'Treatment', y = 'Abnormalities Observed')+  
  scale_fill_brewer(palette="Set3")

#pluteus ####
HE_abncountP<-read.csv("2024_WBML/HE_Pluteus7dpf_abn-counts.csv") #count data
HE_abntypeP<-read.csv("2024_WBML/HE_Pluteus7dpf_abn-type.csv") #categorical data
###count data first####
#remove NA sections
HE_abncountP<- HE_abncountP %>%
  select_if(~ !all(is.na(.))) %>% 
  filter_all(all_vars(!is.na(.)))

#Organize Populations with better labels
HE_abncountP$Population<- factor(x=HE_abncountP$Population, labels=c(
  "Cab" ="Urban, Cabrillo Beach",
  "Trp" = "Nonurban, Treasure Park",
  "Twp" = "Nonurban, Twin Pts",
  "Wh" = "Urban, White Pt"
))

#population differences
p_urb<-
ggplot(HE_abncountP, aes(x=Urban, y=ratio, fill=Treatment)) +
  geom_boxplot()+ theme_box()+
 # facet_wrap(~Population)+
#  geom_line(aes(color= MatePair, group=MatePair))+
  labs(title = 'C. Pluteus', x = '',
       y = 'Proportion of Abnormalities')+
  scale_fill_manual(values = c("NP" = "#C77DFF", "C" = "#C8E9A0"))

ggsave("pluteus_urban.png",p_urb, width=15, height=15, units = "cm")

#test for abn variation with treatment, blocking by population
pstats<-lmer(ratio~Treatment*Urban+ (1|Experiment), data=HE_abncountP) 
anova(pstats) 
#lets test for the amount of variance within vs between populations for survival relative to control
#Levene's test, more flexible if data isn't normally distributed
leveneTest(ratio~Treatment*Urban, data=HE_abncountP) #no difference in homogeneity


##separate dataframe for control and dataframe for treatment
HE_abncountP_split<-HE_abncountP %>% 
  spread("Treatment", "ratio") %>% 
  rename(prop_abnC= C, prop_abnNP=NP)

#Order by population and MP
HE_abncountP_C<-HE_abncountP_split %>% 
  select(BeakerID,Population,MatePair,Experiment,Urban,prop_abnC)%>% 
  unite("MP_Population", MatePair:Population, remove=FALSE)%>% 
  drop_na()

HE_abncountG_NP<-HE_abncountP_split %>% 
  select(BeakerID,Population,MatePair,Experiment,Urban,prop_abnNP)%>% 
  unite("MP_Population", MatePair:Population, remove=FALSE)%>% 
  drop_na()

#now merge back together so we can look at the difference
HE_abncountP_full<-left_join(HE_abncountP_C, abncountG_NP, by="MP_Population")

#remove duplicated columns and rename, and add new column for difference between control and treatment
HE_abncountP_full<-HE_abncountP_full %>% 
  select(-c(BeakerID.y,MatePair.y, Population.y, Urban.y, Experiment.y)) %>% 
  rename(BeakerID=BeakerID.x, Population=Population.x, MatePair=MatePair.x, Urban=Urban.x, Experiment= Experiment.x) %>% 
  mutate (prop_abn2C=(prop_abnNP-prop_abnC))

#Reorder so urban sites are first and nonurban are second
HE_abncountP_full$Urban=factor(HE_abncountP_full$Urban,levels=c("Urban", "Nonurban"))

#rel2CAbnG<-
  ggplot(HE_abncountP_full, aes(x=Urban, y=prop_abn2C, fill=Urban)) +
  geom_boxplot()+ theme_box()+
  # facet_wrap(~Experiment)+
  geom_hline(yintercept = 0, color="gray")+
  labs(title = 'Pluteus', x = 'Site Classification', 
       y = 'Relative Proportion Abnormalities')+
  scale_color_manual(values = c("Urban" = "#5C5649", "Nonurban" = "#99D1EC")) +  scale_fill_manual(values = c("Urban" = "#5C5649", "Nonurban" = "#99D1EC"))

  
ggplot(HE_abncountP_full, aes(x=Population, y=prop_abn2C, fill=Population)) +
    geom_boxplot()+ theme_box()+
    geom_hline(yintercept = 0, color="gray")+
    labs(title = 'Pluteus', x = 'Site Classification', 
         y = 'Relative Proportion Abnormalities')+ 
    scale_color_manual(values = c("Urban, Cabrillo Beach"="#5C5649", "Urban, White Pt"="#A69C8B", "Nonurban, Treasure Park"="#0F6DD2", "Nonurban, Twin Pts"="#99D1EC")) + 
    scale_fill_manual(values = c("Urban, Cabrillo Beach"="#5C5649", "Urban, White Pt"="#A69C8B", "Nonurban, Treasure Park"="#0F6DD2", "Nonurban, Twin Pts"="#99D1EC"))
  
  
  
  
ggplot(HE_abncountP_full, aes(x=Population, y=prop_abn2C, fill=Population)) +
    geom_boxplot()+ theme_box()+
    # facet_wrap(~Experiment)+
    geom_hline(yintercept = 0, color="gray")+
    labs(title = '', x = 'Site Classification', 
         y = "Relative Proportion Abnormalities")

leveneTest(prop_abn2C ~ Urban, data = HE_abncountP_full)

#test for abn variation with treatment, blocking by population
modelz<-glm(prop_abn2C~ Population, data=HE_abncountP_full) 
summary(modelz) #diff across populations & urban vs nonurban

###type####
#now let's do type

#Organize Populations with better labels
HE_abntypeP$Population<- factor(x=HE_abntypeP$Population, labels=c(
  "Cab" ="Urban, Cabrillo Beach",
  "Trp" = "Nonurban, Treasure Park",
  "Twp" = "Nonurban, Twin Pts",
  "Wh" = "Urban, White Pt"
))

#add in total count column from other dataset
HE_abntypeP_full<- HE_abntypeP %>% 
  group_by(BeakerID, Population, Treatment, Experiment, Urban, abnormality_type) %>%
  tally() %>%
  rename(abn_count=n)

#group with data from count information
HE_abntypeP_full<-HE_abntypeP_full %>%
  group_by(BeakerID, Population, Treatment, Urban, abnormality_type) %>% 
  left_join(., HE_abncountP %>% select(BeakerID, total))

#get ratio of abn types
HE_abntypeP_full<-HE_abntypeP_full %>%
  mutate(prop_abn=(abn_count/total)) %>% 
  group_by(BeakerID, Population, Treatment, Urban, Experiment, abnormality_type)%>%
  summarize(totalprop_abn=sum(prop_abn), total=mean(total))


typesgraph<-
ggplot(HE_abntypeP_full, aes(fill=abnormality_type, y=totalprop_abn, x=Treatment)) +
  geom_bar(position="stack", stat="identity")+ 
  facet_wrap("Urban")+
  theme_box()+
  labs(title = 'D. Pluteus', x = 'Treatment', y = 'Proportion of Abnormality Types')+  
  scale_fill_brewer(palette="Blues")

ggsave("plut_abntypes.png",typesgraph, width=20, height=15, units = "cm")

###remove small plut abn####
HE_abntypeP_nosmall<-HE_abntypeP_full %>% 
  filter(!(abnormality_type %in% "small_plut"))


###now I want to replot general abn ratios without small pluts
#add together abn ratios
HE_abntypeP_ratios<-HE_abntypeP_nosmall %>%
  group_by(BeakerID, Population, Experiment, Treatment, Urban, total)%>%
  summarize(  totalprop_abn = sum(totalprop_abn, na.rm = TRUE))

nosmallsplut<-
ggplot(HE_abntypeP_ratios, aes(x=Urban, y=totalprop_abn, fill=Treatment)) +
  geom_boxplot()+ theme_box()+
  labs(title = 'E. Pluteus (subset)', x = '',
       y = 'Proportion of Abnormalities')+
  scale_fill_manual(values = c("NP" = "#C77DFF", "C" = "#C8E9A0"))

ggsave("pluteus_subset_urban.png",nosmallsplut, width=15, height=15, units = "cm")

#stats
#test for abn variation with treatment, blocking by population
pstatsW<-lmer(totalprop_abn~Treatment+Urban+ (1|Experiment), data=HE_abntypeP_ratios) 
anova(pstatsW) 
