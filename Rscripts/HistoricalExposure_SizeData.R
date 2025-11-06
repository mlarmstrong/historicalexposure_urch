##HISTORICALEXPOSURE STUDY SCRIPT CH2###
##SIZE DATA###

setwd("~/Desktop/purp development")

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

###BLASTULA####
#Import the blastula measurements, select and rename useful columns, and
#add a column for experiments 1 and 2
Bsize<-read.csv("HE_Blastula_measurements.csv") 
#view(Bsize)
Bsize<- Bsize %>% 
  group_by(BeakerID, Population, MatePair, Treatment, Experiment, Urban) %>% 
  rename (size=Horizontal.mm.)

Bsize<- Bsize %>% 
  mutate(
    Urban = case_when(
      Population %in% c("Wh", "Cab") ~ "urban",     # replace with your actual "urban" populations
      Population %in% c("Trp", "Twp") ~ "nonurban"))


Bsize<-Bsize %>% select_if(~ !all(is.na(.))) #remove NA columns
Bsize<-Bsize %>% na.omit() ##remove NA rows

#Organize Populations with better labels
Bsize$Population<- factor(x=Bsize$Population, labels=c(
  "Cab" ="Urban, Cabrillo Beach",
  "Trp" = "Nonurban, Treasure Park",
  "Twp" = "Nonurban, Twin Pts",
  "Wh" = "Urban, White Pt"
))

#Reorder so urban sites are first and nonurban are second
Bsize$Population=factor(Bsize$Population,levels=c("Urban, Cabrillo Beach", "Urban, White Pt", "Nonurban, Treasure Park", "Nonurban, Twin Pts"))


#Box plots
bigbois<-
ggplot(Bsize, aes(x=Population, y=size, fill=Treatment)) +
    geom_boxplot()+ theme_classic()+
    labs(title = 'A. Blastula', x = 'Population', y = 'Blastula Size (nm)')+
    scale_fill_brewer(palette = "PRGn")
  
#ggsave("blast_size.png",bigbois, width=32, height=15, units = "cm") 

#subset out giant Cabrillo Blastula
Bsize.sub<- Bsize %>% 
  group_by(BeakerID, Population, MatePair, Treatment, size, Experiment, Urban) %>% 
  filter(size <= 0.13) %>% 
  select(size)
View(Bsize.sub)

treat.colors<-(c("NP" = "#C77DFF", "C" = "#C8E9A0"))
subsetbois<-
ggplot(Bsize.sub, aes(x=Urban, y=size, fill=Treatment)) +
  geom_boxplot()+ theme_box()+
  #facet_wrap("Experiment")+
  labs(title = 'A. Blastula', x = 'subset', y = 'Blastula Size (mm)')+
  scale_fill_manual(values=treat.colors)

ggsave("subblast_size.png",subsetbois, width=20, height=15, units = "cm")

Bsize_means <- Bsize %>%
  group_by(Urban) %>%
  summarise(mean_size = mean(size, na.rm = TRUE))
print(Bsize_means)

#account for experiment as a random variable
subblast<-lmer(size~Treatment*Population+ (1|Experiment), data=Bsize.sub) 
anova(subblast) #diff in treatment and populations!
summary(subblast)

###GASTRULA####
##measurement ####
#figure theme to keep things consistent
#new theme for size
theme_box <- function(base_size = 11, base_family = '') {
  theme_box() %+replace% 
    theme(text = element_text(size = 20)) 
}
####Populations#####
#Import the measurements, select and rename useful columns, and
#add a column for experiments 1 and 2
Gsize<-read.csv("2024_WBML/HE_Gastrula_measurements.csv") 

Gsize<-Gsize %>% select_if(~ !all(is.na(.))) #remove NA columns
Gsize<-Gsize %>% na.omit() ##remove NA rows

#Organize Populations with better labels
Gsize$Population<- factor(x=Gsize$Population, labels=c(
  "Cab" ="Urban, Cabrillo Beach",
  "Trp" = "Nonurban, Treasure Park",
  "Twp" = "Nonurban, Twin Pts",
  "Wh" = "Urban, White Pt"
))

#Reorder so urban sites are first and nonurban are second
Gsize$Population=factor(Gsize$Population,levels=c("Urban, Cabrillo Beach", "Urban, White Pt", "Nonurban, Treasure Park", "Nonurban, Twin Pts"))

#Box plots-Height
ggplot(Gsize, aes(x=Population, y=Height_mm, fill=Treatment)) +
  geom_boxplot()+ theme_box()+
  labs(title = 'Gastrula Height (mm)', x = 'Population', y = 'Gastrula Height (mm)')

#subsetting big boys
Gsize.sub1<- Gsize %>% 
  group_by(Stage, BeakerID, Population, MatePair, Treatment, Urban, Height_mm, Stomach_Length_mm, Experiment, Beaker) %>%
  filter(Height_mm <= 0.75) %>% 
  select(Height_mm)

heightplot<-
ggplot(Gsize.sub1, aes(x=Urban, y=Height_mm, fill=Treatment)) +
  geom_boxplot()+ theme_box()+
  labs(title = 'B. Gastrula', x = '', y = 'Gastrula Height (mm)')+
  scale_fill_manual(values=treat.colors)

stomachplot<-
ggplot(Gsize.sub1, aes(x=Urban, y=Stomach_Length_mm, fill=Treatment)) +
  geom_boxplot()+ theme_box()+
  labs(title = 'A. Gastrula', x = '', y = 'Gastrula Stomach Length (mm)')+
  scale_fill_manual(values=treat.colors)

#save plots as pictures
ggsave("gast_size.png",heightplot, width=20, height=15, units = "cm")
ggsave("gast_stomach.png",stomachplot, width=20, height=15, units = "cm") 

#test for size variation with treatment, blocking by population
Gh<-lmer(Height_mm~Treatment*Urban + (1|Experiment), data=Gsize.sub1) 
anova(Gh) 

Gs<-lmer(Stomach_Length_mm~Treatment*Urban+ (1|Experiment), data=Gsize.sub1)
anova(Gs)

Gsize_means <- Gsize.sub1 %>%
  group_by(Urban) %>%
  summarise(mean.SL = mean(Stomach_Length_mm, na.rm = TRUE))
print(Gsize_means)

##subset out white point to see if pattern still holds
Gsize.subW<- Gsize.sub1 %>% 
  group_by(Stage, BeakerID, Population, MatePair, Treatment, Urban, Height_mm, Stomach_Length_mm, Experiment, Beaker) %>%
  subset(Population!="Urban, White Pt")

#heightplot<-
ggplot(Gsize.subW, aes(x=Urban, y=Height_mm, fill=Treatment)) +
  geom_boxplot()+ theme_box()+
  labs(title = 'Gastrula Height (mm)', x = 'Population', y = 'Gastrula Height (mm)')+
  scale_fill_manual(values=treat.colors)

#stomachplot<-
ggplot(Gsize.subW, aes(x=Urban, y=Stomach_Length_mm, fill=Treatment)) +
  geom_boxplot()+ theme_box()+
  labs(title = 'Gastrula SL w/o Wh', x = 'Population', y = 'Gastrula Stomach Length (mm)')+
  scale_fill_manual(values=treat.colors)

GhW<-lmer(Height_mm~Treatment*Urban + (1|Experiment), data=Gsize.subW) 
anova(GhW) 

GsW<-lmer(Stomach_Length_mm~Treatment*Urban+ (1|Experiment), data=Gsize.subW)
anova(GsW)


Gsize_meansW <- Gsize.subW %>%
  group_by(Treatment) %>%
  summarise(mean.size = mean(Stomach_Length_mm, na.rm = TRUE))
print(Gsize_meansW)
######

###PLUTEUS####
Psize<-read.csv("2024_WBML/HE_Pluteus7dpf_measurements.csv")
#View(Psize)

#Psize<-Psize %>% select_if(~ !all(is.na(.))) #remove NA columns

#Organize Populations with better labels
Psize$Population<- factor(x=Psize$Population, labels=c(
  "Cab" ="Urban, Cabrillo Beach",
  "Trp" = "Nonurban, Treasure Park",
  "Twp" = "Nonurban, Twin Pts",
  "Wh" = "Urban, White Pt"
))

#Reorder so urban sites are first and nonurban are second
Psize$Population=factor(Psize$Population,levels=c("Urban, Cabrillo Beach", "Urban, White Pt", "Nonurban, Treasure Park", "Nonurban, Twin Pts"))

#Box plots-Body Length

sizeplut<-
ggplot(Psize, aes(x=Urban, y=BL_mm, fill=Treatment)) +
  geom_boxplot()+ theme_box()+
  labs(title = 'C. Pluteus', x = '', y = 'Pluteus Body Length (mm)')+
  scale_fill_manual(values = c("NP" = "#C77DFF", "C" = "#C8E9A0"))

ggsave("plut_sizeurban.png",sizeplut, width=20, height=15, units = "cm")

#average arm lengths
Psize <- Psize %>% 
  mutate(avg_AL=((AL1_mm + AL2_mm)/2))

#Box plots-Arm Length
armplut<-ggplot(Psize, aes(x=Urban, y=avg_AL, fill=Treatment)) +
  geom_boxplot()+ theme_box()+
  labs(title = 'C. Pluteus', x = '', y = 'Pluteus Arm Length (mm)')+
  scale_fill_manual(values = c("NP" = "#C77DFF", "C" = "#C8E9A0"))
ggsave("plut_armlenghturban.png",armplut, width=20, height=15, units = "cm")

pbl_u<-lmer(BL_mm~Treatment*Population + (1|Experiment), data=Psize) 
anova(pbl_u)

pal_u<-lmer(avg_AL~Treatment*Urban +(1|Experiment), data=Psize) 
anova(pal_u)

##testing variance in data!!
leveneTest(BL_mm ~ Treatment*Urban, data = Psize)


Psize_mean<- Psize %>%
  group_by(Urban) %>%
  summarise(mean.size = mean(avg_AL, na.rm = TRUE))
print(Psize_mean)
