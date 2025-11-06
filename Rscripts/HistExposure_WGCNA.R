##Historical Exposure
#WGCNA Script
#Madison Armstrong 
#Last Modified 7/16/2025
##tutorial followed: https://bioinformaticsworkbook.org/tutorials/wgcna.html#gsc.tab=0

setwd("~/Desktop/purp dev data/2024_WBML")
#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
#BiocManager::install("WGCNA")
BiocManager::install("impute")

library(tidyverse)     # tidyverse will pull in ggplot2, readr, other useful libraries
library(magrittr) # provides the %>% operator
library(WGCNA)    #no impute package so need to install separately
library(impute)
library(dplyr)
library(DESeq2)
library(tibble)

#figure theme to keep things consistent
theme_box <- function(base_size = 11, base_family = '') {
  theme_classic() %+replace% 
    theme(text = element_text(size = 20),  fill = NA, size = 1)
}


##Read in Metadata####
metadata <- read.table(file.path("HE_RNAdata/RNAhist_metadata.txt"), header=TRUE)
#load in normalized counts from "HistExposure_DEGs.R" script from DESeq2
#blastula, gastrula and pluteus DESeq matrices, "vsd.stage"
#vsd.blast, vsd.gast, vsd.plut


##1. Blastula####

####prep data####
vsd.blast #from deseq
library(genefilter) #filter out 


dim(vsd.blast) ##need to have rows be treatments & columns be genes, currently it is the opposite!
vsd.blast.mat <- assay(vsd.blast)
vsd.blast.flip <- t(vsd.blast.mat)
dim(vsd.blast.flip) ## now it is flipped woo

####WGCNA time####
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))

# Call the network topology analysis function, didn't use threads so will be slow
sft.blast = pickSoftThreshold(
  vsd.blast.flip,             # <= Input data
  powerVector = powers,
  verbose = 5
)

####plots####
##Pick a soft threshold power near the curve of the plot, so use plots to visualize!
par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft.blast$fitIndices[, 1],
     -sign(sft.blast$fitIndices[, 3]) * sft.blast$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)
text(sft.blast$fitIndices[, 1],
     -sign(sft.blast$fitIndices[, 3]) * sft.blast$fitIndices[, 2],
     labels = powers, cex = cex1, col = "purple"
)
abline(h = 0.90, col = "purple")
plot(sft.blast$fitIndices[, 1],
     sft.blast$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft.blast$fitIndices[, 1],
     sft.blast$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "purple")

##looks like best fit is 8!
picked_power.blast = 8
temp_cor.blast <- cor      
cor <- WGCNA::cor        # Force it to use WGCNA cor function (fix a namespace conflict issue)
netwk.blast <- blockwiseModules(vsd.blast.flip,                # <= input here
                          
                          # == Adjacency Function ==
                          power = picked_power.blast,         # <= power here
                          networkType = "signed",
                          
                          # == Tree and Block Options ==
                          deepSplit = 2,
                          pamRespectsDendro = FALSE,
                          # detectCutHeight = 0.75,
                          minModuleSize = 30,
                          maxBlockSize = 4000,
                          
                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          
                          # == TOM == Archive the run results in TOM file (saves time)
                          saveTOMs = TRUE,
                          saveTOMFileBase = "ER",
                          
                          # == Output Options
                          numericLabels = TRUE,
                          verbose = 3)

cor <- temp_cor.blast     # Return cor function to original namespace

####look at modules!####
# Convert labels to colors for plotting
mergedColors.blast = labels2colors(netwk.blast$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(
  netwk.blast$dendrograms[[1]],
  mergedColors.blast[netwk.blast$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05)

#pull out list
module_df.blast <- data.frame(
  gene_id = names(netwk.blast$colors),
  colors = labels2colors(netwk.blast$colors)
)

write_delim(module_df.blast,
            file = "blast.WGCNA.gene_modules.txt",
            delim = "\t") #all 34,000 genes and what module they correlate to

sort(table(module_df.blast$colors), decreasing = TRUE)
unique(module_df.blast$colors)
####eigengenes####
# Get Module Eigengenes per cluster
MEs0.blast <- moduleEigengenes(vsd.blast.flip, mergedColors.blast)$eigengenes

# Reorder modules so similar modules are next to each other
MEs0.blast <- orderMEs(MEs0.blast)
module_order.blast = names(MEs0.blast) %>% gsub("ME","", .)

# Add sample names
MEs0.blast$treatment = row.names(MEs0.blast)

# tidy & plot data
mME.blast = MEs0.blast %>%
  pivot_longer(-treatment) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order.blast)
  )

#read in metadata with sample names 
metadata.blast <- subset(metadata, Stage == "Blastula")

#merge treatment and pop for better visualization
metadata.blast <- metadata.blast %>%
  unite("Treatment_pop", Treatment, Population, sep = "_", remove = FALSE)

#join
mME.blast.joined <- mME.blast %>%
  left_join(metadata.blast, by = c("treatment" = "Library")) %>% 
  rename("treatment"="indiv")

#binarize variables of interest
metadata.blast.np<-metadata.blast %>% 
  mutate(Treatment.bin=ifelse(grepl("NP", Treatment), 1, 0)) %>% 
  select(12)

# Create binary dummy variables, one column per population level
metadata.blast.pop <- model.matrix(~ Population - 1, data = metadata.blast)
colnames(metadata.blast.pop) <- sub("^Population", "", colnames(metadata.blast.pop))

traits.blast<-cbind(metadata.blast.pop,metadata.blast.np)

#define number of genes & samples 

order <- c("C_Cab", "NP_Cab", "C_Wh", "NP_Wh", "C_Trp", "NP_Trp", "C_Twp", "NP_Twp")
mME.blast.joined$Treatment_pop <- factor(mME.blast.joined$Treatment_pop, levels = order)

###modules to NP/pops####
mME.blast.joined %>% ggplot(., aes(x=Treatment, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue3",
    high = "gold2",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-sample Relationships", y = "Modules", fill="corr")

####modules of interest####
###subset specific modules of interest to plot them!
mME.blast.joined.brown <-mME.blast.joined %>% 
  subset(mME.blast.joined$name=="brown")

ggplot(mME.blast.joined.brown, aes(x=urban, y=value, fill=Treatment))+
  geom_boxplot()+ theme_classic()+
  labs(title = 'Blastula: WGCNA brown module', x = 'Population', y = 'eigenvalues of brown module')+
  scale_fill_manual(values = c("NP" = "#C77DFF", "C" = "#C8E9A0"))


mME.blast.joined.yg <-mME.blast.joined %>% 
  subset(mME.blast.joined$name=="yellowgreen")

ggplot(mME.blast.joined.yg, aes(x=urban, y=value, fill=Treatment))+
  geom_boxplot()+ theme_box()+
  labs(title = 'Blastula', x = 'Population', y = 'eigenvalues of yellowgreen module')+
  scale_fill_manual(values = c("NP" = "#C77DFF", "C" = "#C8E9A0"))

yg<-lmer(value~Treatment*urban + (1|Experiment), data=mME.blast.joined.yg) 
anova(yg) 

##accounting for experiment by plotting the residuals
yg.exp<-lmer(value~(1|Experiment), data=mME.blast.joined.yg) 
res.yg.blast<-resid(yg.exp)
mME.blast.joined.yg<-cbind(mME.blast.joined.yg, res.yg.blast)

#replot data
WGCNA_yg<-
  ggplot(mME.blast.joined.yg, aes(x=urban, y=res.yg.blast, fill=Treatment))+
  geom_boxplot()+ theme_box()+
  labs(title = 'A. Blastula', x = 'YellowGreen', y = 'eigenvalues of module A (exp-controlled)')+
  scale_fill_manual(values = c("NP" = "#C77DFF", "C" = "#C8E9A0"))

ggsave("WGCNA_blastula_yg.png",WGCNA_yg, width=20, height=15, units = "cm") 
  
##line graphs#
# pick out a few modules of interest here--not a lot of overlap but let's try darkorange
modules_of_interest = ("blue")

# Pull out list of genes in that module
submods = module_df.blast %>%
  subset(colors %in% modules_of_interest)

row.names(module_df.blast) = module_df.blast$gene_id

# Get normalized expression for those genes
vsd.blast.mat [1:5, 1:10]

subexpr = vsd.blast.mat [submods$gene_id,]

submods_df <- as_tibble(subexpr, .name_repair = "minimal") %>%
  mutate(gene_id = rownames(subexpr)) %>%
  pivot_longer(-gene_id, names_to = "name", values_to = "value") %>%
  mutate(module = module_df.blast[gene_id, "colors"])


##add metadata
submods_df2 <- left_join(submods_df, metadata.blast %>% select(Library, Population),
                       by = c("name" = "Library"))

submods_df2 %>% ggplot(., aes(x=Population.y, y=value, group=gene_id)) +
  geom_line(aes(color = module),
            alpha = 0.2) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90)
  ) +
  facet_grid(rows = vars(module)) +
  labs(x = "treatment",
       y = "normalized expression")


####compare with outlier genes####
DEGs_blast<-read.csv("RNA_data/HE.blast.results.sig.urbtreat.full.genes.csv", header=TRUE) ##outlier genes
head(module_df.blast) #WGCNA modules
module_df.blast$gene_id<- sub("\\..*", "", module_df.blast$gene_id) #remove version

sig_mods.blast<-left_join(DEGs_blast, module_df.blast, by= c("X"= "gene_id"))
sort(table(sig_mods.blast$colors), decreasing = TRUE)
write.csv(sig_mods.blast,"RNA_data/HE.blast.WGCNA.sigmods.csv")

##2. Gastrula####

####prep data####
vsd.gast
library(genefilter) #filter out 


dim(vsd.gast) ##need to have rows be treatments & columns be genes, currently it is the opposite!
vsd.gast.mat <- assay(vsd.gast)
vsd.gast.flip <- t(vsd.gast.mat)
dim(vsd.gast.flip) ## now it is flipped woo

####WGCNA time####
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))

# Call the network topology analysis function, didn't use threads so will be slow
sft.gast = pickSoftThreshold(
  vsd.gast.flip,             # <= Input data
  powerVector = powers,
  verbose = 5
)

####plots####
##Pick a soft threshold power near the curve of the plot, so use plots to visualize!
par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft.gast$fitIndices[, 1],
     -sign(sft.gast$fitIndices[, 3]) * sft.gast$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)
text(sft.gast$fitIndices[, 1],
     -sign(sft.gast$fitIndices[, 3]) * sft.gast$fitIndices[, 2],
     labels = powers, cex = cex1, col = "purple"
)
abline(h = 0.90, col = "purple")
plot(sft.gast$fitIndices[, 1],
     sft.gast$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft.gast$fitIndices[, 1],
     sft.gast$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "purple")

##looks like best fit is 6!
picked_power.gast = 6
temp_cor.gast <- cor      
cor <- WGCNA::cor        # Force it to use WGCNA cor function (fix a namespace conflict issue)
netwk.gast <- blockwiseModules(vsd.gast.flip,                # <= input here
                                
                                # == Adjacency Function ==
                                power = picked_power.gast,         # <= power here
                                networkType = "signed",
                                
                                # == Tree and Block Options ==
                                deepSplit = 2,
                                pamRespectsDendro = FALSE,
                                # detectCutHeight = 0.75,
                                minModuleSize = 30,
                                maxBlockSize = 4000,
                                
                                # == Module Adjustments ==
                                reassignThreshold = 0,
                                mergeCutHeight = 0.25,
                                
                                # == TOM == Archive the run results in TOM file (saves time)
                                saveTOMs = TRUE,
                                saveTOMFileBase = "ER",
                                
                                # == Output Options
                                numericLabels = TRUE,
                                verbose = 3)

cor <- temp_cor.gast     # Return cor function to original namespace

####look at modules!####
# Convert labels to colors for plotting
mergedColors.gast = labels2colors(netwk.gast$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(
  netwk.gast$dendrograms[[1]],
  mergedColors.gast[netwk.gast$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )

#pull out list
module_df.gast <- data.frame(
  gene_id = names(netwk.gast$colors),
  colors = labels2colors(netwk.gast$colors)
)

write_delim(module_df.gast,
            file = "gast.WGCNA.gene_modules.txt",
            delim = "\t") #all 34,000 genes and what module they correlate to

sort(table(module_df.gast$colors), decreasing = TRUE)
unique(module_df.gast$colors)

####eigengenes####
# Get Module Eigengenes per cluster
MEs0.gast <- moduleEigengenes(vsd.gast.flip, mergedColors.gast)$eigengenes

# Reorder modules so similar modules are next to each other
MEs0.gast <- orderMEs(MEs0.gast)
module_order.gast = names(MEs0.gast) %>% gsub("ME","", .)

# Add sample names
MEs0.gast$treatment = row.names(MEs0.gast)


# tidy & plot data
mME.gast = MEs0.gast %>%
  pivot_longer(-treatment) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order.gast)
  )

#read in metadata with sample names 
metadata.gast <- subset(metadata, Stage == "Gastrula")

#merge treatment and pop for better visualization
metadata.gast <- metadata.gast %>%
  unite("Treatment_pop", Treatment, Population, sep = "_", remove = FALSE)

#join
mME.gast.joined <- mME.gast %>%
  left_join(metadata.gast, by = c("treatment" = "Library")) %>% 
  rename("treatment"="indiv")


#binarize variables of interest
metadata.gast.np<-metadata.gast %>% 
  mutate(Treatment.bin=ifelse(grepl("NP", Treatment), 1, 0)) %>% 
  select(12)

# Create binary dummy variables, one column per population level
metadata.gast.pop <- model.matrix(~ Population - 1, data = metadata.gast)
colnames(metadata.gast.pop) <- sub("^Population", "", colnames(metadata.gast.pop))

traits.gast<-cbind(metadata.gast.pop,metadata.gast.np)

#define number of genes & samples 

order <- c("C_Cab", "NP_Cab", "C_Wh", "NP_Wh", "C_Trp", "NP_Trp", "C_Twp", "NP_Twp")
mME.gast.joined$Treatment_pop <- factor(mME.gast.joined$Treatment_pop, levels = order)

###modules to NP/pops####
mME.gast.joined %>% ggplot(., aes(x=Treatment_pop, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "orange",
    high = "green3",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-sample Relationships", y = "Modules", fill="corr")

####modules of interest####
###subset specific modules of interest to plot them!
mME.gast.joined.cyan <-mME.gast.joined %>% 
  subset(mME.gast.joined$name=="cyan")

ggplot(mME.gast.joined.cyan, aes(x=urban, y=value, fill=Treatment))+
  geom_boxplot()+ theme_box()+
  labs(title = 'Gastrula', x = 'Population', y = 'eigenvalues of cyan module')+
  scale_fill_manual(values = c("NP" = "#C77DFF", "C" = "#C8E9A0"))

cy<-lmer(value~Treatment+ urban+ (1|Experiment), data=mME.gast.joined.cyan) 
anova(cy) 

##accounting for experiment by plotting the residuals
cy.exp<-lmer(value~(1|Experiment), data=mME.gast.joined.cyan) 
res.cy.gast<-resid(cy.exp)
mME.gast.joined.cyan<-cbind(mME.gast.joined.cyan, res.cy.gast)

#replot data
WGCNA_cy<-
  ggplot(mME.gast.joined.cyan, aes(x=urban, y=res.cy.plut, fill=Treatment))+
  geom_boxplot()+ theme_box()+
  labs(title = 'B. Gastrula', x = 'Cyan', y = 'eigenvalues of module A (exp-controlled)')+
  scale_fill_manual(values = c("NP" = "#C77DFF", "C" = "#C8E9A0"))

ggsave("WGCNA_gastrula_cy.png",WGCNA_cy, width=20, height=15, units = "cm") 

mME.gast.joined.turq <-mME.gast.joined %>% 
  subset(mME.gast.joined$name=="turquoise")

ggplot(mME.gast.joined.turq, aes(x=urban, y=value, fill=Treatment))+
  geom_boxplot()+ theme_box()+
  labs(title = 'Gastrula', x = 'Population', y = 'eigenvalues of turquoise module')+
  scale_fill_manual(values = c("NP" = "#C77DFF", "C" = "#C8E9A0"))

tq.g<-lmer(value~Treatment+urban + (1|Experiment), data=mME.gast.joined.turq) 
anova(tq.g)  #no sig effect for this module

####compare with outlier genes####
DEGs_gast<-read.csv("RNA_data/HE.gast.results.sig.urbtreat.full.genes.csv", header=TRUE) ##outlier genes
head(module_df.gast) #WGCNA modules
module_df.gast$gene_id<- sub("\\..*", "", module_df.gast$gene_id) #remove version

sig_mods.gast<-left_join(DEGs_gast, module_df.gast, by= c("X"= "gene_id"))
sort(table(sig_mods.gast$colors), decreasing = TRUE)

write.csv(sig_mods.gast,"RNA_data/HE.gast.WGCNA.sigmods.csv")

##3. Pluteus####

####prep data####
vsd.plut
library(genefilter) #filter out 


dim(vsd.plut) ##need to have rows be treatments & columns be genes, currently it is the opposite!
vsd.plut.mat <- assay(vsd.plut)
vsd.plut.flip <- t(vsd.plut.mat)
dim(vsd.plut.flip) ## now it is flipped woo

####WGCNA time####
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))

# Call the network topology analysis function, didn't use threads so will be slow
sft.plut = pickSoftThreshold(
  vsd.plut.flip,             # <= Input data
  powerVector = powers,
  verbose = 5
)

####plots####
##Pick a soft threshold power near the curve of the plot, so use plots to visualize!
par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft.plut$fitIndices[, 1],
     -sign(sft.plut$fitIndices[, 3]) * sft.plut$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)
text(sft.plut$fitIndices[, 1],
     -sign(sft.plut$fitIndices[, 3]) * sft.plut$fitIndices[, 2],
     labels = powers, cex = cex1, col = "blue"
)
abline(h = 0.90, col = "blue")
plot(sft.plut$fitIndices[, 1],
     sft.plut$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft.plut$fitIndices[, 1],
     sft.plut$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "blue")

##looks like best fit is 7!
picked_power.plut = 7
temp_cor.plut <- cor      
cor <- WGCNA::cor        # Force it to use WGCNA cor function (fix a namespace conflict issue)
netwk.plut <- blockwiseModules(vsd.plut.flip,                # <= input here
                               
                               # == Adjacency Function ==
                               power = picked_power.plut,         # <= power here
                               networkType = "signed",
                               
                               # == Tree and Block Options ==
                               deepSplit = 2,
                               pamRespectsDendro = FALSE,
                               # detectCutHeight = 0.75,
                               minModuleSize = 30,
                               maxBlockSize = 4000,
                               
                               # == Module Adjustments ==
                               reassignThreshold = 0,
                               mergeCutHeight = 0.25,
                               
                               # == TOM == Archive the run results in TOM file (saves time)
                               saveTOMs = TRUE,
                               saveTOMFileBase = "ER",
                               
                               # == Output Options
                               numericLabels = TRUE,
                               verbose = 3)

cor <- temp_cor.plut     # Return cor function to original namespace

####look at modules!####
# Convert labels to colors for plotting
mergedColors.plut = labels2colors(netwk.plut$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(
  netwk.plut$dendrograms[[1]],
  mergedColors.plut[netwk.plut$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )

#pull out list
module_df.plut <- data.frame(
  gene_id = names(netwk.plut$colors),
  colors = labels2colors(netwk.plut$colors)
)

write_delim(module_df.plut,
            file = "plut.WGCNA.gene_modules.txt",
            delim = "\t") #all 34,000 genes and what module they correlate to

sort(table(module_df.plut$colors), decreasing = TRUE)
unique(module_df.gast$colors)
####eigengenes####
# Get Module Eigengenes per cluster
MEs0.plut <- moduleEigengenes(vsd.plut.flip, mergedColors.plut)$eigengenes

# Reorder modules so similar modules are next to each other
MEs0.plut <- orderMEs(MEs0.plut)
module_order.plut = names(MEs0.plut) %>% gsub("ME","", .)

# Add sample names
MEs0.plut$treatment = row.names(MEs0.plut)


# tidy & plot data
mME.plut = MEs0.plut %>%
  pivot_longer(-treatment) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order.plut)
  )

#read in metadata with sample names 
metadata.plut <- subset(metadata, Stage == "Pluteus")

#merge treatment and pop for better visualization
metadata.plut <- metadata.plut %>%
  unite("Treatment_pop", Treatment, Population, sep = "_", remove = FALSE)

#join
mME.plut.joined <- mME.plut %>%
  left_join(metadata.plut, by = c("treatment" = "Library")) %>% 
  rename("treatment"="indiv")


#binarize variables of interest
metadata.plut.np<-metadata.plut %>% 
  mutate(Treatment.bin=ifelse(grepl("NP", Treatment), 1, 0)) %>% 
  select(12)

# Create binary dummy variables, one column per population level
metadata.plut.pop <- model.matrix(~ Population - 1, data = metadata.plut)
colnames(metadata.plut.pop) <- sub("^Population", "", colnames(metadata.plut.pop))

traits.plut<-cbind(metadata.plut.pop,metadata.plut.np)

#define number of genes & samples 

order <- c("C_Cab", "NP_Cab", "C_Wh", "NP_Wh", "C_Trp", "NP_Trp", "C_Twp", "NP_Twp")
mME.plut.joined$Treatment_pop <- factor(mME.plut.joined$Treatment_pop, levels = order)

###modules to NP/pops####
mME.plut.joined %>% ggplot(., aes(x=Population, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "maroon",
    high = "grey",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-sample Relationships", y = "Modules", fill="corr")

####modules of interest####
mME.plut.joined.pturq <-mME.plut.joined %>% 
  subset(mME.plut.joined$name=="paleturquoise")

  ggplot(mME.plut.joined.pturq, aes(x=urban, y=value, fill=Treatment))+
  geom_boxplot()+ theme_box()+
  labs(title = 'Pluteus', x = 'Population', y = 'eigenvalues of paleturquoise module')+
  scale_fill_manual(values = c("NP" = "#C77DFF", "C" = "#C8E9A0"))

ptq<-lmer(value~Treatment+ urban + (1|Experiment), data=mME.plut.joined.pturq) 
anova(ptq) 

##accounting for experiment by plotting the residuals
ptq.exp<-lmer(value~(1|Experiment), data=mME.plut.joined.pturq) 
res.ptq.plut<-resid(ptq.exp)
mME.plut.joined.pturq<-cbind(mME.plut.joined.pturq, res.ptq.plut)

#replot data
WGCNA_plut_paleturq<-
ggplot(mME.plut.joined.pturq, aes(x=urban, y=res.ptq.plut, fill=Treatment))+
  geom_boxplot()+ theme_box()+
  labs(title = 'C. Pluteus', x = 'Pale Turq', y = 'eigenvalues of module A (exp-controlled)')+
  scale_fill_manual(values = c("NP" = "#C77DFF", "C" = "#C8E9A0"))

ggsave("WGCNA_plut_paleturq.png",WGCNA_plut_paleturq, width=20, height=15, units = "cm") 


###subset specific modules of interest to plot them!
mME.plut.joined.turq <-mME.plut.joined %>% 
  subset(mME.plut.joined$name=="turquoise")

ggplot(mME.plut.joined.turq, aes(x=urban, y=value, fill=Treatment))+
  geom_boxplot()+ theme_box()+
  labs(title = 'Pluteus', x = 'Population', y = 'eigenvalues of turquoise module')+
  scale_fill_manual(values = c("NP" = "#C77DFF", "C" = "#C8E9A0"))

tq<-lmer(value~Treatment+urban + (1|Experiment), data=mME.plut.joined.turq) 
anova(tq) 

##accounting for experiment by plotting the residuals
tq.exp<-lmer(value~(1|Experiment), data=mME.plut.joined.turq) 
res.tq.plut<-resid(tq.exp)
mME.plut.joined.turq.exp<-cbind(mME.plut.joined.turq, res.tq.plut)

#replot data
WGCNA_plut_turq<-
ggplot(mME.plut.joined.turq.exp, aes(x=urban, y=res.tq.plut, fill=Treatment))+
  geom_boxplot()+ theme_box()+
  labs(title = 'D. Pluteus', x = 'Turquoise', y = 'eigenvalues of module B (exp-controlled)')+
  scale_fill_manual(values = c("NP" = "#C77DFF", "C" = "#C8E9A0"))

ggsave("WGCNA_plut_turq.png",WGCNA_plut_turq, width=20, height=15, units = "cm") 


####compare with outlier genes####
DEGs_plut<-read.csv("RNA_data/HE.plut.results.sig.urbtreat.full.genes.csv", header=TRUE) ##outlier genes
head(module_df.plut) #WGCNA modules
module_df.plut$gene_id<- sub("\\..*", "", module_df.plut$gene_id) #remove version

sig_mods.plut<-left_join(DEGs_plut, module_df.plut, by= c("X"= "gene_id"))
sort(table(sig_mods.plut$colors), decreasing = TRUE)

write.csv(sig_mods.plut,"RNA_data/HE.plut.WGCNA.sigmods.csv")
