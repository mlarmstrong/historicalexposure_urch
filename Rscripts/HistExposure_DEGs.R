#Dose Exposure--Gene Expression Scripts
#website with tutorial information: 
#https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#rich-visualization-and-reporting-of-results

#Madison Armstrong 
#Last Modified 7/14/2025

setwd("~/Desktop/purp dev data/2024_WBML")


#Libraries
library(ggplot2)
library(tidyr)
library(dplyr)
library(tidyverse)
library(RColorBrewer)
library(DESeq2)
library(tximport)
library(variancePartition)
library(pheatmap)
library(ashr)
library(viridis)

#figure theme to keep things consistent
theme_box <- function(base_size = 11, base_family = '') {
  theme_classic() %+replace% 
    theme(text = element_text(size = 20), 
          panel.border = element_rect(color = "black", fill = NA, size = 1)) 
}
##Organize and Read in Data ####
#create a vector of filenames, reading in a table that contains the sample IDs
dir <- "~/Desktop/purp dev data/2024_WBML/Salmon_quant_HE"

#read in metadata
metadata <- read.table(file.path("HE_RNAdata/RNAhist_metadata.txt"), header=TRUE, row.names=1)

#subset data
files_list<-paste(rownames(metadata))
files <- file.path(dir,files_list, "quant.sf")
all(file.exists(files))
names <- gsub("_quant|_USP.*", "", files_list)
names(files) <- names

## Import quantification tables 
#Read in matrix of RSEM expected read counts 
#all_quant.sf is normalized... we don't want that!!
all.quant <- read.delim("Salmon_quant_HE/all_quant.normalized.sf")
tx2gene<-data.frame(all.quant[,1], all.quant[,1])
colnames(tx2gene)<-c("TXNAME", "GENEID")

data <- tximport(files, type="salmon", tx2gene=tx2gene)


#make sure things are in the correct order
all(rownames(metadata) %in% colnames(data))
all(rownames(metadata) == colnames(data))


dds.start <- DESeqDataSetFromTximport(data,
                                colData = metadata,
                                design = ~ Treatment + urban) 

dds.start

#pre-filtering step: using same as Leslie
keep <- rowSums(counts(dds.start)) >= 10
dds.start <- dds.start[keep,]

# identify genes that pass expression cutoff
isexpr <- rowSums(fpm(dds.start) > 1) >= 0.5 * ncol(dds.start)

quantLog.start <- log2(fpm(dds.start)[isexpr, ] + 1)
# Check sample names alignment
all(colnames(quantLog.start) == rownames(metadata))  #should be true!!

#varpar analysis to understand what variation is driving the pcas
#https://www.bioconductor.org/packages/release/bioc/vignettes/variancePartition/inst/doc/variancePartition.html

# Specify variables to consider
# Treatment is a fixed effect, main thing I am investigating + I want to look at Population (urban nested in that) and Stage
form <- ~ Treatment + Population + Stage
form2 <- ~ Treatment + urban + Stage
# Fit model and extract results
# 1) fit linear mixed model on gene expression
# If categorical variables are specified, a linear mixed model is used
# If all variables are modeled as fixed effects, a linear model is used
# each entry in results is a regression model fit on a single gene
# 2) extract variance fractions from each model fit
# for each gene, returns fraction of variation attributable
#       to each variable
# Interpretation: the variance explained by each variables
# after correcting for all other variables
# Note that geneExpr can either be a matrix, and EList output by voom() in the limma package,or an ExpressionSet

# Re-run variance partitioning
varPart <- fitExtractVarPartModel(quantLog.start, form, metadata)
varPart2 <- fitExtractVarPartModel(quantLog.start, form2, metadata)
# sort variables (i.e. columns) by median fraction of variance explained
vp <- sortCols(varPart)
vp2 <- sortCols(varPart2)
# Bar plot of variance fractions for the first 10 genes
plotPercentBars(vp[1:10, ], colors=my_colors)

# violin plot of contribution of each variable to total variance
plotVarPart(vp)
C <- canCorPairs(form2, metadata)
# Plot correlation matrix between all pairs of variables
plotCorrMatrix(C)

##SEPARATE OUT STAGE FIRST####
    ###1. Blastula####
####READ IN####
metadata.blast <- subset(metadata, Stage == "Blastula")

#subset data
files_list_blast<-paste(rownames(metadata.blast))
files.blast <- file.path(dir,files_list_blast, "quant.sf")
all(file.exists(files.blast))
names.blast <- gsub("_quant|_USP.*", "", files_list_blast)
names(files.blast) <- names.blast
data.blast <- tximport(files.blast, type="salmon", tx2gene=tx2gene) #52 blastula

####DeSeq####
#With the count matrix, cts, and the sample information, coldata, we can construct a DESeqDataSet:
dds.blast <- DESeqDataSetFromTximport(data.blast,
                                     colData = metadata.blast,
                                     design = ~ Treatment + urban) 

dds.blast 

#pre-filtering step: using same as Leslie
keep <- rowSums(counts(dds.blast)) >= 10
dds.blast <- dds.blast[keep,]

#variance stabilizing transformation-- extracting transformed values
vsd.blast <- vst(dds.blast, blind=FALSE)

#### DEGs####
hist.blast <- DESeq(dds.blast)
# Report: Number of transcripts
num_transcripts.B <- nrow(hist.blast)
num_transcripts.B #34008 transcripts
resultsNames(hist.blast) #populations compared & treatment vs control compared

#urban contrast
#filter by padj=0.05 and then extract logfoldchange
resblast<-results(hist.blast,(c(contrast="urban", "urban", "nonurban")) ) #gives readout of padj, log2FoldChange and other stuff
summary(resblast)
resblast
#####SAVE FOR TOPGO
#write.csv(resblast,file="RNA_data/HEresults.blast.urb.csv", row.names=TRUE)

##look at treatment comparison
#resblast.treat <- results(hist.blast, (c(contrast="Treatment", "NP", "C")))
#ummary(resblast.treat)
#####SAVE FOR TOPGO
#write.csv(resblast.treat,file="RNA_data/HEresults.blast.urbtreat.csv", row.names=TRUE)

####Heat Map####
# Subset colData to only include samples from the relevant contrast-- pulls outliers
resblast_outliers <- resblast %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < 0.05)
dim(resblast_outliers)

# Remove genes with zero counts across all samples (to prevent vst() errors)
hist.blast <- hist.blast[rowSums(counts(hist.blast)) > 0, ]

# Perform Variance Stabilizing Transformation (VST)
hist.blast.vsd <- vst(hist.blast, blind=FALSE)

# Extract the transformed expression matrix for the top genes
resblast_top_genes_matrix_var <- assay(hist.blast.vsd)[resblast_outliers$gene, ]

df.b <- as.data.frame(colData(hist.blast.vsd)[,c("Treatment","Population")])
df.b$Treatment <- factor(df.b$Treatment, levels = c("C", "NP")) 
df.b$Population <- factor(df.b$Population, levels = c("Cab", "Wh", "Trp", "Twp")) 

##Get things in order!
blast_ordered_samples <- df.b %>%
  arrange(Population, Treatment) %>%
  rownames()

blast_ordered_matrix <- resblast_top_genes_matrix_var[, blast_ordered_samples]
df.b.ordered <- df.b[blast_ordered_samples, ]
                                 
annot_colors=list(
  Population=c("Cab"="#5C5649", "Wh"="#A69C8B", "Trp"="#0F6DD2", "Twp"="#99D1EC"),
  Treatment=c("NP" = "#C77DFF", "C" = "#C8E9A0"))

pheatmap(blast_ordered_matrix, cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df.b.ordered, 
         annotation_colors=annot_colors,
         scale="row",
         color = magma(100))

####Varpar####
# identify genes that pass expression cutoff
isexpr <- rowSums(fpm(dds.blast) > 1) >= 0.5 * ncol(dds.blast)

quantLog.b <- log2(fpm(dds.blast)[isexpr, ] + 1)
# Check sample names alignment
all(colnames(quantLog.b) == rownames(metadata.blast))  #should be true!!

form.b <- ~ Treatment + Population +Experiment
varPartb <- fitExtractVarPartModel(quantLog.b, form.b, metadata.blast)

vpb <- sortCols(varPartb)
# violin plot of contribution of each variable to total variance
plotVarPart(vpb)


####PCAS####
pcaData_blast<-
  plotPCA(vsd.blast, intgroup=c("Treatment", "Population", "urban", "Experiment"),  pcsToUse=1:2, returnData=TRUE) #using ntop=500 top features by variance
percentVar <- round(100 * attr(pcaData_blast, "percentVar"))


ggplot(pcaData_blast, aes(x=PC1, y=PC2, shape=Experiment, color=urban)) +
  geom_point(size=4) +
  #stat_ellipse(aes(x = PC1,y=PC2,color=factor(Experiment))) +
  #facet_wrap(~Population)+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  labs(title= "Blastula")+
  scale_color_brewer(palette="Set1")+
  theme(legend.position="right") + theme_box()

#linear models
##PC5 sig by population (and pop MP) AND experiment
## PC1 sig by experiment...
blast1<-lm(PC1~Population+Treatment+Experiment, data=pcaData_blast) 
anova(blast1)


####2. Gastrula####
###READ IN####
metadata.gast <- subset(metadata, Stage == "Gastrula")

#subset data
files_list_gast<-paste(rownames(metadata.gast))
files.gast <- file.path(dir,files_list_gast, "quant.sf")
all(file.exists(files.gast))
names.gast <- gsub("_quant|_USP.*", "", files_list_gast)
names(files.gast) <- names.gast
data.gast <- tximport(files.gast, type="salmon", tx2gene=tx2gene) #86 gastula

###DeSeq####
#With the count matrix, cts, and the sample information, coldata, we can construct a DESeqDataSet:
dds.gast <- DESeqDataSetFromTximport(data.gast,
                                      colData = metadata.gast,
                                      design = ~Treatment+ urban) 

dds.gast 

#pre-filtering step: using same as Leslie
keep <- rowSums(counts(dds.gast)) >= 10
dds.gast <- dds.gast[keep,]

#variance stabilizing transformation-- extracting transformed values
vsd.gast <- vst(dds.gast, blind=FALSE)

#### DEGs####
hist.gast <- DESeq(dds.gast)
# Report: Number of transcripts
num_transcripts.G <- nrow(hist.gast)
num_transcripts.G 
resultsNames(hist.gast) #populations compared & treatment vs control compared


#filter by padj=0.05 and then extract logfoldchange
resgast<-results(hist.gast, (c(contrast="urban", "urban", "nonurban")) ) #gives readout of padj, log2FoldChange and other stuff
summary(resgast)
resgast

#####SAVE FOR TOPGO
#write.csv(resgast,file="RNA_data/HEresults.gast.urb.csv", row.names=TRUE)


##look at treatment comparison
#resgast.treat <- results(hist.gast, (c(contrast="Treatment", "NP", "C")))
#summary(resgast.treat)
#####SAVE FOR TOPGO
#write.csv(resgast.treat,file="RNA_data/HEresults.gast.urbtreat.csv", row.names=TRUE)

####Heat Map####
# Subset colData to only include samples from the relevant contrast-- pulls outliers
resgast_outliers <- resgast %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < 0.05)
dim(resgast_outliers) #18 outliers woo// for treatment + urb only 5

# Remove genes with zero counts across all samples (to prevent vst() errors)
hist.gast <- hist.gast[rowSums(counts(hist.gast)) > 0, ]

# Perform Variance Stabilizing Transformation (VST)
hist.gast.vsd <- vst(hist.gast, blind=FALSE)

# Extract the transformed expression matrix for the top genes
resgast_top_genes_matrix_var <- assay(hist.gast.vsd)[resgast_outliers$gene, ]

df.g <- as.data.frame(colData(hist.gast.vsd)[,c("Treatment","Population")])
df.g$Treatment <- factor(df.g$Treatment, levels = c("C", "NP")) 
df.g$Population <- factor(df.g$Population, levels = c("Cab", "Wh", "Trp", "Twp")) 

##Get things in order!
gast_ordered_samples <- df.g %>%
  arrange(Population, Treatment) %>%
  rownames()

gast_ordered_matrix <- resgast_top_genes_matrix_var[, gast_ordered_samples]
df.g.ordered <- df.g[gast_ordered_samples, ]

annot_colors=list(
  Population=c("Cab"="#5C5649", "Wh"="#A69C8B", "Trp"="#0F6DD2", "Twp"="#99D1EC"),
  Treatment=c("NP" = "#C77DFF", "C" = "#C8E9A0"))

pheatmap(gast_ordered_matrix, cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df.g.ordered, 
         annotation_colors=annot_colors,
         scale="row",
         color = magma(100))

##Varpar####
# identify genes that pass expression cutoff
G <- rowSums(fpm(dds.gast) > 1) >= 0.5 * ncol(dds.blast)

quantLog.g <- log2(fpm(dds.gast)[G, ] + 1)

form.g <- ~ Treatment + Population +Experiment
varPartg <- fitExtractVarPartModel(quantLog.g, form.g, metadata.gast)

vpg <- sortCols(varPartg)
# violin plot of contribution of each variable to total variance
plotVarPart(vpg)
### PCAS####
pcaData_gast<-
  plotPCA(vsd.gast, intgroup=c("Treatment", "Population", "urban", "Pop_MP", "Experiment"),  pcsToUse=5:6, returnData=TRUE) #using ntop=500 top features by variance
percentVar <- round(100 * attr(pcaData_gast, "percentVar"))


ggplot(pcaData_gast, aes(x=PC5, y=PC6,  shape=Treatment, color=Population)) +
  geom_point(size=4) +
  #stat_ellipse() +
  facet_wrap(~Experiment)+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  labs(title= "Gastrula")+
  scale_color_brewer(palette="Set2")+
  theme(legend.position="right") + theme_box()

#linear models
gast1<-lm(PC5~Population+Treatment+Experiment, data=pcaData_gast) 
anova(gast1)

####3. Pluteus####
###READ IN####
metadata.plut <- subset(metadata, Stage == "Pluteus")

#subset data
files_list_plut<-paste(rownames(metadata.plut))
files.plut <- file.path(dir,files_list_plut, "quant.sf")
all(file.exists(files.plut))
names.plut <- gsub("_quant|_USP.*", "", files_list_plut)
names(files.plut) <- names.plut
data.plut <- tximport(files.plut, type="salmon", tx2gene=tx2gene) #86 pluteus

###DeSeq####
#With the count matrix, cts, and the sample information, coldata, we can construct a DESeqDataSet:
dds.plut <- DESeqDataSetFromTximport(data.plut,
                                     colData = metadata.plut,
                                     design = ~ Treatment + urban) 

dds.plut 

#pre-filtering step: using same as Leslie
keep <- rowSums(counts(dds.plut)) >= 10
dds.plut <- dds.plut[keep,]

#variance stabilizing transformation-- extracting transformed values
vsd.plut <- vst(dds.plut, blind=FALSE)

#### DEGs####
hist.plut <- DESeq(dds.plut)
# Report: Number of transcripts
num_transcripts.P <- nrow(hist.plut)
num_transcripts.P # 37315 transcripts
resultsNames(hist.plut) #populations compared & treatment vs control compared

#filter by padj=0.05 and then extract logfoldchange
##population or urban comparison
resplut<-results(hist.plut) #gives readout of padj, log2FoldChange and other stuff
summary(resplut)
#####SAVE FOR TOPGO
#write.csv(resplut,file="RNA_data/HEresults.plut.urb.csv", row.names=TRUE)

##look at treatment comparison
resplut.treat <- results(hist.plut, (c(contrast="Treatment", "NP", "C")))
summary(resplut.treat)

#####SAVE FOR TOPGO
write.csv(resplut.treat,file="RNA_data/HEresults.plut.urbtreat.csv", row.names=TRUE)

####Heat Map####
# Subset colData to only include samples from the relevant contrast-- pulls outliers
resplut_outliers <- resplut.treat %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < 0.05)
dim(resplut_outliers) #23 outliers ~treat+pop // 25 outliers for treatment + urban

##look at urban comparison
resplut.urban <- results(hist.plut, (c(contrast="urban", "urban", "nonurban")))

resplut_outliers.urban <- resplut.urban %>% #pulling outliers sig for int
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < 0.05)

dim(resplut_outliers.urban) #25 DEGs

# Remove genes with zero counts across all samples (to prevent vst() errors)
hist.plut <- hist.plut[rowSums(counts(hist.plut)) > 0, ]

# Perform Variance Stabilizing Transformation (VST)
hist.plut.vsd <- vst(hist.plut, blind=FALSE)

# Extract the transformed expression matrix for the top genes
resplut_top_genes_matrix_var <- assay(hist.plut.vsd)[resplut_outliers.urban$gene, ]

df.p <- as.data.frame(colData(hist.plut.vsd)[,c("Treatment","Population")])
df.p$Treatment <- factor(df.p$Treatment, levels = c("C", "NP")) 
df.p$Population <- factor(df.p$Population, levels = c("Cab", "Wh", "Trp", "Twp")) 

##Get things in order!
plut_ordered_samples <- df.p %>%
  arrange(Population, Treatment) %>%
  rownames()

plut_ordered_matrix <- resplut_top_genes_matrix_var[, plut_ordered_samples]
df.p.ordered <- df.p[plut_ordered_samples, ]

annot_colors=list(
  Population=c("Cab"="#5C5649", "Wh"="#A69C8B", "Trp"="#0F6DD2", "Twp"="#99D1EC"),
  Treatment=c("NP" = "#C77DFF", "C" = "#C8E9A0"))

pheatmap(plut_ordered_matrix, cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df.p.ordered, 
         annotation_colors=annot_colors,
         scale="row",
         color = magma(100))

##Varpar####
# identify genes that pass expression cutoff
P <- rowSums(fpm(dds.plut) > 1) >= 0.5 * ncol(dds.plut)

quantLog.p <- log2(fpm(dds.plut)[P, ] + 1)

form.p <- ~ Treatment + urban +Experiment
varPartp <- fitExtractVarPartModel(quantLog.p, form.p, metadata.plut)

vp <- sortCols(varPartp)
# violin plot of contribution of each variable to total variance
plotVarPart(vp)

plotPercentBars(vp[1:5, ], colors=my_colors)

### PCAS####

pcaData_plut<-
  plotPCA(vsd.plut, intgroup=c("Treatment", "Population", "urban", "Pop_MP", "Experiment"),  pcsToUse=3:4, returnData=TRUE) #using ntop=500 top features by variance
percentVar <- round(100 * attr(pcaData_plut, "percentVar"))



ggplot(pcaData_plut, aes(x=PC3, y=PC4, shape=urban, color=Treatment)) +
  geom_point(size=4) +
  stat_ellipse(aes(x = PC3,y=PC4, lty=factor(urban),color=factor(Treatment))) +
  #facet_wrap(~Population)+
  xlab(paste0("PC3: ",percentVar[1],"% variance")) +
  ylab(paste0("PC4: ",percentVar[2],"% variance")) + 
  labs(title= "Pluteus")+
  scale_color_brewer(palette="Set2")+
  theme(legend.position="right") + theme_box()

#linear models
plut1<-lm(PC6~Population+Treatment+Experiment, data=pcaData_plut) 
anova(plut1)
