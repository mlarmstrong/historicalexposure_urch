#Dose Exposure--Gene Expression Scripts
#website with tutorial information: 
#https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#rich-visualization-and-reporting-of-results

#Madison Armstrong 
#Last Modified 9/29/25

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

design.int<-as.formula(~Treatment*urban)

dds.start.interaction <- DESeqDataSetFromTximport(data,
                                colData = metadata,
                                design = design.int) 

dds.start.interaction

#pre-filtering step: using same as Leslie
keep <- rowSums(counts(dds.start.interaction)) >= 10
dds.start.interaction<- dds.start.interaction[keep,]

# identify genes that pass expression cutoff
isexpr <- rowSums(fpm(dds.start.interaction) > 1) >= 0.5 * ncol(dds.start.interaction)

quantLog.start <- log2(fpm(dds.start.interaction)[isexpr, ] + 1)
# Check sample names alignment
all(colnames(quantLog.start) == rownames(metadata))  #should be true!!

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
design.int<-as.formula(~Treatment*urban)
#With the count matrix, cts, and the sample information, coldata, we can construct a DESeqDataSet:
dds.blast.interaction <- DESeqDataSetFromTximport(data.blast,
                                     colData = metadata.blast,
                                     design = design.int) 

dds.blast.interaction 

#pre-filtering step: using same as Leslie
keep <- rowSums(counts(dds.blast.interaction)) >= 10
dds.blast.interaction <- dds.blast.interaction[keep,]

#variance stabilizing transformation-- extracting transformed values
vsd.blast.interaction <- vst(dds.blast.interaction, blind=FALSE)

#### DEGs####
hist.blast.interaction <- DESeq(dds.blast.interaction)
# Report: Number of transcripts
num_transcripts.B <- nrow(hist.blast.interaction)
num_transcripts.B #34008
resultsNames(hist.blast.interaction) #interaction

#filter by padj=0.05 and then extract logfoldchange
resblast.int<-results(hist.blast.interaction) #gives readout of padj, log2FoldChange and other stuff for interaction
resblast.int #log2 fold change (MLE): TreatmentNP.urbanurban 

#####SAVE FOR TOPGO
write.csv(resblast.int,file="RNA_data/HEresults.blast.interaction.urb.csv", row.names=TRUE)

####Heat Map####
# Subset colData to only include samples from the relevant contrast-- pulls outliers
resblast_outliers.int <- resblast.int %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < 0.05)
dim(resblast_outliers.int) #22 outliers woo

# Remove genes with zero counts across all samples (to prevent vst() errors)
hist.blast.interaction <- hist.blast.interaction[rowSums(counts(hist.blast.interaction)) > 0, ]

# Perform Variance Stabilizing Transformation (VST)
hist.blast.interaction.vsd <- vst(hist.blast.interaction, blind=FALSE)

# Extract the transformed expression matrix for the top genes
resblast_top_genes_matrix_var <- assay(hist.blast.interaction.vsd)[resblast_outliers.int$gene, ]

df.b <- as.data.frame(colData(hist.blast.interaction.vsd)[,c("Treatment","Population")])
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

####PCAS####
pcaData_blast<-
  plotPCA(vsd.blast.interaction, intgroup=c("Treatment", "Population", "urban", "Experiment"),  pcsToUse=5:6, returnData=TRUE) #using ntop=500 top features by variance
percentVar <- round(100 * attr(pcaData_blast, "percentVar"))


ggplot(pcaData_blast, aes(x=PC5, y=PC6, shape=Experiment, color=urban)) +
  geom_point(size=4) +
  #stat_ellipse(aes(x = PC1,y=PC2,color=factor(Experiment))) +
  #facet_wrap(~Population)+
  xlab(paste0("PC5: ",percentVar[1],"% variance")) +
  ylab(paste0("PC6: ",percentVar[2],"% variance")) + 
  labs(title= "Blastula")+
  scale_color_manual(values=c("nonurban"="#99D1EC","urban"="#5C5649"))+
  theme(legend.position="right") + theme_box()

#linear models
blast1<-lm(PC5~Population+Treatment+ Experiment, data=pcaData_blast) 
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
design.int<-as.formula(~Treatment*urban)
#With the count matrix, cts, and the sample information, coldata, we can construct a DESeqDataSet:
dds.gast.interaction <- DESeqDataSetFromTximport(data.gast,
                                      colData = metadata.gast,
                                      design = design.int) 

dds.gast.interaction

#pre-filtering step: using same as Leslie
keep <- rowSums(counts(dds.gast.interaction)) >= 10
dds.gast.interaction <- dds.gast.interaction[keep,]

#variance stabilizing transformation-- extracting transformed values
vsd.gast.interaction <- vst(dds.gast.interaction, blind=FALSE)

#### DEGs####
hist.gast.interaction <- DESeq(dds.gast.interaction)
# Report: Number of transcripts
num_transcripts.G..interaction <- nrow(hist.gast.interaction)
resultsNames(hist.gast.interaction) #populations compared & treatment vs control compared


#filter by padj=0.05 and then extract logfoldchange
resgast.int<-results(hist.gast.interaction) #gives readout of padj, log2FoldChange and other stuff
resgast.int
#####SAVE FOR TOPGO
write.csv(resgast.int,file="RNA_data/HEresults.gast.interaction.urb.csv", row.names=TRUE)

####Heat Map####
# Subset colData to only include samples from the relevant contrast-- pulls outliers
resgast_outliers.int <- resgast.int %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < 0.05)
dim(resgast_outliers.int) #10 DEGs

# Remove genes with zero counts across all samples (to prevent vst() errors)
hist.gast.interaction <- hist.gast.interaction[rowSums(counts(hist.gast.interaction)) > 0, ]

# Perform Variance Stabilizing Transformation (VST)
hist.gast.interaction.vsd <- vst(hist.gast.interaction, blind=FALSE)

# Extract the transformed expression matrix for the top genes
resgast_top_genes_matrix_var.int <- assay(hist.gast.interaction.vsd)[resgast_outliers.int$gene, ]

df.g <- as.data.frame(colData(hist.gast.interaction.vsd)[,c("Treatment","Population")])
df.g$Treatment <- factor(df.g$Treatment, levels = c("C", "NP")) 
df.g$Population <- factor(df.g$Population, levels = c("Cab", "Wh", "Trp", "Twp")) 

##Get things in order!
gast_ordered_samples <- df.g %>%
  arrange(Population, Treatment) %>%
  rownames()

gast_ordered_matrix <- resgast_top_genes_matrix_var.int[, gast_ordered_samples]
df.g.ordered <- df.g[gast_ordered_samples, ]

annot_colors=list(
  Population=c("Cab"="#5C5649", "Wh"="#A69C8B", "Trp"="#0F6DD2", "Twp"="#99D1EC"),
  Treatment=c("NP" = "#C77DFF", "C" = "#C8E9A0"))

pheatmap(gast_ordered_matrix, cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df.g.ordered, 
         annotation_colors=annot_colors,
         scale="row",
         color = magma(100))

### PCAS####
pcaData_gast<-
  plotPCA(vsd.gast, intgroup=c("Treatment", "Population", "urban", "Pop_MP", "Experiment"),  pcsToUse=3:4, returnData=TRUE) #using ntop=500 top features by variance
percentVar <- round(100 * attr(pcaData_gast, "percentVar"))


ggplot(pcaData_gast, aes(x=PC3, y=PC4,  shape=Experiment, color=urban)) +
  geom_point(size=4) +
  #stat_ellipse() +
 # facet_wrap(~Experiment)+
  xlab(paste0("PC3: ",percentVar[1],"% variance")) +
  ylab(paste0("PC4: ",percentVar[2],"% variance")) + 
  labs(title= "Gastrula")+
  scale_color_manual(values=c("nonurban"="#99D1EC","urban"="#5C5649"))+
  theme(legend.position="right") + theme_box()

#linear models
gast1<-lm(PC4~urban+Treatment+Experiment, data=pcaData_gast) 
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
design.int<-as.formula(~Treatment*urban)
#With the count matrix, cts, and the sample information, coldata, we can construct a DESeqDataSet:
dds.plut.int <- DESeqDataSetFromTximport(data.plut,
                                     colData = metadata.plut,
                                     design =design.int) 

dds.plut.int 

#pre-filtering step: using same as Leslie
keep <- rowSums(counts(dds.plut.int)) >= 10
dds.plut.int <- dds.plut.int[keep,]

#variance stabilizing transformation-- extracting transformed values
vsd.plut.int <- vst(dds.plut.int, blind=FALSE)

#### DEGs####
hist.plut.int <- DESeq(dds.plut.int)
# Report: Number of transcripts
num_transcripts.P <- nrow(hist.plut.int)
num_transcripts.P # 37315 transcripts
resultsNames(hist.plut.int) #populations compared & treatment vs control compared

#filter by padj=0.05 and then extract logfoldchange
##population or urban comparison
resplut.int<-results(hist.plut.int) #gives readout of padj, log2FoldChange and other stuff
resplut.int
#####SAVE FOR TOPGO
write.csv(resplut.int,file="RNA_data/HEresults.plut.interaction.urb.csv", row.names=TRUE)

####Heat Map####
# Subset colData to only include samples from the relevant contrast-- pulls outliers
resplut_outliers.int <- resplut.int %>% #pulling outliers sig for int
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < 0.05)
dim(resplut_outliers.int) #12

# Remove genes with zero counts across all samples (to prevent vst() errors)
hist.plut.int <- hist.plut.int[rowSums(counts(hist.plut.int)) > 0, ]

# Perform Variance Stabilizing Transformation (VST)
hist.plut.vsd.int <- vst(hist.plut.int, blind=FALSE)

##look at urban comparison
resplut.urb <- results(hist.plut.int, (c(contrast="urban", "urban", "nonurban")))

resplut_outliers.urb <- resplut.urb %>% #pulling outliers sig for int
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < 0.05)

dim(resplut_outliers.urb) #29 DEGs

# Extract the transformed expression matrix for the top genes
resplut_top_genes_matrix_var.urb <- assay(hist.plut.vsd.int)[resplut_outliers.urb$gene, ]

df.p <- as.data.frame(colData(hist.plut.vsd.int)[,c("Treatment","Population")])
df.p$Treatment <- factor(df.p$Treatment, levels = c("C", "NP")) 
df.p$Population <- factor(df.p$Population, levels = c("Cab", "Wh", "Trp", "Twp")) 

##Get things in order!
plut_ordered_samples <- df.p %>%
  arrange(Population, Treatment) %>%
  rownames()

plut_ordered_matrix <- resplut_top_genes_matrix_var.urb[, plut_ordered_samples]
df.p.ordered <- df.p[plut_ordered_samples, ]

annot_colors=list(
  Population=c("Cab"="#5C5649", "Wh"="#A69C8B", "Trp"="#0F6DD2", "Twp"="#99D1EC"),
  Treatment=c("NP" = "#C77DFF", "C" = "#C8E9A0"))

pheatmap(plut_ordered_matrix, cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df.p.ordered, 
         annotation_colors=annot_colors,
         scale="row",
         color = magma(100))
