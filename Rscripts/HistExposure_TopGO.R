##Historical Exposure
#TopGO script
#Madison Armstrong 
#Last Modified 7/14/2025


setwd("~/Desktop/purp dev data/2024_WBML")

#Libraries
#if(!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#BiocManager::install(c("biomaRt","topGO", "stringr"))

library(biomaRt)
##https://www.ensembl.org/info/data/biomart/index.html
library(topGO)
library(tidyverse)
library(stringr)
library(dplyr)

#read in files####
#read in .csv from DEG script 

##change ".pop" to ".urb". for urban/nonurban comparisons
##look at treatment effects for each model by including .poptreat or .urbtreat

###blastula####
blast.results<-read.csv("RNA_data/HEresults.blast.urbtreat.csv") #for DEGs

parname <- "padj"  # Corrected column name
# Define significance threshold (typically 0.05 for DEGs)
significance_threshold <- 0.05  
# Create a new column to indicate significant genes (1) and non-significant genes (0)
blast.results$significant <- ifelse(blast.results[[parname]] < significance_threshold, 1, 0)
parsig.blast = which(blast.results$significant==1)
length(parsig.blast)
blast.results.sig <- blast.results[!is.na(blast.results$significant) & blast.results$significant != 0, ] #get rid of NAs & 0s
#write.csv(blast.results.sig,"RNA_data/HE.blast.results.sig.urbtreat.full.csv")

###gastrula####
gast.results<-read.csv("RNA_data/HEresults.gast.urbtreat.csv")
parname <- "padj"  # Corrected column name
# Define significance threshold (typically 0.05 for DEGs)
significance_threshold <- 0.05  
# Create a new column to indicate significant genes (1) and non-significant genes (0)
gast.results$significant <- ifelse(gast.results[[parname]] < significance_threshold, 1, 0)
parsig.gast = which(gast.results$significant==1)
length(parsig.gast)
gast.results.sig <- gast.results[!is.na(gast.results$significant) & gast.results$significant != 0, ] #get rid of NAs & 0s
#write.csv(gast.results.sig,"RNA_data/HE.gast.results.sig.urbtreat.full.csv")

###pluteus####
plut.results<-read.csv("RNA_data/HEresults.plut.urbtreat.csv")
parname <- "padj"  # Corrected column name
# Define significance threshold (typically 0.05 for DEGs)
significance_threshold <- 0.05  
# Create a new column to indicate significant genes (1) and non-significant genes (0)
plut.results$significant <- ifelse(plut.results[[parname]] < significance_threshold, 1, 0)
parsig.plut = which(plut.results$significant==1)
length(parsig.plut)
plut.results.sig <- plut.results[!is.na(plut.results$significant) & plut.results$significant != 0, ] #get rid of NAs & 0s
#write.csv(plut.results.sig,"RNA_data/HE.plut.results.sig.urbtreat.full.csv")

#Get annotations from biomaRt####
metazoa_mart <- useMart(biomart = 'metazoa_mart',
                        host = "https://metazoa.ensembl.org/")
Spurp_mart <- useDataset(mart = metazoa_mart, dataset = "spurpuratus_eg_gene")
listAttributes(Spurp_mart)
Spurp_genes <- getBM(attributes = c("ensembl_gene_id", "description", "chromosome_name", "start_position", "end_position",
                                    "external_gene_name", "gene_biotype", "go_id", "name_1006", "definition_1006", 
                                    "namespace_1003","ensembl_transcript_id"),
                     filters = c("with_geneid"),
                     values = TRUE,
                     mart = Spurp_mart)
head(Spurp_genes)
#write.csv(Spurp_genes,"RNA_data/Spurp_genes.csv")

### Create GO Key for all transcripts in the whole genome
uniqTranscripts <- unique(Spurp_genes$ensembl_transcript_id)  # Get unique transcript IDs
gos <- data.frame(ensembl_id = uniqTranscripts, GO=NA)  # Initialize data frame

for (t in uniqTranscripts) {
  sub <- subset(Spurp_genes, ensembl_transcript_id == t)  # Subset data for each transcript
  gostring <- paste(sub$go_id, collapse=",")  # Concatenate GO terms
  gos$GO[match(t, gos$ensembl_id)] <- gostring # Assign GO terms to the transcript
}
#write.table(gos,"RNA_data/Spurp_GOmap.txt",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")

#use S.purp dataset to pull gene_info to fill in gene info
#ensembl_gene_id + description
Spurp_genes_sub<- Spurp_genes[c("ensembl_gene_id", "description", "ensembl_transcript_id")]


#pull out versions 
blast.results.sig.sep<-read.csv("RNA_data/HE.blast.results.sig.urbtreat.full.csv")
#pull out version numbers
blast.results.sig.sep$X<- sub("\\..*", "", blast.results.sig.sep$X)

blast.genes<-left_join(blast.results.sig.sep, Spurp_genes_sub, by=join_by(X==ensembl_transcript_id))
blast.genes<- blast.genes %>%
  distinct(X, .keep_all = TRUE)
#write.csv(blast.genes,"RNA_data/HE.blast.results.sig.urbtreat.full.genes.csv")

#pull out versions
gast.results.sig.sep<-read.csv("RNA_data/HE.gast.results.sig.urbtreat.full.csv")
#pull out version numbers
gast.results.sig.sep$X<- sub("\\..*", "", gast.results.sig.sep$X)

gast.genes<-left_join(gast.results.sig.sep, Spurp_genes_sub, by=join_by(X==ensembl_transcript_id))
gast.genes <- gast.genes %>%
  distinct(X, .keep_all = TRUE)
#write.csv(gast.genes,"RNA_data/HE.gast.results.sig.urbtreat.full.genes.csv")

#pull out versions 
plut.results.sig.sep<-read.csv("RNA_data/HE.plut.results.sig.urbtreat.full.csv")
#pull out version numbers
plut.results.sig.sep$X<- sub("\\..*", "", plut.results.sig.sep$X)

plut.genes<-left_join(plut.results.sig.sep, Spurp_genes_sub, by=join_by(X==ensembl_transcript_id))
plut.genes <- plut.genes %>%
  distinct(X, .keep_all = TRUE)

#write.csv(plut.genes,"RNA_data/HE.plut.results.sig.urbtreat.full.genes.csv")

#compare datasets
shared.genes.poptreat <- blast.genes %>%
  inner_join(gast.genes, by = "X") %>%
  inner_join(plut.genes, by = "X")

#remove unwanted columns, like duplicate columns
#shared.genes.poptreat<-shared.genes.poptreat %>% select(-c(X.1.x, significant.x, significant.y, ensembl_gene_id.y, description.y, ensembl_gene_id.x, description.x))
#write.csv(shared.genes.poptreat, "RNA_data/HE.shared.genes.poptreat.csv")

#GO enrichment for DEGs####
geneID2GO <- readMappings(file="RNA_data/Spurp_GOmap.txt", sep="\t", IDsep=",")
# Read in set of all genes and make 'background' gene set
all.genes.blast <- blast.results$X #all blastula genes
all.genes.gast <- gast.results$X #all gastrula genes
all.genes.plut <- plut.results$X #all pluteus genes

####Blastula####
candgenes.blast <- blast.results$X[parsig.blast] #for regular DEG input
myIG.blast <-factor(as.numeric(all.genes.blast%in%candgenes.blast))
myIG.blast <- as.numeric(as.character(myIG.blast)) #make sure it is 0s & 1s
names(myIG.blast) <- all.genes.blast
names(myIG.blast) <- sub("\\..*", "", names(myIG.blast))  # Remove everything after the first dot
topGO_selectFun <- function(x) return(x == 1)  # Select genes labeled as "1"

parname = "Blastula"
##now separate out by category
##BP
GOdata.blast <- new("topGOdata",ontology="BP",allGenes=myIG.blast, nodeSize=10,
                   annotationFun=annFUN.gene2GO,gene2GO=geneID2GO,
                   geneSelectionFun=topGO_selectFun)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
resultFisher <- getSigGroups(GOdata.blast,test.stat)
res <- GenTable(GOdata.blast,classic=resultFisher,topNodes=length(resultFisher@score),numChar=100)
filt <- res[res$classic<0.05 & res$Significant>2,]
write.csv(filt,paste("RNA_data/",parname,".GO.urbtreat_BP.HE.csv",sep=""))

##MF
GOdata <- new("topGOdata",ontology="MF",allGenes=myIG.blast, nodeSize=10,
              annotationFun=annFUN.gene2GO,gene2GO=geneID2GO,  geneSelectionFun=topGO_selectFun)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
resultFisher <- getSigGroups(GOdata,test.stat)
res <- GenTable(GOdata,classic=resultFisher,topNodes=length(resultFisher@score),numChar=100)
filt <- res[res$classic<0.05 & res$Significant>2,]
write.csv(filt,paste("RNA_data/",parname,".GO.urbtreat_MF.HE.csv",sep=""))

##CC
GOdata <- new("topGOdata",ontology="CC",allGenes=myIG.blast, nodeSize=10,
              annotationFun=annFUN.gene2GO,gene2GO=geneID2GO,  geneSelectionFun=topGO_selectFun)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
resultFisher <- getSigGroups(GOdata,test.stat)
res <- GenTable(GOdata,classic=resultFisher,topNodes=length(resultFisher@score),numChar=100)
filt <- res[res$classic<0.05 & res$Significant>2,]
write.csv(filt,paste("RNA_data/",parname,".GO.urbtreat_CC.HE.csv",sep=""))


####Gastrula####
candgenes.gast <- gast.results$X[parsig.gast] #
myIG.gast <-factor(as.numeric(all.genes.gast%in%candgenes.gast))
myIG.gast <- as.numeric(as.character(myIG.gast)) #make sure it is 0s & 1s
names(myIG.gast) <- all.genes.gast
names(myIG.gast) <- sub("\\..*", "", names(myIG.gast))  # Remove everything after the first dot
topGO_selectFun <- function(x) return(x == 1)  # Select genes labeled as "1"

parname = "Gastrula"

##now separate out by category
##BP
GOdata.gast <- new("topGOdata",ontology="BP",allGenes=myIG.gast, nodeSize=10,
                    annotationFun=annFUN.gene2GO,gene2GO=geneID2GO,
                    geneSelectionFun=topGO_selectFun)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
resultFisher <- getSigGroups(GOdata.gast,test.stat)
res <- GenTable(GOdata.gast,classic=resultFisher,topNodes=length(resultFisher@score),numChar=100)
filt <- res[res$classic<0.05 & res$Significant>2,]
write.csv(filt,paste("RNA_data/",parname,".GO.poptreat_BP.HE.csv",sep=""))

##MF
GOdata <- new("topGOdata",ontology="MF",allGenes=myIG.gast, nodeSize=10,
              annotationFun=annFUN.gene2GO,gene2GO=geneID2GO,  geneSelectionFun=topGO_selectFun)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
resultFisher <- getSigGroups(GOdata,test.stat)
res <- GenTable(GOdata,classic=resultFisher,topNodes=length(resultFisher@score),numChar=100)
filt <- res[res$classic<0.05 & res$Significant>2,]
write.csv(filt,paste("RNA_data/",parname,".GO.poptreat_MF.HE.csv",sep=""))

##CC
GOdata <- new("topGOdata",ontology="CC",allGenes=myIG.gast, nodeSize=10,
              annotationFun=annFUN.gene2GO,gene2GO=geneID2GO,  geneSelectionFun=topGO_selectFun)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
resultFisher <- getSigGroups(GOdata,test.stat)
res <- GenTable(GOdata,classic=resultFisher,topNodes=length(resultFisher@score),numChar=100)
filt <- res[res$classic<0.05 & res$Significant>2,]
write.csv(filt,paste("RNA_data/",parname,".GO.poptreat_CC.HE.csv",sep=""))

####Pluteus####
candgenes.plut <- plut.results$X[parsig.plut] #
myIG.plut <-factor(as.numeric(all.genes.plut%in%candgenes.plut))
myIG.plut <- as.numeric(as.character(myIG.plut)) #make sure it is 0s & 1s
names(myIG.plut) <- all.genes.plut
names(myIG.plut) <- sub("\\..*", "", names(myIG.plut))  # Remove everything after the first dot
topGO_selectFun <- function(x) return(x == 1)  # Select genes labeled as "1"

parname = "Pluteus"

##now separate out by category
##BP
GOdata.plut <- new("topGOdata",ontology="BP",allGenes=myIG.plut, nodeSize=10,
                   annotationFun=annFUN.gene2GO,gene2GO=geneID2GO,
                   geneSelectionFun=topGO_selectFun)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
resultFisher <- getSigGroups(GOdata.plut,test.stat)
res <- GenTable(GOdata.plut,classic=resultFisher,topNodes=length(resultFisher@score),numChar=100)
filt <- res[res$classic<0.05 & res$Significant>2,]
write.csv(filt,paste("RNA_data/",parname,".GO.poptreat_BP.HE.csv",sep=""))

##MF
GOdata <- new("topGOdata",ontology="MF",allGenes=myIG.plut, nodeSize=10,
              annotationFun=annFUN.gene2GO,gene2GO=geneID2GO,  geneSelectionFun=topGO_selectFun)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
resultFisher <- getSigGroups(GOdata,test.stat)
res <- GenTable(GOdata,classic=resultFisher,topNodes=length(resultFisher@score),numChar=100)
filt <- res[res$classic<0.05 & res$Significant>2,]
write.csv(filt,paste("RNA_data/",parname,".GO.poptreat_MF.HE.csv",sep=""))

##CC
GOdata <- new("topGOdata",ontology="CC",allGenes=myIG.plut, nodeSize=10,
              annotationFun=annFUN.gene2GO,gene2GO=geneID2GO,  geneSelectionFun=topGO_selectFun)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
resultFisher <- getSigGroups(GOdata,test.stat)
res <- GenTable(GOdata,classic=resultFisher,topNodes=length(resultFisher@score),numChar=100)
filt <- res[res$classic<0.05 & res$Significant>2,]
write.csv(filt,paste("RNA_data/",parname,".GO.poptreat_CC.HE.csv",sep=""))

##for WGCNA 

# Read in set of all genes and make 'background' gene set but remove versions for WGCNA analyses
all.genes.blast.w <- blast.results$X #all blastula genes
all.genes.blast.w<- sub("\\..*", "", all.genes.blast.w)  # Remove everything after the first dot

all.genes.gast.w <- gast.results$X #all gastrula genes
all.genes.gast.w<- sub("\\..*", "", all.genes.gast.w)

all.genes.plut.w <- plut.results$X #all pluteus genes
all.genes.plut.w<- sub("\\..*", "", all.genes.plut.w)

blast.WGCNA<-read.csv("RNA_data/HE.blast.WGCNA.sigmods.csv") #for DEGs with WGCNA modules
gast.WGCNA<-read.csv("RNA_data/HE.gast.WGCNA.sigmods.csv")
plut.WGCNA<-read.csv("RNA_data/HE.plut.WGCNA.sigmods.csv")


####WGCNA Blastula####
Sg = which(blast.WGCNA$colors=="yellowgreen")
candgenes.blast<-blast.WGCNA$X[Sg] #for WGCNA input
myIG.blast <-factor(as.numeric(all.genes.blast.w%in%candgenes.blast))
myIG.blast <- as.numeric(as.character(myIG.blast)) #make sure it is 0s & 1s
names(myIG.blast) <- all.genes.blast.w
topGO_selectFun <- function(x) return(x == 1)  # Select genes labeled as "1"

parname = "Blastula-yg"
##now separate out by category
##BP
GOdata.blast <- new("topGOdata",ontology="BP",allGenes=myIG.blast, nodeSize=10,
                    annotationFun=annFUN.gene2GO,gene2GO=geneID2GO,
                    geneSelectionFun=topGO_selectFun)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
resultFisher <- getSigGroups(GOdata.blast,test.stat)
res <- GenTable(GOdata.blast,classic=resultFisher,topNodes=length(resultFisher@score),numChar=100)
filt <- res[res$classic<0.05 & res$Significant>2,]
write.csv(filt,paste("RNA_data/",parname,".GO.WGCNAurbtreat_BP.HE.csv",sep=""))

##MF
GOdata <- new("topGOdata",ontology="MF",allGenes=myIG.blast, nodeSize=10,
              annotationFun=annFUN.gene2GO,gene2GO=geneID2GO,  geneSelectionFun=topGO_selectFun)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
resultFisher <- getSigGroups(GOdata,test.stat)
res <- GenTable(GOdata,classic=resultFisher,topNodes=length(resultFisher@score),numChar=100)
filt <- res[res$classic<0.05 & res$Significant>2,]
write.csv(filt,paste("RNA_data/",parname,".GO.WGCNAurbtreat_MF.HE.csv",sep=""))

##CC
GOdata <- new("topGOdata",ontology="CC",allGenes=myIG.blast, nodeSize=10,
              annotationFun=annFUN.gene2GO,gene2GO=geneID2GO,  geneSelectionFun=topGO_selectFun)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
resultFisher <- getSigGroups(GOdata,test.stat)
res <- GenTable(GOdata,classic=resultFisher,topNodes=length(resultFisher@score),numChar=100)
filt <- res[res$classic<0.05 & res$Significant>2,]
write.csv(filt,paste("RNA_data/",parname,".GO.WGCNAurbtreat_CC.HE.csv",sep=""))


####WGCNA Gastrula####
Sg = which(gast.WGCNA$colors=="cyan")
candgenes.gast <- gast.WGCNA$X[Sg] #
myIG.gast <-factor(as.numeric(all.genes.gast.w%in%candgenes.gast))
myIG.gast <- as.numeric(as.character(myIG.gast)) #make sure it is 0s & 1s
names(myIG.gast) <- all.genes.gast.w
topGO_selectFun <- function(x) return(x == 1)  # Select genes labeled as "1"

parname = "Gastrula-cy"

##now separate out by category
##BP
GOdata.gast <- new("topGOdata",ontology="BP",allGenes=myIG.gast, nodeSize=10,
                   annotationFun=annFUN.gene2GO,gene2GO=geneID2GO,
                   geneSelectionFun=topGO_selectFun)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
resultFisher <- getSigGroups(GOdata.gast,test.stat)
res <- GenTable(GOdata.gast,classic=resultFisher,topNodes=length(resultFisher@score),numChar=100)
filt <- res[res$classic<0.05 & res$Significant>2,]
write.csv(filt,paste("RNA_data/",parname,".GO.WGCNAurbtreat_BP.HE.csv",sep=""))

##MF
GOdata <- new("topGOdata",ontology="MF",allGenes=myIG.gast, nodeSize=10,
              annotationFun=annFUN.gene2GO,gene2GO=geneID2GO,  geneSelectionFun=topGO_selectFun)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
resultFisher <- getSigGroups(GOdata,test.stat)
res <- GenTable(GOdata,classic=resultFisher,topNodes=length(resultFisher@score),numChar=100)
filt <- res[res$classic<0.05 & res$Significant>2,]
write.csv(filt,paste("RNA_data/",parname,".GO.WGCNAurbtreat_MF.HE.csv",sep=""))

##CC
GOdata <- new("topGOdata",ontology="CC",allGenes=myIG.gast, nodeSize=10,
              annotationFun=annFUN.gene2GO,gene2GO=geneID2GO,  geneSelectionFun=topGO_selectFun)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
resultFisher <- getSigGroups(GOdata,test.stat)
res <- GenTable(GOdata,classic=resultFisher,topNodes=length(resultFisher@score),numChar=100)
filt <- res[res$classic<0.05 & res$Significant>2,]
write.csv(filt,paste("RNA_data/",parname,".GO.WGCNAurbtreat_CC.HE.csv",sep=""))

#####2nd module####
Sg = which(gast.WGCNA$colors=="turquoise")
candgenes.gast <- gast.WGCNA$X[Sg] #
myIG.gast <-factor(as.numeric(all.genes.gast.w%in%candgenes.gast))
myIG.gast <- as.numeric(as.character(myIG.gast)) #make sure it is 0s & 1s
names(myIG.gast) <- all.genes.gast.w
topGO_selectFun <- function(x) return(x == 1)  # Select genes labeled as "1"

parname = "Gastrula-tq"

##now separate out by category
##BP
GOdata.gast <- new("topGOdata",ontology="BP",allGenes=myIG.gast, nodeSize=10,
                   annotationFun=annFUN.gene2GO,gene2GO=geneID2GO,
                   geneSelectionFun=topGO_selectFun)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
resultFisher <- getSigGroups(GOdata.gast,test.stat)
res <- GenTable(GOdata.gast,classic=resultFisher,topNodes=length(resultFisher@score),numChar=100)
filt <- res[res$classic<0.05 & res$Significant>2,]
write.csv(filt,paste("RNA_data/",parname,".GO.WGCNAurbtreat_BP.HE.csv",sep=""))

##MF
GOdata <- new("topGOdata",ontology="MF",allGenes=myIG.gast, nodeSize=10,
              annotationFun=annFUN.gene2GO,gene2GO=geneID2GO,  geneSelectionFun=topGO_selectFun)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
resultFisher <- getSigGroups(GOdata,test.stat)
res <- GenTable(GOdata,classic=resultFisher,topNodes=length(resultFisher@score),numChar=100)
filt <- res[res$classic<0.05 & res$Significant>2,]
write.csv(filt,paste("RNA_data/",parname,".GO.WGCNAurbtreat_MF.HE.csv",sep=""))

##CC
GOdata <- new("topGOdata",ontology="CC",allGenes=myIG.gast, nodeSize=10,
              annotationFun=annFUN.gene2GO,gene2GO=geneID2GO,  geneSelectionFun=topGO_selectFun)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
resultFisher <- getSigGroups(GOdata,test.stat)
res <- GenTable(GOdata,classic=resultFisher,topNodes=length(resultFisher@score),numChar=100)
filt <- res[res$classic<0.05 & res$Significant>2,]
write.csv(filt,paste("RNA_data/",parname,".GO.WGCNAurbtreat_CC.HE.csv",sep=""))

####WGCNA Pluteus####
Sg = which(plut.WGCNA$colors=="paleturquoise")
candgenes.plut<-plut.WGCNA$X[Sg] #for WGCNA input
myIG.plut <-factor(as.numeric(all.genes.plut.w%in%candgenes.plut))
myIG.plut <- as.numeric(as.character(myIG.plut)) #make sure it is 0s & 1s
names(myIG.plut) <- all.genes.plut.w
topGO_selectFun <- function(x) return(x == 1)  # Select genes labeled as "1"

parname = "Pluteus-ptq"

##now separate out by category
##BP
GOdata.plut <- new("topGOdata",ontology="BP",allGenes=myIG.plut, nodeSize=10,
                   annotationFun=annFUN.gene2GO,gene2GO=geneID2GO,
                   geneSelectionFun=topGO_selectFun)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
resultFisher <- getSigGroups(GOdata.plut,test.stat)
res <- GenTable(GOdata.plut,classic=resultFisher,topNodes=length(resultFisher@score),numChar=100)
filt <- res[res$classic<0.05 & res$Significant>2,]
write.csv(filt,paste("RNA_data/",parname,".GO.WGCNAurbtreat_BP.HE.csv",sep=""))

##MF
GOdata <- new("topGOdata",ontology="MF",allGenes=myIG.plut, nodeSize=10,
              annotationFun=annFUN.gene2GO,gene2GO=geneID2GO,  geneSelectionFun=topGO_selectFun)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
resultFisher <- getSigGroups(GOdata,test.stat)
res <- GenTable(GOdata,classic=resultFisher,topNodes=length(resultFisher@score),numChar=100)
filt <- res[res$classic<0.05 & res$Significant>2,]
write.csv(filt,paste("RNA_data/",parname,".GO.WGCNAurbtreat_MF.HE.csv",sep=""))

##CC
GOdata <- new("topGOdata",ontology="CC",allGenes=myIG.plut, nodeSize=10,
              annotationFun=annFUN.gene2GO,gene2GO=geneID2GO,  geneSelectionFun=topGO_selectFun)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
resultFisher <- getSigGroups(GOdata,test.stat)
res <- GenTable(GOdata,classic=resultFisher,topNodes=length(resultFisher@score),numChar=100)
filt <- res[res$classic<0.05 & res$Significant>2,]
write.csv(filt,paste("RNA_data/",parname,".GO.WGCNAurbtreat_CC.HE.csv",sep=""))

#####2nd module####
Sg = which(plut.WGCNA$colors=="turquoise")
candgenes.plut<-plut.WGCNA$X[Sg] #for WGCNA input
myIG.plut <-factor(as.numeric(all.genes.plut.w%in%candgenes.plut))
myIG.plut <- as.numeric(as.character(myIG.plut)) #make sure it is 0s & 1s
names(myIG.plut) <- all.genes.plut.w
topGO_selectFun <- function(x) return(x == 1)  # Select genes labeled as "1"

parname = "Pluteus-tq"

##now separate out by category
##BP
GOdata.plut <- new("topGOdata",ontology="BP",allGenes=myIG.plut, nodeSize=10,
                   annotationFun=annFUN.gene2GO,gene2GO=geneID2GO,
                   geneSelectionFun=topGO_selectFun)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
resultFisher <- getSigGroups(GOdata.plut,test.stat)
res <- GenTable(GOdata.plut,classic=resultFisher,topNodes=length(resultFisher@score),numChar=100)
filt <- res[res$classic<0.05 & res$Significant>2,]
write.csv(filt,paste("RNA_data/",parname,".GO.WGCNAurbtreat_BP.HE.csv",sep=""))

##MF
GOdata <- new("topGOdata",ontology="MF",allGenes=myIG.plut, nodeSize=10,
              annotationFun=annFUN.gene2GO,gene2GO=geneID2GO,  geneSelectionFun=topGO_selectFun)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
resultFisher <- getSigGroups(GOdata,test.stat)
res <- GenTable(GOdata,classic=resultFisher,topNodes=length(resultFisher@score),numChar=100)
filt <- res[res$classic<0.05 & res$Significant>2,]
write.csv(filt,paste("RNA_data/",parname,".GO.WGCNAurbtreat_MF.HE.csv",sep=""))

##CC
GOdata <- new("topGOdata",ontology="CC",allGenes=myIG.plut, nodeSize=10,
              annotationFun=annFUN.gene2GO,gene2GO=geneID2GO,  geneSelectionFun=topGO_selectFun)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
resultFisher <- getSigGroups(GOdata,test.stat)
res <- GenTable(GOdata,classic=resultFisher,topNodes=length(resultFisher@score),numChar=100)
filt <- res[res$classic<0.05 & res$Significant>2,]
write.csv(filt,paste("RNA_data/",parname,".GO.WGCNAurbtreat_CC.HE.csv",sep=""))