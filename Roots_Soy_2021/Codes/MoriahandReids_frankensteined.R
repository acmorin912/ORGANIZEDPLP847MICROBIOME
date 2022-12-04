---
  title: "PLB 847 Group Project - LTER Soybean Roots and their Fungal Community"
author: "Moriah Young" "Ashlynn Morin" "Nico"
date: "2022-12-4"
output: pdf_document
---
  
  ```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{r}
rm(list = ls()) # clear working environment
# set working directory - don't need this, working directory is git repository
#setwd("/Users/moriahyoung/Desktop/PLB 847/PLB 847 Group Project/")
# load packages
library(tidyverse)
library(stringr)
library(plotrix) #std.error
library(phyloseq) 
library(tibble)
library(Biostrings)
#library(microViz)
# Three tables are needed for the phyloseq object
# OTU
# Taxonomy
# Samples
# Sequences (as far as I can tell, this isn't necessary???)
#Ashlynn Added sequences here just in case
otu_mat <- read.delim("otu_table_ITS_UPARSE_R1.txt") # OTUs - fyi, mat = matrix
tax_mat <- read.delim("constax_taxonomy.txt") # taxonomy
samples_df <- read.delim("root_fun_map_soybean_2021_1.txt") # metadata
sample_sequences<- readDNAStringSet("otus_R1.fasta", format = "fasta", seek.first.rec = TRUE, use.names = TRUE)
# change column name in "otu_mat" dataframe to match the "taxa" dataframe column
colnames(otu_mat)[colnames(otu_mat) == "X.OTU.ID"] <- "OTU_ID" 
# make it so the first column becomes the row names
# this is setting it up for the phyloseq object
otu_mat <- otu_mat %>%
  tibble::column_to_rownames("OTU_ID") 
tax_mat <- tax_mat %>% 
  tibble::column_to_rownames("OTU_ID") 
# remove the _1 that comes after every taxa name
tax_mat[] <- lapply(tax_mat, function(x) sub("_1", "", x, fixed = TRUE))
colnames(samples_df)[colnames(samples_df) == "X.SampleID"] <- "SampleID" # get rid of that X. that shows up
samples_df <- samples_df %>% 
  tibble::column_to_rownames("SampleID") 
#samples_df2 <- samples_df[-c("NC1root", "PC1root","NC1root"),]
# turn data frames into matrices
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)

physeq_object_roots <- phyloseq(OTU, TAX, samples, sample_sequences)
physeq_object_roots
# check out the phyloseq object
sample_names(physeq_object_roots)
rank_names(physeq_object_roots)
sample_variables(physeq_object_roots)
physeq_object_roots <- subset_samples(physeq_object_roots, Sample_or_Control =="True Sample")
sort(unique(as.data.frame(sample_data(physeq_object_roots))$Sample_or_Control)) # get rid of non true samples
sample_names(physeq_object_roots) # success
tax_table(physeq_object_roots)[tax_table(physeq_object_roots)==""]<- NA
# checking the phyloseq object
sort(unique(as.data.frame(tax_table(physeq_object_roots))$Kingdom)) # not everything is Fungi
# remove non fungi 
physeq_object_roots <- subset_taxa(physeq_object_roots, Kingdom!="Anthophyta")
physeq_object_roots <- subset_taxa(physeq_object_roots, Kingdom!="Alveolata")
physeq_object_roots <- subset_taxa(physeq_object_roots, Kingdom!="Ichthyosporia")
physeq_object_roots <- subset_taxa(physeq_object_roots, Kingdom!="Protista")
physeq_object_roots <- subset_taxa(physeq_object_roots, Kingdom!="Metazoa")
physeq_object_roots <- subset_taxa(physeq_object_roots, Kingdom!="Rhizaria")
physeq_object_roots <- subset_taxa(physeq_object_roots, Kingdom!="Viridiplantae")
sort(unique(as.data.frame(tax_table(physeq_object_roots))$Kingdom)) # now we're good
```

```{r}

##decontam from Nico code
count(as.data.frame(as.matrix(sample_data(
  physeq_object_roots))), bar_label, BarcodePlate)

table(physeq_object_roots@sam_data$BarcodePlate)
# adding column for controls from nico code
physeq_object_roots@sam_data$is.neg <- physeq@sam_data$Description == "Control"
head(sample_data(physeq))

contaminants <-
  isContaminant(
    physeq,
    method = "prevalence",
    batch = "BarcodePlate",
    neg = "is.neg",
    threshold = 0.4)

table(contaminants$contaminant) 


#Over 1000 reads after decontam Moriahs
samplesover1000_all <- subset_samples(physeq_object_roots, sample_sums(physeq_object_roots) > 1000)


##rarify, needs updated for Moriah's name but with Nico's code
set.seed(2018)
sort(colSums(otu_table(samplesover1000_all)))

physeq_rare = rarefy_even_depth(samplesover1000_all, sample.size = 8000, 
                                rngseed = FALSE, replace = TRUE, 
                                trimOTUs = TRUE, verbose = TRUE)

otu_table(physeq_rare) <- otu_table(physeq_rare)[which(rowSums(otu_table(physeq_rare)) >= 1),]
physeq_rare

colSums(otu_table(physeq_rare))
any(taxa_sums(physeq_rare) == 0)

physeq_rare

otu_table(physeq_rare)[otu_table(physeq_rare) <= 4] <- 0
physeq_rare
sample_data(physeq_rare)
# removes any OTUs that has less than 10 total reads across all samples
otu_table(physeq_rare) <- otu_table(physeq_rare)[which(rowSums(otu_table(physeq_rare)) >= 10),]
sample_data(physeq_rare)
###barplots by fungal genus---------------------------
library(data.table)
library(dplyr)
library(ggplot2)
#roots
root_barplots <- merge_samples(physeq_rare, "bar_label")
sample_data(root_barplots)
sample_data(root_barplots)$Sample <- factor(sample_data(root_barplots)$bar_label,
                                            levels = c("Conventional V2", "Conventional R2", "Conventional R6", "No-Till V2", "No-Till R2", "No-Till R6","Organic V2", "Organic R2", "Organic R6"))
```

#Phylum Level
```{r}
root_barplots <- root_barplots %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate at Genus level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format                       # Filter out low abundance taxa
  arrange(Genus)           # Sort data frame alphabetically by Genus
root_barplots
```

#Genus Level
```{r}
root_barplots <- root_barplots %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate at Genus level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format                       # Filter out low abundance taxa
  arrange(Genus)           # Sort data frame alphabetically by Genus
root_barplots
dat_root_bp <- data.table(root_barplots)
dat_root_bp[(Abundance <= 0.04), Genus := "Other"]
```


```{r}
bar_ITS_root= ggplot(dat_root_bp, aes(x = Sample, y = Abundance, fill = Genus)) + 
  #facet_wrap(~Stage, strip.position = "bottom") +
  theme(axis.text.x = element_text(angle = 90))+
  geom_bar(stat = "identity")
plot(bar_ITS_root)