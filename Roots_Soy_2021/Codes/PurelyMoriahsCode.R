---
  title: "PLB 847 Group Project - LTER Soybean Roots and their Fungal Community"
author: "Moriah Young"
date: "2022-11-29"
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
#library(microViz)
# Three tables are needed for the phyloseq object
# OTU
# Taxonomy
# Samples
# Sequences (as far as I can tell, this isn't necessary???)
otu_mat <- read.delim("otu_table_ITS_UPARSE_R1.txt") # OTUs - fyi, mat = matrix
tax_mat <- read.delim("constax_taxonomy.txt") # taxonomy
samples_df <- read.delim("root_fun_map_soybean_2021_1.txt") # metadata
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

physeq_object_roots <- phyloseq(OTU, TAX, samples)
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
samplesover1000_all <- subset_samples(physeq_object_roots, sample_sums(physeq_object_roots) > 1000)
#any(taxa_sums(samplesover1000_all) == 0)
set.seed(81)
# rarefy samples
rarefy_samplesover1000_all <- rarefy_even_depth(samplesover1000_all, rngseed= 81, sample.size = min(sample_sums(samplesover1000_all)))
roots <- rarefy_samplesover1000_all
```

# this chunk of code is the closest I got to what we want. I need to figure out how to get an "Other" category.
# https://github.com/joey711/phyloseq/issues/901 - this might help
```{r}
# Genus Level
# this works
# I got this code from someone online so I'm not 100% sure what all the lines of code do
top30 <- names(sort(taxa_sums(roots), decreasing=TRUE)[1:25])
top30 #shows 30 results
dat.aglo = tax_glom(roots, taxrank = "Genus")
dat.trans = transform_sample_counts(dat.aglo, function(x) x/sum(x)) #calculating abundances
prune.dat.two = prune_taxa(top30, dat.trans)
dat.dataframe = psmelt(prune.dat.two)
dat.agr = aggregate(Abundance~bar_label+Genus+Management+Description, data=dat.dataframe, FUN=mean) # calculating the mean abundances by management type and growth stage (which is bar_label column from the "root_fun_map_soybean_2021_1.txt" file)
ggplot(dat.agr, 
       aes(x=Description, y=Abundance, fill=Genus)) + 
  geom_bar(stat="identity",  position="fill") + 
  facet_grid(~Management, scale="free") # this groups the bars by management type