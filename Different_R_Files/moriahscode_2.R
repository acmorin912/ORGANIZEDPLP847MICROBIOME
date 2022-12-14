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
```

```{r}
# this didn't do what I wanted it to do - another attempt from code I found online
top20otus = names(sort(taxa_sums(roots), TRUE)[1:20])
taxtab20 = cbind(tax_table(roots), genus_20 = NA)
taxtab20[top20otus, "genus_20"] <- as(tax_table(roots)[top20otus, "Genus"],
                                      "character")
taxtab20
tax_table(roots) <- tax_table(taxtab20)
roots <- transform_sample_counts(roots, function(x) 100 * x/sum(x))
title = "Relative Abundance"
genus_plot2 <- plot_bar(roots, "bar_label", fill = "genus_20", title = title)
print(genus_plot2)
```

# Below is code that works but still doesn't do exactly what I want it to do. 

# Genus Level
```{r}
genus_abundance <- roots %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate at genus level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  arrange(Genus) 
head(genus_abundance)
genus_abundance <- data.table(genus_abundance)
genus_abundance[(Abundance <= 0.04), Genus := "Other"]
all_genus <- genus_abundance %>%
  select(Genus, Sample, Abundance, bar_label) %>%
  group_by(Genus, bar_label) %>%
  summarize(
    avg_abundance = mean(Abundance)
  ) %>%
  filter(avg_abundance > 0.04)
head(all_genus)
#phylum_colors <- c(
#  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
#  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
#  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861")
ggplot(all_genus) +
  geom_col(mapping = aes(x = bar_label, y = avg_abundance, fill = Genus), position = "fill", show.legend = TRUE)+
  ylab("Relative Abundance") +
  #scale_fill_manual(values = phylum_colors) +
  xlab(NULL)+
  theme_minimal()+
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle = 90, vjust = 0.5, hjust = 1),
        axis.title.y = element_text(size = 20),
        title = element_text(size = 25))
ggplot(genus_abundance) +
  geom_col(mapping = aes(x = bar_label, y = Abundance, fill = Genus), position = "fill", show.legend = TRUE)+
  ylab("Relative Abundance") +
  #scale_fill_manual(values = phylum_colors) +
  xlab(NULL)+
  theme_minimal()+
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle = 90, vjust = 0.5, hjust = 1),
        axis.title.y = element_text(size = 20),
        title = element_text(size = 25))
```


# Phylum Level
```{r}
phylum_abundance <- roots %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  arrange(Phylum) 
head(phylum_abundance)
phylum <- data.table(phylum_abundance)
phylum[(Abundance <= 0.01), Phylum := "Other"]
all_phylum <- phylum %>%
  select(Phylum, Sample, Abundance, bar_label) %>%
  group_by(Phylum, bar_label) %>%
  summarize(
    avg_abundance = mean(Abundance)
  ) %>%
  filter(avg_abundance > 0.01)
head(all_phylum)
phylum_colors <- c(
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861")
ggplot(all_phylum)+
  geom_col(mapping = aes(x = bar_label, y = avg_abundance, fill = Phylum), position = "fill", show.legend = TRUE)+
  ylab("Relative Abundance") +
  scale_fill_manual(values = phylum_colors) +
  xlab(NULL)+
  theme_minimal()+
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle = 90, vjust = 0.5, hjust = 1),
        axis.title.y = element_text(size = 20),
        title = element_text(size = 25))
```

# below is code copied from Reid's script but it wasn't working well for me
```{r}
# filtering otus---------------------------------------------------------
# will filter now before creating barplots
# any sample with less than 5 reads for a particular otu will be placed to 0
otu_table(physeq_object_roots)[otu_table(physeq_object_roots) <= 4] <- 0
physeq_object_roots
sample_data(physeq_object_roots)
# removes any OTUs that has less than 10 total reads across all samples
otu_table(physeq_object_roots) <- otu_table(physeq_object_roots)[which(rowSums(otu_table(physeq_object_roots)) >= 10),]
sample_data(physeq_object_roots)
###barplots by fungal genus---------------------------
library(data.table)
library(dplyr)
library(ggplot2)
#roots
root_barplots <- merge_samples(physeq_object_roots, "bar_label")
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
  geom_bar(stat = "identity") +
  scale_x_discrete(limits= c("Conventional V2", "Conventional R2", "Conventional R6", "No-Till V2", "No-Till R2", "No-Till R6","Organic V2", "Organic R2", "Organic R6"))+
  scale_fill_manual(values = c("Albifimbria" ="#652926",
                               "Alternaria_1" = "steelblue", 
                               "Arnium" = "#C84248",
                               "Articulospora_1" = "darksalmon",
                               "Cladosporium_1" = "green",
                               "Clonostachys_1" = "#CD9BCD",
                               "Conlarium aquaticum" = "#AD6F3B",
                               "Corynespora_1" = "#673770", 
                               "Devriesia_1" = "#D14285",
                               "Dictyochaeta_1" = "#652926",
                               "Other" = "Blue",
                               "Didymella_1" ="pink",
                               "Fusarium_1" ="#673770",
                               "Epicoccum_1" = "#AD6F3B",
                               "Exophiala_1" ="#CBD588",
                               "Funneliformis_1" = "#5F7FC7", 
                               "Glomus" = "orange",
                               "Lachnum" = "#DA5724",
                               "Macrophomina_1" = "#508578",
                               "Monocillium_1" = "#CD9BCD",
                               "Mortierella_1" = "tan",
                               "Mrakia_1" = "magenta1",
                               "Neoascochyta_1" = "gray52",
                               "Neosetophoma_1" = "darkorange4",
                               "Occultifur_1" = "lightsalmon1"
  ))+
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  theme(legend.key.height = unit(0.15, "cm"), legend.key.width = unit(0.25, "cm")) +
  theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.text.x = element_text(size = 8, face = "bold")) +
  theme_classic()+
  ggtitle("Roots")+
  theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1)) +
  theme(plot.title = element_text(size = 14, hjust = 0.5)) +
  theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1)) +
  #theme(axis.ticks.x = element_blank()) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) +
  theme(legend.position="right")+
  
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Genera > 4%) \n") +
  xlab("")
plot(bar_ITS_root)
```



```{r}
roots %>%                                                              #phyloseq object
  plot_richness(
    x = "sampletype",                                                #compare diversity of datatype
    measures = c("Observed", "Shannon")) +                           #choose diversity measures
  geom_violin(aes(fill = sampletype), show.legend = FALSE)+          #make violin plot, set fill aes to sampletype
  geom_boxplot(width=0.1) +                                          #add boxplot, set width
  theme_classic()+                                                   #change theme to classic
  xlab(NULL)+                                                        #no label on x-axis
  theme(axis.text.y.left = element_text(size = 20),                  #adjust y-axis text
        axis.text.x = element_text(size = 20, hjust = 0.5, vjust = 1, angle = 0),           #adjust x-axis label position
        axis.title.y = element_text(size = 20))+                     #adjust y-axis title
  theme(strip.text = element_text(face = "bold", size = 25))+        #adjust headings
  #scale_fill_manual(values = sample_colors)+   #set fill colors
  scale_x_discrete(                                                  #change x-axis labels
    breaks = sample_types, 
    labels = sample_labels)+                   
  ggtitle("Alpha Diversity") +                                       #add title
  theme(plot.title=element_text(size = 25, face = "bold", hjust = 0.5)) #change title size, face and position
```