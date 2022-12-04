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
bar_ITS_root <- 
  ggplot(dat_root_bp, aes(x = Sample, y = Abundance, fill = Genus)) + 
  #facet_wrap(~Stage, strip.position = "bottom") +
  theme(axis.text.x = element_text(angle = 90))+
  geom_bar(stat = "identity") +
  scale_x_discrete(limits= c("Conventional V2", "Conventional R2", "Conventional R6", "No-Till V2", "No-Till R2", "No-Till R6","Organic V2", "Organic R2", "Organic R6"))+
  scale_fill_manual(values = c("Albifimbria_1" ="#652926",
                               "Alternaria_1" = "steelblue", 
                               "Arnium_1" = "#C84248",
                               "Articulospora_1" = "darksalmon",
                               "Cladosporium_1" = "green",
                               "Clonostachys_1" = "#CD9BCD",
                               "Conlarium aquaticum_1" = "#AD6F3B",
                               "Corynespora_1" = "#673770", 
                               "Devriesia_1" = "#D14285",
                               "Dictyochaeta_1" = "#652926",
                               "Other_1" = "Blue",
                               "Didymella_1" ="pink",
                               "Fusarium_1" ="#673770",
                               "Epicoccum_1" = "#AD6F3B",
                               "Exophiala_1" ="#CBD588",
                               "Funneliformis_1" = "#5F7FC7", 
                               "Glomus_1" = "orange",
                               "Lachnum_1" = "#DA5724",
                               "Macrophomina_1" = "#508578",
                               "Monocillium_1" = "#CD9BCD",
                               "Mortierella_1" = "tan",
                               "Mrakia_1" = "magenta1",
                               "Neoascochyta_1" = "gray52",
                               "Neosetophoma_1" = "darkorange4",
                               "Occultifur_1" = "lightsalmon1"
                               "Metacordyceps_1" = "darkolivegreen1"
                               "Paraphoma_1" = "brown"
                               "Plectosphaerella_1" = "darkviolet"
                               "Setophoma_1" = "darkgoldenrod2"
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