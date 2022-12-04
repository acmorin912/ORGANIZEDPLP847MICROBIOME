library(phyloseq)
library(Biostrings)
library(tidyverse)
library(vegan)
library(decontam)
library(ggpubr)
#import
mapping <- read.delim("FIXED_FIXED_root_fun_map_soybean_2021.txt", row.names = 1, header= TRUE, sep = "\t")
head(mapping)

otutable<- read.delim("otu_table_ITS_UPARSE_R1.txt", row.names = 1, header= TRUE, sep = "\t")
head(otutable)

taxonomy<- read.delim("constax_taxonomy.txt", row.names = 1, header= TRUE, sep= "\t")
head(taxonomy)

sequences<- readDNAStringSet("otus_R1.fasta", format = "fasta", seek.first.rec = TRUE, use.names = TRUE)
sequences
taxonomy[, "Kingdom"] <- as.factor(gsub("_1", "", taxonomy[, "Kingdom"]))
taxonomy[, "Phylum"] <- as.factor(gsub("_1", "", taxonomy[, "Phylum"]))
taxonomy[, "Class"] <- as.factor(gsub("_1", "", taxonomy[, "Class"]))
taxonomy[, "Order"] <- as.factor(gsub("_1", "",taxonomy[, "Order"]))
taxonomy[, "Family"] <- as.factor(gsub("_1", "", taxonomy[, "Family"]))
taxonomy[, "Genus"] <- as.factor(gsub("_1", "", taxonomy[, "Genus"]))
taxonomy[, "Species"] <- as.factor(gsub("_1", "", taxonomy[, "Species"]))

table(taxonomy$High_level_taxonomy)
table(taxonomy$Kingdom)

untarget <- c("Ichthyosporia","Metazoa", "Protista", 
              "Rhizaria", "Viridiplantae", "Stramenopila")

apply(taxonomy, 2, function(x) which(x %in% untarget))

###ISSUES HERE
taxonomy_bad <- 
  subset(
    taxonomy,
    High_level_taxonomy%in%untarget |
      Kingdom%in%untarget  &
      HL_hit_query_cover<60 |
      HL_hit_percent_id<60  )

dim(taxonomy_bad)
head(taxonomy_bad)

taxonomy_filt <-
  taxonomy[!(rownames(taxonomy) %in% rownames(taxonomy_bad)), ]
dim(taxonomy_filt)
### Proceeding WITHOUT the above. come back if needed!

#Phyloseq object
physeq<- phyloseq(otu_table(otutable, taxa_are_rows = TRUE), sample_data(mapping), tax_table(as.matrix(taxonomy)), sequences)

physeq
head(sample_data(physeq))
tax_table(physeq)[tax_table(physeq)==""]<- NA
head(tax_table(physeq))

#decontamination, revist? Possible error?
count(as.data.frame(as.matrix(sample_data(
  physeq))), Description, BarcodePlate)

table(physeq@sam_data$BarcodePlate)
# adding column for controls
physeq@sam_data$is.neg <- physeq@sam_data$Description == "Control"
head(sample_data(physeq))

contaminants <-
  isContaminant(
    physeq,
    method = "prevalence",
    batch = "BarcodePlate",
    neg = "is.neg",
    threshold = 0.4)

table(contaminants$contaminant) 


# remove contaminants ---------------------------------------------------
remove_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  myTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(myTaxa, physeq))
}
  
physeq_filt <-
    remove_taxa(physeq, rownames(subset(
      contaminants, contaminant%in%c("TRUE"))))
  
physeq_filt
  
  
  # plot depth --------------------------------------------------
physeq_filt@sam_data$LibrarySize <- as.numeric(as.character(sample_sums(physeq_filt)))
physeq_filt@sam_data$Index <- as.numeric(as.character(seq(nrow(physeq_filt@sam_data))))
head(physeq_filt@sam_data)
  
  
  # function to plot depth
PlotDepth <- function(physeq){
df <-
      as(sample_data(physeq), "matrix")
df <- 
      as.data.frame(df)
    # reconvert to numeric
df$LibrarySize <- as.numeric(as.character(df$LibrarySize))
df$Index <- as.numeric(as.character(df$Index))
    # order
df <- df[order(df$LibrarySize), ]
df$Index <- seq(nrow(df))
    # inspect
ggplot(data=df, aes(x=Index, y=LibrarySize, color=is.neg)) +
      geom_point(alpha =0.7) -> plot_dist
    return(plot_dist)  
  }

PlotDepth(physeq_filt) +
    labs(title="OTU", 
         subtitle = "Samples read depth", 
         x="Sample index", 
         y="Read number")
  # remove control samp[les  
physeq_clean <-
    subset_samples(physeq_filt, is.neg%in%c("FALSE"))
  
otu_table(physeq_clean) <-
    otu_table(physeq_clean)[which(rowSums(otu_table(physeq_clean)) > 0), ]
  
physeq_clean
  # rarefying dataset at even depth -----------------------------------------------------------------------
set.seed(2018)
sort(colSums(otu_table(physeq_clean)))
  
physeq_rare = rarefy_even_depth(physeq_clean, sample.size = 8000, 
                                  rngseed = FALSE, replace = TRUE, 
                                  trimOTUs = TRUE, verbose = TRUE)
  
otu_table(physeq_rare) <- otu_table(physeq_rare)[which(rowSums(otu_table(physeq_rare)) >= 1),]
physeq_rare
  
colSums(otu_table(physeq_rare))
any(taxa_sums(physeq_rare) == 0)
  
physeq_rare
  
  
  # Diversity analysis ------------------------------------------
  
  # Barplot -----------------------------------------------------
ExtractBar <- function(physeq, last_gen, Rank){
    print(head(tax_table(physeq)))
    top30 <- 
    names(sort(taxa_sums(physeq), TRUE)[1:last_gen]) # get 30 top genera
 physeq_gen <- 
    prune_taxa(top30, physeq)
    print(head(physeq_gen@sam_data))
    df_bar <- 
    physeq_gen %>%
    tax_glom(taxrank = Rank) %>%                     
    transform_sample_counts(function(x) {x/sum(x)} ) %>% 
    psmelt() %>%                                         
    #filter(Abundance > 0.01) %>%                         
     arrange(get(Rank))                                     
    print(levels(df_bar[,Rank]))
    return(df_bar)
  }
  
PlotBar <- function(df, Rank){
  barplot <-
  ggplot(df, aes(x = Description, # you can also change x to another variable 
                     y = Abundance, 
                     fill = get(Rank))) + 
  geom_bar(stat = "identity") +
    #facet_grid(~Niche, scales = "free_x", space="free_x") +
  facet_wrap(~ Description, scales = "free_x") + #This is the facet, change this to the variable you want to see
  guides(fill = guide_legend(ncol = 1, title = "Taxa") )
  return(barplot)
}
  
head(physeq_rare@sam_data)
physeq_rare@sam_data$Description
PlotBar(ExtractBar(physeq_rare, 30, "Genus"), "Genus") 
  
  # with your own color palette
palette_new30 <- c("#560d0d","#a35151", "#dba4a4", "#cc1c1c","#111b77","#283dff","#636bb7",
                     "#bfc5ff","#014443","#195637","#117744","#60ffaf","#b7ffdb","#825121",
                     "#ea7f17","#fcb067","#ffe8d3","#d8d6d4","#82807f", "#3f3e3d","#000000",
                     "#5b5b19","#fcfc00","#ffff9e","#ffb7ef","#fa7efc","#ae09ea","#521899",
                     "#a0fffc","#1e0047")
  # Ordination ----------------------------------------------------------------
  
nmds <- ordinate(physeq_rare, method ="NMDS", distance="bray", trymax=100) # this doesnt converge, need more test
pcoa <- ordinate(physeq_rare, method ="PCoA", distance="bray") 
  
pcoa 
  
  # taxa
p1 <- plot_ordination(physeq_rare, pcoa, type="taxa", color="Phylum", title="Taxa")
p1
  
  # samples
p2 <- plot_ordination(physeq_rare, pcoa, type="samples", color="Target", title="Target")
p2
  
  # try faceting, you have to play with this yourself
head(physeq_rare@sam_data)
p2 + facet_wrap(~Host, scales = "free_x")
p2 + facet_wrap(~Treatment, scales = "free_x")
  
  
  # biplot 
p3 <- plot_ordination(physeq_rare, pcoa, type="biplot", color="Phylum", shape="Host", title="Biplot")
p3
  
  