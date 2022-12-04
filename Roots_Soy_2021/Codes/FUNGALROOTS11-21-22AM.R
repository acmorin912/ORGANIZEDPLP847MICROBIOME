#Ashlynn Morin, HPCC R, Fungal Roots Analysis
#New code for 2021 data
#to get here, you will need to set your working directory as one of mine for all the data.
#PATHWAY FOR WORKING DIRECTORY: /mnt/research/bonito_lab/Morin/fungal_roots_2021/Needed data now
#hopefully accessible? If not I will add to my github as well :)

#Get Phyloseq (since install ____ isnt working)
source("https://raw.githubusercontent.com/joey711/phyloseq/master/inst/scripts/installer.R",
       local = TRUE)

#Getting Biostrings (since install ______ isn't working)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Biostrings")

#Getting "SpiecEasi"
install.packages(devtools)
library(devtools)
install_github("zdk123/SpiecEasi")
library(SpiecEasi)

#make sure packages there
library("Biostrings")
library("phyloseq")
library("ggplot2")
library("SpiecEasi")
library("igraph")
library("vegan")
library("stringi")
library("rhdf5")
library("zlibbioc")
library("S4Vectors")
library("yaml")
library("colorspace")
library("indicspecies")

#Time to get to work
#Create phyloseq ITS_OTU table
ITS_otus<- read.delim("otu_table_ITS_UPARSE_R1.txt",
                      row.names=1) 
ITS_otus

#ITS_OTU_PHY
ITS_otus_phy <-otu_table(ITS_otus,
                         taxa_are_rows = TRUE)
ITS_otus_phy

#ITS_Metadata (WILL NEED TO UPDATE THIS NAME!!!!)
ITS_metadata <-read.delim("root_fun_map_soybean_2021.txt",
                          row.names=1)
ITS_metadata
ITS_metadata_phy <-sample_data(ITS_metadata)

#ITS taxonomy
  ITS_taxonomy<- read.delim("constax_taxonomy.txt",
                            header= TRUE,
                            row.names =1)
  ITS_taxonomy
ITS_taxonomy_phy <- tax_table(as.matrix(ITS_taxonomy))
ITS_taxonomy_phy

#Sequences
ITS_sequences<- readDNAStringSet("otus_R1.fasta", format="fasta", seek.first.rec = TRUE, use.names = TRUE)
ITS_sequences

physeq_object_Fungi <- phyloseq(ITS_otus_phy, ITS_metadata_phy, ITS_taxonomy_phy, ITS_sequences)
physeq_object_Fungi
tax_table(physeq_object_Fungi)
sample_data(physeq_object_Fungi)

colnames(tax_table(physeq_object_Fungi)) = c(k = "Kingdom", p = "Phylum", c = "Class", o = "Order", f = "Family", g = "Genus", s = "Species")
tax_table(physeq_object_Fungi)


tax_table(physeq_object_Fungi)[, "Kingdom"] <- gsub("d:", "", tax_table(physeq_object_Fungi)[, "Kingdom"])
tax_table(physeq_object_Fungi)[, "Phylum"] <- gsub("p:", "", tax_table(physeq_object_Fungi)[, "Phylum"])
tax_table(physeq_object_Fungi)[, "Class"] <- gsub("c:", "", tax_table(physeq_object_Fungi)[, "Class"])
tax_table(physeq_object_Fungi)[, "Order"] <- gsub("o:", "", tax_table(physeq_object_Fungi)[, "Order"])
tax_table(physeq_object_Fungi)[, "Family"] <- gsub("f:", "", tax_table(physeq_object_Fungi)[, "Family"])
tax_table(physeq_object_Fungi)[, "Genus"] <- gsub("g:", "", tax_table(physeq_object_Fungi)[, "Genus"])
tax_table(physeq_object_Fungi)[, "Species"] <- gsub("s:", "", tax_table(physeq_object_Fungi)[, "Species"])
tax_table(physeq_object_Fungi)

#remove chloropast,mitochondria,cyanobacteria
physeq_object_Fungi <- subset_taxa(physeq_object_Fungi, Phylum!="Chloroplast")
physeq_object_Fungi <- subset_taxa(physeq_object_Fungi, Class!="Chloroplast")
physeq_object_Fungi <- subset_taxa(physeq_object_Fungi, Order!="Chloroplast")
physeq_object_Fungi <- subset_taxa(physeq_object_Fungi, Family!="Chloroplast")
physeq_object_Fungi <- subset_taxa(physeq_object_Fungi, Genus!="Chloroplast")
tax_table(physeq_object_Fungi)

physeq_object_Fungi <- subset_taxa(physeq_object_Fungi, Phylum!="Mitochondria")
physeq_object_Fungi <- subset_taxa(physeq_object_Fungi, Class!="Mitochondria")
physeq_object_Fungi <- subset_taxa(physeq_object_Fungi, Order!="Mitochondria")

physeq_object_Fungi <- subset_taxa(physeq_object_Fungi, Family!="Mitochondria")
physeq_object_Fungi <- subset_taxa(physeq_object_Fungi, Genus!="Mitochondria")
tax_table(physeq_object_Fungi)

physeq_object_Fungi_roots <- subset_samples(physeq_object_Fungi,origin%in%c("root"))

sample_data(physeq_object_Fungi_roots)

library(devtools)
library(processx)
devtools::install_github("benjjneb/decontam", force = TRUE)
library(decontam)


write.csv(sample_data(physeq_object_Fungi_roots), file = "sample_check1_roots.csv")
df_roots <- as.data.frame(sample_data(physeq_object_Fungi_roots)) # Put sample_data into a ggplot-friendly data.frame
df_roots$LibrarySize_roots <- sample_sums(physeq_object_Fungi_roots)
df_roots <- df_roots[order(df_roots$LibrarySize_roots),]
df_roots$Index <- seq(nrow(df_roots))
write.csv(df_roots, file = "rank_sums_roots.csv")
ggplot(data=df_roots, aes(x=Index, y=LibrarySize_roots)) + geom_point() #TOOK OUT COLOR CUZ IT DIDNT WORK WITH IT AS REID HAD IT

####HERE UNTIL THE NEXT 4 HASTTAGS WILL BE PROBLEMATIC ERRORS!
# filter by prevelance 
sample_data(physeq_object_Fungi_roots)$is.neg <- sample_data(physeq_object_Fungi_roots)$Sample_or_Control == "Control Sample"
contamdf.prev_roots <- isContaminant(physeq_object_roots, method="prevalence", neg="is.neg")
table(contamdf.prev_roots$contaminant)



# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa_roots <- transform_sample_counts(physeq_object_Fungi_roots, function(abund) 1*(abund>0))
ps.pa.neg_roots <- prune_samples(sample_data(ps.pa_roots)$Sample_or_Control == "Control Sample", ps.pa_roots)
ps.pa.pos_roots <- prune_samples(sample_data(ps.pa_roots)$Sample_or_Control == "True Sample", ps.pa_roots)
# Make data.frame of prevalence in positive and negative samples
df.pa_roots <- data.frame(pa.pos_roots=taxa_sums(ps.pa.pos_roots), pa.neg_roots=taxa_sums(ps.pa.neg_roots),
                          contaminant=contamdf.prev_roots$contaminant)
ggplot(data=df.pa_roots, aes(x=pa.neg_roots, y=pa.pos_roots, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

# remove contaminants
ps.noncontam_roots <- prune_taxa(!contamdf.prev_roots$contaminant, physeq_object_roots)
# with contaminants removed
otu_table(ps.noncontam_roots)
# remove negative controls
#AH SHIT I HAVE TO FIGURE OUT NEGATIVE CONTROLS?
ps.noncontam_roots <- subset_samples(ps.noncontam_roots, experiment%in%c("obj_1","obj_2"))
otu_table(ps.noncontam_roots) <- otu_table(ps.noncontam_roots)[which(rowSums(otu_table(ps.noncontam_roots)) >= 1),]
ps.noncontam_roots



write.csv(otu_table(ps.noncontam_roots), file = "filtering_low_roots.csv")


otu_table(ps.noncontam_roots) <- subset(otu_table(ps.noncontam_roots),
                                        select = -c(T1R5CBR3R,T2R1CBR6R,T4R5CR2R,T4R2BR2R,T4R1AR2R,T1R6FCR4R,T2R5AR2R,T1R1BR6R,T1R2FCR4R,T1R2AR2R,T1R2FAR4R,T1R6CCR4R,T1R2AR2R,T1R5FAR3R,T1R1FAR3R,T1R6CR2R,T1R2FAR3R,T4R6BR2R))
ps.noncontam_roots

#### END PROBLEMATIC?


