
######## Phyloseq object assembly code from Jenn 

#Load libraries
library(phyloseq)
library(Biostrings)
library(ggplot2)
library(tidyverse)

#load the phyloseq package and import CSV to matrix for combined (normal and preEclampsia) Vaginal microbiome samples

OTU_VC_mat <- (as.matrix(read_csv("OTU_C.csv")))

Taxonomy_C_mat <- (as.matrix(read_csv("Taxonomy_C.csv")))

#import the OTU and taxonomy matrix into phyloseq with both N and PE with new species from specialteID

OTU_C = otu_table(OTU_VC_mat, taxa_are_rows = TRUE)

TAX_C = tax_table(Taxonomy_C_mat)

OTU_C

TAX_C

physeq_C = phyloseq(OTU_C, TAX_C)

physeq_C

sample_names(physeq_C)

#load metadata file

library(tidyverse)

meta_data <- read_csv("Metadata_table_C.csv")

sdataC <- meta_data

head(sdataC)

#creates a dataframe that can be read by phyloseq for the metadata file

sampledataC = sample_data(data.frame(
  
  sdataC, row.names=sample_names(physeq_C), stringsAsFactors=FALSE))

sampledataC

#merges metadata with phyloseq data

physeq4 = merge_phyloseq(physeq_C, sampledataC)

physeq4



#plots the genus level for each sample plus code to put the names in order on the x axis

plot_bar(physeq4, x= "sample_name", fill="Genus") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

#plots the species level for each sample plus code to put the names in order on the x axis

plot_bar(physeq4, x= "sample_name", fill="Species") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

#groups by Treatment "SampleType" the genus level abidance

plot_bar(physeq4, x= "SampleType", fill = "Genus", title = "title")

#plots the alpha diversity for each sample

plot_richness(physeq4, x= "sample_name", measures="Shannon", color= "SampleType")

#normalize data for all 20 samples plus code to put the names in order on the x axis

ps2 <- transform_sample_counts(physeq4, function(x) x / sum(x))

ps2

plot_bar(ps2, x= "sample_name", fill="Genus") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

plot_bar(ps2, x= "sample_name", fill="Species") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

plot_bar(ps2,"SampleType", fill = "Genus", title = "Vaginal Micrrobiome")

plot_bar(ps2, "Family", fill="Genus", facet_grid=~SampleType)

plot_bar(ps2, "Family", fill="SampleType", facet_grid=~sample)

#shannon and Observed alpha diversity measures combined

plot_richness(physeq4, x = "SampleType", measures = c("Observed", "Shannon")) +
  
  geom_boxplot(aes(fill = SampleType), show.legend = FALSE)

#shannon and Observed alpha diversity measures independent samples

plot_richness(physeq4, x = "sample_name", measures = c("Observed", "Shannon")) +
  
  geom_boxplot(aes(fill = SampleType), show.legend = FALSE)

#alpha diversity measures independent samples

plot_richness(physeq4, x = "sample_name")

#ordination plotting

ORD <- ordinate(physeq4, method=‘MDS’,distance=‘bray’)

plot_ordination(physeq4,ORD,color=‘SampleType’) + geom_point(size=5)

#Kriskall-Wallis test for stat difference across all samples combined

alphadiv <- estimate_richness(physeq4, measures = c("Observed", "Shannon")) %>%
  
  rownames_to_column(var = "SampleID") %>%
  
  left_join(as.data.frame(sample_data(physeq4)), by = "SampleID")

alphadiv

#Observed

Kruskal_Observed<-kruskal.test(Observed ~ SampleType, data = alphadiv)

#Shannon

Kruskal_Shannon<-kruskal.test(Shannon ~ SampleType, data = alphadiv)

#
write.csv(alphadiv, "PreEclampsia_alphadiv.csv")