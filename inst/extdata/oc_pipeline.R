
######################################################
#DATA SOURCE: 
# https://github.com/fernanda-t/microbiome-oc
# Terrazas, Fernanda, et al. "Influence of menstrual cycle and oral contraception on taxonomic composition and gas production in the gut microbiome." Journal of Medical Microbiology 74.3 (2025): 001987.
######################################################


######################################################
#Assemble phyloseq objects
######################################################

#Files were downloaded from GitHub and manually rearranged in a spreadsheet in order to generate the phyloseq object
library(phyloseq)
library(ggplot2)
library(DESeq2)

# Assemble phyloseq object
taxa <- read.csv("Terrazas_taxtab.csv", sep=",", row.names=1)
samdf <- read.csv("Terrazas_samdf.csv", row.names=1)
seqtab <- read.csv("Terrazas_seqtab.csv", sep=",", row.names=1)


# Correct factors to create categorical variables out of numerical data
samdf$subject <- as.factor(samdf$subject)
samdf$sample_id <- as.factor(samdf$sample_id)
samdf$OCE2Dosage <- as.factor(samdf$OCE2Dosage)

# Construct phyloseq object (straightforward from dada2 outputs)
#Don't forget as matrix for tax_table!!! For some reason, will throw an error that the taxa/OTU names don't match otherwise
oc_counts <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), 
                          sample_data(samdf), 
                          tax_table(as.matrix(taxa)))

# Remove any samples or ASVs that have 0 reads across all samples 
oc_counts <- prune_samples(sample_sums(oc_counts) >= 1, oc_counts)
oc_counts <- prune_taxa(taxa_sums(oc_counts) >= 1, oc_counts)


# Create normalized data using proportions
oc <- transform_sample_counts(oc_counts, function(otu) otu/sum(otu))


# Rename the samples for convenience's sake if desired
#newnames <- scan(file.choose(), what = "", quiet=TRUE)
#sample_names(ps)<-newnames


######################################################
#Remove outlier 948N2
######################################################
oc_no948N2 = subset_samples(oc, sample_id != "948N2")
oc_no948N2_counts = subset_samples(oc_counts, sample_id != "948N2")


######################################################
################# Creating object including only subjects with recorded dosage data
######################################################

# Assemble phyloseq object
taxa <- read.csv("Terrazas_taxtab.csv", sep=",", row.names=1)
samdf <- read.csv("Terrazas_dosage_samdf.csv", row.names=1)
seqtab <- read.csv("Terrazas_dosage_seqtab.csv", sep=",", row.names=1)

# Correct factors
samdf$subject <- as.factor(samdf$subject)
samdf$OCE2Dosage <- as.factor(samdf$OCE2Dosage)
samdf$sample_id <- as.factor(samdf$sample_id)

# Construct phyloseq object (straightforward from dada2 outputs)
#Don't forget as matrix for tax_table!!! For some reason, will throw an error that the taxa/OTU names don't match otherwise
oc_dosage_counts <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), 
                             sample_data(samdf), 
                             tax_table(as.matrix(taxa)))


# Remove any samples and ASVs that have 0 reads across all samples (In this case 0 ASVs)
oc_dosage_counts <- prune_samples(sample_sums(oc_dosage_counts) >= 1, oc_dosage_counts)
oc_dosage_counts <- prune_taxa(taxa_sums(oc_dosage_counts) >= 1, oc_dosage_counts)


# Create normalized data using proportions
oc_dosage <- transform_sample_counts(oc_dosage_counts, function(otu) otu/sum(otu))



######################################################
################# Creating object including only subjects with recorded agefirstOCuse and totalyearsonOC data
######################################################

# Assemble phyloseq object
taxa <- read.csv("Terrazas_taxtab.csv", sep=",", row.names=1)
samdf <- read.csv("Terrazas_ageOCuse_samdf.csv", row.names=1)
seqtab <- read.csv("Terrazas_ageOCuse_seqtab.csv", sep=",", row.names=1)


# Correct factors
samdf$subject <- as.factor(samdf$subject)
samdf$OCE2Dosage <- as.factor(samdf$OCE2Dosage)
samdf$sample_id <- as.factor(samdf$sample_id)

# Construct phyloseq object (straightforward from dada2 outputs)
#Don't forget as matrix for tax_table!!! For some reason, will throw an error that the taxa/OTU names don't match otherwise
oc_agefirstOCuse_counts <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), 
                                    sample_data(samdf), 
                                    tax_table(as.matrix(taxa)))


# Remove any ASVs that have 0 reads across all samples and any samples that have 0 reads across all ASVs
oc_agefirstOCuse_counts <- prune_samples(sample_sums(oc_agefirstOCuse_counts) >= 1, oc_agefirstOCuse_counts)
oc_agefirstOCuse_counts <- prune_taxa(taxa_sums(oc_agefirstOCuse_counts) >= 1, oc_agefirstOCuse_counts)


# Create normalized data using proportions
oc_agefirstOCuse <- transform_sample_counts(oc_agefirstOCuse_counts, function(otu) otu/sum(otu))


#Export phyloseq objects
saveRDS(oc, "oc.RDS")
saveRDS(oc_counts, "oc_counts.RDS")
saveRDS(oc_no948N2, "oc_no948N2.RDS")
saveRDS(oc_no948N2_counts, "oc_no948N2_counts.RDS")
saveRDS(oc_dosage, "oc_dosage.RDS")
saveRDS(oc_dosage_counts, "oc_dosage_counts.RDS")
saveRDS(oc_agefirstOCuse, "oc_agefirstOCuse.RDS")
saveRDS(oc_agefirstOCuse_counts, "oc_agefirstOCuse_counts.RDS")



# Test barplots
plot_bar(oc, x="OCE2Dosage", fill = "Class", title = "Test") +
  geom_bar(aes(color = Class, fill = Class), stat = "identity", position = "fill")


plot_bar(oc, x="OCE2Dosage", fill = "Class", title = "Terrazas") +
  geom_bar(aes(color = Class, fill = Class), stat = "identity", position = "stack")


merge <- merge_samples(oc, "OCE2Dosage")
merge <- transform_sample_counts(merge, function(x) 100 * x/sum(x))
sample_data(merge)$OCE2Dosage = sample_names(merge)
plot_bar(merge, x="OCE2Dosage", fill = "Class", title = "Dosage") +
  geom_bar(aes(color = Class, fill = Class), stat = "identity", position = "stack")




# Test alpha diversity plot
alpha <-plot_richness(oc_counts, x="subject", 
                             color="ethnicity", 
                             shape="treatment", 
                             measures="Shannon")
alpha

# Test NMDS plot
# In this study a single outlier exists
ord.bray <- ordinate(oc, method="NMDS", distance="bray")
p <- plot_ordination(oc, ord.bray, color="ethnicity", shape = "treatment", title="NMD")
p + geom_point(size=3, position="jitter") + 
  theme_classic()
  #geom_text(aes(label=sample_id), nudge_x=0.45, nudge_y=0.1)

#Differential abundance
# STEP 1: Convert the phyloseq object to a DESeq2 object and tell R the experimental design
my_DESeq2 <- phyloseq_to_deseq2(oc_counts, design = ~ age)

# STEP 2: Select the groups to compare
name_for_comparision <-"age"

# STEP 3: Run the differential abundance analysis at the chosen p-value
Significant_DEseq2_ASVs<-Differential_Abundance_Continuous(my_DESeq2, name_for_comparision, 0.05)

# STEP 4: Retrieve the list of ASVs with a significant difference in abundance between the chosen groups
Significant_DEseq2_ASVs

# STEP 5: Plot the results with your chosen x axis and legend
ggplot(Significant_DEseq2_ASVs, aes(x = Phylum, y=log2FoldChange, color= Class, shape = Order)) + geom_point(size=6) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))



# Common core microbiome
#https://microbiome.github.io/tutorials/Core.html

library(microbiome)

#Subset phyloseq object into the categories you want to look for the core microbiome in (ex. healthy vs disease state samples)
#Using count data
user_samples <- subset_samples(oc_counts, treatment=="User")
control_samples <- subset_samples(oc_counts, treatment=="Control")


# Calculate compositional version of the data
# (relative abundances)
rel.user_samples <- microbiome::transform(user_samples, "compositional")
rel.control_samples <- microbiome::transform(control_samples, "compositional")

#names of core microbiome:
#setting to found in 70%
core.taxa.standard.user_samples <- core_members(rel.user_samples, detection = 0, prevalence = 70/100)
core.taxa.standard.control_samples <- core_members(rel.control_samples, detection = 0, prevalence = 70/100)
core.taxa.standard.user_samples
core.taxa.standard.control_samples




#Rearranging samples in barplot
desiredOrder <- c("701C", "702C", "705C", "706C", "916C", "917C", "918C", "919C", "920C", 
                  "921C", "922C", "923C", "924C", "925C", "926C", "927C", "928C", "929C", 
                  "930C", "701C2", "702C2", "705C2", "706C2", "918C2", "919C2", "921C2", 
                  "922C2", "923C2", "924C2", "925C2", "926C2", "927C2", "928C2", "929C2", 
                  "732N", "734N", "735P", "736N", "737P", "739N", "740P", "741P", "742P", 
                  "743P", "744P", "947N", "948N", "949N", "950N", "951N", "952P", "953P", 
                  "954P", "955N", "956N", "957P", "958P", "960N", "962P", "963P", "966N",
                  "967N", "970P", "971P", "972P", "973P", "975N", "978P", "980N", "981P",
                  "982N", "983N", "984N", "985P", "986N", "987P", "988P", "989N", "990P",
                  "732N2", "735P2", "736N2", "737P2", "739N2", "741P2", "742P2", "947N2",
                  "948N2", "949N2", "950N2", "951N2", "952P2", "953P2", "954P2", "955N2", 
                  "956N2", "957P2", "958P2", "960N2", "962P2", "963P2", "970P2", "971P2", 
                  "972P2", "973P2", "975N2", "978P2", "979N2", "981P2", "984N2", "985P2", 
                  "986N2", "988P2", "989N2", "990P2")


glom <- tax_glom(oc, taxrank = 'Class')
glom # Squish all ASVs of the same taxonomic rank together
melt <- psmelt(glom) # create dataframe from phyloseq object
melt$sample_id <- factor(melt$sample_id, levels = desiredOrder) #Reorder samples (by reordering factors) of the melt dataframe

p <- ggplot(data=melt, aes(x=sample_id, y=Abundance, fill=Class))
p + geom_bar(aes(), stat="identity", position="stack") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))





#######################################################
# Loading libraries and RDS
#######################################################

library(phyloseq)
library(ggplot2)
library(DESeq2)


oc<-readRDS("oc.RDS")
oc_counts<-readRDS("oc_counts.RDS")
oc_no948N2<-readRDS("oc_no948N2.RDS")
oc_no948N2_counts<-readRDS("oc_no948N2_counts.RDS")
oc_dosage <- readRDS("oc_dosage.RDS")
oc_dosage_counts <- readRDS("oc_dosage_counts.RDS")
oc_agefirstOCuse <- readRDS("oc_agefirstOCuse.RDS")
oc_agefirstOCuse_counts <- readRDS("oc_agefirstOCuse_counts.RDS")

