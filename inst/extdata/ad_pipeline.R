
#PIPELINE HELPFUL LINKS
# https://astrobiomike.github.io/amplicon/dada2_workflow_ex#generating-a-count-table
# https://benjjneb.github.io/dada2/tutorial.html
# https://www.bioconductor.org/packages/release/bioc/vignettes/dada2/inst/doc/dada2-intro.html#filter-and-trim

#DATA SOURCE: 
# BIOPROJECT: PRJNA554111
# Zhuang, Zhen-Qian, et al. "Gut microbiota is altered in patients with Alzheimer’s disease." Journal of Alzheimer’s Disease 63.4 (2018): 1337-1346.

####################################################################
#Processing raw data, assigning taxonomy, creating phyloseq objects
####################################################################

#load requisite packages
library(dada2)
library(ggplot2)
library(phyloseq)
library(DESeq2)

# Set path to unzipped, renamed fastq files downloaded to SRA
path <- "./Raw data"

# Get a list of all files in the directory at the end of path
fns <- list.files(path)

# Only select the .fastq files.
fastqs <- fns[grepl(".fastq", fns)]
# Sort the fastqs so they are all in order
fastqs <- sort(fastqs)
# Designate which files represent the forward and which represent the reverse reads based on the filename
fnFs <- fastqs[grepl("_1", fastqs)] 
fnRs <- fastqs[grepl("_2", fastqs)] 

# Get sample names, assuming files named as so: SAMPLENAME_XXX.fastq; modify as neccessary.
# It is okay if the designations for forward and reverse are lost at this point.
sample.names.F <- sapply(strsplit(fnFs, "_"), `[`, 2)
sample.names.R <- sapply(strsplit(fnRs, "_"), `[`, 2) 

# Get the file paths for the forward and reverse reads
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)

####Inspect Read Quality Profiles####
#Using the package "dada2", lets look at quality profile of forward reads
plotQualityProfile(fnFs[c(1,2,3,4,5,6,7,8,9)])
plotQualityProfile(fnFs[c(10,11,12,13,14,15,16,17,18)])
plotQualityProfile(fnFs[c(19,20,21,22,23,24,25,26,27)])
plotQualityProfile(fnFs[c(28,29,30,31,32,33,34,35,36)])
plotQualityProfile(fnFs[c(37,38,39,40,41,42,43,44,45)])
plotQualityProfile(fnFs[c(46,47,48,49,55,51,52,53,54)])
plotQualityProfile(fnFs[c(55,56,57,58,59,60,61,62,63)])
plotQualityProfile(fnFs[c(64,65,66,67,68,69,70,71,72)])
plotQualityProfile(fnFs[c(73,74,75,76,77,78,79,80,81)])
plotQualityProfile(fnFs[c(82,83,84,85,86,87,88,89,90)])

#Now, lets look at quality profile of reverse reads
plotQualityProfile(fnRs[c(1,2,3,4,5,6,7,8,9)])
plotQualityProfile(fnRs[c(10,11,12,13,14,15,16,17,18)])
plotQualityProfile(fnRs[c(19,20,21,22,23,24,25,26,27)])
plotQualityProfile(fnRs[c(28,29,30,31,32,33,34,35,36)])
plotQualityProfile(fnRs[c(37,38,39,40,41,42,43,44,45)])
plotQualityProfile(fnRs[c(46,47,48,49,55,51,52,53,54)])
plotQualityProfile(fnRs[c(55,56,57,58,59,60,61,62,63)])
plotQualityProfile(fnRs[c(64,65,66,67,68,69,70,71,72)])
plotQualityProfile(fnRs[c(73,74,75,76,77,78,79,80,81)])
plotQualityProfile(fnRs[c(82,83,84,85,86,87,88,89,90)])


#Make directory and file names for the trimmed/filtered fastqs
filt_path <- file.path(path, "trimmed")
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sample.names.F, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names.R, "_R_filt.fastq.gz"))


### Filter and Trim ####
# Filter using the "filterandTrim" function in dada2
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                    truncLen=c(250,250), 
                    maxN=0, #DADA does not allow Ns
                    maxEE=c(1,1), #allow 1 expected errors, where EE = sum(10^(-Q/10)); more conservative, model converges
                    truncQ=2,
                    rm.phix=TRUE, #remove reads matching phiX genom
                    compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE

# Learn Error Rates
# Note that this steps fails if there are no differences in the quality scores of reads. This can happen sometimes where the quality score is replaced on databases such as SRA (although I have had some luck using SRA explorer  without encourntering this issue)
setDadaOpt(MAX_CONSIST=16) #increase number of cycles to allow convergence
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE) 
plotErrors(errR, nominalQ=TRUE) 


##### Dereplicate reads 
#Dereplication combines all identical sequencing reads into into "unique sequences" with a corresponding "abundance": the number of reads with that unique sequence. 
#Dereplication substantially reduces computation time by eliminating redundant comparisons.
#DADA2 retains a summary of the quality information associated with each unique sequence. The consensus quality profile of a unique sequence is the average of the positional qualities from the dereplicated reads. These quality profiles inform the error model of the subsequent denoising step, significantly increasing DADA2’s accuracy.
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names.F
names(derepRs) <- sample.names.R


##### Infer Sequence Variants 

dadaFs <- dada(derepFs, err=errF, multithread=FALSE)
dadaRs <- dada(derepRs, err=errR, multithread=FALSE)

#now, look at the dada class objects by sample
#will tell how many 'real' variants in unique input seqs
#By default, the dada function processes each sample independently, but pooled processing is available with pool=TRUE and that may give better results for low sampling depths at the cost of increased computation time. See our discussion about pooling samples for sample inference. 
dadaFs[[1]]
dadaRs[[1]]


##### Merge paired reads 

#To further cull spurious sequence variants
#Merge the denoised forward and reverse reads
#Paired reads that do not exactly overlap are removed

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

summary((mergers[[1]]))


##### Construct sequence table 
#a higher-resolution version of the "OTU table" produced by classical methods

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

# Look at the length of sequences. Do they seem reasonable for the expected size of your barcode?
# V3-V4 case we can expect anywhere from 460-480bp depending on primers...
# I wasn't that conservative here as some variation in size can also occur naturally...
plot(table(nchar(getSequences(seqtab)))) 

#The sequence table is a matrix with rows corresponding to (and named by) the samples, and 
#columns corresponding to (and named by) the sequence variants. 
#Sequences that are much longer or shorter than expected may be the result of non-specific priming, and may be worth removing


##### Remove chimeras 
#The core dada method removes substitution and indel errors, but chimeras remain. 
#Fortunately, the accuracy of the sequences after denoising makes identifying chimeras easier 
#than it is when dealing with fuzzy OTUs: all sequences which can be exactly reconstructed as 
#a bimera (two-parent chimera) from more abundant sequences.

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE, verbose=TRUE)
dim(seqtab.nochim) #SEQTAB is a count table

sum(seqtab.nochim)/sum(seqtab)
#The fraction of chimeras varies based on factors including experimental procedures and sample complexity, 
#but can be substantial. 

###### Get read counts for every step in the process and save them in a nice .csv file

getN <- function(x) sum(getUniques(x))
track_F <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
colnames(track_F) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track_F) <- sample.names.F
head(track_F)
tail(track_F)

#Generates a new file for forward reads
write.csv(track_F,file="readstats_F.csv",row.names=TRUE,quote=FALSE)

getN <- function(x) sum(getUniques(x))
track_R <- cbind(out, sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
colnames(track_R) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track_R) <- sample.names.R
head(track_R)
tail(track_R)


#Generates a new file for reverse reads
write.csv(track_R,file="readstats_R.csv",row.names=TRUE,quote=FALSE)

####### Assign taxonomy using Silva database formatted for dada2
# I usually save the environment or objects and reload the entire RStudio environment to clear the memory
# Assigning taxonomy is a very computationally heavy task. You can split the data and do it in parts or run it over night if needed.

#Assign Taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "../silva_nr99_v138.2_toSpecies_trainset.fa.gz",tryRC=TRUE)

# I believe the assignment of the Species level has already been done by the previous line, but in previous versions of the database, we had to assign the species separately
#unname(head(taxa))
#taxa.plus <- addSpecies(taxa, "../silva_v138.2_assignSpecies.fa.gz",tryRC=TRUE,verbose=TRUE)

# Save as .RDS 
# I restart the RStudio environment to clean it out and get memory back.
seqtab.nochim<-readRDS("seqtab.nochim.RDS")
taxa<-readRDS("taxa.RDS")

# Import dataframe holding sample information
# Sample information must be made manually
samdf<-read.csv("PRJNA554111_samdf.csv")
head(samdf)
rownames(samdf) <- samdf$sample

# Construct phyloseq object (straightforward from dada2 outputs)
ad_counts <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))

# Remove any ASVs that have 0 reads across all samples (in this case, no ASVs)
ad_counts <- prune_samples(sample_sums(ad_counts) >= 1, ad_counts)
oc_agefirstOCuse_counts <- prune_taxa(taxa_sums(oc_agefirstOCuse_counts) >= 1, oc_agefirstOCuse_counts)


# Rename the samples for convenience's sake. Not a critical step.
newnames <- scan("PRJNA554111_sample_names.txt", what = "", quiet=TRUE)
sample_names(ad_counts)<-newnames

# Normlaize data by proportions
ad <- transform_sample_counts(ad_counts, function(otu) otu/sum(otu))

# Save the phyloseq objects as .RDS files
saveRDS(ad, "ad.RDS")
saveRDS(ad_counts, "ad_counts.RDS")


####################################################################
#Testing
####################################################################

library(ggplot2)
library(phyloseq)
library(DESeq2)

#Load RDS if needed here before testing
ad<-readRDS("ad.RDS")
ad_counts<-readRDS("ad_counts.RDS")

# Test alpha diversity plot
plot_richness(ad_counts, x="group", measures=c("Shannon", "Simpson"), color="sex")

# Test multidimensional analysis plot
# In this study, PRJNA554111, it seems like some outliers exist.
# Error is occuring showing the geom_text labels...not sure why.
ord.bray <- ordinate(ad, method="PCoA", distance="bray")
p <- plot_ordination(ad, ord.bray, color="group", shape = "sex", title="PCoA")
p + geom_point(size=3) + 
  theme_classic() +
  geom_text(aes(label=sample_id, alpha=0.6), check_overlap = TRUE) 

ord.bray <- ordinate(ad, method="NMDS", distance="bray")
p <- plot_ordination(ad, ord.bray, color="group", shape = "sex", title="NMDS")
p + geom_point(size=3) + 
  theme_classic() +
  geom_text(aes(label=sample_id, alpha=0.6), check_overlap = TRUE) 

# Test barplot
plot_bar(ad, x="sample_id", fill = "Class", title = "PRJNA554111") +
  geom_bar(aes(color = Class, fill = Class), stat = "identity", position = "stack")

####################################################################
#Common core microbiome
####################################################################

library(microbiome)

#Subset phyloseq object into the categories you want to look for the core microbiome in (ex. healthy vs disease state samples)
#Using count data
ad_samples <- subset_samples(ad_counts, group=="control")
control_samples <- subset_samples(ad_counts, group=="AD")

ad_samples
control_samples


# Calculate compositional version of the data
# (relative abundances)
rel.ad_samples <- microbiome::transform(ad_samples, "compositional")
rel.control_samples <- microbiome::transform(control_samples, "compositional")


#names of core microbiome:
#setting to found in 70%
core.taxa.standard.ad_samples <- core_members(rel.ad_samples, detection = 0, prevalence = 70/100)
core.taxa.standard.control_samples <- core_members(rel.control_samples, detection = 0, prevalence = 70/100)

core.taxa.standard.ad_samples
core.taxa.standard.control_samples


####################################################################
#### Alpha diversity Kruskal Wallis test
####################################################################
alpha_ad<- estimate_richness(ad_counts, measures = c("Shannon", "Simpson", "Observed"))
alpha_ad$group <- samdf$group
alpha_ad$age <- samdf$age
alpha_ad$sex <- samdf$sex

write.csv(alpha_ad, file="alpha_div_ad.csv")

var.test(alpha_ad$Observed ~ alpha_ad$group) #F-test to check for difference in variances. If variances are not significantly different, use argument "var.equal=TRUE"
t.test(alpha_ad$Observed ~ alpha_ad$group, var.equal=TRUE)

Kruskal_Simpson<-kruskal.test(Simpson ~ group, data = alpha_ad)
Kruskal_Simpson

Kruskal_Shannon<-kruskal.test(Shannon ~ group, data = alpha_ad)
Kruskal_Shannon

Kruskal_Observed<-kruskal.test(Observed ~ group, data = alpha_ad)
Kruskal_Observed


