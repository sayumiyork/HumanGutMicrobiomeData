
#load libraries
library(phyloseq)

# Load Data #################

#FeFiFo
# Wastyk, Hannah C., et al. "Gut-microbiota-targeted diets modulate human immune status." Cell 184.16 (2021): 4137-4153.
#Extracted from files on GitHub
#https://github.com/SonnenburgLab/fiber-fermented-study
fefifo <- readRDS("data-raw/fefifo.RDS") #Normalized by proportions
fefifo_counts <- readRDS("data-raw/fefifo_counts.RDS") #Raw counts

#Oral contraceptive
# Terrazas, Fernanda, et al. "Influence of menstrual cycle and oral contraception on taxonomic composition and gas production in the gut microbiome." Journal of Medical Microbiology 74.3 (2025): 001987.
oc <- readRDS("data-raw/oc.RDS") #Normalized by proportions
oc_counts <- readRDS("data-raw/oc_counts.RDS") #Raw counts
oc_dosage <- readRDS("data-raw/oc_dosage.RDS") #Normalized by proportions, only samples with dosage data
oc_dosage_counts <- readRDS("data-raw/oc_dosage_counts.RDS") #Raw counts, only samples with dosage data
oc_no948N2  <- readRDS("data-raw/oc_no948N2.RDS") #Normalized by proportions, no sample 948N2, identified as an outlier in NMDS plot
oc_no948N2_counts <- readRDS("data-raw/oc_no948N2_counts.RDS") #Raw counts,  no sample 948N2, identified as an outlier in NMDS plot
oc_agefirstOCuse <- readRDS("data-raw/oc_agefirstOCuse.RDS") #Normalized by proportions, only samples with agefirstOCuse data
oc_agefirstOCuse_counts <- readRDS("data-raw/oc_agefirstOCuse_counts.RDS") #Raw counts, only samples with agefirstOCuse data

#Alzheimer's disease
#Zhuang, Zhen-Qian, et al. "Gut microbiota is altered in patients with Alzheimer’s disease." Journal of Alzheimer’s Disease 63.4 (2018): 1337-1346.
ad <- readRDS("data-raw/ad.RDS") #Normalized by proportions
ad_counts <- readRDS("data-raw/ad_counts.RDS") #Raw counts


#Use this ############################

usethis::use_data(fefifo, overwrite = TRUE)
usethis::use_data(fefifo_counts, overwrite = TRUE)

usethis::use_data(oc, overwrite = TRUE)
usethis::use_data(oc_counts, overwrite = TRUE)
usethis::use_data(oc_dosage, overwrite = TRUE)
usethis::use_data(oc_dosage_counts, overwrite = TRUE)
usethis::use_data(oc_no948N2, overwrite = TRUE)
usethis::use_data(oc_no948N2_counts, overwrite = TRUE)
usethis::use_data(oc_agefirstOCuse, overwrite = TRUE)
usethis::use_data(oc_agefirstOCuse_counts, overwrite = TRUE)

usethis::use_data(ad, overwrite = TRUE)
usethis::use_data(ad_counts, overwrite = TRUE)
