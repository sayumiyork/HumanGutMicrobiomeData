#' Data of Sonnenburg Lab's FeFiFo study in normalized proportions
#'
#'  Relative abundance data for fecal samples representing the gut microbiome of subjects on fiber-rich and fermented food diets
#'
#' @format A phyloseq object:
#' \describe{
#'   \item{otu_table}{Relative abundance of ASVs per each sample. 2134 ASVs are represented across 311 samples}
#'   \item{sample_data}{Metadata table connecting metadata to each sample. Contains 311 samples by 27 sample variables}
#'   \item{tax_table}{A table listing the taxonomic assignment of each ASV. Contains 2134 taxa across 7 taxonomic ranks}
#' }
#' @source \url{https://github.com/SonnenburgLab/fiber-fermented-study}
"fefifo"


#' Data of Sonnenburg Lab's FeFiFo study in counts
#'
#'  Count data for fecal samples representing the gut microbiome of subjects on fiber-rich and fermented food diets
#'
#' @format A phyloseq object:
#' \describe{
#'   \item{otu_table}{Count data of ASVs per each sample. 2134 ASVs are represented across 311 samples}
#'   \item{sample_data}{Metadata table connecting metadata to each sample. Contains 311 samples by 27 sample variables}
#'   \item{tax_table}{A table listing the taxonomic assignment of each ASV. Contains 2134 taxa across 7 taxonomic ranks}
#' }
#' @source \url{https://github.com/SonnenburgLab/fiber-fermented-study}
"fefifo_counts"


#' Data of Terrazas et al.'s 2025 study on the effect of oral contraceptives on the gut microbiome in normalized proportions
#'
#'  Relative abundance data for fecal samples representing the gut microbiome of subjects on various dosages including none of oral contraceptives
#'
#' @format A phyloseq object:
#' \describe{
#'   \item{otu_table}{Relative abundance of ASVs per each sample. 5457 ASVs are represented across 115 samples}
#'   \item{sample_data}{Metadata table connecting metadata to each sample. Contains 115 samples by 25 sample variables}
#'   \item{tax_table}{A table listing the taxonomic assignment of each ASV. Contains 5457 taxa across 7 taxonomic ranks}
#' }
#' @source \url{https://github.com/fernanda-t/microbiome-oc}
"oc"

#' Data of Terrazas et al.'s 2025 study on the effect of oral contraceptives on the gut microbiome in counts
#'
#'  Count data for fecal samples representing the gut microbiome of subjects on various dosages including none of oral contraceptives
#'
#' @format A phyloseq object:
#' \describe{
#'   \item{otu_table}{Couts of ASVs per each sample. 5457 ASVs are represented across 115 samples}
#'   \item{sample_data}{Metadata table connecting metadata to each sample. Contains 115 samples by 25 sample variables}
#'   \item{tax_table}{A table listing the taxonomic assignment of each ASV. Contains 5457 taxa across 7 taxonomic ranks}
#' }
#' @source \url{https://github.com/fernanda-t/microbiome-oc}
"oc_counts"

#' Data of Terrazas et al.'s 2025 study on the effect of oral contraceptives on the gut microbiome in normalized proportions. Only participants with dosage data are included.
#'
#'  Relative abundance data for fecal samples representing the gut microbiome of subjects on various dosages including none of oral contraceptives. Only participants with dosage data are included.
#'
#' @format A phyloseq object:
#' \describe{
#'   \item{otu_table}{Relative abundance of ASVs per each sample. 5457 ASVs are represented across 73 samples}
#'   \item{sample_data}{Metadata table connecting metadata to each sample. Contains 73 samples by 25 sample variables}
#'   \item{tax_table}{A table listing the taxonomic assignment of each ASV. Contains 5457 taxa across 7 taxonomic ranks}
#' }
#' @source \url{https://github.com/fernanda-t/microbiome-oc}
"oc_dosage"


#' Data of Terrazas et al.'s 2025 study on the effect of oral contraceptives on the gut microbiome in counts. Only participants with dosage data are included.
#'
#'  Count data for fecal samples representing the gut microbiome of subjects on various dosages including none of oral contraceptives. Only participants with dosage data are included.
#'
#' @format A phyloseq object:
#' \describe{
#'   \item{otu_table}{Counts of ASVs per each sample. 5457 ASVs are represented across 73 samples}
#'   \item{sample_data}{Metadata table connecting metadata to each sample. Contains 73 samples by 25 sample variables}
#'   \item{tax_table}{A table listing the taxonomic assignment of each ASV. Contains 5457 taxa across 7 taxonomic ranks}
#' }
#' @source \url{https://github.com/fernanda-t/microbiome-oc}
"oc_dosage_counts"


#' Data of Terrazas et al.'s 2025 study on the effect of oral contraceptives on the gut microbiome in normalized proportions, excluding sample 948N2
#'
#'  Relative abundance data for fecal samples representing the gut microbiome of subjects on various dosages including none of oral contraceptives, excluding sample 948N2
#'
#' @format A phyloseq object:
#' \describe{
#'   \item{otu_table}{Relative abundance of ASVs per each sample. 5457 ASVs are represented across 114 samples}
#'   \item{sample_data}{Metadata table connecting metadata to each sample. Contains 114 samples by 25 sample variables}
#'   \item{tax_table}{A table listing the taxonomic assignment of each ASV. Contains 5457 taxa across 7 taxonomic ranks}
#' }
#' @source \url{https://github.com/fernanda-t/microbiome-oc}
"oc_no948N2"

#' Data of Terrazas et al.'s 2025 study on the effect of oral contraceptives on the gut microbiome in count data excluding sample 948N2
#'
#'  Count data for fecal samples representing the gut microbiome of subjects on various dosages including none of oral contraceptives, excluding sample 948N2
#'
#' @format A phyloseq object:
#' \describe{
#'   \item{otu_table}{Counts of ASVs per each sample. 5457 ASVs are represented across 114 samples}
#'   \item{sample_data}{Metadata table connecting metadata to each sample. Contains 114 samples by 25 sample variables}
#'   \item{tax_table}{A table listing the taxonomic assignment of each ASV. Contains 5457 taxa across 7 taxonomic ranks}
#' }
#' @source \url{https://github.com/fernanda-t/microbiome-oc}
"oc_no948N2_counts"

#' Data of Zhuang et al's 2018 study examining fecal samples representing the gut microbiome of patients with a clinical diagnosis of Alzheimer's disease and age-sex matched controls in normalized proportions. Data were reprocessed in dada2 using Silva database 138.2.
#'
#'  Relative abundance data for fecal samples representing the gut microbiome of subjects with a clinical diagnosis of Alzheimer's disease and age-sex matched controls.
#'
#' @format A phyloseq object:
#' \describe{
#'   \item{otu_table}{Relative abundance of ASVs per each sample. 16359 ASVs are represented across 86 samples}
#'   \item{sample_data}{Metadata table connecting metadata to each sample. Contains 86 samples by 27 sample variables}
#'   \item{tax_table}{A table listing the taxonomic assignment of each ASV. Contains 16359 taxa across 7 taxonomic ranks}
#' }
#' @source \url{https://doi.org/10.3233/JAD-180176}
"ad"

#' Data of Zhuang et al's 2018 study examining fecal samples representing the gut microbiome of patients with a clinical diagnosis of Alzheimer's disease and age-sex matched controls in counts. Data were reprocessed in dada2 using Silva database 138.2.
#'
#'  Count data for fecal samples representing the gut microbiome of subjects with a clinical diagnosis of Alzheimer's disease and age-sex matched controls.
#'
#' @format A phyloseq object:
#' \describe{
#'   \item{otu_table}{Counts of ASVs per each sample. 16359 ASVs are represented across 86 samples}
#'   \item{sample_data}{Metadata table connecting metadata to each sample. Contains 86 samples by 27 sample variables}
#'   \item{tax_table}{A table listing the taxonomic assignment of each ASV. Contains 16359 taxa across 7 taxonomic ranks}
#' }
#' @source \url{https://doi.org/10.3233/JAD-180176}
"ad_counts"



