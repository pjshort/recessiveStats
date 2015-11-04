# get autozygous rates for CNEs and analyse and write out to file

library(Cairo)
library(dplyr)

AUTOZYGOSITY_DIR = "/lustre/scratch113/projects/ddd/users/jm33/autozygosity"
KINSHIP_PATH = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2014-11-04/kinship_and_pca_trios.txt"
FAMILIES_PATH = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2014-11-04/family_relationships.txt"
SIGNIFICANT_PATH = "/nfs/users/nfs_j/jm33/apps/recessiveStats/results/recessive.combined_tests.ver2.txt"

get_consanguinous_probands <- function(path) {
  cohort = read.table(path, sep="\t", header=TRUE, stringsAsFactors=FALSE)
  
  probands = cohort[cohort$king_kinship > 0, "proband_stable_id"]
  
  return(probands)
}

#' read the autozygosity output (which are individual files for each proband)
get_autozygous_regions <- function(path) {
  paths = Sys.glob(file.path(path, "*"))
  
  regions = vector("list", length(paths))
  
  for (pos in 1:length(paths)) {
    temp = read.table(paths[pos], sep="\t", header=TRUE,
                      colClasses=c("character", "character", "numeric", "numeric"))
    regions[[pos]] = temp
  }
  
  regions = data.frame(rbind_all(regions))
  
  return(regions)
}

#' load a dataframe of family IDs and individual IDs
get_families <- function(path) {
  families = read.table(path, sep="\t", header=TRUE, stringsAsFactors=FALSE)
  families = families[, c("family_id", "individual_id")]
  
  return(families)
}

#' for a list of probands, make sure we only count one proband from each family
get_independent_probands <- function(families, probands) {
  
  # find the family IDs for all the probands, then get the subset of probands
  # where the family ID is not a duplicate.
  family_ids = families$family_id[families$individual_id %in% probands]
  independent_probands = probands[!duplicated(family_ids)]
  
  return(independent_probands)
}

get_autozygosity_per_CNE <- function(path, regions, families, CNEs, subset=NULL) {
  
  # CNEs should have chr, start, and end as column names
  
  # Define the cohort size (we can't rely upon the number of unique probands
  # in the regions dataframe, since many probands will not have autozygous
  # regions).
  # probands_n = length(Sys.glob(file.path(path, "*")))
  probands_n = 3071
  
  # if we want to estimate the rates for consanguinous probands only, we take
  # a subset of the autozygous regions
  if (!is.null(subset)) {
    probands_n = length(subset)
    regions = regions[regions$sample_id %in% subset, ]
  }
  
  CNE_ids = paste(CNEs$chr, CNEs$start, CNEs$end, sep = ".")
  
  rates = data.frame(region_id = CNE_ids, chrom=CNEs$chr, count=NA, rate=NA,
                     setNames(replicate(length(unique(regions$sample_id)), FALSE, simplify=F), sort(unique(regions$sample_id))))
  
  # sapply(sort(unique(regions$sample_id)), function(x) rates[[x]] = NA)
  for (pos in 1:nrow(CNEs)) {
    # get the coordinates for the gene
    chrom = CNEs$chr[pos]
    start_pos = CNEs$start[pos]
    end_pos = CNEs$end[pos]
    
    # find the autozygous regions that overlap the gene region
    in_chrom = regions[regions$chrom == chrom, ]
    overlapping = in_chrom[in_chrom$end_pos > start_pos
                           & in_chrom$start_pos < end_pos, ]
    
    probands = unique(overlapping$sample_id)
    probands = get_independent_probands(families, probands)
    rates$count[pos] = length(probands)
    rates$rate[pos] = rates$count[pos]/probands_n
    rates[pos, probands] = TRUE
  }
  
  # remove the chrX rates, since all of the male probands appear autozygous
  rates = rates[rates$chrom != "X", ]
  
  return(rates)
}



# run from some sub-directory within ~/software/CNE
CNEs = read.table("../data/noncoding_regions.txt", sep = "\t", header = TRUE)
CNEs = CNEs[,c("chr", "start", "stop")]
colnames(CNEs) = c("chr", "start", "end")

proband_ancestry = read.table("../data/proband_ethnicity_predictions.txt", sep = "\t", header = TRUE)
eur_probands = proband_ancestry$person_id[proband_ancestry$group == "EUR"]

families = get_families(FAMILIES_PATH)
regions_per_proband = get_autozygous_regions(AUTOZYGOSITY_DIR)

#restrict only to probands with EUR ancestry
regions_per_proband = regions_per_proband[regions_per_proband$sample_id %in% eur_probands]

all_autozygous = get_autozygosity_per_CNE(AUTOZYGOSITY_DIR, regions_per_proband, families, CNEs)

all_autozygous = all_autozygous[, 1:4]
write.table(all_autozygous, file="../data/autozygosity.all_CNEs.EUR_probands.tsv", sep="\t", row.names=FALSE, quote=FALSE)

