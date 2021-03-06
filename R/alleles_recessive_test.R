

#' test for enrichment of inherited variants in the DDD and ExAC datasets
#'
#' @param hgnc HGNC symbol for a gene.
#' @param chrom chromosome that the gene is on.
#' @param biallelic_lof number of probands with inherited biallelic LoF variants
#'        in the gene.
#' @param biallelic_func number of probands with inherited biallelic functional
#'        variants in the gene.
#' @param lof_func number of probands with inherited Lof/Func variants in the gene.
#' @param start start position of region to be investigated (if checking for
#'        gene region defined by chromosome coordinates rather than a HGNC-based
#'        gene region).
#' @param end end position of region to be investigated (if checking for
#'        gene region defined by chromosome coordinates rather than a HGNC-based
#'        gene region).
#' @param probands vector of probands who have inherited recessive variants in
#'     the gene, or NULL.
#' @param cohort_n number of probands in population.
#' @param check_last_base whether to correct missense or synonymous G alleles at
#'     the last base of exons to a LoF consequence.
#' @param autozygous_rate rate of autozygosity within the gene in the probands.
#' @export
#'
#' @return a list of P values from tests using the DDD population, the ExAC
#'     population, under LoF and functional tests.
analyse_inherited_enrichment <- function(hgnc, chrom, biallelic_lof, biallelic_func, lof_func, start=NULL, end=NULL, probands=NULL, cohort_n=3072, check_last_base=FALSE, autozygous_rate=0) {
    
    cat("extracting ddd frequencies\n")
    ddd = try(get_ddd_variants_for_gene(hgnc, chrom, probands, start=start,
        end=end, check_last_base=check_last_base), silent=TRUE)
    if (class(ddd) != "try-error") {
        ddd = get_cumulative_frequencies(ddd)
        ddd = test_enrichment(ddd, biallelic_lof, biallelic_func, lof_func, sum(unlist(cohort_n)), autozygous_rate)
    } else {
        ddd=list(lof=NA, func=NA, biallelic_lof_p=NA, lof_func_p=NA, biallelic_func_p=NA)
    }
    
    cat("extracting ExAC frequencies\n")
    exac = get_exac_variants_for_gene(hgnc, chrom, start=start, end=end, check_last_base=check_last_base)
    exac = get_cumulative_frequencies(exac)
    
    if (!is.list(cohort_n)) {
        exac = test_enrichment(exac[["NFE"]], biallelic_lof, biallelic_func, lof_func, cohort_n, autozygous_rate)
    } else {
        exac = test_enrichment_across_multiple_populations(exac, biallelic_lof, biallelic_func, lof_func, cohort_n, autozygous_rate)
    }
    
    p_values = list(ddd=ddd, exac=exac)
    
    return(p_values)
}

#' test for enrichment of inherited variants in the DDD and ExAC datasets
#'
#' @param hgnc set to NULL
#' @param chrom chromosome that the gene is on.
#' @param rare_homozygotes number of probands homozygous for rare noncoding variants
#' @param rare_compound_hets number of probands with two (possibly compound het) rare vars in region
#' @param start start position of region to be investigated (if checking for
#'        gene region defined by chromosome coordinates rather than a HGNC-based
#'        gene region).
#' @param end end position of region to be investigated (if checking for
#'        gene region defined by chromosome coordinates rather than a HGNC-based
#'        gene region).
#' @param probands vector of probands who have inherited recessive variants in
#'     the gene, or NULL.
#' @param cohort_n number of probands in population.
#' @param check_last_base whether to correct missense or synonymous G alleles at
#'     the last base of exons to a LoF consequence.
#' @param autozygous_rate rate of autozygosity within the gene in the probands.
#' @export
#'
#' @return a list of P values from tests using the DDD population, the ExAC
#'     population, under LoF and functional tests.
analyse_inherited_noncoding_enrichment <- function(chrom, rare_homozygotes, start=NULL, end=NULL, probands=NULL, cohort_n=3072, check_last_base=FALSE, autozygous_rate=0) {
  
  cat("extracting ddd frequencies\n")
  ddd = try(get_ddd_variants_for_gene(hgnc = NULL, chrom, probands, start=start,
                                      end=end, check_last_base=check_last_base), silent=TRUE)
  if (class(ddd) != "try-error") {
    ddd = get_cumulative_frequencies(ddd)
    p_values = biallelic_lof_enrichment(ddd, rare_homozygotes, sum(unlist(cohort_n)), autozygous_rate)
  } else {
    p_values=list(lof=NA, func=NA, biallelic_lof_p=NA, lof_func_p=NA, biallelic_func_p=NA)
  }
  
  return(p_values)
}

#' test for enrichment of inherited variants in the DDD and ExAC datasets
#'
#' @param chrom chromosome that the gene target/CNEs are on
#' @param rare_homozygotes number of probands homozygous for rare noncoding variants
#' @param start vector of start positions of regions to be investigated 
#' @param end vector of end positions of region to be investigated 
#' @param probands vector of probands who have inherited recessive variants in
#'     the gene, or NULL.
#' @param cohort_n number of probands in population.
#' @param check_last_base whether to correct missense or synonymous G alleles at
#'     the last base of exons to a LoF consequence.
#' @param autozygous_rate weighted rate of autozygosity across the set of variants
#' @export
#'
#' @return a list of P values from tests using the DDD population, the ExAC
#'     population, under LoF and functional tests.
analyse_gene_target_CNE_enrichment <- function(chrom, rare_homozygotes, start=NULL, end=NULL, probands=NULL, cohort_n=3072, check_last_base=FALSE, autozygous_rate=0) {
  # re-factoring of analyse_inherited_noncoding_enrichment to group together several CNEs targetting the same gene
  # chrom, start, and end are assumed to be a vector of multiple regions
  cat("extracting ddd frequencies\n")
  
  # probands should be LIST of character vectors with proband ids for probands with rare homs in each region
  
  ddd = get_ddd_variants_for_CNEs(chrom, start, end, probands, check_last_base=check_last_base)
  
  # we can feed the full list of variants pulled from get_ddd_variants_for_gene to get_cumulative frequencies
  # this will give us cumulative frequency of rare vars over all of the CNEs combined
  
  # for small regions, there is a chance that DDD will have no rare variants - in these cases,
  # we may need to add a pseudocount?
  
  if (class(ddd) != "try-error") {
    ddd = get_cumulative_frequencies(ddd)
    p_values = biallelic_lof_enrichment(ddd, rare_homozygotes, sum(unlist(cohort_n)), autozygous_rate)
  } else {
    p_values=list(lof=NA, func=NA, biallelic_lof_p=NA, lof_func_p=NA, biallelic_func_p=NA)
  }
  
  return(p_values)
}

#' test for enrichment of inherited variants
#'
#' @param freq list of cumulative frequencies of variation in a population for
#'        rare LoF variants, and rare functional variants.
#' @param biallelic_lof number of probands with inherited Lof/LoF variants in
#'        the gene.
#' @param biallelic_func number of probands with inherited func/func variants in
#'        the gene.
#' @param lof_func number of probands with inherited Lof/Func variants in the gene.
#' @param cohort_n number of probands in population.
#' @param autozygous_rate rate of autozygosity in the cohort being investigated.
#' @export
#'
#' @return a list of P values from tests, under LoF and functional tests.
#' 
test_enrichment <- function(freq, biallelic_lof, biallelic_func, lof_func, cohort_n, autozygous_rate=0) {
    
    # get the probability of getting more than or equal to the number of
    # observed inherited events
    freq$biallelic_lof_p = biallelic_lof_enrichment(freq, biallelic_lof, cohort_n, autozygous_rate)
    freq$lof_func_p = lof_func_enrichment(freq, biallelic_lof + lof_func, cohort_n, autozygous_rate)
    freq$biallelic_func_p = biallelic_func_enrichment(freq, biallelic_func, cohort_n, autozygous_rate)
    
    return(freq)
}

#' test for enrichment of biallelic LoF inherited variants
#'
#' @param freq list of cumulative frequencies of variation in a population for
#'        rare LoF variants, and rare functional variants.
#' @param count number of probands with inherited variants in the gene.
#' @param cohort_n number of probands in population.
#' @param autozygous_rate rate of autozygosity in the cohort being investigated.
#' @export
#'
#' @return P-value from testing for biallelic LoF variants.
biallelic_lof_enrichment <- function(freq, count, cohort_n, autozygous_rate=0) {
    rate = (freq$lof ** 2) * (1 - autozygous_rate) + freq$lof * autozygous_rate
    return(pbinom(count - 1, cohort_n, prob=rate, lower.tail=FALSE))
}

#' test enrichment of inherited biallelic LoF and compound heterozygous LoF/Func
#'
#' @param freq list of cumulative frequencies of variation in a population for
#'        rare LoF variants, and rare functional variants.
#' @param count number of probands with inherited variants in the gene.
#' @param cohort_n number of probands in population.
#' @param autozygous_rate rate of autozygosity in the cohort being investigated.
#' @export
#'
#' @return P-value from testing for biallelic LoF and Lof/Func variants.
lof_func_enrichment <- function(freq, count, cohort_n, autozygous_rate=0) {
    rate = (freq$lof ** 2) * (1 - autozygous_rate) +
        (freq$lof) * autozygous_rate +
        (2 * freq$lof * (1 - freq$lof) * freq$functional)
    return(pbinom(count - 1, cohort_n, prob=rate, lower.tail=FALSE))
}

#' test enrichment of inherited biallelic functional variants
#'
#' @param freq list of cumulative frequencies of variation in a population for
#'        rare LoF variants, and rare functional variants.
#' @param count number of probands with inherited variants in the gene.
#' @param cohort_n number of probands in population.
#' @param autozygous_rate rate of autozygosity in the cohort being investigated.
#' @export
#'
#' @return P-value from testing for biallelic functional variants.
biallelic_func_enrichment <- function(freq, count, cohort_n, autozygous_rate=0) {
    rate = (freq$functional ** 2) * (1 - autozygous_rate) + freq$functional * autozygous_rate
    return(pbinom(count - 1, cohort_n, prob=rate, lower.tail=FALSE))
}

#' test for enrichment of inherited variants in multiple populations
#'
#' @param exac list of frequency estimates for each ExAC population.
#' @param biallelic_lof number of probands with inherited biallelic LoF variants in the gene.
#' @param biallelic_func number of probands with inherited biallelic functional variants in the gene.
#' @param lof_func number of probands with inherited Lof/Func variants in the gene.
#' @param cohort_n list of number of probands in each population.
#' @param autozygous_rate rate of autozygosity in the cohort being investigated.
#' @export
#'
#' @return a list of P values from tests across the populations, under LoF and
#'         functional tests.
test_enrichment_across_multiple_populations <- function(exac, biallelic_lof, biallelic_func, lof_func, cohort_n, autozygous_rate) {
    # define the ExAC different populations (AFR="African/African American",
    # EAS="East Asian", NFE="Non-Finnish European", SAS="South Asian")
    populations = names(cohort_n)
    biallelic_lof_combos = get_count_combinations(populations, biallelic_lof)
    biallelic_func_combos = get_count_combinations(populations, biallelic_func)
    lof_func_combos = get_count_combinations(populations, biallelic_lof + lof_func)
    
    p_values = list(lof=NA, functional=NA)
    p_values$biallelic_lof_p = sum_combo_tests(exac, cohort_n,
        biallelic_lof_combos, biallelic_lof_enrichment, autozygous_rate)
    p_values$lof_func_p = sum_combo_tests(exac, cohort_n,
        lof_func_combos, lof_func_enrichment, autozygous_rate)
    p_values$biallelic_func_p = sum_combo_tests(exac, cohort_n,
        biallelic_func_combos, biallelic_func_enrichment, autozygous_rate)
    
    return(p_values)
}

#' get a p-value that sums across different population possibilities
#'
#' @param exac list of frequency estimates for each ExAC population.
#' @param cohort_n list of number of probands in each population.
#' @param combos a dataframe of the possible count combinations for a functional
#'        type.
#' @param enrich_function function to test enrichment.
#' @param autozygous_rate rate of autozygosity in the cohort being investigated.
#' @export
#'
#' @return a p-value from testing for
sum_combo_tests <- function(exac, cohort_n, combos, enrich_function, autozygous_rate=0) {
    summed_p_value = 0
    for (pos in 1:nrow(combos)) {
        row_p_value = 1
        for (pop in names(cohort_n)) {
            # get the functional variants for the population from the row
            count = combos[[pop]][pos]
            
            # For populations where the count is not the highest in the row, we
            # want the probability of the population having that count families.
            # Otherwise we want the probability of having that count or greater.
            if (count == max(combos[pos, ])) {
                p_value = enrich_function(exac[[pop]], count, cohort_n[[pop]], autozygous_rate)
            } else {
                inclusive = enrich_function(exac[[pop]], count, cohort_n[[pop]], autozygous_rate)
                exclusive = enrich_function(exac[[pop]], count + 1, cohort_n[[pop]], autozygous_rate)
                p_value = inclusive - exclusive
            }
            
            # the p-value for the row is the product of the p-values for every
            # population
            row_p_value = row_p_value * p_value
        }
        # the overall p-value is the sum of p-values for each row
        summed_p_value = summed_p_value + row_p_value
    }
    
    return(summed_p_value)
}

#' Get all the combinations of spreading the probands across populations.
#'
#' @param populations a vector of population names to be tested
#' @param count number of probands with inherited variants in the gene.
#' @export
#'
#' @return a list of count dataframes for each functional type
get_count_combinations <- function(populations, count) {
    
    stopifnot(count >= 0)
    
    # get a matrix of count combinations
    combos = expand.grid(rep(list(seq(0, count)), each=length(populations)))
    
    n_parity = ceiling(count/length(populations))
    # make sure each of the rows sums to the correct value, so that we only use
    # rows where the counts are dispersed correctly amongst the populations
    combos = combos[rowSums(combos) == count |
        (rowSums(combos) > count & apply(combos, 1, max) <= n_parity), ]
    
    # tidy up the dataframe so that it is of a standard format
    combos = data.frame(combos)
    row.names(combos) = 1:nrow(combos)
    
    # name the combinations by the population names
    names(combos) = populations
    
    return(combos)
}
