library(dplyr)
library(tidyr)

#' read in a haplot RDS and call genotypes
#' 
#' @param path the path to the RDS file from microhaplot
#' @param min_depth1 required minimum read depth for the most prevalent read.  This is the required depth of the rank 1 allele.
#' @param min_depth2 required minimum read depth for the rank-2 allele.
#' @param min_balance required minimum allele balance
#' @details Here is how this works.  It first filters things down to just rank1 and rank2 (if present) haplotypes.
#' (The rank1 haplotype is the highest read depth haplotype).  Then it tosses out any rank2 haplotypes if their
#' read depth divided by the read depth of the rank1 allele is less than min_balance.  After that, loci are not
#' called if, either: 1) the rank1 allele has read depth less than min_depth1, or 2) the rank2 allele (if present)
#' has a read_depth less than min_depth2.  For example,  if you had a rank1 allele read depth of 12 and a rank1 allele
#' depth of 7, you probably want to call that as a heterozygote, but you couldn't get that with a single
#' read depth criterion of 10, say.
#' @return This returns a data frame with columns: id, locus, gene_copy, allele, depth, allele.balance.  It should be noted that 
#' homozygotes are reported to have the same depth for gene_copy1 and gene_copy2.  Missing genotypes are explicitly NAs.
call_genos_from_haplotRDS <- function(path, min_depth1 = 10, min_depth2 = 6, min_balance = 0.4) {
  
  rds <- readRDS(path) %>%
    tbl_df() %>%
    select(-sum.Phred.C, -max.Phred.C, -group) %>%
    filter(rank <= 2) %>%
    arrange(id, locus, rank) %>%
    filter(allele.balance >= min_balance)      # filter on allele balance.
    
  
  # now, some stuff to toss all individual x locus combos that
  # have a rank2 or a rank1 read depth of less than 10.  So, I don't want to
  # call it if has 10 reads of one allele and 3 of another.  Better to just leave
  # it uncalled, I think
  rds2 <- rds %>%
    group_by(id, locus) %>%
    filter((n()==2 && depth[2] < min_depth2) + (n() == 1 && depth[1] < min_depth1) == 0) %>%
    ungroup()
  
  # now we want to fill that out so everyone has a rank2 row, but it is explicitly NA for our homozygotes.
  # I should be able to do this with tidyr::complete, but that is always baffling to me.
  rds3 <- expand.grid(id = unique(rds$id), locus = unique(rds$locus), rank = 1:2, stringsAsFactors = FALSE) %>%
    tbl_df() %>%
    arrange(id, locus, rank) %>%
    left_join(., rds2)

  # and now we assign gene_copy1 and gene_copy2
  rds4 <- rds3 %>%
    group_by(id, locus) %>%
    mutate(allele = ifelse(is.na(haplo), haplo[1], haplo),
           depth = ifelse(is.na(depth), depth[1], depth)) %>%
    rename(gene_copy = rank) %>%
    select(id, locus, gene_copy, allele, depth, allele.balance) %>%
    ungroup()
  
  rds4
}
