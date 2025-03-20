



#' Generate well x gene level data for a pooled screen
#'
#' Given assay parameters run the following simulation
#'    1) Select genotypes from the library as hits/non-hits
#'    2) For each well and genotype pair
#'        - sample the count of cells with that genotype that are in that well
#'        - and among those cells the number are classified as positive.
#'
#' The count for genoptype i in well j is sampled from the following binomial
#' distribution
#'
#'      count_ij ~ binomial(
#'          size = gene_well_dispersion,                           # check this
#'          prob = (n_cells_per_well / n_genes_per_library) *
#'               gene_abundance_factor * well_abundance_factor)
#'
#' where the gene and well abundance factors control scale expected abundance
#' based on the geneotype or well level effects, the gene abundance factor for
#' geneotype i is sampled from the following gamma distribution:
#'
#'      gene_abundance_factor_i ~ gamma(
#'        shape = gene_abundance_factor_mu^2 / gene_abundance_factor_var,
#'        rate = gene_abundance_factor_mu / gene_abundance_factor_var)
#'
#' where gene_abundance_factor_mu and gene_abundance_factor_var control the variation
#' across different geneoptyes. Similarly, the abundance factor for well j is
#' sampled from the following gamma distribution:
#'
#'      well_abundance_factor_i ~ gamma(
#'        shape = well_abundance_factor_mu^2 / well_abundance_factor_var,
#'        rate = well_abundance_factor_mu / well_abundance_factor_var)
#'
#' The number of positive cells for geneotype i in well j is sampled from a binomial
#' distribution is defined by
#'
#'    positive_ij ~ binom(
#'      size = count_ij,
#'      prob = ifelse(hit, pos_hit_rate, neg_hit_rate))
#'
#' where the pos_hit_rate and neg_hit_rate model the accuracy of classifying a
#' cell as positive or negative given if it has a hit geneotype or not, these are
#' defined are sampled from beta distributions defined by the mean and standard
#' deviation of the beta distirbution and then translated into shape parameters
#'
#'    pos_hit_rate_ij ~ beta(
#'          shape1 = class_pos_mu * ((class_pos_mu * (1 - class_pos_mu) / class_pos_var) - 1),
#'          shape2 = (1 - class_pos_mu) * ((class_pos_mu * (1 - class_pos_mu) / class_pos_var) - 1))
#'
#'    neg_hit_rate_ij ~ beta(
#'          n = dplyr::n(),
#'          shape1 = class_neg_mu * ((class_neg_mu * (1 - class_neg_mu) / class_neg_var) - 1),
#'          shape2 = (1 - class_neg_mu) * ((class_neg_mu * (1 - class_neg_mu) / class_neg_var) - 1)))))
#'
#'
#' @param n_genes_per_library integer number of genes in the library
#' @param gene_abundance_factor_mu numeric
#' @param gene_abundance_factor_var numeric
#' @param gene_hit_rate numeric
#' @param n_wells_per_screen integer
#' @param well_abundance_factor_mu numeric
#' @param well_abundance_factor_var numeric
#' @param n_genes_per_well integer
#' @param n_cells_per_well_mu numeric
#' @param n_cells_per_well_var numeric
#' @param class_pos_mu numeric
#' @param class_pos_var numeric
#' @param class_neg_mu numeric
#' @param class_neg_var numeric
#'
#' @export
simulate_screen <- function(
    n_genes_per_library,
    gene_abundance_factor_mu,
    gene_abundance_factor_var,
    gene_hit_rate,
    n_wells_per_screen,
    well_abundance_factor_mu,
    well_abundance_factor_var,
    n_genes_per_well,
    n_cells_per_well_mu,
    n_cells_per_well_var,
    class_pos_mu,
    class_pos_var,
    class_neg_mu,
    class_neg_var) {

  # sanity_checks
  assertthat::assert_that(
    n_genes_per_well <= n_genes_per_library,
    msg = "At most as genes in a well in than there are genes in the library")

  assertthat::assert_that(
    n_genes_per_library > 0,
    msg = "Positive number of genes in library")

  assertthat::assert_that(
    gene_hit_rate >= 0 && gene_hit_rate <= 1,
    msg = "gene_hit_rate is in [0,1]")

  assertthat::assert_that(
    n_genes_per_well > 0,
    msg = "Target number of genes per well is positive")


  gene_covariates <- data.frame(
    gene = 1:n_genes_per_library,

    # mean = shape / rate
    # variance = shape / rate^2
    gene_abundance_factor = rgamma(
      n = n_genes_per_library,
      shape = gene_abundance_factor_mu^2 / gene_abundance_factor_var,
      rate = gene_abundance_factor_mu / gene_abundance_factor_var),
    hit = rbinom(
      n = n_genes_per_library,
      size = 1,
      prob = gene_hit_rate))

  well_covariates <- data.frame(
    well = 1:n_wells_per_screen,
    well_abundance_factor = rgamma(
      n = n_wells_per_screen,
      shape = well_abundance_factor_mu^2 / well_abundance_factor_var,
      rate = well_abundance_factor_mu / well_abundance_factor_var))

  well_gene <- tidyr::expand_grid(
    gene = 1:n_genes_per_library,
    well = 1:n_wells_per_screen) |>
    dplyr::mutate(
      # add parameters
      n_genes_per_library = n_genes_per_library,
      gene_abundance_factor_mu = gene_abundance_factor_mu,
      gene_abundance_factor_var = gene_abundance_factor_var,
      gene_hit_rate = gene_hit_rate,
      n_wells_per_screen = n_wells_per_screen,
      well_abundance_factor_mu = well_abundance_factor_mu,
      well_abundance_factor_var = well_abundance_factor_var,
      n_genes_per_well = n_genes_per_well,
      n_cells_per_well_mu = n_cells_per_well_mu,
      n_cells_per_well_var = n_cells_per_well_var,
      class_pos_mu = class_pos_mu,
      class_pos_var = class_pos_var,
      class_neg_mu = class_neg_mu,
      class_neg_var = class_neg_var) |>
    dplyr::left_join(gene_covariates, by = "gene") |>
    dplyr::left_join(well_covariates, by = "well") |>
    dplyr::mutate(
      gene_in_well = rbinom(
        n = dplyr::n(),
        size = 1,
        prob = pmin((n_genes_per_well / n_genes_per_library) * gene_abundance_factor, 1)),
      n_cells_per_gene_per_well_mu = well_abundance_factor *
        n_cells_per_well_mu / n_genes_per_well,
      n_cells_per_gene_per_well_var = well_abundance_factor *
        n_cells_per_well_var / n_genes_per_well,
      count = gene_in_well *
        rnbinom(
          n = dplyr::n(),
          size = (n_cells_per_gene_per_well_mu)^2 /
            (n_cells_per_gene_per_well_var - n_cells_per_gene_per_well_mu),
          prob = n_cells_per_gene_per_well_mu / n_cells_per_gene_per_well_var),

      positive = rbinom(
        n = dplyr::n(),
        size = count,
        prob = ifelse(
          hit,
          rbeta(
            n = dplyr::n(),
            shape1 = class_pos_mu * ((class_pos_mu * (1 - class_pos_mu) / class_pos_var) - 1),
            shape2 = (1 - class_pos_mu) * ((class_pos_mu * (1 - class_pos_mu) / class_pos_var) - 1)),
          rbeta(
            n = dplyr::n(),
            shape1 = class_neg_mu * ((class_neg_mu * (1 - class_neg_mu) / class_neg_var) - 1),
            shape2 = (1 - class_neg_mu) * ((class_neg_mu * (1 - class_neg_mu) / class_neg_var) - 1)))))
}
