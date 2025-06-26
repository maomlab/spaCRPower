


#' Parameterize the gamma distribution by the mean and variance
#'
#' For a gamma distribution parameterized by the shape and rate is
#'
#'     mean = shape / rate
#'     variance = shape / rate^2
#'
#' solving both for alpha and setting them equal gives
#'
#'    shape = mean*rate = variance*rate^2
#'
#' solving for rate in terms of mean and variance gives
#'
#'    rate = mean / variance
#'
#' substituting that back in to solve for shape gives
#'
#'    shape = mean^2 / variance
#'
rgamma_mean_variance <- function(n, mean, var) {
  rgamma(
    n = n,
    shape = mean^2 / var,
    rate = mean / var)
}

rnbinom_mean_variance <- function(n, mean, var) {
  assertthat::assert_that(
    all(0 < mean),
    msg = paste0(
      "The mean must be positive, but the requested mean is ", mean))

  assertthat::assert_that(
    all(var >= mean),
    msg = paste0(
      "The variance must be between greater than the  mean, ",
      "but the requested variance is ", var))

  rnbinom(
    n = n,
    size = mean^2 / (var - mean),
    prob = mean / var)
}

rbeta_mean_variance <- function(n, mean, var) {
  assertthat::assert_that(
    all(0 <= mean & mean <= 1),
    msg = paste0(
      "The mean must be in [0, 1], but the requested mean is ", mean))

  assertthat::assert_that(
    all(0 < var & var < mean*(1-mean)),
    msg = paste0(
      "The variance must be [0, mean*(1-mean)] = [0, ", mean*(1-mean), "], ",
      "but the requested variance is ", var))

  rbeta(
    n = n,
    shape1 = mean * ((mean * (1 - mean) / var) - 1),
    shape2 = (1 - mean) * ((mean * (1 - mean) / var) - 1))
}

anscombe_mean <- function(n, mean) {
  rpois(
    n = n,
    lambda = (mean / 2)^2 - 3/8)
}


#' Simulate a library of gene perturbations
#'
#' Given the library size, sample the abundance and hit status for each gene in
#' the library:
#'
#'    abundance ~ dirichlet(gene_abundance_alpha*n)
#'    for i in 1..n:
#'      hit_i ~ binom(gene_hit_rate)
#'
#' For the gene_abundance_alpha parameter larger values have abundances that are
#' are more evenly distributed
#'   alpha -> Infinity: abundance_i -> 1/n
#'   alpha -> 0:
#'     k ~ sample(1, n)
#'     abundance_i = Indicator(i = k)
#' When alpha is 1, then the Gini index = 0.5.
#'
#'
#' @param n_genes_in_library integer number of genes in the library
#' @param gene_abundance_alpha numeric positive real number giving the
#'   how biased the gene abundance is across the library
#' @param gene_hit_rate numeric probability of each gene being a hit
#' @returns data.frame with `n_genes_in_library` rows and columns
#'   `\[gene, gene_abundance_alpha, gene_hit_rate, gene_abundance, hit\]`, where
#'   `gene` is 1-based index for the gene, `gene_abundance_alpha`, and
#'   `gene_hit_rate` record the simulation parameters, `gene_abundance` for gene
#'   i is the fraction of the library represented by the gene i, such that the
#'   sum of the `gene_abundance` values sum to `1`, and for gene i `hit` is
#'   whether or not gene i is a hit.
#'
#' @examples
#' \dontrun{
#'   scores <- spaCRPower::simulate_library(
#'     n_genes_in_library = 200,
#'     gene_abundance_alpha = 100,
#'     gene_hit_rate = .05)
#' }
#'
#' @export
simulate_library <- function(
  n_genes_in_library,
  gene_abundance_alpha,
  gene_hit_rate) {

  assertthat::assert_that(
    n_genes_in_library > 0,
    msg = "Positive number of genes in library")

  assertthat::assert_that(
    gene_abundance_alpha > 0,
    msg = "gene_abundance_alpha should be positive")

  assertthat::assert_that(
    gene_hit_rate >= 0 && gene_hit_rate <= 1,
    msg = "gene_hit_rate is in [0,1]")

  data.frame(
    gene = 1:n_genes_in_library,
    gene_abundance_alpha = gene_abundance_alpha,
    gene_hit_rate = gene_hit_rate,

    gene_abundance = extraDistr::rdirichlet(
      n = 1,
      alpha = rep_len(
        x = gene_abundance_alpha,
        length.out = n_genes_in_library))[1,],
    hit = rbinom(
      n = n_genes_in_library,
      size = 1,
      prob = gene_hit_rate))
}


#' Simulate a spot plate
#'
#' Given the number of wells and parameters for the well abundance, simulate
#' whether or not gene i is in well j for all i and j.
#'
#'    well_abundance_j ~ gamma(
#'      mean = well_abundance_factor_mu,
#'      var = well_abundance_factor_var)
#'    gene_in_well_ij ~ binom(gene_abundance_i * well_abundance_j)
#'
#' Note we use `rgamma_mean_variance()` to get the gamma distribution
#' parameterized by the mean and variance.
#'
#'
#' @param gene_library data.frame result of calling `simulate_library()`
#' @param n_wells_per_screen integer number of wells in the screen
#' @param well_abundance_facotr_mu numeric mean of the per-well abundance
#'  factor (must be greater than 0).
#' @param well_abundance_factor_var numeric variance of the per-well abundance
#'  factor (must be greater than 0).
#' @returns data.frame with one row for each gene, well pair and columns
#'  `\[gene, well, <gene_covariates>, <well_covariates>, gene_in_well\]`, where
#'  the `gene` and `well` columns are 1-based indices for each gene and well,
#'  `<gene_covariates>` are columns from the inpute `gene_library`,
#'  `<well_covariates>` are the `well_abundance_factor_mu` parameters
#'  `well_abundance_factor_var`, and the sampled `well_abundance` value.
#'
#'
#' @export
simulate_spot_plate <- function(
  gene_library,
  n_wells_per_screen,
  well_abundance_factor_mu,
  well_abundance_factor_var) {

  assertthat::assert_that(
    well_abundance_factor_mu > 0,
    msg = "well_abundance_factor_mu should be positive")

  assertthat::assert_that(
    well_abundance_factor_var > 0,
    msg = "well_abundance_factor_var should be positive")

  well_covariates <- data.frame(
    well = 1:n_wells_per_screen,
    well_abundance_factor_mu = well_abundance_factor_mu,
    well_abundance_factor_var = well_abundance_factor_var,
    well_abundance = rgamma_mean_variance(
      n = n_wells_per_screen,
      mean = well_abundance_factor_mu,
      var = well_abundance_factor_var))

  tidyr::expand_grid(
    gene = 1:nrow(gene_library),
    well = 1:nrow(well_covariates)) |>
    dplyr::left_join(gene_library, by = c("gene")) |>
    dplyr::left_join(well_covariates, by = c("well")) |>
    dplyr::mutate(
      gene_in_well = rbinom(
        n = dplyr::n(),
        size = 1,
        prob = gene_abundance * well_abundance))
}

#' Simulate an imaging plate
#'
#'
#'
#'
#' @param spotplate data.frame results of calling `simulate_spot_plate()`
#' @param n_cells_per_well_lambda numeric mean of the number of cells for each
#'   gene in each well (must be greater than 0)
#' @param n_cells_per_well_nu numeric dispersion of the number of wells for each
#'   per-well number of cells
#' @param class_pos_mu numeric mean of the probability of being called positive
#'   for cells that are hits
#' @param class_pos_var numeric variance of the probability of being called
#'   positive for cells that are hits
#' @param class_neg_mu numeric mean of the probability of being called positive
#'   for cells that are not hits
#' @param class_neg_var numeric variance of the probability of being called
#'   positive for cells that are not hits
#' @returns data.frame
#'
#'
#' @export
simulate_imaging_plate <- function(
  spot_plate,
  imaging_n_cells_per_well_mu,
  imaging_n_cells_per_well_var,
  class_pos_mu,
  class_pos_var,
  class_neg_mu,
  class_neg_var) {

  assertthat::assert_that(
    imaging_n_cells_per_well_mu > 0,
    msg = "imaging_n_cells_per_well_lambda should be positive")

  assertthat::assert_that(
    imaging_n_cells_per_well_var > imaging_n_cells_per_well_mu,
    msg = paste0(
      "imaging_n_cells_per_well_var should be greater than ",
      "imaging_n_cells_per_well_mu"))

  assertthat::assert_that(
    class_pos_mu >= 0 && class_pos_mu <= 1,
    msg = "class_pos_mu should be in [0, 1]")

  assertthat::assert_that(
    class_neg_mu >= 0 && class_neg_mu <= 1,
    msg = "class_neg_mu should be in [0, 1]")

  imaging_plate <- spot_plate |>
    dplyr::group_by(well) |>
    dplyr::do({
      spot_well <- .
      if (sum(spot_well$gene_in_well) == 0){
        spot_well |>
          dplyr::mutate(imaging_n_cells_per_gene_per_well = 0)
      } else {
        spot_well |>
          dplyr::mutate(
            imaging_n_cells_per_gene_per_well =
              stats::rmultinom(
                n = 1,
                size = rnbinom_mean_variance(
                  n = 1,
                  mean = imaging_n_cells_per_well_mu,
                  var = imaging_n_cells_per_well_var),
                prob = gene_in_well) |>
              as.numeric())
      }
    }) |>
    dplyr::ungroup()

  imaging_plate |>
    dplyr::transmute(
      gene,
      well,
      imaging_n_cells_per_well_mu = imaging_n_cells_per_well_mu,
      imaging_n_cells_per_well_var = imaging_n_cells_per_well_var,
      imaging_n_cells_per_gene_per_well = imaging_n_cells_per_gene_per_well,
      class_pos_mu = class_pos_mu,
      class_pos_var = class_pos_var,
      class_neg_mu = class_neg_mu,
      class_neg_var = class_neg_var,
      positive = rbinom(
        n = dplyr::n(),
        size = imaging_n_cells_per_gene_per_well,
        prob = ifelse(
          hit,
          rbeta_mean_variance(
            n = dplyr::n(), mean = class_pos_mu, var = class_pos_var),
          rbeta_mean_variance(
            n = dplyr::n(), mean = class_neg_mu, var = class_neg_var))))
}


#' @export
simulate_sequencing_plate <- function(
  spot_plate,
  sequencing_n_cells_per_well_lambda,
  sequencing_n_cells_per_well_nu,
  pcr_factor_mu,
  pcr_factor_var,
  n_reads_total) {

  assertthat::assert_that(
    sequencing_n_cells_per_well_lambda > 0,
    msg = "sequencing_n_cells_per_well_lambda should be positive")

  assertthat::assert_that(
    sequencing_n_cells_per_well_nu > 0,
    msg = "sequencing_n_cells_per_well_nu should be positive")

  well_summary <- spot_plate |>
    dplyr::group_by(well) |>
    dplyr::summarize(
      n_genes_in_well = sum(gene_in_well),
      .groups = "drop") |>
    dplyr::mutate(
      pcr_factor = rlnorm(
        n = dplyr::n(),
        meanlog = pcr_factor_mu,
        sdlog = sqrt(pcr_factor_var)))

  sequencing_plate <- spot_plate |>
    dplyr::left_join(well_summary, by = "well")

  # if the dispersion parameter is 1, then this is the poisson distribution
  # so sample from rpois, because it is fast, and sample from the
  # COM-Poisson distribution otherwise
  if (sequencing_n_cells_per_well_nu == 1) {
    sequencing_plate <- sequencing_plate |>
      dplyr::mutate(
        sequencing_n_cells_per_gene_per_well = gene_in_well *
          rpois(
            n = dplyr::n(),
            lambda = sequencing_n_cells_per_well_lambda))
  } else {
    sequencing_plate <- sequencing_plate |>
      dplyr::mutate(
        sequencing_n_cells_per_gene_per_well = gene_in_well *
          COMPoissonReg::rcmp(
            n = dplyr::n(),
            lambda = sequencing_n_cells_per_well_lambda,
            nu = sequencing_n_cells_per_well_nu))
  }

  sequencing_plate <- sequencing_plate |>
    dplyr::group_by(well) |>
    dplyr::do({
      well_data <- .
      n_cells_in_well <- sum(well_data$sequencing_n_cells_per_gene_per_well)

      n_reads_per_well <- round(n_reads_total / nrow(well_data))

      if (sum(well_data$sequencing_n_cells_per_gene_per_well) > 0) {
        n_barcodes_per_genes_per_well <- round(
          well_data$sequencing_n_cells_per_gene_per_well * well_data$pcr_factor)
        n_reads_per_gene_per_well <- extraDistr::rmvhyper(
          nn = 1,
          n = n_barcodes_per_genes_per_well,
          k = min(n_cells_in_well * well_data$pcr_factor[1], n_reads_per_well)) |>
          as.numeric()
      } else {
        n_reads_per_gene_per_well <- rep(0, length.out = nrow(well_data))
      }

      tibble::tibble(
        well = well_data$well,
        gene = well_data$gene,
        sequencing_n_cells_per_well_lambda = sequencing_n_cells_per_well_lambda,
        sequencing_n_cells_per_well_nu = sequencing_n_cells_per_well_nu,
        n_reads_total,
        pcr_factor = well_data$pcr_factor,
        sequencing_n_cells_per_gene_per_well =
          well_data$sequencing_n_cells_per_gene_per_well,
        n_barcodes_per_genes_per_well = n_barcodes_per_genes_per_well,
        n_reads_per_gene_per_well = n_reads_per_gene_per_well)
    }) |>
    dplyr::ungroup()

}

#' Generate well x gene level data for a pooled screen
#'
#' Given assay parameters run the following simulation
#'    1) Select genotypes from the library as hits/non-hits
#'    2) For each well and genotype pair sample the count of cells with that
#'       genotype that are in that well
#'    3) For each cell, based on its genotype, determine if it is classified
#'       as positive or negative
#'    4) For each well, determine the read counts of each geneotype
#'
#' The count for genoptype i in well j is sampled from the following binomial
#' distribution
#'
#'    gene_well_count_ij ~ binomial(
#'        size = gene_well_dispersion,                           # check this
#'        prob = (n_cells_per_well / n_genes_in_library) *
#'             gene_abundance_factor * well_abundance_factor)
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
#' @param n_genes_in_library integer number of genes in the library
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
#'
#' @examples
#' \dontrun{
#'   scores <- spaCRPower::scan_parameters(
#'     n_genes_in_library = 200,
#'     gene_abundance_factor_mu = 1,
#'     gene_abundance_factor_var = .01,
#'     gene_hit_rate = 0.05,
#'     n_wells_per_screen = seq(30, 200, length.out = 35),
#'     well_abundance_factor_mu = 1,
#'     well_abundance_factor_var = .01,
#'     n_genes_per_well = c(1, 3, 10, 30),
#'     n_cells_per_well_mu = 850,
#'     n_cells_per_well_var = 1000,
#'     class_pos_mu = 0.99,
#'     class_pos_var = 0.0001,
#'     class_neg_mu = 0.01,
#'     class_neg_var = 0.0001)
#' }
#'
#'
#' @export
simulate_screen <- function(
  n_genes_in_library,
  gene_abundance_alpha,
  gene_hit_rate,
  n_wells_per_screen,
  well_abundance_factor_mu,
  well_abundance_factor_var,
  imaging_n_cells_per_well_mu,
  imaging_n_cells_per_well_var,
  class_pos_mu,
  class_pos_var,
  class_neg_mu,
  class_neg_var,
  sequencing_n_cells_per_well_lambda,
  sequencing_n_cells_per_well_nu,
  pcr_factor_mu,
  pcr_factor_var,
  n_cells_per_gene_per_well,
  n_barcodes_per_gene_per_well,
  n_reads_total) {

  gene_library <- simulate_library(
    n_genes_in_library = n_genes_in_library,
    gene_abundance_alpha = gene_abundance_alpha,
    gene_hit_rate = gene_hit_rate)

  spot_plate <- simulate_spot_plate(
    gene_library,
    n_wells_per_screen,
    well_abundance_factor_mu,
    well_abundance_factor_var)

  imaging_plate <- simulate_imaging_plate(
    spot_plate,
    imaging_n_cells_per_well_mu = imaging_n_cells_per_well_mu,
    imaging_n_cells_per_well_var = imaging_n_cells_per_well_var,
    class_pos_mu = class_pos_mu,
    class_pos_var = class_pos_var,
    class_neg_mu = class_neg_mu,
    class_neg_var = class_neg_var)

  sequencing_plate <- simulate_sequencing_plate(
    spot_plate,
    sequencing_n_cells_per_well_lambda = sequencing_n_cells_per_well_lambda,
    sequencing_n_cells_per_well_nu = sequencing_n_cells_per_well_nu,
    pcr_factor_mu = pcr_factor_mu,
    pcr_factor_var = pcr_factor_var,
    n_reads_total = n_reads_total)

  screen_data <- spot_plate |>
    dplyr::left_join(imaging_plate, by = c("well", "gene")) |>
    dplyr::left_join(sequencing_plate, by = c("well", "gene"))
  screen_data
}
