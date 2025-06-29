
#' Prepare data for modeling
#'
#' Given a dataset, compute for each (well, gene) pair, the number of positive
#' and log10 expression, where expression is the the fraction of reads for that
#' gene in the well + .0001
#'
#' @param data data.frame with columns \[`well`, `gene`, `positive`,
#'   `n_reads_per_gene_per_well`\]
#'
#'
#' @return `tibble::tibble` with columns \[`well`, `Npositive`,
#'   `log10expression`\] where the types are `factor`, `numeric`, `matrix` and
#'   the matrix has `1` row and `n_gene` columns.
#'
#' @examples
#' \dontrun{
#'   data <- simulate_screen(...)
#'   model_data <- prepare_model_data(data)
#'   model <- compile_model(model_data)
#'   model_fit <- fit_model(model_data, model)
#'   model_estimate <- gather_model_estimate(model_fit)
#'   model_evaluation <- evaluate_model_fit(data, model_estimate)
#' }
#'
#'
#'
#' @export
prepare_model_data <- function(data) {
  data_model <- data |>
    dplyr::group_by(well) |>
    dplyr::do({
      well_data <- .
      total_reads <- sum(well_data$n_reads_per_gene_per_well)

      tibble::tibble(
        well = as.factor(well_data$well[1]),
        Npositive = sum(well_data$positive),
        Ntotal = sum(well_data$imaging_n_cells_per_gene_per_well),
        log10expression = matrix(
          log10(well_data$n_reads_per_gene_per_well/total_reads + .0001),
          nrow = 1))
    }) |>
    dplyr::ungroup() |>
    dplyr::filter(Ntotal > 0)
}

#' Compile a brms model
#'
#' Given a dataset that has been prepared for model to predict the number of
#' positive, cells given a matrix of cell counts for each gene in each well,
#' but don't run it.
#'
#' @param model_data `tibble::tibble` result of running `prepare_model_data`
#'
#' @return `brms::brmsfit` object
#'
#'
#' @examples
#' \dontrun{
#'.   data <- simulate_screen(...)
#'.   model_data <- prepare_model_data(data)
#'.   model <- compile_model(model_data)
#'.   model_fit <- fit_model(model_data, model)
#'.   model_estimate <- gather_model_estimate(model_fit)
#'.   model_evaluation <- evaluate_model_fit(data, model_estimate)
#' }
#'
#' @export
compile_model <- function(
    model_data,
    backend = "cmdstanr") {

  brms::brm(
    formula = Npositive ~ 1 + log10expression + offset(log(Ntotal)),
    data = model_data,
    family = poisson,
    prior = brms::prior(horseshoe(df = 10), class = "b"),
    iter = 6000,
    chains = 0,
    control = list(
      max_treedepth = 12),
    backend = backend,
    cores = 4)
}


#' Fit a brms model
#'
#' Given a model prepared by `compile_model` fit it with the given data
#'
#' Note by separating out the compilation step it is one model can be
#' fit multiple times.
#'
#' @param model_data `tibble::tibble` data prepared for the model with
#'   `model_prepare_data`
#' @param model `brms::brmfit` a compiled model the result of `compile_model`
#' @param chains int number of parellel chains to run
#' @param iter int number of iterations for each chain
#' @param cores int number of cores to use (e.g. the same number of chains)
#' @param control list parameters to control the mcmc sampling
#'
#' @returns `brms::brmsfit` object
#'
#' @examples
#' \dontrun{
#'   data <- simulate_screen(...)
#'   model_data <- prepare_model_data(data)
#'   model <- compile_model(model_data)
#'   model_fit <- fit_model(model_data, model)
#'   model_estimate <- gather_model_estimate(model_fit)
#'   model_evaluation <- evaluate_model_fit(data, model_estimate)
#' }
#'
#' @export
fit_model <- function(
    model_data,
    model,
    chains = 4,
    iter = 8000,
    cores = 4,
    control = list(
      max_treedepth = 12),
    backend = "cmdstanr",
    ...) {

  model |>
    stats::update(
      newdata = model_data,
      chains = chains,
      iter = iter,
      cores = cores,
      control = control,
      backend = backend,
      recompile = FALSE,
      ...)
}

#' Gather model parameter estimates
#'
#' Gather the model parameter estimates into table with one row per gene
#' and parameter estimate where higher means more likely to be a hit
#'
#' @param model `brms::brmsfit` fit using `fit_model`.
#'
#' @returns `data.frame` with columns
#'
#' @examples
#' \dontrun{
#'   data <- simulate_screen(...)
#'   model_data <- prepare_model_data(data)
#'   model <- compile_model(model_data)
#'   model_fit <- fit_model(model_data, model)
#'   model_estimate <- gather_model_estimate(model_fit)
#'   model_evaluation <- evaluate_model_fit(data, model_estimate)
#' }
#'
#'
#' @export
gather_model_estimate <- function(model_fit) {
  model_fit |>
    posterior::summarize_draws() |>
    dplyr::filter(variable |> stringr::str_detect("expression")) |>
    dplyr::transmute(
      gene = variable |> stringr::str_extract("[0-9]+$") |> as.integer(),
      variable_label = paste0("Gene ", gene) |> forcats::fct_inorder(),
      mean, q5, q95)
}

#' Evaluate the model performance against ground truth hit/no-hit status
#'
#' Evaluate the average precision and area under the ROC curve as quality
#' metrics for predicting the hit-status of each gene.
#'
#' @param data `data.frame` e.g. generated by `simulate_data` having columns
#'   \[`gene`, hit`\]
#' @param model_estimate `data.frame` e.g. produced by `gather_model_estimate`
#'   having columns \[`gene`, `mean`\]
#'
#' @return `data.frame` with columns \[`model_ap`, `model_auroc`\]
#'
#'
#' @examples
#' \dontrun{
#'   data <- simulate_screen(...)
#'   model_data <- prepare_model_data(data)
#'   model <- compile_model(model_data)
#'   model_fit <- fit_model(model_data, model)
#'   model_estimate <- gather_model_estimate(model_fit)
#'   model_evaluation <- evaluate_model_fit(data, model_estimate)
#' }
#'
#'
#' @export
evaluate_model_fit <- function(
  data,
  model_estimate) {

  model_estimate <- model_estimate |>
    dplyr::left_join(
      data |> dplyr::distinct(gene, hit),
      by = c("gene")) |>
    dplyr::mutate(
      mean_inv = -mean,
      hit = hit |> factor(
        levels = c(0, 1),
        labels = c("no", "yes")))

  model_ap <- model_estimate |>
    yardstick::average_precision(mean_inv, truth = hit)
  model_auroc <- model_estimate |>
    yardstick::roc_auc(mean_inv, truth = hit)

  tibble::tibble(
    model_ap = model_ap$.estimate[1],
    model_auroc = model_auroc$.estimate[1])
}


