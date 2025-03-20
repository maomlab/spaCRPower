

#' Plot model estimate
#'
#' For each gene (vertical) plot the model estimate (horizontal) as
#' the mean and 95% credible interval
#'
#' @param model_estimate data.frame result of calling `gather_model_estimate`
#'
#' @returns `ggplot2::ggplot` object
#'
#'
#' @examples
#' \dontrun{
#'    data <- simulate_data(...)
#'    model_data <- prepare_model_data(data)
#'    model <- compile_model(model_data)
#'    model_fit <- fit_model(model_data, model)
#'    model_estimate <- gather_model_estimate(model_fit)
#'    plot <- plot_model_estimate(model_estimate)
#'    ggplot2::ggsave(
#'      filename = "model_estimate.png",
#'      plot = plot,
#'      width = 5,
#'      height = 4)
#' }
#'
#' @export
plot_model_estimate <- function(model_estimate) {

  plot_data <- model_estimate |>
    dplyr::filter(variable |> stringr::str_detect("count")) |>
    dplyr::transmute(
      variable =
        paste0("Gene ", variable |> stringr::str_extract("[0-9]+$")) |>
        forcats::fct_inorder(),
      mean,
      q5,
      q95)


  ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous("Estimated Gene Score") +
    ggplot2::scale_y_discrete("Gene ID") +
    ggplot2::geom_segment(
      data = plot_data,
      mapping = ggplot2::aes(
        x = q5,
        y = variable,
        xend = q95,
        yend = variable),
      color = "grey50",
      linewidth = 1.5) +
    ggplot2::geom_point(
      data = plot_data,
      mapping = ggplot2::aes(
        x = q5,
        y = variable),
      color = "grey20") +
    ggplot2::geom_point(
      data = plot_data,
      mapping = ggplot2::aes(
        x = q95,
        y = variable),
      color = "grey20") +
    ggplot2::geom_point(
      data = plot_data,
      mapping = ggplot2::aes(
        x = mean,
        y = variable),
      color = "orange",
      size = 1)
}
