

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
    dplyr::arrange(dplyr::desc(mean)) |>
    dplyr::mutate(variable_label = variable_label |> forcats::fct_inorder())

  ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous("Estimated Gene Score") +
    ggplot2::scale_y_discrete("Gene ID") +
    ggplot2::geom_segment(
      data = plot_data,
      mapping = ggplot2::aes(
        x = q5,
        y = variable_label,
        xend = q95,
        yend = variable_label),
      color = "grey50",
      linewidth = 1.5) +
    ggplot2::geom_point(
      data = plot_data,
      mapping = ggplot2::aes(
        x = q5,
        y = variable_label),
      color = "grey20") +
    ggplot2::geom_point(
      data = plot_data,
      mapping = ggplot2::aes(
        x = q95,
        y = variable_label),
      color = "grey20") +
    ggplot2::geom_point(
      data = plot_data,
      mapping = ggplot2::aes(
        x = mean,
        y = variable_label),
      color = "orange",
      size = 1)  
}
