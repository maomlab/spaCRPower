

#' Plot positivity rate by well
#'
#' Given a dataset generate a scatter plot comparing number of positive
#' cells by the number of total cells for each well.
#'
#' @param data data.frame with columns \[`well`, `n_cells_per_gene_per_well`,
#'   `positive`\]
#'
#' @returns `ggplot2::ggplot` object
#'
#' @examples
#' \dontrun{
#'    data <- simulate_data(...)
#'    plot <- plot_positivity_rate_by_well(data)
#'    ggplot2::ggsave(
#'      filename = "positivity_rate_by_well.png",
#'      plot = plot,
#'      width = 5,
#'      height = 4)
#' }
#'
#' @export
plot_positivity_rate_by_well <- function(data) {

  plot_data <- data |>
    dplyr::group_by(well) |>
    dplyr::summarize(
      n_cells = sum(n_cells_per_gene_per_well),
      n_positive = sum(positive),
      .groups = "drop")

  ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::ggtitle(
      "Distribution of positivity rate by well") +
    ggplot2::geom_point(
      data = plot_data,
      mapping = ggplot2::aes(
        x = n_cells,
        y = n_positive)) +
    ggplot2::scale_x_continuous(
      "Cells per well") +
    ggplot2::scale_y_continuous(
      "Predicted positive cells per well")


}
