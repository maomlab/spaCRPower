
#' Plot cells per well histogram
#'
#' From the given dataset, plot the distribution of cells per well
#' as a histogram
#'
#' @param data data.frame with columns \[`well`,
#'   `n_cells_per_gene_per_well`\] and each row is a well x gene pair
#'
#' @returns `ggplot2::ggplot` object
#'
#' @examples
#' \dontrun{
#'    data <- simulate_data(...)
#'    plot <- plot_cells_per_well(data)
#'    ggplot2::ggsave(
#'      filename = "cells_per_well.png",
#'      plot = plot,
#'      width = 5,
#'      height = 4)
#' }
#'
#' @export
plot_cells_per_well <- function(data) {

  plot_data <- data |>
    dplyr::group_by(well) |>
    dplyr::summarize(
      n_cells = sum(n_cells_per_gene_per_well),
      .groups = "drop")

  plot_data_mean <- plot_data |>
    dplyr::summarize(
      mean_n_cells = mean(n_cells),
      .groups = "drop")

  max_height <- plot_data |>
    dplyr::count(n_cells) |>
    purrr::pluck("n") |>
    max()

  ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::geom_bar(
      data = plot_data,
      mapping = ggplot2::aes(
        x = n_cells)) +
    ggplot2::geom_vline(
      data = plot_data_mean,
      mapping = ggplot2::aes(
        xintercept = mean_n_cells)) +
    ggrepel::geom_text_repel(
      data = plot_data_mean,
      mapping = ggplot2::aes(
        x = mean_n_cells,
        label = signif(mean_n_cells, digits = 2),
        y = 1.1 * max_height),
      direction = "x",
      nudge_x = .1) +
    ggplot2::scale_y_continuous("Well Count") +
    ggplot2::scale_x_continuous("Cells per Well Count")
}
