
#' Plot wells per gene histogram
#'
#' From the given dataset, plot the distribution of wells per gene
#' as a histogram
#'
#' @param data data.frame with columns \[`well`, `gene`\] and each row is a well
#'   well x gene pair
#'
#' @returns `ggplot2::ggplot` object
#'
#'
#' @examples
#' \dontrun{
#'    data <- simulate_data(...)
#'    plot <- plot_wells_per_gene(data)
#'    ggplot2::ggsave(
#'      filename = "wells_per_gene.png",
#'      plot = plot,
#'      width = 5,
#'      height = 4)
#' }
#'
#' @export
plot_wells_per_gene <- function(data) {

  plot_data <- data |>
    dplyr::group_by(gene) |>
    dplyr::summarize(
      n_wells = sum(count > 0),
      .groups = "drop")

  plot_data_mean <- plot_data |>
    dplyr::summarize(
      mean_n_wells = mean(n_wells),
      .groups = "drop")

  ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::geom_bar(
      data = plot_data,
      mapping = ggplot2::aes(
        x = n_wells)) +
    ggplot2::geom_vline(
      data = plot_data_mean,
      mapping = ggplot2::aes(
        xintercept = mean_n_wells)) +
    ggrepel::geom_text_repel(
      data = plot_data_mean,
      mapping = ggplot2::aes(
        x = mean_n_wells,
        label = round(mean_n_wells),
        y = 50),
      direction = "x",
      nudge_x = 1) +
    ggplot2::scale_y_continuous("Gene Count") +
    ggplot2::scale_x_continuous("Wells per Gene Count")
}
