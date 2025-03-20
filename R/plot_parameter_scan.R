

#' Plot results of a parameter scan
#'
#' Given a dataframe with parameters and scores, plot the scores
#' across specified parameters
#'
#' @param parameter_scan data.frame with columns specified by the parameters
#'   \[`<param_x>`, `<param_y>`, `<param_color>`].
#' @param param_x character column name of parameter to be mapped to the x-axis
#' @param param_color character column name of parameter to be mapped to the
#'   color dimension
#' @param param_y character column name of parameter to be mapped to the y-axis
#'
#' @returns `ggplot2::ggplot` object
#'
#' @examples
#' \dontrun{
#'   parameter_scan <- readr::read_tsv(file = "parameter_scan.tsv")
#'   plot <- plot_parameter_scan(parameter_scan)
#'   ggplot2::ggsave(
#'     filename = "parameter_scan_plot.png",
#'     width = 6,
#'     height = 5)
#' }
#'
#'
#' @export
plot_parameter_scan <- function(
    parameter_scan,
    param_x = "n_wells_per_screen",
    param_color = "n_genes_per_well",
    param_y = "model_auroc",
    label_x = "N Wells Per Screen",
    label_y = "Model AUROC (larger is better)",
    label_color = "N Genes Per Well",
    title = "Parameter Scan") {

  # ggplot2::ggplot() +
  #   ggplot2::theme_bw() +
  #   ggplot2::ggtitle(title) +
  #   ggplot2::geom_line(
  #     data = parameter_scan,
  #     mapping = ggplot2::aes(
  #       x = .data[[param_x]],
  #       color = .data[[param_color]],
  #       group = .data[[param_color]],
  #       y = .data[[param_y]]),
  #     size = 1.3) +
  #   ggplot2::scale_x_continuous(label_x) +
  #   ggplot2::scale_y_continuous(label_y) +
  #   ggplot2::scale_color_continuous(label_color) +
  #   ggplot2::theme(legend.position = "bottom")

  ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::ggtitle(title) +
    ggplot2::geom_point(
      data = parameter_scan,
      mapping = ggplot2::aes(
        x = .data[[param_x]],
        color = .data[[param_color]],
        group = .data[[param_color]],
        y = .data[[param_y]]),
      size = 1.3) +
    ggplot2::geom_smooth(
      data = parameter_scan,
      mapping = ggplot2::aes(
        x = .data[[param_x]],
        color = .data[[param_color]],
        group = .data[[param_color]],
        y = .data[[param_y]]),
      method = "loess",
      formula = y ~ x,
      size = 1.3) +
    ggplot2::scale_x_continuous(label_x) +
    ggplot2::scale_y_continuous(label_y) +
    viridis::scale_color_viridis(
      label_color,
      option = "C",
      begin = 0.2,
      end = 0.8) +
    ggplot2::theme(legend.position = "bottom")

}
