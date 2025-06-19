
#' Plot Classification Scores
#'
#' Given parameters for the distribution of scores for negative and positively
#' scored objects, generate a plot of the distribution of scores for each class
#'
#' The functional for each class is a beta distribution, where the parameters
#' are given in the mean and variance form.
#'
#' @param class_pos_mu numeric mean classification score of the positive class
#' @param class_pos_var numeric variance of classification score of the positive
#'        class
#' @param class_pos_mu numeric mean classification score of the negative class
#' @param class_pos_var numeric variance of classification score of the negative
#'        class
#' @returns `ggplot2::ggplot` object.
#'
#' @examples
#' \dontrun{
#'    plot <- plot_classification_scores(
#'      class_pos_mu = 0.99,
#'      class_pos_var = 0.0001,
#'      class_neg_mu = 0.01,
#'      class_neg_var = 0.0001)
#'    ggplot2::ggsave(
#'      filename = "classification_scores.png",
#'      plot = plot,
#'      width = 5,
#'      height = 4)
#' }
#'
#' @export
plot_classification_scores <- function(
  class_pos_mu,
  class_pos_var,
  class_neg_mu,
  class_neg_var) {

  plot_data <- data.frame(
    x = seq(0, 1, length.out = 1000)) |>
    dplyr::mutate(
      class_pos = dbeta(
        x = x,
        shape1 = class_pos_mu *
          ((class_pos_mu * (1 - class_pos_mu) / class_pos_var) - 1),
        shape2 = (1 - class_pos_mu) *
          ((class_pos_mu * (1 - class_pos_mu) / class_pos_var) - 1)),
      class_neg = dbeta(
        x = x,
        shape1 = class_neg_mu *
          ((class_neg_mu * (1 - class_neg_mu) / class_neg_var) - 1),
        shape2 = (1 - class_neg_mu) *
          ((class_neg_mu * (1 - class_neg_mu) / class_neg_var) - 1)))

  ggplot2::ggplot(data = plot_data) +
    ggplot2::theme_bw() +
    ggplot2::geom_area(
      mapping = ggplot2::aes(
        x = x,
        y = log(class_neg+1)),
      fill = "blue",
      alpha = 0.6) +
    ggplot2::geom_area(
      mapping = ggplot2::aes(
        x = x,
        y = log10(class_pos+1)),
      fill = "orange",
    alpha = 0.6) +
    ggplot2::scale_x_continuous(
      "Cell classification score",
      limits = c(0, 1)) +
    ggplot2::scale_y_log10(
      "Density score")
}
