% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_well_gene_heatmap.R
\name{plot_well_gene_heatmap}
\alias{plot_well_gene_heatmap}
\title{Plot well x gene heatmap}
\usage{
plot_well_gene_heatmap(data)
}
\arguments{
\item{data}{data.frame with columns [\code{gene}, \code{well}, \code{count}]}
}
\value{
\code{ggplot2::ggplot} object
}
\description{
Given a dataset plot heatmap where each box corresponds to a well
(row) by gene (column) pair. The values are the \code{log10(count + 1)} on
a scale from blue to yellow.
}
\examples{
\dontrun{
   data <- simulate_data(...)
   plot <- plot_well_gene_heatmap(data)
   ggplot2::ggsave(
     filename = "well_gene_heatmap.png",
     plot = plot,
     width = 5,
     height = 4)
}

}
