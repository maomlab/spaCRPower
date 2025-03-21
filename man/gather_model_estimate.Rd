% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_model.R
\name{gather_model_estimate}
\alias{gather_model_estimate}
\title{Gather model parameter estimates}
\usage{
gather_model_estimate(model_fit)
}
\arguments{
\item{model}{\code{brms::brmsfit} fit using \code{fit_model}.}
}
\value{
\code{data.frame} with columns
}
\description{
Gather the model parameter estimates into table with one row per gene
and parameter estimate where higher means more likely to be a hit
}
\examples{
\dontrun{
  data <- simulate_screen(...)
  model_data <- prepare_model_data(data)
  model <- compile_model(model_data)
  model_fit <- fit_model(model_data, model)
  model_estimate <- gather_model_estimate(model_fit)
  model_evaluation <- evaluate_model_fit(data, model_estimate)
}


}
