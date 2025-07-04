% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_screen.R
\name{rgamma_mean_variance}
\alias{rgamma_mean_variance}
\title{Parameterize the gamma distribution by the mean and variance}
\usage{
rgamma_mean_variance(n, mean, var)
}
\description{
For a gamma distribution parameterized by the shape and rate is
}
\details{
\if{html}{\out{<div class="sourceCode">}}\preformatted{mean = shape / rate
variance = shape / rate^2
}\if{html}{\out{</div>}}

solving both for alpha and setting them equal gives

shape = mean\emph{rate = variance}rate^2

solving for rate in terms of mean and variance gives

rate = mean / variance

substituting that back in to solve for shape gives

shape = mean^2 / variance
}
