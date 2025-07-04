% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_screen.R
\name{simulate_library}
\alias{simulate_library}
\title{Simulate a library of gene perturbations}
\usage{
simulate_library(n_genes_in_library, gene_abundance_alpha, gene_hit_rate)
}
\arguments{
\item{n_genes_in_library}{integer number of genes in the library}

\item{gene_abundance_alpha}{numeric positive real number giving the
how biased the gene abundance is across the library}

\item{gene_hit_rate}{numeric probability of each gene being a hit}
}
\value{
data.frame with \code{n_genes_in_library} rows and columns
\verb{\[gene, gene_abundance_alpha, gene_hit_rate, gene_abundance, hit\]}, where
\code{gene} is 1-based index for the gene, \code{gene_abundance_alpha}, and
\code{gene_hit_rate} record the simulation parameters, \code{gene_abundance} for gene
i is the fraction of the library represented by the gene i, such that the
sum of the \code{gene_abundance} values sum to \code{1}, and for gene i \code{hit} is
whether or not gene i is a hit.
}
\description{
Given the library size, sample the abundance and hit status for each gene in
the library:
}
\details{
abundance ~ dirichlet(gene_abundance_alpha*n)
for i in 1..n:
hit_i ~ binom(gene_hit_rate)

For the gene_abundance_alpha parameter larger values have abundances that are
are more evenly distributed
alpha -> Infinity: abundance_i -> 1/n
alpha -> 0:
k ~ sample(1, n)
abundance_i = Indicator(i = k)
When alpha is 1, then the Gini index = 0.5.
}
\examples{
\dontrun{
  scores <- spaCRPower::simulate_library(
    n_genes_in_library = 200,
    gene_abundance_alpha = 100,
    gene_hit_rate = .05)
}

}
