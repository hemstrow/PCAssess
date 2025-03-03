% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/permutation.R
\name{plot_permutation_res}
\alias{plot_permutation_res}
\title{Plot PCA permutation results}
\usage{
plot_permutation_res(permute_res, n_perm_pcas = 1, plot_observed_pcas = TRUE)
}
\arguments{
\item{permute_res}{a list containing the null distribution, observed values,
p-vaules from PCA permutation testing, and possibly PCAs, like that
produced by \code{\link{run_permutation}}.}

\item{n_perm_pcas}{numeric, default 1. The number of permuted PCAs to plot.
Only used if PCA data is included in \code{permute_res}.}

\item{plot_observed_pcas}{logical, default TRUE. If TRUE and PCA data is
present in \code{permute_res}, observed PCAs will be plotted alongside
bootstraps and/or p-value distribution.}
}
\value{
a ggplot of the null distribution of F-statistic and observed
  F-statistic from PCA permutation testing. If PCAs are present, also
  optionally includes plots of the observed PCA and n number of permuted
  PCAs. PCAs for observed data will have a gray background, PCAs for
  permutations will not.
}
\description{
Plots results from PCA permutation testing.
}
\examples{
# run some permutations
res <- run_permutation(mon_sn[1:100,],
                        facet = sample(LETTERS[1:4], ncol(mon_sn), TRUE),
                        n = 10, fst_cut = 0.95, store_pca = TRUE)

# plot observed PCA and all 10 permutations
plot_permutation_res(res, 10)

# plot only the observed PCA
plot_permutation_res(res, 0)

# plot no PCAs, only the distribution of p-values
plot_permutation_res(res, 0, plot_observed_pcas = FALSE)

# plot only the real PCAs and the distribution of observed p-values
plot_permutation_res(res, 0, plot_observed_pcas = TRUE)
}
\author{
William Hemstrom

Andy Lee
}
