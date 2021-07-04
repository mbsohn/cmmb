\name{cmmb}
\alias{cmmb}
\title{
Compositional Mediation Model for Binary Outcomes
}
\description{
Estimate direct and indirect effects of treatment on binary outcomes transmitted through compositional mediators
}
\usage{
cmmb(Y, M, tr, X, n.cores=NULL, n.boot=2000, ci.method="empirical",
     p.value=FALSE, ForSA=FALSE, max.rho=0.5, sig.level=0.05, FWER=FALSE,
     w=rep(1,length(Y)), prec=1e-4, max.iter=1000)
}
\arguments{
  \item{Y}{a vector of binary outcomes}
  \item{M}{a matrix of compositional data}
  \item{tr}{a vector of continuous or binary treatments}
  \item{X}{a matrix of covariates}
  \item{n.cores}{a number of CPU cores for parallel processing}
  \item{n.boot}{a number of bootstrap samples}
  \item{ci.method}{options for bootstrap confidence interval. It can be either "empirical" (default) or "percentile".}
  \item{p.value}{a logical value for calculating the p value. It is inactive when \emph{ci.method="percentile"}.}
  \item{ForSA}{a logical value for sensitivity analysis}
  \item{max.rho}{a maximum correlation allowed between mediators and an outcome}
  \item{sig.level}{a significance level to estimate bootstrap confidence intervals for direct and indirect effects of treatment}
  \item{FWER}{a logical value for family-wise error rate for direct and total indirect effects. If \emph{FWER=TRUE}, the Bonferroni correct will be applied.}
  \item{w}{a vector of weights on samples. If measurements in a sample is more reliable than others, this argument can be used to take that information into the model.}
  \item{prec}{an error tolerance or a stopping criterion for the debiasd procedure}
  \item{max.iter}{a maximum number of iteration in the debias procedure}

  Note: the range of rho is not from -1 to 1 when the number of components is more than two because the correlation between them is not zero, and the range gets smaller as the number of components increases.
}
\value{
If \emph{ForSA=FALSE},
  \item{total}{contains estimated direct and total indirect effects with their confidence limits}
  \item{cwprod}{contains component-wise products of path coefficients with their confidence limits}
If \emph{ForSA=TRUE},
  \item{total}{contains estimated direct and total indirect effects with their confidence limits}
  \item{cwprod}{contains component-wise products of path coefficients with their confidence limits}
  \item{cide.rho}{contains estimated indirect effects and corresponding pointwise 95\% confidence intervals, given correlations between mediators and an outcome}
}
\references{
  Sohn, M.B., Lu, J. and Li, H. (2021). \emph{A Compositional Mediation Model for Binary Outcome: Application to Microbiome Studies} (Submitted)
}
\author{
  Michael B. Sohn

  Maintainer: Michael B. Sohn <michael_sohn@urmc.rochester.edu>
}
\examples{
\dontrun{
# Load a simulated dataset
data(cmmb_deme_data)
# Run CMM for binary outcomes
rslt <- cmmb(Y=cmmb_deme_data$Y, M=cmmb_deme_data$M,
             tr=cmmb_deme_data$tr, X=cmmb_deme_data$X)
rslt
# Plot products of component-wise path coefficients
plot_cw_ide(rslt)
}
}