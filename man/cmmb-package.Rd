\name{cmmb-package}
\alias{cmmb-package}
\docType{package}
\title{
Compositional Mediation Model for Binary Outcomes
}
\description{
Estimate direct and indirect effects of treatment on binary outcomes transmitted through compositional mediators
}
\author{
Michael B. Sohn

Maintainer: Michael B. Sohn <michael_sohn@urmc.rochester.edu>
}
\references{
Sohn, M.B., Lu, J. and Li, H. (2021). A Compositional Mediation Model for Binary Outcome: Application to Microbiome Studies (Submitted)
}
\examples{
\dontrun{
# Load a simulated dataset
data(cmmb_demo_data)
# Run CMM for binary outcomes
rslt <- cmmb(Y=cmmb_demo_data$Y, M=cmmb_demo_data$M,
             tr=cmmb_demo_data$tr, X=cmmb_demo_data$X, ForSA=TRUE)
rslt$total
# Plot products of component-wise path coefficients
plot_cw_ide(rslt)
# Plot sensitivity of the estimated total indirect effect
plot_cmmb_sa(rslt)
}
}