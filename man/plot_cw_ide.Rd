\name{plot_cw_ide}
\alias{plot_cw_ide}
\title{
Plot sensitivity of the estimated total indirect effect
}
\description{
Plot sensitivity of the estimated total indirect effect between mediators and an outcome
}
\usage{
plot_cw_ide(cmmb.output)
}
\arguments{
  \item{cmmb.output}{an output of the \emph{cmmb} function}
}
\value{It plots estimated component-wise products of path coefficients}
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
data(cmmb_demo_data)
# Run CMM for binary outcomes
rslt <- cmmb(Y=cmmb_demo_data$Y, M=cmmb_demo_data$M,
             tr=cmmb_demo_data$tr, X=cmmb_demo_data$X)
rslt
# Plot products of component-wise path coefficients
plot_cw_ide(rslt)
}
}