% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calibration.R
\name{calibration}
\alias{calibration}
\title{Cross Validation}
\usage{
calibration(train, test = train, FittingFunction, covariates, time,
  by = 0.1, ...)
}
\arguments{
\item{train}{a data.frame which includes follow up time, event, and covariates. In a k-fold cross validation this is a sample proportion (k-1)/k of the original data set.}

\item{test}{a data.frame which includes follow up time, event, and covariates. In a k-fold cross validation this is a sample proportion 1/k of the original data set.}

\item{FittingFunction}{a function which returns list of survival probabilities of individuals where each element in the list represent a specific time point.}

\item{covariates}{character vector specifying the names of covariates}

\item{time}{a numeric value specifying at what time survival probability is to be calculated.}

\item{by}{a numeric value specifying the spacing between values in the proportion vector ranging from 0 to 1.}
}
\value{
The function returns a data frame of the proportion that had an event for a given survival probability including confidence intervals
}
\description{
Cross Validation
}
\examples{
calibration(train,test)
}

