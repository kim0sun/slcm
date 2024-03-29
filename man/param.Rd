% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{param}
\alias{param}
\alias{param.slcm}
\title{Printing Estimated Parameters of \code{slcm} Object}
\usage{
param(object, ...)

\method{param}{slcm}(
   object, type = c("probs", "logit"),
   se = FALSE, index = FALSE ..
)
}
\arguments{
\item{object}{an object of class \code{slcm} and \code{estimated}.}

\item{what}{specifies which parameter types to display. Valid options are \code{"pi"}, \code{"tau"}, and \code{"rho"}. \code{pi} represents membership probabilities of root variable, \code{tau} denotes conditional probabilities between latent class variables, and \code{rho} corresponds to item response probabilities for each measurement model. Multiple types can be chosen to display.}

\item{type}{a character string indicating the format in which the estimated parameters should be returned. If set to \code{"probs"}, estimates are returned in probability form, while \code{"logit"} returns them in log-odds (logit) form. Default is \code{"probs"}.}

\item{index}{a logical value. If \code{TRUE}, the indices of the estimated parameters are included in the output, enclosed in parentheses. If \code{FALSE}, they are omitted.}

\item{digits}{an integer that sets the number of decimal places for rounding the output.}
}
\value{
A list containing the specified estimated parameters (or corresponding standard errors for \code{se} function):
\item{pi}{Membership probabilities of the root variable.}
\item{tau}{Conditional probabilities between latent class variables, represented with uppercase alphabets for considering measurement invariance.}
\item{rho}{Item response probabilities for each measurement model, represented with lowercase alphabets for considering measurement invariance.}
}
\description{
By passing \code{estimated} \code{slcm} object, the function prints estimated parameters of the slcm model.
}
