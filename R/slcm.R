#' Construct structural latent class model
#'
#' Function for constructing multivariate latent class model structure
#'
#'
#' @param formula a formula for specifying latent structure. The details of model specification are given under 'Details'.
#' @param ... additional formulae of the model
#' @param constraints a list containing constraints for the measurement invariance assumption. The details of posing constraints are given under 'Details'
#'
#' @details
#' 1. Define latent class variable with manifest indicators
#' \preformatted{LC1[k] ~ x1 + x2 + x3
#' LC2[k] ~ y1 + y2 + y3
#' LC3[k] ~ z1 + z2 + z3}
#'
#' 2. Relate latent class variables each other
#' \preformatted{LC1 ~ LC2}
#'
#' 3. Define higher-level latent class variable
#' \preformatted{P[k] ~ LC1 + LC2 + LC3}
#'
#'
#' 1. Measurement invariance for measurement model
#' \preformatted{c("LC1", "LC2", "LC3")}
#'
#' 2. Measurement invariance for structural model
#' \preformatted{c("P ~ LC1", "P -> LC2")}
#'
#' @example examples/slcm.R
#'
#'
#' @export
slcm = function(
   formula = NULL, ..., constraints = NULL
) {
   if (!is.list(formula)) formulae <- list(formula, ...)
   else formulae <- formula
   formulae <- formulae[sapply(formulae, inherits, "formula")]

   res <- list()
   res$model <- proc_formula(formulae, constraints)
   class(res) <- "slcm"
   res
}
