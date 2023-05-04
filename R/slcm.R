#' Construct structural latent class model
#'
#' Function for constructing multivariate latent class model structure
#'
#' @param formula a formula for specifying latent structure. The details of model specification are below.
#' @param ... additional formulae of the model
#' @param constraints a list containing constraints of the measurement invariance.
#'
#' @details formula
#'
#'
#' @export
slcm = function(
   x = NULL, ..., constraints = NULL
) {
   if (!is.list(x)) formulae <- list(x, ...)
   else formulae <- x
   formulae <- formulae[sapply(formulae, inherits, "formula")]
   model <- proc_formula(formulae, constraints)

   res = list()
   res$model <- model
   res$args <- args_return(model)
   res$fitted <- FALSE

   class(res) <- "slcm"
   res
}