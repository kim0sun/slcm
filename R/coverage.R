#' Coverage Probabilities for Parameter Estimates for `slcm` model
#'
#' @param object a \code{slcm} object
#' @param nsim number of simulation sample
#' @param parm parameter to be returned
#' @param nobs number of observations for
#' @param level significnace level
#' @param verbose print
#' @param control slcm control
#'
#' @export
coverage <- function(
   object, nsim = 100, parm,
   level = 0.95, verbose = TRUE, control, ...
) {
   if (!object$fitted) return(NULL)

   method <- object$method
   ctrl <- object$control
   data <- object$data
   args <- object$args
   parm <- object$estimates$param
   restr <- object$restriction
   if (!missing(control))
      ctrl[names(control)] <- control
   control <- ctrl
   control$init.param <- parm
   control$verbose = FALSE

   cover <- par <- unlist(parm)
   cover[] <- 0

   for (i in 1:nsim) {
      if (verbose && i %% 50 == 0) cat(".")
      sim <- object %>% simulate(nsim = object$data$nobs, params = parm)
      fit <- estimate(method, control, data, args, control$init.param, restr)
      ci <- confint.slcm(fit, level = level)
      cover <- cover + as.numeric(par >= ci[, 1] & par <= ci[, 2])
   }
   if (verbose) cat("DONE. \n")
   cover / nsim
}
