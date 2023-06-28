#' Estimation for Parameters of \code{slcm} Object
#' @aliases estimate estimate.slcm
#' @usage
#'
#' estimate(object, ...)
#'
#' \method{estimate}{slcm}(object, data,
#'     method = c("em", "hybrid", "nlm"),
#'     restriction = NULL,
#'     control = slcmControl(), ...)
#'
#' @param object a \code{slcm} object which defines latent structure to be estimated.
#' @param data a \code{data.frame} object which contains observed categorical variables beloning to the latent structure.
#' @param formula a \code{formula} object introducing exogenous covariates.
#' @param nrep number of trial
#' @param method estimating method for slcm parameters.
#' @param restriction a \code{list} of parameters to be restricted to zero.
#' @param control slcm control.
#'
#' @export
estimate <- function(object, ...) UseMethod("estimate")

#' @export
estimate.slcm <- function(
      object, data, formula,
      method = c("em", "hybrid", "nlm"),
      fix2zero = NULL,
      control = slcmControl(), ...
) {
   method <- match.arg(method)
   if (inherits(object, "estimated")) mf <- object$mf
   else {
      if (missing(data)) data = parent.frame()
      mf <- proc_data(data, object$model, object$num)
   }

   if (!inherits(control, "slcmControl")) {
      ctrl <- slcmControl()
      index <- match(names(control), names(ctrl), 0L)
      ctrl[index] <- control[!is.na(index)]
      control <- ctrl
   }

   arg <- arguments(object$model, mf, fix2zero)

   if (inherits(object, "estimated")) {
      par <- object$par
   } else if (!is.null(control$init.param)) {
      init.param <- unlist(control$init.param)
      if (all(init.param >= 0))
         par <- unlist(tapply(init.param, arg$id, norm1), use.names = FALSE)
      else
         par <- unlist(tapply(init.param, arg$id, norm2), use.names = FALSE)
   } else {
      for (i in 1:control$nrep) {
         init.param <- runif(length(arg$id))
         init.param[arg$fix0] <- -Inf
         par <- unlist(tapply(init.param, arg$id, norm2), use.names = FALSE)
      }
   }
   par[arg$fix0] <- -Inf
   est <- estModel(method, control, par, mf, arg)
   par <- est$par
   conv <- est$conv
   logit <- par - par[arg$ref_idx[arg$id]]

   etc <- calcModel(
      attr(mf, "y"), arg$nobs, arg$nvar, unlist(arg$nlev),
      par, arg$fix0, arg$ref - 1, arg$nlv, arg$nrl, arg$nlf,
      arg$npi, arg$ntau, arg$nrho, arg$ul, arg$vl,
      arg$lf, arg$tr, arg$rt, arg$eqrl, arg$eqlf,
      arg$nc, arg$nk, arg$nl, arg$ncl,
      arg$nc_pi, arg$nk_tau, arg$nl_tau, arg$nc_rho, arg$nr_rho
   )
   score <- relist(etc$score, arg$skeleton$score)
   score <- t(do.call(rbind, score))
   dimnames(score) <- list(dimnames(mf)[[1]], unlist(arg$par_index))

   post <- relist(exp(etc$post), arg$skeleton$post)
   joint <- relist(exp(etc$joint), arg$skeleton$joint)

   object$method = method
   object$arg <- arg
   object$mf <- mf
   object$par <- par
   object$logit <- logit
   object$fix2zero <- fix2zero
   object$score <- score
   object$posterior <- list(
      marginal = lapply(post, t), joint = joint
   )
   object$convergence <- conv
   object$loglikelihood <- etc$ll
   object$control <- control

   class(object) <- c("slcm", "estimated")
   object
}


