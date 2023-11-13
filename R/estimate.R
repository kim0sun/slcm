#' Estimating Parameters of `slcm` Object
#'
#' Estimate the parameters of model constructed using the `slcm` function.
#'
#' @aliases estimate estimate.slcm
#' @usage
#' estimate(object, ...)
#'
#' \method{estimate}{slcm}(object, data,
#'     method = c("em", "hybrid", "nlm"),
#'     fix2zero = NULL,
#'     control = slcmControl(), ...)
#'
#' @param object an `slcm` object defining SLCM model to be estimated.
#' @param data a `data.frame` object containing observed categorical variables incorporated in the model.
#' @param method estimation method for SLCM parameters. The default is `"em"`, which employs expecation-maximization (EM) algorithm for estimation; the alternative `"nlm"`, utilizes `nlm` function for Newton-Raphson algorithm. The `"hybrid"` method begins with the EM algorithm and concludes with the `nlm` function for refined estimation.
#' @param fix2zero a `vector` of parameters to be restricted to zero. The details of restriction is given under 'Details'
#' @param na.rm a logical value whether to remove NA.
#' @param control a `list` of control for the estimation procedure. Used to modify default values in [slcmControl].
#'
#' @details
#' To constrain certain parameters to zero, use the `fix2zero` argument. Each parameter is associated with a unique index. You can identify the index of a specific parameter by invoking the \link[slcm]{param()} function with the `index = TRUE` argument. To apply these constraints, include the relevant parameter indices in the `fix2zero` argument.
#'
#' @returns
#' An object of class `slcm` and `estimated` with an following elements:
#' \item{model}{a `list` describing of the model.}
#' \item{method}{the method used for estimation}
#' \item{arg}{the brief model description used during the estimation.}
#' \item{mf}{the data.frame used for estimation.}
#' \item{par}{the log of the estimated paramters.}
#' \item{logit}{the log-odds of the estimated parameters.}
#' \item{score}{the score function for the estimated parameters.}
#' \item{posterior}{the `list` of posterior probablities for each latent class variable.}
#' \item{convergence}{a logical indicator of whether convergence was achieved.}
#' \item{loglikelihood}{the loglikelihood of the estimated model.}
#' \item{control}{the control values used during the estimation process.}
#'
#' This returned object can be further processed using the \link[slcm]{param()} or \link[slcm]{se()} functions to extract the estimated parameters or their respective standard errors. Additionally, the \link[slcm]{regress()} function enables logistic regression analysis using three-step approach to evaluate the effect of external variables on latent class variables.
#'
#' @seealso \link[slcm]{slcm} \link[slcm]{param} \link[slcm]{regress} \link[slcm]{slcmControl}
#'
#' @export
estimate <- function(object, ...) UseMethod("estimate")

#' @exportS3Method slcm::estimate slcm
estimate.slcm <- function(
      object, data, formula,
      method = c("em", "hybrid", "nlm"),
      fix2zero = NULL,
      control = slcmControl(), ...
) {
   method <- match.arg(method)
   if (!inherits(control, "slcmControl")) {
      ctrl <- slcmControl()
      index <- match(names(control), names(ctrl), 0L)
      ctrl[index] <- control[!is.na(index)]
      control <- ctrl
   }
   na.rm <- control$na.rm
   if (inherits(object, "estimated")) mf <- object$mf
   else {
      if (missing(data)) data = parent.frame()
      mf <- proc_data(data, object$model, na.rm)
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


