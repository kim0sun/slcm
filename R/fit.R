#' Estimation for Parameters of \code{slcm} Object
#' @aliases fit fit.slcm
#' @usage
#'
#' fit(object, ...)
#'
#' \method{fit}{slcm}(object, data,
#'     method = c("em", "hybrid", "nlm"),
#'     restriction = NULL,
#'     control = slcm.control(), ...)
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
fit <- function(object, ...) UseMethod("fit")

#' @export
fit.slcm <- function(
      object, data, formula,
      method = c("em", "hybrid", "nlm"),
      restriction = NULL,
      control = slcm.control(), ...
) {
   method <- match.arg(method)
   if (is.null(object$data)) {
      if (missing(data)) {
         data = parent.frame()
      }
      data <- proc_data(data, object$model)
      args <- update_args(object$args, data)
   } else if (missing(data)) {
      data <- object$data
      args <- object$args
   } else {
      data <- proc_data(data, object$model)
      args <- update_args(object$args, data)
   }
   object$data <- data
   object$args <- args

   if (!inherits(control, "slcm.control")) {
      ctrl <- slcm.control()
      index <- match(names(control), names(ctrl), 0L)
      ctrl[index] <- control[index > 0]
      control <- ctrl
   }

   init.param <- control$init.param
   if (object$fitted)
      init.param <- rapply(object$estimates$param, log, how = "list")
   else {
      is.init <- init_validate(init.param, args)
      init.param <- par_gnr(
         args$nobs, args$nvar, args$ncat, args$nroot,
         args$nlink_unique, args$nleaf_unique,
         args$root - 1, args$u - 1, args$v - 1,
         args$nclass, args$nclass_leaf,
         args$nclass_u, args$nclass_v, is.init, init.param
      )
   }

   if (!is.null(object$restriction$target))
      restriction <- union(object$restriction$target, restriction)
   restr <- proc_restr(restriction, init.param, args)
   if (any(is.na(unlist(restr$ref)))) stop("Invalid restriction.")
   if (length(restriction) > 0 && control$verbose)
      cat("Restricted to zero:", restriction, "\n")
   args$npar[2:3] <- args$npar[2:3] -
      sapply(restr$restr0, function(x) sum(unlist(x)))
   init.param <- restr$param0

   est <- estimate(method, control, data, args, init.param, restr)
   log_par <- est$par
   convergence <- est$conv

   # GIBBS
   etc <- calcModel(
      log_par, data$y, args$nobs, args$nvar, args$ncat,
      args$nlv, args$nroot, args$nlink, args$nleaf,
      args$nlink_unique, args$nleaf_unique,
      args$root - 1, args$tree_index - 1,
      args$u - 1, args$v - 1, args$leaf - 1,
      args$cstr_link - 1, args$cstr_leaf - 1,
      args$nclass, args$nclass_leaf,
      args$nclass_u, args$nclass_v,
      restr$ref
   )

   logits <- log2logit(
      log_par, length(unlist(log_par)), args$ncat,
      args$nroot, args$nlink_unique, args$nleaf_unique,
      args$nclass_root, args$nclass_leaf,
      args$nclass_u, args$nclass_v,
      unlist(restr$restr0), unlist(restr$ref) - 1
   )

   logit_par <- param2list(
      logits, args$ncat, args$nroot,
      args$nlink_unique, args$nleaf_unique,
      args$root, args$u, args$v,
      args$nclass[args$root],
      args$nclass_u, args$nclass_v,
      args$nclass_leaf
   )

   logit_cov <- vcov_logit_fi(etc$score, args, unlist(restr$restr0))
   logit_se  <- param2list(
      get_se(logit_cov), args$ncat, args$nroot,
      args$nlink_unique, args$nleaf_unique, args$root,
      args$u, args$v, args$nclass[args$root],
      args$nclass_u, args$nclass_v, args$nclass_leaf
   )
   param_cov <- vcov_transform(log_par, logit_cov, args)
   param_se  <- param2list(
      get_se(param_cov), args$ncat, args$nroot,
      args$nlink_unique, args$nleaf_unique, args$root,
      args$u, args$v, args$nclass[args$root],
      args$nclass_u, args$nclass_v, args$nclass_leaf
   )

   object$estimates = list(
      param = output_param(log_par, object$model, args),
      logit = output_param(logit_par, object$model, args, FALSE),
      vcov = list(param = param_cov, logit = logit_cov),
      se = list(param = output_param(param_se, object$model, args, FALSE),
                logit = output_param(logit_se, object$model, args, FALSE))
   )

   object$method = method
   object$args <- args
   object$args$log_par <- log_par
   object$llik <- etc$ll
   object$posterior <- output_posterior(etc$post, object$model, data)
   object$joint <- etc$joint
   object$score <- unlist(rapply(etc$score, rowSums))
   object$convergence <- convergence
   object$fitted <- TRUE
   object$control <- control

   class(object) <- "slcm"
   object
}


