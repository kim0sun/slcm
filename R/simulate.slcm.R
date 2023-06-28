#' Simulate \code{slcm} Object
#' Simulate data from \code{slcm} model.
#'
#' @param object a \code{slcm} object
#' @param nsim number of response data to simulate. Defaults to 100.
#' @param params parameters to be designated
#' @param seed random seed
#'
#' @export
simulate.slcm <- function(
   object, nobs = 500, nlevel = 2,
   what = c("response", "class", "probs"),
   params = NULL, seed, ...
) {
   cl <- match.call()
   model <- object$model

   if (inherits(object, "estimated")) {
      level <- levels(object$mf)
      arg <- arguments_nmf(model, nobs, level)
      par <- object$par
   } else {
      level <- sapply(unlist(model$latent$children), function(x)
         seq_len(nlevel), simplify = FALSE)
      arg <- arguments_nmf(model, nobs, level)
      if (missing(params)) {
         params <- runif(length(arg$id))
         par <- unlist(tapply(params, arg$id, norm1), use.names = FALSE)
      } else {
         if (is.list(params)) params <- unlist(params)
         if (all(params >= 0))
         par <- unlist(tapply(params, arg$id, norm1), use.names = FALSE)
         else
         par <- unlist(tapply(params, arg$id, norm2), use.names = FALSE)
      }
   }
   sim <- simModel(
      nobs, arg$nvar, arg$nlev, object$par,
      arg$nlv, arg$nrl, arg$nlf, arg$npi, arg$ntau, arg$nrho,
      arg$ul, arg$vl, arg$lf, arg$rt, arg$eqrl, arg$eqlf,
      arg$nc, arg$nk, arg$nl, arg$ncl,
      arg$nc_pi, arg$nk_tau, arg$nl_tau, arg$nc_rho, arg$nr_rho
   )
   class <- data.frame(sim$class)
   names(class) <- model$latent$label

   # data.name
   y <- data.frame(do.call(cbind, sim$y))
   colnames(y) <- unlist(model$latent[model$latent$leaf, "children"])
   mf <- proc_data(y, model)

   res <- list(response = mf, class = sim$class,
        probs = relist(exp(par), arg$skeleton$par))

   if (!("what" %in% names(cl)))
      return(res$response)
   res[match.arg(what, several.ok = TRUE)]
}
