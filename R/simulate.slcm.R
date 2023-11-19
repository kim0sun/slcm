#' Simulate \code{slcm} Object
#' Simulate data from \code{slcm} model.
#'
#' @param object a \code{slcm} object
#' @param nsim number of response data to simulate. Defaults to 100.
#' @param nlevel number of categories for each manifest item.
#' @param what objects to be printed among response, class, probs (params)
#' @param params parameters to be designated
#' @param seed random seed
#'
#' @exportS3Method stats::simulate slcm
simulate.slcm <- function(
   object, nsim = 500, nlevel = 2,
   random.param = TRUE, seed = NULL, ...
) {
   model <- object$model
   arg <- object$arg

   if (inherits(object, "estimated")) {
      level <- levels(object$mf)
      par <- object$par
   } else {
      level <- sapply(unlist(arg$vars, use.names = FALSE), function(x)
         seq_len(nlevel), simplify = FALSE)
      arg <- arg_sim(arg, level)
      if (random.param) {
         par <- runif(length(arg$id))
      } else {
         par <- numeric(length(arg$id))
         cat("Insert PI\n")
         pos <- 0
         for (r in seq_len(arg$npi)) {
            cat(rownames(model$latent)[arg$rt[r] + 1], ":")
            input <- abs(scan(n = arg$nc_pi[r]))
            if (length(input == arg$nc_pi[r])) {
               par[(pos + 1):(pos + arg$nc_pi[r])] <- input
            } else {
               par[(pos + 1):(pos + arg$nc_pi[r])] <- runif(arg$nc_pi[r])
            }
            pos <- pos + arg$nc_pi[r]
         }
         cat("Insert TAU\n")
         for (k in seq_len(arg$ntau)) {
            cat(LETTERS[k], ":")
            input <- abs(scan(n = arg$nk_tau[k] * arg$nl_tau[k]))
            if (length(input) == arg$nk_tau[k] * arg$nl_tau[k]) {
               par[(pos + 1):(pos + arg$nk_tau[k] * arg$nl_tau[k])] <-
                  input
            } else {
               par[(pos + 1):(pos + arg$nk_tau[k] * arg$nl_tau[k])] <-
                  runif(arg$nk_tau[k] * arg$nl_tau[k])
            }
            pos <- pos + arg$nk_tau[k] * arg$nl_tau[k]
         }
         cat("Insert RHO\n")
         for (r in seq_len(arg$nrho)) {
            cat(letters[r], ":")
            input <- abs(scan(n = arg$nc_rho[r] * arg$nr_rho[r]))
            if (length(input) == arg$nc_rho[r] * arg$nr_rho[r]) {
               par[(pos + 1):(pos + arg$nc_rho[r] * arg$nr_rho[r])] <-
                  input
            } else {
               par[(pos + 1):(pos + arg$nc_rho[r] * arg$nr_rho[r])] <-
                  runif(arg$nc_rho[r] * arg$nr_rho[r])
            }
            pos <- pos + arg$nc_rho[r] * arg$nr_rho[r]
         }
      }
      par <- unlist(tapply(par, arg$id, norm1), use.names = FALSE)
   }
   sim <- simModel(
      nsim, arg$nvar, arg$nlev, par,
      arg$nlv, arg$nrl, arg$nlf, arg$npi, arg$ntau, arg$nrho,
      arg$ul, arg$vl, arg$lf, arg$rt, arg$eqrl, arg$eqlf,
      arg$nc, arg$nk, arg$nl, arg$ncl,
      arg$nc_pi, arg$nk_tau, arg$nl_tau, arg$nc_rho, arg$nr_rho
   )
   class <- data.frame(sim$class)
   names(class) <- model$latent$label

   # data.name
   y <- data.frame(do.call(cbind, sim$y))
   items <- unlist(model$latent[model$latent$leaf, "children"],
                   use.names = FALSE)
   colnames(y) <- items
   y[] <- lapply(items, function(x)
      factor(y[[x]], labels = level[[x]]))
   mf <- proc_data(y, model, FALSE)

   class <- sim$class + 1
   colnames(class) <- row.names(model$latent)
   rownames(class) <- row.names(mf)

   list(response = mf, class = class)
}
