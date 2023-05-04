output_param <- function(param, model, args, log = TRUE) {
   if (log) {
      pi <- lapply(param$pi, exp)
      tau <- lapply(param$tau, exp)
      rho <- lapply(param$rho, exp)
   } else {
      pi <- param$pi
      tau <- param$tau
      rho <- param$rho
   }

   names(pi) <- model$root
   for (r in seq_len(args$nroot)) {
      pi[[r]] <- matrix(pi[[r]], ncol = length(pi[[r]]))
      dimnames(pi[[r]]) <- list(
         "", class = seq_len(length(pi[[r]]))
      )
   }
   names(tau) <- LETTERS[seq_len(args$nlink_unique)]
   for (d in seq_len(args$nlink_unique)) {
      dimnames(tau[[d]]) <- list(
         seq_len(args$nclass_u[d]),
         seq_len(args$nclass_v[d])
      )
      names(dimnames(tau[[d]])) <- c("child", "parent")
   }

   names(rho) <- letters[seq_len(args$nleaf_unique)]
   var <- split(model$leaf, args$cstr_leaf)
   item <- split(model$vars$manifest, args$cstr_leaf)
   for (v in seq_len(args$nleaf_unique)) {
      rho[[v]] <- matrix(rho[[v]], ncol = args$nclass_leaf[v])
      dimnames(rho[[v]]) <- list(
         reponse = sapply(args$ncat[[v]], seq_len),
         class = 1:ncol(rho[[v]])
      )

      var_index <- cumsum(c(1, args$ncat[[v]][-args$nvar[v]]))
      rownames(rho[[v]])[var_index] <-
         paste0("1 (item ", seq(var_index), ")")
   }

   return(list(pi = pi, tau = tau, rho = rho))
}


output_posterior <- function(post, model, data) {
   names(post) = model$label
   lapply(post, function(x) {
      dimnames(x) <- list(data$dimnames[[1]], class = seq(ncol(x)))
      exp(x)
   })
}

