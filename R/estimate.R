estimate <- function(method, control, data, args, init.param, restr) {
   if (method == "em") {
      if (control$verbose) cat("EM iteration begin.\n")
      em_fit <- emFit(
         data$y, args$nobs, args$nvar, args$ncat,
         args$nlv, args$nroot, args$nlink, args$nleaf,
         args$nlink_unique, args$nleaf_unique,
         args$tree_index - 1, args$root - 1,
         args$u - 1, args$v - 1, args$leaf - 1,
         args$cstr_link - 1, args$cstr_leaf - 1,
         args$nclass, args$nclass_leaf,
         args$nclass_u, args$nclass_v,
         init.param, unlist(restr$restr0),
         control$em.iterlim, control$em.tol,
         control$verbose, control$new.iter
      )
      log_par <- em_fit$params
      em.conv <- em_fit$converged
      nlm.conv <- NA
      if (control$verbose) cat(".. done.\n")
   } else if (method == "nlm") {
      logits <- log2logit(
         init.param, length(unlist(init.param)), args$ncat,
         args$nroot, args$nlink_unique, args$nleaf_unique,
         args$nclass_root, args$nclass_leaf,
         args$nclass_u, args$nclass_v,
         unlist(restr$restr0), unlist(restr$ref) - 1
      )
      if (control$verbose) cat("nlm iteration begin.\n")
      find <- is.finite(logits)
      nlm_fit <- nlm(
         floglik, logits[find], y = data$y,
         nobs = args$nobs, nvar = args$nvar, ncat = args$ncat,
         nlv = args$nlv, nroot = args$nroot,
         nlink = args$nlink, nlink_unique = args$nlink_unique,
         nleaf = args$nleaf, nleaf_unique = args$nleaf_unique,
         root = args$root - 1, ulv = args$u - 1,
         vlv = args$v - 1, cstr_link = args$cstr_link - 1,
         leaf = args$leaf - 1, cstr_leaf = args$cstr_leaf - 1,
         nclass = args$nclass, nclass_leaf = args$nclass_leaf,
         nclass_u = args$nclass_u, nclass_v = args$nclass_v,
         restr0 = unlist(restr$restr0), ref = unlist(restr$ref) - 1,
         iterlim = control$nlm.iterlim, steptol = control$nlm.tol
      )
      logits[find] <- nlm_fit$estimate
      log_par <- logit2log(
         logits, args$ncat, args$nroot,
         args$nlink_unique, args$nleaf_unique,
         args$root - 1, args$u - 1, args$v - 1,
         args$nclass_root, args$nclass_leaf,
         args$nclass_u, args$nclass_v,
         unlist(restr$restr0), unlist(restr$ref) - 1
      )
      em.conv <- NA
      nlm.conv <- nlm_fit$code < 3
      if (control$verbose) cat(".. done.\n")
   } else if (method == "hybrid") {
      if (control$verbose) cat("EM iteration begin.\n")
      em_fit <- emFit(
         data$y, args$nobs, args$nvar, args$ncat,
         args$nlv, args$nroot, args$nlink, args$nleaf,
         args$nlink_unique, args$nleaf_unique,
         args$tree_index - 1, args$root - 1,
         args$u - 1, args$v - 1, args$leaf - 1,
         args$cstr_link - 1, args$cstr_leaf - 1,
         args$nclass, args$nclass_leaf,
         args$nclass_u, args$nclass_v,
         init.param, unlist(restr$restr0),
         control$em.iterlim, control$em.tol,
         control$verbose, control$new.iter
      )
      log_par <- em_fit$params
      logits <- log2logit(
         log_par, length(unlist(log_par)), args$ncat,
         args$nroot, args$nlink_unique, args$nleaf_unique,
         args$nclass_root, args$nclass_leaf,
         args$nclass_u, args$nclass_v,
         unlist(restr$restr0), unlist(restr$ref) - 1
      )
      if (control$verbose) cat(".. done. \nnlm iteration begin.\n")
      find <- is.finite(logits)
      nlm_fit <- nlm(
         floglik, logits[find], y = data$y,
         nobs = args$nobs, nvar = args$nvar, ncat = args$ncat,
         nlv = args$nlv, nroot = args$nroot,
         nlink = args$nlink, nlink_unique = args$nlink_unique,
         nleaf = args$nleaf, nleaf_unique = args$nleaf_unique,
         root = args$root - 1, ulv = args$u - 1,
         vlv = args$v - 1, cstr_link = args$cstr_link - 1,
         leaf = args$leaf - 1, cstr_leaf = args$cstr_leaf - 1,
         nclass = args$nclass, nclass_leaf = args$nclass_leaf,
         nclass_u = args$nclass_u, nclass_v = args$nclass_v,
         restr0 = unlist(restr$restr0), ref = unlist(restr$ref) - 1,
         iterlim = control$nlm.iterlim, steptol = control$nlm.tol
      )
      logits[find] <- nlm_fit$estimate
      log_par <- logit2log(
         logits, args$ncat, args$nroot,
         args$nlink_unique, args$nleaf_unique,
         args$root - 1, args$u - 1, args$v - 1,
         args$nclass_root, args$nclass_leaf,
         args$nclass_u, args$nclass_v,
         unlist(restr$restr0), unlist(restr$ref) - 1
      )

      em.conv  <- em_fit$converged
      nlm.conv <- nlm_fit$code < 3
      if (control$verbose) cat(".. done.\n")
   }

   list(par = log_par, conv = c(EM = em.conv, nlm = nlm.conv))
}
