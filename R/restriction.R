proc_restr <- function(restriction, param, args) {
   restr <- list()
   num_param <- number_par(param, restriction)
   numbers <- num_param$numbers
   param0  <- num_param$param0

   tau_ind <- as.numeric(gsub("tau", "", restriction[grep("tau", restriction)]))
   rho_ind <- as.numeric(gsub("rho", "", restriction[grep("rho", restriction)]))

   restr_tau <- lapply(numbers$tau, apply, 1:2, "%in%", tau_ind)
   restr_rho <- lapply(numbers$rho, apply, 1:2, "%in%", rho_ind)
   restr0 <- list(restr_tau, restr_rho)
   ref <- list(
      args$nclass_root,
      lapply(restr_tau, apply, 2, function(x)
         if (length(which(!x)) == 0) NA
         else max(which(!x))),
      mapply(function(x, y, z) {
         lev <- rep(y, z)
         ind <- rep.int(seq_along(lev), lev)
         tapply(x, ind, function(i)
            if (length(which(!i)) == 0) NA
            else max(which(!i)))
      }, restr_rho, args$ncat, args$nclass_leaf, SIMPLIFY = FALSE)
   )

   list(target = restriction, restr0 = restr0, ref = ref, param0 = param0)
}
