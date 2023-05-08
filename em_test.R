num <- list(nobs = args$nobs, nvar = args$nvar, nresp = args$ncat,
           nlv = args$nlv, nr = args$nroot, nrl = args$nlink, nlf = args$nleaf,
           ntau = args$nlink_unique, nrho = args$nleaf_unique)
idx <- list(ul = args$u - 1, vl = args$v - 1, lf = args$leaf - 1, tr = args$tree_index - 1,
            rt = args$root - 1, eqrl = args$cstr_link - 1, eqlf = args$cstr_leaf - 1)

nclass <- list(nc = args$nclass, nk = args$nclass[args$u], nl = args$nclass[args$v],
               ncl = args$nclass[args$leaf])
dim <- list(nc_pi = args$nclass_root,
            nk_tau = args$nclass_u, nl_tau = args$nclass_v, nc_rho = args$nclass_leaf)

reg = c(FALSE, FALSE, FALSE)
restr0 <- unlist(restr$restr0)
control = list(max_iter = control$em.iterlim,
               tol = control$em.tol,
               verbose = control$verbose,
               newiter = control$new.iter)


emFit1(data$y, num, idx, nclass, dim, unlist(init.param), reg, restr0, control)
