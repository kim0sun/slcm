arguments <- function(model, mf, fix2zero) {
   lt <- model$latent
   mr <- model$measure
   st <- model$struct
   label <- row.names(lt)
   leaf <- row.names(mr)

   lf <- factor(leaf, levels = label)
   rt <- factor(label[lt$root], levels = label)
   tr <- factor(st$root, levels = unique(st$root))

   # index vector for constraints
   eqrl <- st$constraint
   eqlf <- mr$constraint

   child <- mr$indicator
   names(child) <- leaf
   var <- split(child, eqlf)
   levs <- lapply(var, lapply, function(x)
      unname(levels(mf)[x]))
   lev <- lapply(levs, unique)
   if (!all(lengths(lev) == 1)) {
      stop("levels for constrained to equal are not equal")
   }
   nlev <- lapply(lev, function(x) lengths(x[[1]]))

   vars <- lapply(var, function(x) {
      res <- do.call(rbind, x)
      dimnames(res) <- list(
         names(x), paste0("V", seq_len(ncol(res)))
      )
      res
   })
   rel <- split(st, eqrl)

   # number of node, edge, leaf
   nlv <- nrow(lt)
   nlf <- nrow(mr)
   nrl <- nrow(st)

   # numbers for data
   nobs <- nrow(mf)
   nvar <- lengths(nlev)

   # number of parameter list
   npi <- length(rt)
   ntau <- length(unique(eqrl))
   nrho <- length(unique(eqlf))

   # index vector for model
   ul <- as.numeric(st$child)
   vl <- as.numeric(st$parent)
   lf <- as.numeric(lf)
   tr <- as.numeric(tr)
   rt <- as.numeric(rt)

   # number of classes for index vector
   nc <- lt$nclass
   nk <- lt$nclass[st$child]
   nl <- lt$nclass[st$parent]
   ncl <- mr$nclass

   # dimension of parameter list
   nc_pi <- lt[rt, "nclass"]
   nk_tau <- lt[st[!duplicated(eqrl), "child"], "nclass"]
   nl_tau <- lt[st[!duplicated(eqrl), "parent"], "nclass"]
   nc_rho <- mr[!duplicated(eqlf), "nclass"]
   nr_rho <- sapply(nlev, sum)

   # parameters
   ref <- c(nc_pi, rep(nk_tau, nl_tau),
            unlist(mapply(rep, nlev, nc_rho)))
   id <- rep(seq_along(ref), ref)

   pi <- lapply(seq_len(npi), function(x)
      matrix(ncol = nc_pi[x], dimnames = list(
         "", class = seq_len(nc_pi[x]))))
   names(pi) <- label[rt]
   tau <- lapply(seq_len(ntau), function(x) {
      res <- matrix(nrow = nk_tau[x], ncol = nl_tau[x])
      dimnames(res) <- list(child = seq_len(nk_tau[x]),
                            parent = seq_len(nl_tau[x]))
      rl <- rel[[x]][,1:3]
      row.names(rl) <- NULL
      attr(res, "vars") <- rl
      res
   })
   names(tau) <- LETTERS[seq_along(tau)]
   rho <- lapply(seq_len(nrho), function(x) {
      res <- matrix(nrow = nr_rho[x], ncol = nc_rho[x])
      rn <- unlist(lev[[x]])
      idx <- cumsum(c(1, nlev[[x]][-nvar[x]]))
      rn[idx] <- paste0(rn[idx], "(V", seq_len(nvar[x]), ")")
      dimnames(res) <- list(
         response = rn,
         class = seq_len(nc_rho[x])
      )
      attr(res, "vars") <- vars[[x]]
      res
   })
   names(rho) <- letters[seq_along(rho)]

   skeleton_par <- list(pi = pi, tau = tau, rho = rho)
   skeleton_score <- lapply(c(nc_pi, nk_tau * nl_tau, nr_rho * nc_rho),
                            function(x) matrix(nrow = x, ncol = nobs))
   skeleton_post <- lapply(nc, function(x)
      matrix(nrow = x, ncol = nobs, dimnames = list(
         class = seq_len(x), dimnames(mf)[[1]])))
   names(skeleton_post) <- label
   skeleton_joint <- lapply(seq_len(nrl), function(x)
      array(dim = c(nk[x], nl[x], nobs), dimnames = list(
         child = seq_len(nk[x]), parent = seq_len(nl[x]),
         dimnames(mf)[[1]]
      )))
   if (length(skeleton_joint))
      names(skeleton_joint) <- paste(st$parent, "->", st$child)
   par_index <- relist(paste0("(", seq_along(id), ")"), skeleton_par)

   fix0 <- logical(length(id))
   fix0[fix2zero] <- TRUE

   non0 <- tapply(!fix0, id, sum)
   par_n <- tapply(seq_along(id), id, paste0, collapse = ", ")
   if (!all(non0 > 0)) {
      stop("all parameters in the following parameter set(s) cannot be fixed to zero:\n",
           paste0("(", par_n[non0 == 0], ")", collapse = ", "))
   }

   ref_idx <- cumsum(ref)
   while (any(cond <- ref_idx %in% fix2zero)) {
      ref[cond] <- ref[cond] - 1
      ref_idx[cond] <- ref_idx[cond] - 1
   }
   df <- length(id) - length(ref)

   list(nlv = nlv, nrl = nrl, nlf = nlf,
        npi = npi, ntau = ntau, nrho = nrho,
        ul = ul - 1, vl = vl - 1, lf = lf - 1,
        tr = tr - 1, rt = rt - 1,
        eqrl = as.numeric(factor(eqrl)) - 1,
        eqlf = as.numeric(factor(eqlf)) - 1,
        nc = nc, nk = nk, nl = nl, ncl = ncl,
        nobs = nobs, nvar = nvar, nlev = nlev,
        nc_pi = nc_pi, nk_tau = nk_tau, nl_tau = nl_tau,
        nc_rho = nc_rho, nr_rho = nr_rho,
        par_index = par_index, id = id, df = df,
        fix0 = fix0, ref = ref, ref_idx = ref_idx,
        skeleton = list(par = skeleton_par,
                        score = skeleton_score,
                        post = skeleton_post,
                        joint = skeleton_joint))
}


arguments_nmf <- function(model, nobs, nlevel) {
   lt <- model$latent
   mr <- model$measure
   st <- model$struct
   label <- row.names(lt)

   lf <- factor(leaf, levels = label)
   rt <- factor(label[lt$root], levels = label)
   tr <- factor(st$root, levels = unique(st$root))

   # index vector for constraints
   eqrl <- st$constraint
   eqlf <- lt$constraint[lt$leaf]

   child <- mr$indicator
   names(child) <- leaf
   var <- split(child, eqlf)
   vars <- lapply(var, function(x) {
      res <- do.call(rbind, x)
      dimnames(res) <- list(
         names(x), paste0("V", seq_len(ncol(res)))
      )
      res
   })
   levs <- lapply(var, lapply, function(x) unname(level[x]))
   lev <- lapply(levs, unique)
   if (!all(lengths(lev) == 1)) {
      stop("levels for constrained to equal are not equal")
   }
   nlev <- lapply(lev, function(x) lengths(x[[1]]))

   # number of node, edge, leaf
   nlv <- nrow(lt)
   nrl <- nrow(st)
   nlf <- length(lf)

   # numbers for data
   nvar <- lengths(nlev)

   # number of parameter list
   npi <- length(rt)
   ntau <- length(unique(eqrl))
   nrho <- length(unique(eqlf))

   # index vector for model
   ul <- as.numeric(st$child)
   vl <- as.numeric(st$parent)
   lf <- as.numeric(lf)
   tr <- as.numeric(tr)
   rt <- as.numeric(rt)

   # number of classes for index vector
   nc <- lt$nclass
   nk <- lt$nclass[st$child]
   nl <- lt$nclass[st$parent]
   ncl <- lt$nclass[lf]

   # dimension of parameter list
   nc_pi <- lt[rt, "nclass"]
   nk_tau <- lt[st[!duplicated(eqrl), "child"], "nclass"]
   nl_tau <- lt[st[!duplicated(eqrl), "parent"], "nclass"]
   nc_rho <- lt[lf[!duplicated(eqlf)], "nclass"]
   nr_rho <- sapply(nlev, sum)

   # parameters
   pi <- lapply(seq_len(npi), function(x)
      matrix(ncol = nc_pi[x], dimnames = list("", class = seq_len(nc_pi[x]))))
   names(pi) <- label[rt]
   tau <- lapply(seq_len(ntau), function(x)
      matrix(nrow = nk_tau[x], ncol = nl_tau[x],
             dimnames = list(child = seq_len(nk_tau[x]),
                             parent = seq_len(nl_tau[x]))))
   names(tau) <- LETTERS[seq_along(tau)]
   rho <- lapply(seq_len(nrho), function(x) {
      res <- matrix(nrow = nr_rho[x], ncol = nc_rho[x])
      rn <- unlist(lev[[x]])
      idx <- cumsum(c(1, nlev[[x]][-nvar[x]]))
      rn[idx] <- paste0(rn[idx], "(V", seq_len(nvar[x]), ")")
      dimnames(res) <- list(
         response = rn,
         class = seq_len(nc_rho[x])
      )
      attr(res, "vars") <- vars[[x]]
      res
   })
   names(rho) <- letters[seq_along(rho)]
   skeleton_par <- list(pi = pi, tau = tau, rho = rho)

   list(nlv = nlv, nrl = nrl, nlf = nlf,
        npi = npi, ntau = ntau, nrho = nrho,
        ul = ul - 1, vl = vl - 1, lf = lf - 1,
        tr = tr - 1, rt = rt - 1,
        eqrl = as.numeric(factor(eqrl)) - 1,
        eqlf = as.numeric(factor(eqlf)) - 1,
        nc = nc, nk = nk, nl = nl, ncl = ncl,
        nobs = nobs, nvar = nvar, nlev = nlev,
        nc_pi = nc_pi, nk_tau = nk_tau, nl_tau = nl_tau,
        nc_rho = nc_rho, nr_rho = nr_rho,
        skeleton = list(par = skeleton_par))
}
