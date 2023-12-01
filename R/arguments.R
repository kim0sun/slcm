arg_mf <- function(model, arg, mf, fix2zero) {
   nobs <- nrow(mf)
   levs <- lapply(arg$vars, apply, 1, function(x)
      unname(levels(mf)[x]), simplify = FALSE)
   lev <- lapply(levs, unique)
   if (!all(lengths(lev) == 1)) {
      stop("levels for variables constrained to equal are different")
   }
   nlev <- lapply(lev, function(x) lengths(x[[1]]))

   # dimension of parameter list
   nc_pi <- arg$nc[arg$rt + 1]
   nk_tau <- arg$nk[!duplicated(arg$eqrl)]
   nl_tau <- arg$nl[!duplicated(arg$eqrl)]
   nc_rho <- arg$ncl[!duplicated(arg$eqlf)]
   nr_rho <- sapply(nlev, sum)

   # parameters
   ref <- c(nc_pi, rep(nk_tau, nl_tau),
            unlist(mapply(rep, nlev, nc_rho)))
   id <- rep(seq_along(ref), ref)

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
   df <- sum(non0) - length(non0)

   arg$nobs <- nobs
   arg$nlev <- nlev
   arg$nc_pi <- nc_pi
   arg$nk_tau <- nk_tau
   arg$nl_tau <- nl_tau
   arg$nc_rho <- nc_rho
   arg$nr_rho <- nr_rho
   arg$id <- id
   arg$df <- df
   arg$fix0 <- fix0
   arg$ref <- ref
   arg$ref_idx <- ref_idx

   arg
}

arg_sim <- function(arg, levels) {
   levs <- lapply(arg$vars, apply, 1, function(x)
      unname(levels[x]), simplify = FALSE)
   lev <- lapply(levs, unique)
   nlev <- lapply(lev, function(x) lengths(x[[1]]))

   # dimension of parameter list
   nc_pi <- arg$nc[arg$rt + 1]
   nk_tau <- arg$nk[!duplicated(arg$eqrl)]
   nl_tau <- arg$nl[!duplicated(arg$eqrl)]
   nc_rho <- arg$ncl[!duplicated(arg$eqlf)]
   nr_rho <- sapply(nlev, sum)

   # parameters
   ref <- c(nc_pi, rep(nk_tau, nl_tau),
            unlist(mapply(rep, nlev, nc_rho)))
   id <- rep(seq_along(ref), ref)

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
   df <- sum(non0) - length(non0)

   arg$nobs <- nobs
   arg$nlev <- nlev
   arg$nc_pi <- nc_pi
   arg$nk_tau <- nk_tau
   arg$nl_tau <- nl_tau
   arg$nc_rho <- nc_rho
   arg$nr_rho <- nr_rho
   arg$id <- id
   arg$df <- df
   arg$fix0 <- fix0
   arg$ref <- ref
   arg$ref_idx <- ref_idx

   arg
}

arguments <- function(model) {
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

   # numbers for data
   nvar <- sapply(vars, ncol)

   list(nlv = nlv, nrl = nrl, nlf = nlf,
        npi = npi, ntau = ntau, nrho = nrho,
        ul = ul - 1, vl = vl - 1, lf = lf - 1,
        tr = tr - 1, rt = rt - 1,
        eqrl = as.numeric(factor(eqrl)) - 1,
        eqlf = as.numeric(factor(eqlf)) - 1,
        nc = nc, nk = nk, nl = nl, ncl = ncl,
        vars = vars, nvar = nvar)
}

get_frame <- function(model, arg, mf) {
   rel <- split(model$struct, arg$eqrl)
   pi <- lapply(seq_len(arg$npi), function(x)
      matrix(ncol = arg$nc_pi[x], dimnames = list(
         "", class = seq_len(arg$nc_pi[x]))))
   names(pi) <- row.names(model$latent)[model$latent$root]
   tau <- lapply(seq_len(arg$ntau), function(x) {
      res <- matrix(nrow = arg$nk_tau[x], ncol = arg$nl_tau[x])
      dimnames(res) <- list(child = seq_len(arg$nk_tau[x]),
                            parent = seq_len(arg$nl_tau[x]))
      rl <- t(rel[[x]][, c(1, 3)])
      colnames(rl) <- rep("", ncol(rl))
      attr(res, "vars") <- rl
      res
   })
   names(tau) <- LETTERS[seq_along(tau)]
   rho <- lapply(seq_len(arg$nrho), function(x) {
      res <- matrix(nrow = arg$nr_rho[x], ncol = arg$nc_rho[x])
      rn <- unlist(sapply(arg$nlev[[x]], seq_len))
      idx <- cumsum(c(1, arg$nlev[[x]][-arg$nvar[x]]))
      rn[idx] <- paste0(rn[idx], "(V", seq_len(arg$nvar[x]), ")")
      dimnames(res) <- list(
         response = rn,
         class = seq_len(arg$nc_rho[x])
      )
      attr(res, "vars") <- arg$vars[[x]]
      res
   })
   names(rho) <- letters[seq_along(rho)]

   skeleton_par <- list(pi = pi, tau = tau, rho = rho)
   skeleton_score <- lapply(
      c(arg$nc_pi, arg$nk_tau * arg$nl_tau, arg$nr_rho * arg$nc_rho),
      function(x) matrix(nrow = x, ncol = arg$nobs)
   )
   skeleton_post <- lapply(arg$nc, function(x)
      matrix(nrow = x, ncol = arg$nobs, dimnames = list(
         class = seq_len(x), dimnames(mf)[[1]])))
   names(skeleton_post) <- row.names(model$latent)
   skeleton_joint <- lapply(seq_len(arg$nrl), function(x)
      array(dim = c(arg$nk[x], arg$nl[x], arg$nobs), dimnames = list(
         child = seq_len(arg$nk[x]), parent = seq_len(arg$nl[x]),
         dimnames(mf)[[1]]
      )))
   if (length(skeleton_joint))
      names(skeleton_joint) <-
      paste(model$struct$parent, "->", model$struct$child)

   list(
      par = skeleton_par, score = skeleton_score,
      post = skeleton_post, joint = skeleton_joint
   )
}
