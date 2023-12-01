estModel <- function(method, control, par, mf, arg) {
   llf <- function(
      logit, fix0, ref, id, y, nobs, nvar, nlev, nlv, nrl, nlf,
      npi, ntau, nrho, ul, vl, lf, tr, rt, eqrl, eqlf,
      nc, nk, nl, ncl, nc_pi, nk_tau, nl_tau, nc_rho, nr_rho
   ) {
      logits <- numeric(length(id))
      logits[-c(fix0, ref)] <- logit
      logits[fix0] <- -Inf
      logits[ref] <- 0
      par <- unlist(tapply(logits, id, norm2))

      fll(y, par, nobs, nvar, nlev, nlv, nrl, nlf,
          npi, ntau, nrho, ul, vl, lf, tr, rt, eqrl, eqlf,
          nc, nk, nl, ncl, nc_pi, nk_tau, nl_tau, nc_rho, nr_rho)
   }

   if (method == "em") {
      if (control$verbose) cat("EM iteration begin.\n")
      em <- em_est(
         attr(mf, "y"), arg$nobs, arg$nvar, unlist(arg$nlev), par, arg$fix0,
         arg$nlv, arg$nrl, arg$nlf, arg$npi, arg$ntau, arg$nrho,
         arg$ul, arg$vl, arg$lf, arg$tr, arg$rt, arg$eqrl, arg$eqlf,
         arg$nc, arg$nk, arg$nl, arg$ncl,
         arg$nc_pi, arg$nk_tau, arg$nl_tau, arg$nc_rho, arg$nr_rho,
         control$em.iterlim, control$em.tol, control$verbose, 100
      )
      par <- em$param
      logit <- par - par[arg$ref_idx[arg$id]]
      em.conv <- em$converged
      em.niter <- em$niter
      nlm.conv <- NA
      if (control$verbose) cat(".. done.\n")
   } else if (method == "nlm") {
      if (control$verbose) cat("nlm iteration begin.\n")
      fix0 <- union(which(arg$fix0), which(is.infinite(par)))
      logit <- par - par[arg$ref_idx[arg$id]]
      nonlm <- nlm(
         llf, logit[-c(which(arg$fix0), arg$ref_idx)],
         which(arg$fix0), arg$ref_idx, arg$id, y = attr(mf, "y"),
         nobs = arg$nobs, nvar = arg$nvar, nlev = unlist(arg$nlev),
         nlv = arg$nlv, nrl = arg$nrl, nlf = arg$nlf,
         npi = arg$npi, ntau = arg$ntau, nrho = arg$nrho,
         ul = arg$ul, vl = arg$vl, lf = arg$lf, tr = arg$tr, rt = arg$rt,
         eqrl = arg$eqrl, eqlf = arg$eqlf,
         nc = arg$nc, nk = arg$nk, nl = arg$nl, ncl = arg$ncl,
         nc_pi = arg$nc_pi, nk_tau = arg$nk_tau, nl_tau = arg$nl_tau,
         nc_rho = arg$nc_rho, nr_rho = arg$nr_rho,
         iterlim = control$nlm.iterlim,
         gradtol = control$nlm.tol, steptol = control$nlm.tol
      )
      logit[-c(fix0, arg$ref_idx)] <- nonlm$estimate
      par <- unlist(tapply(logit, arg$id, norm2))
      em.conv <- NA
      nlm.conv <- nlm_fit$code < 3
      if (control$verbose) cat(".. done.\n")
   } else if (method == "hybrid") {
      if (control$verbose) cat("EM iteration begin.\n")
      em <- em_est(
         attr(mf, "y"), arg$nobs, arg$nvar, unlist(arg$nlev), par, arg$fix0,
         arg$nlv, arg$nrl, arg$nlf, arg$npi, arg$ntau, arg$nrho,
         arg$ul, arg$vl, arg$lf, arg$tr, arg$rt, arg$eqrl, arg$eqlf,
         arg$nc, arg$nk, arg$nl, arg$ncl,
         arg$nc_pi, arg$nk_tau, arg$nl_tau, arg$nc_rho, arg$nr_rho,
         control$em.iterlim, control$em.tol, control$verbose, 100
      )
      par <- em$param
      if (control$verbose) cat(".. done. \nnlm iteration begin.\n")
      fix0 <- union(which(arg$fix0), which(is.infinite(par)))
      logit <- par - par[arg$ref_idx[arg$id]]
      nonlm <- nlm(
         llf, logit[-c(fix0, arg$ref_idx)],
         fix0, arg$ref_idx, arg$id, y = attr(mf, "y"),
         nobs = arg$nobs, nvar = arg$nvar, nlev = unlist(arg$nlev),
         nlv = arg$nlv, nrl = arg$nrl, nlf = arg$nlf,
         npi = arg$npi, ntau = arg$ntau, nrho = arg$nrho,
         ul = arg$ul, vl = arg$vl, lf = arg$lf, tr = arg$tr, rt = arg$rt,
         eqrl = arg$eqrl, eqlf = arg$eqlf,
         nc = arg$nc, nk = arg$nk, nl = arg$nl, ncl = arg$ncl,
         nc_pi = arg$nc_pi, nk_tau = arg$nk_tau, nl_tau = arg$nl_tau,
         nc_rho = arg$nc_rho, nr_rho = arg$nr_rho,
         iterlim = control$nlm.iterlim,
         gradtol = control$nlm.tol, steptol = control$nlm.tol
      )
      logit[-c(fix0, arg$ref_idx)] <- nonlm$estimate
      par <- unlist(tapply(logit, arg$id, norm2))

      em.conv  <- em$converged
      nlm.conv <- nonlm$code < 3
      if (control$verbose) cat(".. done.\n")
   }

   list(par = par, logit = logit, conv = c(EM = em.conv, nlm = nlm.conv))
}
