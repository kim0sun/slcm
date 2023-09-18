#' @export
print.slcm <- function(x, ...) {
   lt <- x$model$latent
   mr <- x$model$measure
   st <- x$model$struct
   cat("Structural Latent Variable Model\n")
   cat("\nLatent variables (Root*):")
   label <- row.names(lt)
   label[lt$root] <- paste0(label[lt$root], "*")
   mat <- rbind(label, lt$nclass)
   dimnames(mat) <- list(c(" Label:", "nclass:"), rep("", ncol(mat)))
   print(mat, quote = FALSE)

   cat("\nMeasurement model:")
   mat <- mr[c("indicator", "constraint")]
   mapping <- function(x) paste("-> {", sapply(x, paste, collapse = ", "), "}")
   mat$indicator <- mapping(lapply(mat$indicator, function(x)
      x[!(x %in% label)]))
   mat$constraint <- paste("", mat$constraint)
   mat <- as.matrix(mat)
   colnames(mat) = rep("", ncol(mat))
   print(mat, quote = FALSE)

   if (nrow(st) > 0) {
      cat("\nStructural model:")
      vars <- tapply(st$child, st$parent, paste, collapse = ", ")
      vars <- vars[!is.na(vars)]
      mat <- cbind("->", paste("{", vars, "}"))
      dimnames(mat) = list(paste0(" ", names(vars)),  rep("", ncol(mat)))
      print(mat, quote = FALSE)

      cat("\nDependency constraints:\n")
      const <- st$constraint
      parent <- split(st$parent, const)
      child <- split(st$child, const)
      mapping <- function(x, y) paste(x, "->", y)
      maps <- mapply(mapping, parent, child, SIMPLIFY = FALSE)
      mat <- matrix("", nrow = max(lengths(maps)), ncol = length(maps))
      for (i in seq_len(length(maps))) {
         mat[seq_len(length(maps[[i]])), i] = maps[[i]]
      }
      dimnames(mat) = list(rep("", nrow(mat)), names(

         maps))
      print(mat, quote = FALSE)
   }
}

#' @export
summary.slcm = function(
      object, ...
) {
   estimated <- inherits(object, "estimated")
   cat("Structural Latent Class Model\n")

   cat("Summary of model:\n")
   lt <- object$model$latent
   mr <- object$model$measure
   st <- object$model$struct
   tr <- object$model$tree

   nvar <- length(setdiff(tr$child, tr$parent))
   nlv <- nrow(lt)
   if (estimated) {
      nobs <- object$arg$nobs
      mat <- rbind(nobs, nvar, nlv)
      dimnames(mat) <- list(
         c(" Number of observations",
           " Number of manifest variables",
           " Number of latent class variables"), ""
      )
   } else {
      mat <- rbind(nvar, nlv)
      dimnames(mat) <- list(
         c(" Number of manifest variables",
           " Number of latent class variables"), ""
      )
   }
   print(mat)

   if (nrow(st) > 0) {
      cat("\nSummary of dependent structure:\n")
      root <- st[st$rank == 1, c(1, 3)]
      names(root) <- c("root", "child1")
      node <- st[, c(1, 3)]
      for (i in seq_len(max(lt$rank) - 1)) {
         names(node) <- paste0("child", c(i, i + 1))
         root <- merge(root, node, by = paste0("child", i),
                       no.dups = TRUE, all.x = TRUE)
      }
      nlev <- row.names(lt)[order(lt$rank, -lt$leaf, lt$nvar)]
      inv_tree <- root[-ncol(root)]
      tree <- inv_tree[rev(seq_len(ncol(inv_tree)))]
      tree[] <- lapply(tree, factor, levels = nlev)
      tree <- tree[do.call(order, tree),]
      tree[] <- lapply(tree, function(x) {
         x[duplicated(x)] <- NA
         x
      })
      mat <- as.matrix(tree)
      dimnames(mat) <- list(rep("", nrow(tree)), rep("", ncol(tree)))
      mat[, -1] <- ifelse(is.na(mat[, -1]), NA, paste0(" -> ", mat[, -1]))
      mat[, 1] <- ifelse(is.na(mat[, 1]), NA, paste0(" ", mat[, 1]))
      print(mat, quote = FALSE, na.print = "")
   }

   if (estimated) {
      cat("\n\nSummary of manifest variables:\n")
      freq <- lapply(object$mf, table)
      lst <- Map(function(x, y) paste0(x, ": ", y),
                 lapply(freq, names), freq)
      mat <- matrix(nrow = nvar, ncol = max(sapply(freq, length)))
      for (i in seq_len(nvar)) {
         mat[i, ] <- lst[[i]]
      }
      dimnames(mat) <- list(paste0(" ", names(freq)), c("", ""))
      print(mat, quote = FALSE, print.gap = 2)

      cat("\n\nSummary of model fit:\n")
      npar <- object$arg$df
      llik <- logLik(object)
      aic <- AIC(object)
      bic <- BIC(object)

      uy <- attr(object$mf, "y_unique")
      mf <- proc_data(uy, object$model, FALSE)
      arg <- arguments(object$model, mf, object$fix2zero)
      ll <- calcModel(
         attr(mf, "y"), arg$nobs, arg$nvar, unlist(arg$nlev),
         object$par, arg$fix0, arg$ref - 1, arg$nlv, arg$nrl, arg$nlf,
         arg$npi, arg$ntau, arg$nrho, arg$ul, arg$vl,
         arg$lf, arg$tr, arg$rt, arg$eqrl, arg$eqlf,
         arg$nc, arg$nk, arg$nl, arg$ncl,
         arg$nc_pi, arg$nk_tau, arg$nl_tau, arg$nc_rho, arg$nr_rho
      )$ll
      mfreq <- exp(rowSums(ll) + log(object$arg$nobs))
      dfreq <- attr(object$mf, "freq")
      chisq <- sum((dfreq - mfreq)^2 / mfreq) + (nobs - sum(mfreq))

      gsq <- attr(object$mf, "loglik") - logLik(object)
      resdf <- attr(object$mf, "df") - npar
      mat <- rbind(npar, llik, NA, aic, bic, NA, resdf,
                   chisq, pchisq(chisq, resdf, lower.tail = FALSE),
                   gsq, pchisq(gsq, resdf, lower.tail = FALSE))
      dimnames(mat) <- list(
         c(" Number of free parameters",
           " Log-likelihood",
           " Information criteria",
           "   Akaike (AIC)",
           "   Bayesian (BIC)",
           " Chi-squared Tests",
           "   Residual degree of freedom (df)",
           "   Pearson Chi-squared (X-squared)",
           "     P(>Chi)",
           "   Likelihood Ratio (G-squared)",
           "     P(>Chi)"
         ), ""
      )
      print(round(mat, 3), na.print = "")
   }
}

#' @export
param <- function(object, ...) UseMethod("param")
#' @export
param.slcm <- function(
      object, what = c("pi", "tau", "rho"),
      type = c("probs", "logit"),
      digits = max(3L, getOption("digits") - 3L),
      index = FALSE, ...
) {
   if (!inherits(object, "estimated")) return(NA)
   what <- match.arg(what, several.ok = TRUE)
   type <- match.arg(type)

   est <- switch(
      type,
      probs = exp(object$par),
      logit = object$logit
   )

   idx <- seq_along(est)
   skeleton <- object$arg$skeleton$par

   val <- sprintf(paste0(" %.", digits, "f"), est)
   if (index) {
      ans <- relist(paste0(val, " (", idx, ")"), skeleton)
   } else ans <- relist(val, skeleton)

   noquote(ans[what], right = TRUE)
}

#' @export
logLik.slcm <- function(object, ...) {
   res <- if (inherits(object, "estimated"))
      structure(sum(object$loglikelihood),
                df = object$arg$df,
                nobs = object$arg$nobs
      ) else NA

   class(res) <- "logLik"
   res
}


#' @method stats::vcov slcm
vcov.slcm <- function(object, type = c("probs", "logit")) {
   if (!inherits(object, "estimated")) return(NA)
   type <- match.arg(type)
   id <- object$arg$id
   score <- object$score
   fi <- crossprod(score)
   vcov <- matrix(0, nrow(fi), ncol(fi), dimnames = dimnames(fi))
   nan <- apply(score, 2, anyNA)
   vcov[!nan, !nan] <- ginv(fi[!nan, !nan])

   if (type == "probs") {
      jac <- bdiag(tapply(object$par, id, jmat))
      vcov <- jac %*% vcov %*% t(jac)
   }

   vcov
}


#' @export
se <- function(object, ...) UseMethod("se")
#' @export
se.slcm <- function(
      object, what = c("pi", "tau", "rho"),
      type = c("probs", "logit"),
      digits = max(3L, getOption("digits") - 3L),
      index = FALSE, ...
) {
   if (!inherits(object, "estimated")) return(NA)
   what <- match.arg(what, several.ok = TRUE)
   type <- match.arg(type)

   vcov <- vcov.slcm(object, type)
   var <- diag(vcov)
   se <- numeric(length(var))
   se[var > 0] <- sqrt(var[var >= 0])
   idx <- seq_along(se)
   skeleton <- object$arg$skeleton$par

   val <- sprintf(paste0(" %.", digits, "f"), se)
   if (index) {
      ans <- relist(paste0(val, " (", idx, ")"), skeleton)
   } else ans <- relist(val, skeleton)

   noquote(ans[what], right = TRUE)
}

#' @export
reorder.slcm = function(x, what, order, ...) {
   latent <- x$model$latent
   const <- x$model$constraint
   leaf <- latent[latent$label == what, "leaf"]
   if (leaf) {
      idx <- const$eqlf[what]
   } else {
      idx <- const$eqlf[what]
   }
}


#' @export
gof <- function(
      object, ...,  test = c("none", "chisq", "boot"), nboot = 100,
      maxiter = 100, tol = 1e-6, verbose = FALSE
) {
   cl <- match.call(expand.dots = FALSE)
   mn <- sapply(c(cl[["object"]], cl[["..."]]), deparse)
   nmodel <- length(mn)
   objects <- list(object, ...)
   test <- match.arg(test)

   df <- sapply(objects, function(x) x$arg$df)
   aic <- sapply(objects, AIC)
   bic <- sapply(objects, BIC)
   ll <- sapply(objects, logLik)
   gsq <- sapply(objects, function(x)
      attr(x$mf, "loglik") - logLik(x))
   resdf <- sapply(objects, function(x)
      attr(x$mf, "df") - x$arg$df)

   tab <- cbind(df, ll, aic, bic, gsq, resdf)
   rownames(tab) <- mn
   dt <- data.frame(tab)
   names(dt) <- c("Df", "logLik", "AIC", "BIC", "Gsq", "Res. Df")

   if (test == "chisq") {
      dt["Pr(>Chi)"] <-
         pchisq(gsq, resdf, lower.tail = FALSE)
   }
   if (test == "boot") {
      pb <- numeric(nmodel)
      if (verbose) cat("START Bootstrap Sampling.\n")
      for (i in seq_len(nmodel)) {
         obj <- objects[[i]]
         lt <- obj$model$latent
         mr <- obj$model$measure
         arg <- obj$arg
         par <- unlist(obj$par)
         gb <- numeric(nboot)

         trace <- paste0(mn[i], "(", i, "/", nmodel, ") : ")
         for (b in 1:nboot) {
            if (verbose) cat("\r", trace, b, "/", nboot, sep = "")
            sim <- simModel(
               arg$nobs, arg$nvar, arg$nlev, par,
               arg$nlv, arg$nrl, arg$nlf, arg$npi, arg$ntau, arg$nrho,
               arg$ul, arg$vl, arg$lf, arg$rt, arg$eqrl, arg$eqlf,
               arg$nc, arg$nk, arg$nl, arg$ncl,
               arg$nc_pi, arg$nk_tau, arg$nl_tau, arg$nc_rho, arg$nr_rho
            )
            y <- data.frame(do.call(cbind, sim$y))
            colnames(y) <- unlist(mr[["indicator"]])
            y[] <- lapply(names(y), function(x) {
               res <- factor(y[[x]])
               levels(res) <- attr(obj$mf, "levels")[[x]]
               res
            })
            mf <- proc_data(y, obj$model)
            con <- obj$control
            con$verbose <- FALSE
            con$em.iterlim <- maxiter; con$nlm.iterlim <- maxiter
            con$em.tol <- tol; con$nlm.tol <- tol
            est <- estModel(obj$method, con, par, mf, arg)
            etc <- calcModel(
               attr(mf, "y"), arg$nobs, arg$nvar, unlist(arg$nlev),
               est$par, arg$fix0, arg$ref - 1, arg$nlv, arg$nrl, arg$nlf,
               arg$npi, arg$ntau, arg$nrho, arg$ul, arg$vl,
               arg$lf, arg$tr, arg$rt, arg$eqrl, arg$eqlf,
               arg$nc, arg$nk, arg$nl, arg$ncl,
               arg$nc_pi, arg$nk_tau, arg$nl_tau, arg$nc_rho, arg$nr_rho
            )
            gb[b] <- attr(mf, "loglik") - sum(etc$ll)
         }
         if (verbose) cat("\n")
         pb[i] <- mean(gb >= gsq[i])
      }
      if (verbose) cat("END.\n")
      dt["Pr(Boot)"] <-pb
   }
   structure(dt, heading = "Analysis of Goodness of Fit Table\n",
             class = c("anova", "data.frame"))
}


#' @export
compare <- function(
      model1, model2, test = c("chisq", "boot"),
      nboot = 100, maxiter = 500, tol = 1e-5,
      verbose = FALSE
) {
   models <- list(model1, model2)
   cl <- match.call()
   test <- match.arg(test)
   name <- c(cl[["model1"]], cl[["model2"]])
   if (any(!sapply(models, inherits, "estimated")))
      stop("both model should be estimated")
   mf1 <- model1$mf
   mf2 <- model2$mf
   attr(mf1, "y") <- NULL
   attr(mf2, "y") <- NULL
   if (!all.equal(mf1, mf2))
      stop("datasets used for models are different")

   df <- sapply(models, function(x) x$arg$df)
   models <- models[order(df)]
   h0 <- models[[1]]; h1 <- models[[2]]
   cat("Model H0:", deparse(name[[which.min(df)]]),
       "\nModel H1:", deparse(name[[which.max(df)]]), "\n\n")

   aic <- sapply(models, AIC)
   bic <- sapply(models, BIC)
   ll <- sapply(models, logLik)
   resdf <- diff(df)

   gsq <- logLik(h1) - logLik(h0)
   tab <- cbind(df, ll, aic, bic, c(NA, gsq), c(NA, resdf))
   rownames(tab) <- unlist(name[order(df)])
   dt <- data.frame(tab)
   names(dt) <- c("Df", "logLik", "AIC", "BIC", "Gsq", "Res. Df")

   if (test == "chisq") {
      dt$`Pr(>Chi)` <- c(NA, pchisq(gsq, resdf, lower.tail = FALSE))
   } else if (test == "boot") {
      gb <- numeric(nboot)
      if (verbose) cat("START Bootstrap Sampling. \n")
      blank <- rep(" ", nchar(nboot))
      for (b in seq_len(nboot)) {
         if (verbose) cat("\r", b, "/", nboot, blank, sep = "")
         arg0 <- h0$arg
         sim <- simModel(
            arg0$nobs, arg0$nvar, arg0$nlev, h0$par,
            arg0$nlv, arg0$nrl, arg0$nlf, arg0$npi, arg0$ntau, arg0$nrho,
            arg0$ul, arg0$vl, arg0$lf, arg0$rt, arg0$eqrl, arg0$eqlf,
            arg0$nc, arg0$nk, arg0$nl, arg0$ncl,
            arg0$nc_pi, arg0$nk_tau, arg0$nl_tau, arg0$nc_rho, arg0$nr_rho
         )
         y <- data.frame(do.call(cbind, sim$y))
         colnames(y) <- unlist(h0$model$measure[["indicator"]])
         y[] <- lapply(names(y), function(x) {
            res <- factor(y[[x]])
            levels(res) <- attr(h0$mf, "levels")[[x]]
            res
         })
         mf0 <- proc_data(y, h0$model)
         mf1 <- proc_data(y, h1$model)

         con <- slcmControl()
         con$verbose <- FALSE
         con$em.iterlim <- maxiter; con$nlm.iterlim <- maxiter
         con$em.tol <- tol; con$nlm.tol <- tol

         arg1 <- h1$arg
         est0 <- estModel(h0$method, con, h0$par, mf0, arg0)
         est1 <- estModel(h1$method, con, h1$par, mf1, arg1)
         etc0 <- calcModel(
            attr(mf0, "y"), arg0$nobs, arg0$nvar, unlist(arg0$nlev),
            est0$par, arg0$fix0, arg0$ref - 1, arg0$nlv, arg0$nrl, arg0$nlf,
            arg0$npi, arg0$ntau, arg0$nrho, arg0$ul, arg0$vl,
            arg0$lf, arg0$tr, arg0$rt, arg0$eqrl, arg0$eqlf,
            arg0$nc, arg0$nk, arg0$nl, arg0$ncl,
            arg0$nc_pi, arg0$nk_tau, arg0$nl_tau, arg0$nc_rho, arg0$nr_rho
         )
         etc1 <- calcModel(
            attr(mf1, "y"), arg1$nobs, arg1$nvar, unlist(arg1$nlev),
            est1$par, arg1$fix0, arg1$ref - 1, arg1$nlv, arg1$nrl, arg1$nlf,
            arg1$npi, arg1$ntau, arg1$nrho, arg1$ul, arg1$vl,
            arg1$lf, arg1$tr, arg1$rt, arg1$eqrl, arg1$eqlf,
            arg1$nc, arg1$nk, arg1$nl, arg1$ncl,
            arg1$nc_pi, arg1$nk_tau, arg1$nl_tau, arg1$nc_rho, arg1$nr_rho
         )
         gb[b] <- sum(etc1$ll) - sum(etc0$ll)
      }
      dt$`Pr(Boot)` <- c(NA, mean(gb >= gsq))
   }

   structure(dt, heading = "Analysis of Relative Model Fit\n",
             class = c("anova", "data.frame"))
}


#' @export
predict.slcm <- function(
      object, newdata, label, type = c("class", "posterior"), ...
) {
   type <- match.arg(type)
   if (!inherits(object, "estimated")) stop("Latent variable model should be estimated.")
   if (missing(label)) label <- row.names(object$model$latent)
   if (missing(newdata)) post <- object$posterior
   else {
      mf <- proc_data(newdata, object$model, FALSE)
      arg <- arguments(object$model, mf, object$fix2zero)

      post <- calcModel(
         attr(mf, "y"), arg$nobs, arg$nvar, unlist(arg$nlev),
         par, arg$fix0, arg$ref - 1, arg$nlv, arg$nrl, arg$nlf,
         arg$npi, arg$ntau, arg$nrho, arg$ul, arg$vl,
         arg$lf, arg$tr, arg$rt, arg$eqrl, arg$eqlf,
         arg$nc, arg$nk, arg$nl, arg$ncl,
         arg$nc_pi, arg$nk_tau, arg$nl_tau, arg$nc_rho, arg$nr_rho
      )$post
      post <- relist(post, arg$skeleton$post)
   }
   impute <- function(x) apply(x, 1, which.max)

   switch(
      type,
      class = if (length(label) == 1)
         impute(post$marginal[[label]])
      else sapply(post$marginal[label], impute),
      posterior = if (length(label) == 1)
         post$marginal[[label]]
      else post$marginal[label]
   )
}


#' Confidence Intervals for Model Parameters
#'
#' Computes confidence intervals for one or more parameters of fitted model. Package \pkg{slcm} adds methods for \code{slcm} fits.
#'
#' @param object a fitted slcm object.
#' @param parm a specification of which parameters are to be given confidence intervals, either a vector of numbers or a vector of names. If missing, all parameters are considered.
#' @param level the confidence level required.
#' @param method calculating method for standard error
#' @param out the type of output
#'
#' @export
confint.slcm <- function(
      object, parm = c("pi", "tau", "rho"), level = 0.95,
      method = c("asymp", "logit"), out = c("param", "logit")
) {
   parm <- match.arg(parm, several.ok = TRUE)
   out <- match.arg(out)

   logit <- object$logit
   vcov <- vcov(object, "logit")
   se <-

      lower <- (1 - level) / 2
   upper <- 1 - lower
   cn <- format.pc(c(lower, upper), 3)

   ci_logit <- logit + se %o% qnorm(c(lower, upper))
   args <- object$args
   restr <- object$restriction

   switch(out, param = {
      ci_par <- sapply(seq_len(ncol(ci_logit)), function(i) {
         cii <- logit2log(
            ci_logit[, i], args$ncat, args$nroot,
            args$nlink_unique, args$nleaf_unique,
            args$root - 1, args$u - 1, args$v - 1,
            args$nclass_root, args$nclass_leaf,
            args$nclass_u, args$nclass_v,
            unlist(restr$restr0), unlist(restr$ref) - 1
         )
         unlist(output_param(cii, object$model, args))
      })
      ci <- t(apply(ci_par, 1, sort))
      colnames(ci) <- cn
      ci
   }, logit = {
      ci <- sapply(seq_len(ncol(ci_logit)), function(i) {
         cii <- param2list(
            ci_logit[, i], args$ncat, args$nroot,
            args$nlink_unique, args$nleaf_unique,
            args$root, args$u, args$v, args$nclass_root,
            args$nclass_u, args$nclass_v, args$nclass_leaf
         )
         unlist(output_param(cii, object$model, args, FALSE))
      })
      colnames(ci) <- cn
      ci
   })
}

format.pc <- function(perc, digits)
   paste(format(100 * perc, trim = TRUE, scientific = FALSE, digits = digits), "%")

