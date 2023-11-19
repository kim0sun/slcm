#' @exportS3Method base::print slcm
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

#' @exportS3Method base::summary slcm
summary.slcm = function(
      object, ...
) {
   estimated <- inherits(object, "estimated")
   cat("Structural Latent Class Model\n")

   cat("\nSummary of model:\n")
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

   cat("\nSummary of model structure\n")
   cat("\n Latent variables (Root*):")
   label <- row.names(lt)
   label[lt$root] <- paste0(label[lt$root], "*")
   mat <- rbind(label, lt$nclass)
   dimnames(mat) <- list(c("  Label:", " nclass:"), rep("", ncol(mat)))
   print(mat, quote = FALSE)

   cat("\n Measurement model:")
   mat <- mr[c("indicator", "constraint")]
   mapping <- function(x) paste("-> {", sapply(x, paste, collapse = ", "), "}")
   mat$indicator <- mapping(lapply(mat$indicator, function(x)
      x[!(x %in% label)]))
   mat$constraint <- paste("", mat$constraint)
   mat <- as.matrix(mat)
   colnames(mat) = rep("", ncol(mat))
   rownames(mat) = paste0("  ", rownames(mat))
   print(mat, quote = FALSE)

   if (nrow(st) > 0) {
      cat("\n Structural model:")
      vars <- tapply(st$child, st$parent, paste, collapse = ", ")
      vars <- vars[!is.na(vars)]
      mat <- cbind("->", paste("{", vars, "}"))
      dimnames(mat) = list(paste0("  ", names(vars)),  rep("", ncol(mat)))
      print(mat, quote = FALSE)

      cat("\n Dependency constraints:\n")
      const <- st$constraint
      parent <- split(st$parent, const)
      child <- split(st$child, const)
      mapping <- function(x, y) paste(x, "->", y)
      maps <- mapply(mapping, parent, child, SIMPLIFY = FALSE)
      mat <- matrix("", nrow = max(lengths(maps)), ncol = length(maps))
      for (i in seq_len(length(maps))) {
         mat[seq_len(length(maps[[i]])), i] = maps[[i]]
      }
      dimnames(mat) = list(rep(" ", nrow(mat)), names(

         maps))
      print(mat, quote = FALSE)
   }

   if (nrow(st) > 0) {
      cat("\n Tree of structural model:")
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
      cat("\nSummary of manifest variables:\n")
      freq <- lapply(object$mf, table)
      lst <- Map(function(x, y) paste0(x, ": ", y),
                 lapply(freq, names), freq)
      mat <- matrix(nrow = nvar, ncol = max(sapply(freq, length)))
      for (i in seq_len(nvar)) {
         mat[i, ] <- lst[[i]]
      }
      dimnames(mat) <- list(paste0(" ", names(freq)), c("", ""))
      print(mat, quote = FALSE, print.gap = 2)

      cat("\nSummary of model fit:\n")
      npar <- object$arg$df
      llik <- logLik(object)
      aic <- AIC(object)
      bic <- BIC(object)

      arg <- object$arg
      mf <- object$mf
      ll <- calcModel(
         attr(mf, "yu"), nrow(attr(mf, "y_unique")),
         arg$nvar, unlist(arg$nlev), object$par, arg$fix0,
         arg$ref - 1, arg$nlv, arg$nrl, arg$nlf,
         arg$npi, arg$ntau, arg$nrho, arg$ul, arg$vl,
         arg$lf, arg$tr, arg$rt, arg$eqrl, arg$eqlf,
         arg$nc, arg$nk, arg$nl, arg$ncl,
         arg$nc_pi, arg$nk_tau, arg$nl_tau, arg$nc_rho, arg$nr_rho
      )$ll
      mfreq <- exp(rowSums(ll) + log(arg$nobs))
      dfreq <- attr(object$mf, "freq")
      chisq <- sum((dfreq - mfreq)^2 / mfreq) +
         (arg$nobs - sum(mfreq))

      gsq <- 2 * (attr(object$mf, "loglik") - logLik(object))
      resdf <- attr(object$mf, "df") - npar
      sprintf("%.0f", npar)
      mat <- rbind(sprintf("%.0f", npar),
                   sprintf("%.3f", llik), NA,
                   sprintf("%.3f", aic),
                   sprintf("%.3f", bic), NA,
                   sprintf("%.0f", resdf),
                   sprintf("%.3f", chisq),
                   sprintf("%.3f", pchisq(chisq, resdf, lower.tail = FALSE)),
                   sprintf("%.3f", gsq),
                   sprintf("%.3f", pchisq(gsq, resdf, lower.tail = FALSE)))
      format(pchisq(gsq, resdf, lower.tail = FALSE), digits = 3)
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
      print(mat, na.print = "", quote = FALSE, right = TRUE)
   }
}

#' Printing Estimated Parameters of `slcm` Object
#'
#' By passing `estimated` `slcm` object, the function prints estimated parameters of the slcm model.
#'
#' @aliases param param.slcm se se.slcm
#'
#' @usage
#' param(object, ...)
#'
#' \method{param}{slcm}(
#'    object, type = c("probs", "logit"),
#'    se = FALSE, index = FALSE ..
#' )
#'
#' @param object an object of class `slcm` and `estimated`.
#' @param what specifies which parameter types to display. Valid options are `"pi"`, `"tau"`, and `"rho"`. `pi` represents membership probabilities of root variable, `tau` denotes conditional probabilities between latent class variables, and `rho` corresponds to item response probabilities for each measurement model. Multiple types can be chosen to display.
#' @param type a character string indicating the format in which the estimated parameters should be returned. If set to `"probs"`, estimates are returned in probability form, while `"logit"` returns them in log-odds (logit) form. Default is `"probs"`.
#' @param index a logical value. If `TRUE`, the indices of the estimated parameters are included in the output, enclosed in parentheses. If `FALSE`, they are omitted.
#' @param digits an integer that sets the number of decimal places for rounding the output.
#'
#' @returns
#' A list containing the specified estimated parameters (or corresponding standard errors for `se` function):
#' \item{pi}{Membership probabilities of the root variable.}
#' \item{tau}{Conditional probabilities between latent class variables, represented with uppercase alphabets for considering measurement invariance.}
#' \item{rho}{Item response probabilities for each measurement model, represented with lowercase alphabets for considering measurement invariance.}
#'
#' @export
param <- function(object, ...) UseMethod("param")
#' @exportS3Method slcm::param slcm
param.slcm <- function(
   object, type = c("probs", "logit"),
   se = FALSE, index = FALSE, ...
) {
   if (!inherits(object, "estimated")) return(NA)
   type <- match.arg(type)

   if (se) {
      vcov <- vcov.slcm(object, type)
      var <- diag(vcov)
      est <- numeric(length(var))
      est[var >= 0] <- sqrt(var[var >= 0])
   } else {
      est <- switch(
         type,
         probs = exp(object$par),
         logit = object$logit,
      )
   }

   skeleton <- object$skeleton$par

   res <- relist(est, skeleton)
   attr(res, "index") <- index
   if (index) {
      attr(res, "idx") <- seq_along(est)
   }
   class(res) <- c("param.slcm", "list")
   res
}

#' @exportS3Method base::print param.slcm
print.param.slcm <- function(
   x, digits = max(3L, getOption("digits") - 3L), ...
) {
   class(x) <- "list"
   index <- attr(x, "index")
   val <- unlist(x)
   val <- sprintf(paste0(" %.", digits, "f"), val)
   if (index) {
      idx <- attr(x, "idx")
      ans <- relist(paste0(val, " (", idx, ")"), x)
   } else ans <- relist(val, x)

   cat("PI :\n")
   for (i in names(ans[["pi"]])) {
      cat(paste0("(", i, ")\n"))
      print(ans[["pi"]][[i]], right = TRUE, quote = FALSE)
   }
   cat("\nTAU :\n")
   for (i in names(ans[["tau"]])) {
      cat(paste0("(", i, ")\n"))
      print(ans[["tau"]][[i]], right = TRUE, quote = FALSE)
      print(attr(x[["tau"]][[i]], "vars"), quote = FALSE)
   }
   cat("\nRHO :\n")
   for (i in names(ans[["rho"]])) {
      cat(paste0("(", i, ")\n"))
      print(ans[["rho"]][[i]], right = TRUE, quote = FALSE)
      cat("\n")
      print(attr(x[["rho"]][[i]], "vars"), quote = FALSE)
   }
}

#' @exportS3Method stats::logLik slcm
logLik.slcm <- function(object, ...) {
   res <- if (inherits(object, "estimated"))
      structure(sum(object$loglikelihood),
                df = object$arg$df,
                nobs = object$arg$nobs
      ) else NA

   class(res) <- "logLik"
   res
}


#' @exportS3Method stats::vcov slcm
vcov.slcm <- function(object, type = c("probs", "logit")) {
   if (!inherits(object, "estimated")) return(NA)
   type <- match.arg(type)
   id <- object$arg$id
   score <- object$score
   fi <- crossprod(score)
   vcov <- matrix(0, nrow(fi), ncol(fi), dimnames = dimnames(fi))
   nan <- apply(score, 2, anyNA)
   vcov[!nan, !nan] <- MASS::ginv(fi[!nan, !nan])

   if (type == "probs") {
      jac <- bdiag(tapply(object$par, id, jmat))
      vcov <- jac %*% vcov %*% t(jac)
   }

   vcov
}

#' Goodness of Fit Tests for Estimated `slcm` Model
#'
#' Provides AIC, BIC and deviance statistic (G-squared) for goodness of fit test for the fitted model. Absolute model fit can be tested with deviance statistics, if `test` argument is specified.
#'
#' @aliases gof gof.slcm
#'
#' @usage
#' gof(object, ...)
#'
#' \method{gof}{slcm}(
#'    object, ..., test = c("none", "chisq", "boot"), nboot = 100,
#'    maxiter = 100, tol = 1e-6, verbose = FALSE
#' )
#'
#' @param object an object of class `slcm` and `estimated`.
#' @param ... additional objects of class `slcm` and `estimated`.
#' @param test a character string specifying the type of test to be conducted. If "none", no test is conducted. If "chisq", a chi-squared test is conducted. If "boot", a bootstrap test is conducted.
#' information on the bootstrap procedure.
#' @param nboot an integer specifying the number of bootstrap rounds to be performed.
#' @param maxiter an integer specifying maximum number of iterations allowed for the estimation process of each bootstrapping round.
#' @param tol a numeric value setting tolerance for the convergence of each bootstrapping round.
#' @param verbose a logical value indicating whether to print progress updates on the number of bootstrapping rounds completed.
#'
#' @returns
#' A `data.frame` containing the number of parameters (Df), loglikelihood, AIC, BIC, G-squared statistics, and the residual degree of freedom for each object.
#' Depending on the `test` argument, the p-value for the corresponding statistical test may also be included.
#'
#' @seealso \link[slcm]{compare}
#'
#' @export
gof <- function(object, ...) UseMethod("gof")

#' @exportS3Method slcm::gof slcm
gof.slcm <- function(
   object, ...,  test = c("none", "chisq", "boot"), nboot = 100,
   maxiter = 100, tol = 1e-6, verbose = FALSE
) {
   cl <- match.call(expand.dots = FALSE)
   objects <- list(object, ...)
   est <- sapply(objects, inherits, "estimated")
   mn <- sapply(c(cl[["object"]], cl[["..."]]), deparse)
   if (all(!est)) stop("at least 1 model should be estimated")
   objects <- objects[est]
   mn <- mn[est]
   nmodel <- length(mn)

   test <- match.arg(test)

   df <- sapply(objects, function(x) x$arg$df)
   aic <- sapply(objects, AIC)
   bic <- sapply(objects, BIC)
   ll <- sapply(objects, logLik)
   gsq <- sapply(objects, function(x)
      2 * (attr(x$mf, "loglik") - logLik(x)))
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
            mf <- proc_data(y, obj$model, obj$control$na.rm)
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
            gb[b] <- 2 * (attr(mf, "loglik") - sum(etc$ll))
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

#' Comparing Two Estimated `slcm` Models
#'
#'
#'
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

   gsq <- 2 * (logLik(h1) - logLik(h0))
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
         mf0 <- proc_data(y, h0$model, h0$control$na.rm)
         mf1 <- proc_data(y, h1$model, h1$control$na.rm)

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
         gb[b] <- 2 * (sum(etc1$ll) - sum(etc0$ll))
      }
      dt$`Pr(Boot)` <- c(NA, mean(gb >= gsq))
   }

   structure(dt, heading = "Analysis of Relative Model Fit\n",
             class = c("anova", "data.frame"))
}

#' @exportS3Method stats::predict slcm
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
      post <- relist(post, skeleton$post)
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
#' @param object an object of class `slcm` and `estimated`.
#' @param level level a numeric value representing the desired confidence level for the intervals. Default is 0.95
#' @param type a character string specifying the format in which the results should be returned. Options are `"probs"` for probability format and `"logit"` for log-odds (logit) format. The default is `"probs"`.
#'
#'
#' @exportS3Method stats::confint slcm
confint.slcm <- function(
   object, level = 0.95,
   type = c("param", "logit"), ...
) {
   if (!inherits(object, "estimated")) return(NA)
   type <- match.arg(type)

   logit <- object$logit
   vcov <- vcov(object, "logit")
   vars <- diag(vcov)
   se <- vars
   se[] <- 0
   se[vars >= 0] <- sqrt(vars[vars >= 0])

   lower <- (1 - level) / 2
   upper <- 1 - lower
   cn <- format.pc(c(lower, upper), 3)

   ci_logit <- logit + se %o% qnorm(c(lower, upper))
   ci_par <- ci_logit
   ci_par[] <- exp(unlist(apply(ci_logit, 2, tapply, object$arg$id, norm2)))

   id <- unlist(object$arg$par_index)
   ci <- switch(type, param = ci_par[id,], logit = ci_logit[id,])
   colnames(ci) <- cn
   ci
}

format.pc <- function(perc, digits)
   paste(format(100 * perc, trim = TRUE, scientific = FALSE, digits = digits), "%")

