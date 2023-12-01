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
   dimnames(mat) = list(
      paste0(" ", row.names(mat)), rep("", ncol(mat))
   )
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
      dimnames(mat) = list(rep("", nrow(mat)), names(maps))
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

