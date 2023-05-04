#' @export
print.slcm <- function(x, digits = 5, ...) {
   cat("Structural Latent Variable Model\n")
   cat("\nLatent variables (Root*) :")
   label <- x$model$label
   label[x$args$root] <- paste0(label[x$args$root], "*")
   mat <- rbind(label, x$model$nclass)
   dimnames(mat) <- list(c(" Label:", "nclass:"), rep("", ncol(mat)))
   print(mat, quote = FALSE)

   cat("\nMeasurement model:\n")
   vars <- x$model$vars$manifest
   collapsed <- sapply(vars, paste, collapse = ", ")
   formula <- paste(names(vars), "-> {", collapsed, "}")
   constr <- letters[x$model$constr$cstr_leaf]
   mat <- rbind(cbind(formula, " ", constr), "")
   dimnames(mat) <- list(rep("", nrow(mat)), mat[1,])
   print(mat[-1, , drop = FALSE], quote = FALSE)

   if (x$args$nlink > 0) {
      cat("Latent dependent structure:\n")
      vars <- x$model$vars$latent
      formula <- sapply(vars, paste, collapse = ", ")
      fmlae <- paste(names(vars), "-> {", formula, "}")
      cat("", paste(fmlae, collapse = "\n "), "\n")

      cstr <- LETTERS[x$model$constr$cstr_link]
      parent <- split(x$model$links$parent, cstr)
      child <- split(x$model$links$child, cstr)
      mapping <- function(x, y) paste(x, "->", y)
      maps <- mapply(mapping, parent, child, SIMPLIFY = FALSE)
      mat <- matrix("", nrow = max(sapply(maps, length)), ncol = length(maps))
      for (i in seq_len(length(maps))) {
         mat[seq_len(length(maps[[i]])), i] = maps[[i]]
      }
      cat("\nDependency constraints:\n")
      dimnames(mat) = list(rep("", nrow(mat)), names(maps))
      print(mat, quote = FALSE)
   }
}

#' @export
summary.slcm = function(object, type = c("data", "model", "posterior", "parameter"), ...) {
   cat("Summary of Structural Latent Class Model\n")

   if ("data" %in% type) {
      cat("Summary of manifest items\n")
   }
   if ("model" %in% type) {
      cat("Summary of model\n")
   }
   if ("posterior" %in% type) {
      cat("Summary of posterior\n")
   }
   if ("parameter" %in% type) {
      cat("Summary of parameter estimation\n")
   }
}

#' @export
logLik.slcm <- function(object, ...) {
   res <- if (is.null(object$estimate)) NA
   else structure(
      sum(object$llik),
      df = sum(object$args$npar),
      nobs = object$data$nobs
   )
   class(res) <- "logLik"
   res
}

# #' @export
# reorder.slcm = function(x, ...) {
#
# }

# #' @export
# anova.slcm = function(x, ...) {
#
# }


#' @export
predict.slcm <- function(object, newdata, label, type = c("class", "posterior"), ...) {
   type <- match.arg(type)
   if (!object$fitted) stop("Latent variable model should be estimated.")

   if (missing(label))
      label <- object$model$label

   if (missing(newdata))
      posterior <- object$posterior
   else {
      data <- proc_data(newdata, object$model, FALSE)
      args <- object$args

      post <- calcPost(
         args$log_par, data$y, data$nobs, args$nvar, args$ncat,
         args$nlv, args$nroot, args$nlink, args$nleaf,
         args$nlink_unique, args$nleaf_unique,
         args$root - 1, args$tree_index - 1,
         args$u - 1, args$v - 1, args$leaf - 1,
         args$cstr_link - 1, args$cstr_leaf - 1,
         args$nclass, args$nclass_leaf,
         args$nclass_u, args$nclass_v
      )
      posterior <- output_posterior(post[label], object$model, data)
   }

   impute <- function(x) apply(x, 1, which.max)

   switch(
      type,
      class = if (length(label) == 1) impute(posterior[[label]])
      else sapply(posterior[label], impute),
      posterior = if (length(label) == 1) posterior[[label]]
      else posterior[label])
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
      object, parm, level = 0.95,
      method = c("asymp", "logit"),
      out = c("param", "logit")
) {
   if (missing(parm)) parm <- c("pi", "tau", "rho")
   out <- match.arg(out)

   logit <- unlist(object$estimates$logit[parm])
   se <- unlist(object$estimates$se$logit[parm])

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

