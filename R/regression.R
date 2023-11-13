#' Regress Exogenous Variables on Latent Variables
#'
#' This function allows you to perform regression analysis to understand the impact of exogenous (external) variables on the latent class variables in the estimated slcm model. Logistic regression and three-step approach is used for this purpose.
#'
#' @param object an object of class `slcm` and `estimated`
#' @param formula a formula defining the regression model. It should include both latent class variables from the estimated model and any exogenous (external) variables
#' @param data an optional data frame containing the exogenous variables in interest.
#' @param imputation imputation method for latent class variables.
#' @param method bias-adjusting method for three-step approach
#'
#'
#' @export
regress <- function(object, ...) UseMethod("regress")

#' @exportS3Method slcm::regress slcm
regress.slcm <- function(
      object, formula, data = parent.frame(),
      imputation = c("modal", "prob"),
      method = c("naive", "BCH", "ML"), ...
) {
   if (!inherits(object, "estimated"))
      stop("Latent variable model should be estimated.")

   # Import
   labels <- all.vars(formula)
   latent <- labels[labels %in% row.names(object$model$latent)]
   imputation <- match.arg(imputation)
   method <- match.arg(method)

   # Imputation
   impute <- function(x, imputation) {
      if (imputation == "modal")
         imputed <- as.factor(apply(x, 1, which.max))
      else
         imputed <- as.factor(apply(x, 1, function(y)
            sample(seq_len(ncol(x)), 1, prob = y)))
      return(imputed)
   }

   if (missing(data)) data <- object$mf
   else data <- data.frame(object$mf, data)
   imputed <- lapply(object$posterior$marginal[latent],
                     impute, imputation)
   data <- data.frame(data, imputed)
   mf <- model.frame(formula, data)

   # Functions
   dn <- function(X, b, y, k)
      -sum(dnorm(X, b[y], b[k], log = TRUE))

   d <- function(X, b, y, k)
      -sum(dnorm(X, b[y], b[k], log = TRUE))

   cprobs <- function(X, b, ref) {
      beta <- matrix(nrow = nrow(b), ncol = ncol(b) + 1)
      beta[, ref] = 0
      beta[, -ref] = b
      xb <- X %*% beta
      exb <- exp(xb)
      denom <- rowSums(exb)
      xb - log(denom)
   }

   y <- model.response(mf)
   X <- model.matrix(formula, mf)
   nr <- nlevels(y) - 1
   nc <- ncol(X)
   init <- numeric(nc * nr)

   if (method == "naive") {
      # naive (biased)
      naive_ll <- function(par, X, y, ref) {
         b <- matrix(par, nrow = ncol(X))
         prob <- cprobs(X, b, ref)
         - sum(prob[cbind(1:nrow(prob), y)])
      }
      fit <- nlm(naive_ll, init, X = X, y = y,
                 ref = nlevels(y), hessian = TRUE)
   } else {
      # bias_adjusted
      p <- object$posterior$marginal[latent][rownames(mf), ]
      w <- switch(
         imputation,
         modal = lapply(p, apply, 1, function(x) x == max(x)),
         prob  = p
      )
      d <- lapply(p, function(pp) (pp %*% w) / rowSums(pp))

      if (method == "BCH") {
         # BCH
         w_ <- Map(function(x, y) x %*% ginv(y), w, d)

         bch_ll <- function(par, X, w_, ref) {
            b <- matrix(par, ncol(X))
            prob <- cprobs(X, b, ref)
            - sum(w_ * prob)
         }
         fit <- nlm(bch_ll, init, X = X, w_ = w_,
                    ref = nlevels(y), hessian = TRUE)
      } else if (method == "ML") {
         # ML
         w_ <- t(sapply(y, function(x) d[,x]))

         ml_ll <- function(par, X, w_, ref) {
            b <- matrix(par, ncol(X))
            prob <- probs(X, b, ref)
            ll <- rowSums(exp(prob + log(w_)))
            - sum(log(ll))
         }
         fit <-  nlm(ml_ll, init, X = X, w_ = w_,
                     ref = nlevels(y), hessian = TRUE)
      }
   }

   rn <- paste0(seq_len(nr), "/", nr + 1)
   cn <- colnames(X)
   coef <- matrix(
      fit$estimate, nr, nc,
      dimnames = list(class = rn, cn)
   )

   dn <- paste0(rep(cn, nr), "|", rep(rn, each = nc))
   vcov <- matrix(
      ginv(fit$hessian), nc * nr, nc * nr,
      dimnames = list(dn, dn)
   )

   se <- matrix(
      diag(vcov), nr, nc, byrow = TRUE,
      dimnames = list(class = rn, cn)
   )

   res <- list()
   res$coefficients <- coef
   res$std.err <- se
   res$vcov <- vcov
   res$dim <- c(nr, nc)
   class(res) <- "reg.slcm"
   res$ll <- -fit$minimum

   return(res)
}


#' @exportS3Method base::print reg.slcm
print.reg.slcm <- function(
      x, digits = 3, wald = TRUE, pval = TRUE
) {
   cat("Coefficients:")
   print.default(format(x$coefficients, digits = digits),
                 print.gap = 2L, quote = FALSE)
   invisible(x)
}

#' @exportS3Method base::summary reg.slcm
summary.reg.slcm <- function(
   x, digits = 3, odds.ratio = FALSE, wald = TRUE
) {
   cat("Coefficients:")
   print.default(format(x$coefficients, digits = digits), print.gap = 2L,
                 quote = FALSE)
   cat("\nStd. Errors:")
   print.default(format(x$std.err, digits = digits), print.gap = 2L,
                 quote = FALSE)
   if (odds.ratio) {
      cat("Odds Ratio:")
      print.default(format(exp(x$coefficients), digits = digits),
                    print.gap = 2L, quote = FALSE)
   }
   if (wald) {
      wald <- x$coefficients / x$std.err
      pval <- pnorm(abs(wald), 1, lower.tail = FALSE)
      cat("\nValue/SE (Wald statistics):")
      print.default(format(wald, digits = digits),
                    print.gap = 2L, quote = FALSE)
      cat("\nPr(>|W|):")
      print.default(format(pval, digits = digits),
                    print.gap = 2L, quote = FALSE)
   }
   invisible(x)
}


#' @exportS3Method stats::confint reg.slcm
confint.reg.slcm <- function(
   object, parm, odds.ratio = FALSE, level = 0.95, ...
) {
   fci <- function(cf, se) {
      a <- (1 - level)/2
      a <- c(a, 1 - a)
      pct <- format.pc(a, 3)
      fac <- qnorm(a)
      ci <- array(NA, dim = c(length(parm), 2L),
                  dimnames = list(parm, pct))
      ci[] <- cf[parm] + se %o% fac
      ci
   }
   cf <- object$coefficients
   se <- object$std.err
   cn <- colnames(cf)
   rn <- rownames(cf)
   if (missing(parm))
      parm <- cn
   else if (is.numeric(parm))
      parm <- cn[parm]
   ci <- lapply(seq_len(nrow(cf)), function(i)
      fci(cf[i, ], se[i,]))
   names(ci) <- rn

   for (i in seq_len(nrow(cf))) {
      cat(rn[i], ":\n")
      print.default(ci[[i]])
   }
   invisible(ci)
}
