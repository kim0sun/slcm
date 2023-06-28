rep_row <- function(x, lev) {
   sb <- as.matrix(expand.grid(lev[x == 0]))
   m <- t(replicate(nrow(sb), x))
   m[, x == 0] <- sb
   m
}

split_by_name <- function(x, f) lapply(f, function(i) x[i])

stretch_data <- function(child, mf) {
   m <- split_by_name(mf, child)
   y <- lapply(m, function(x) t(sapply(x, as.numeric)))
   unlist(unname(y))
}

proc_saturated <- function(mfn, nlev, mis) {
   lev <- lapply(nlev, seq_len)
   if (!mis) {
      yobs <- mfn
      y_aggr <- aggregate(numeric(nrow(yobs)), yobs, length)
      y_unique <- y_aggr[, -ncol(y_aggr)]
      freq <- y_aggr[, ncol(y_aggr)]
      loglik <- c(freq %*% (log(freq) - log(sum(freq))))
   } else {
      na_ind <- rowSums(mfn == 0) > 0
      yobs <- mfn[!na_ind, , drop = FALSE]
      ymis <- mfn[ na_ind, , drop = FALSE]

      yobs0 <- aggregate(numeric(nrow(yobs)), yobs, length)
      ymis0 <- aggregate(numeric(nrow(ymis)), ymis, length)
      y_expand <- do.call(rbind, apply(
         ymis0[, -ncol(ymis0)], 1, rep_row, lev, simplify = FALSE))
      y0 <- rbind(yobs0, cbind(y_expand, x = 0))

      y_aggr <- aggregate(y0[[ncol(y0)]], y0[-ncol(y0)], sum)
      y_unique <- y_aggr[, -ncol(y0)]
      freq <- y_aggr[, ncol(y0)]

      nrep <- apply(ymis0[, -ncol(ymis0)], 1,
                    function(x) prod(nlev[x == 0]))
      miss <- match(apply(y_expand, 1, paste, collapse = ""),
                    apply(y_unique, 1, paste, collapse = "")) - 1

      calc_mis <- calcfreq(
         miss, nrep, nrow(ymis0), ymis0[, ncol(ymis0)],
         freq, nrow(y_aggr), nrow(mfn), 1e-6, 500
      )

      freq <- calc_mis$freq
      loglik <- sum(calc_mis$freq * log(calc_mis$theta), na.rm = TRUE) + calc_mis$loglik
   }

   list(y_unique = y_unique,
        freq = freq,
        loglik = loglik)
}


proc_data <- function(data, model, saturate = TRUE) {
   child <- model$measure$indicator
   f <- paste("~", paste(unlist(child), collapse = "+"))
   mf <- mf_raw <- model.frame(formula(f), data, na.action = NULL)
   mf[] <- lapply(mf, function(x)
      if (is.factor(x)) x else factor(x))
   lev <- lapply(mf, levels)
   nlev <- sapply(mf, nlevels)
   if (any(nlev < 2))
      stop("wrong indicator variable(s) (< 2 level):\n",
           paste0("`", unlist(child)[nlev < 2], "`", collapse = " "))

   attr(mf, "levels") <- lev

   mfn <- mf
   mfn[] <- lapply(mf, as.numeric)
   mfn[is.na(mfn)] <- 0
   y <- lapply(child, function(x) t(mfn[x]))
   attr(mf, "y") <- unlist(y)

   misInd <- which(is.na(mf), arr.ind = TRUE)
   sat <- proc_saturated(mfn, nlev, any(is.na(mf)))
   attr(mf, "misInd") <- misInd
   attr(mf, "y_unique") <- sat$y_unique
   attr(mf, "freq") <- sat$freq
   attr(mf, "df") <- length(sat$freq) - 1
   attr(mf, "loglik") <- sat$loglik

   mf
}

