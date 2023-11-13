rep_row <- function(x, lev) {
   sb <- as.matrix(expand.grid(lev[is.na(x)]))
   m <- t(replicate(nrow(sb), x))
   m[, is.na(x)] <- sb
   m
}

match_row <- function(x, y) {
   nna <- !is.na(x)
   which(!colSums(t(y[,nna, drop = FALSE]) != x[nna]))
}

count_row <- function(mf, na.rm) {
   y_sort <- mf[do.call(order, mf), ]
   if (!sum(is.na(mf))) na.rm = TRUE
   if (na.rm) {
      uniq <- !duplicated(y_sort)
      y_uniq <- y_sort[uniq,]
      freq <- tabulate(cumsum(uniq))
      theta <- freq / sum(freq)
      loglik <- sum(freq * log(theta))
      row.names(y_uniq) <- NULL
      return(list(y_unique = y_uniq, freq = freq,
                  theta = theta, loglik = loglik))
   } else {
      miss <- rowSums(is.na(y_sort)) >  0
      yobs <- y_sort[!miss, ]
      uobs <- !duplicated(yobs)
      obs_uniq <- yobs[uobs, ]
      obs_freq <- tabulate(cumsum(uobs))
      ncell <- length(obs_freq)
      row.names(obs_uniq) <- NULL

      ymis <- y_sort[miss, ]
      umis <- !duplicated(ymis)
      mis_uniq <- ymis[umis, ]
      misobs <- apply(mis_uniq, 1, match_row, obs_uniq)
      rmis <- sapply(misobs, length)
      nmis <- length(misobs)
      imis <- unlist(misobs) - 1
      fmis <- tabulate(cumsum(umis))

      calc_mis <- calcfreq2(
         imis, rmis, nmis, fmis, obs_freq, ncell, nrow(mf), 1e-6, 500
      )
      return(c(list(y_unique = obs_uniq), calc_mis))
   }
}

proc_data <- function(data, model, na.rm) {
   child <- model$measure$indicator
   f <- paste("~", paste(unlist(child), collapse = "+"))
   if (na.rm)
      mf <- model.frame(formula(f), data)
   else
      mf <- model.frame(formula(f), data, na.action = NULL)

   mf[] <- lapply(mf, function(x)
      if (is.factor(x)) x else factor(x))
   lev <- lapply(mf, levels)
   nlev <- sapply(mf, nlevels)
   if (any(nlev < 2))
      stop("wrong indicator variable(s) (< 2 level):\n",
           paste0("`", unlist(child)[nlev < 2], "`", collapse = " "))

   res <- count_row(mf, na.rm)

   nmf <- mf
   nmf[] <- lapply(mf, as.numeric)
   nmf[is.na(nmf)] <- 0
   attr(mf, "y") <- unlist(lapply(child, function(x) t(nmf[x])),
                           use.names = FALSE)
   attr(mf, "levels") <- lev
   attr(mf, "y_unique") <- res$y_unique
   nyu <- res$y_unique
   nyu[] <- lapply(nyu, as.numeric)
   attr(mf, "yu") <- unlist(lapply(child, function(x)
      t(nyu[x])), use.names = FALSE)
   attr(mf, "freq") <- res$freq
   attr(mf, "theta") <- res$theta
   attr(mf, "loglik") <- res$loglik
   attr(mf, "df") <- prod(nlev) - 1

   mf
}

