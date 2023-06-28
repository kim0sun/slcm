norm1 <- function(x) log(x) - log(sum(x))

norm2 <- function(x) x - log_sum_exp(x)

ginv <- function(x, tol = sqrt(.Machine$double.eps)) {
   xsvd <- svd(x)
   pos <- xsvd$d > tol * xsvd$d[1]

   if (all(pos))
      inv <- xsvd$v %*% (t(xsvd$u) / xsvd$d)
   else
      inv <- xsvd$v[, pos, drop = FALSE] %*%
      (t(xsvd$u[, pos, drop = FALSE]) / xsvd$d[pos])
   structure(inv, dimnames = dimnames(x)[2:1])
}

get_se <- function(vcov) {
   var <- diag(vcov)
   var[var < 0] = 0
   sqrt(var)
}

bdiag <- function(x) {
   nr <- sum(sapply(x, nrow))
   nc <- sum(sapply(x, ncol))
   ans <- matrix(0, nr, nc)
   ibegin <- 1; iend <- 0
   jbegin <- 1; jend <- 0

   for (m in x){
      iend <- iend + nrow(m)
      jend <- jend + ncol(m)
      ans[ibegin:iend, jbegin:jend] <- m
      ibegin <- ibegin + nrow(m)
      jbegin <- jbegin + ncol(m)
   }

   ans
}

jmat <- function(x, simplify = TRUE) {
   ex <- exp(x)
   a <- diag(ex)
   b <- ex %o% ex
   a - b
}

