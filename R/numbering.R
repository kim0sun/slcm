number_par <- function(param, restriction) {
   numbers <- param[-1]
   pos = 0
   r <- as.numeric(gsub("tau", "", restriction[grep("tau", restriction)]))
   for (i in seq_along(param$tau)) {
      numbers$tau[[i]][] <- pos + seq(length(numbers$tau[[i]]))
      param$tau[[i]][which(numbers$tau[[i]] %in% r)] <- -Inf
      pos = pos + length(numbers$tau[[i]])
   }

   pos = 0
   r <- as.numeric(gsub("rho", "", restriction[grep("rho", restriction)]))
   for (i in seq_along(param$rho)) {
      numbers$rho[[i]][] <- pos + seq(length(numbers$rho[[i]]))
      param$rho[[i]][which(numbers$rho[[i]] %in% r)] <- -Inf
      pos = pos + length(numbers$rho[[i]])
   }

   list(numbers = numbers, param0 = param)
}

#' Indexing \code{slcm} Parameters
#'
#' For restriction of parameters,
#'
#' @param object a \code{slcm} object
#'
#' @export
parameter_numbers <- function(object) {
   number_par(object$estimates$param, NULL)$numbers
}
