#' @export
estimates <- function(object, ...) UseMethod("estimates")

#' @export
estimates.slcm <- function(
   object, which = c("pi", "tau", "rho"),
   type = c("param", "logit"),
   digits = max(3L, getOption("digits") - 3L),
   round = NULL, ...
) {
   type <- match.arg(type)
   which <- match.arg(which, several.ok = TRUE)
   if (!object$fitted) return(NULL)
   if (type == "param")
      est <- object$estimates$param[which]
   if (type == "logit")
      est <- object$estimates$logit[which]

   if (is.numeric(round))
      est <- rapply(est, base::round, how = "list", digits = round)
   print(est, digits = digits)
   invisible(est)
}

#' @export
std.err <- function(object, ...) UseMethod("std.err")

#' @export
std.err.slcm <- function(
   object, which = c("pi", "tau", "rho"),
   type = c("param", "logit"),
   digits = max(3L, getOption("digits") - 3L), ...
) {
   type <- match.arg(type)
   which <- match.arg(which, several.ok = TRUE)
   if (!object$fitted) return(NULL)

   if (type == "param")
      se <- object$estimates$se$param[which]
   if (type == "logit")
      se <- object$estimates$se$logit[which]

   print(se, digits = digits)

   invisible(se)
}
