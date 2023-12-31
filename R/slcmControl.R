#' Control Parameters for slcm Estimation
#'
#' @param em.iterlim an integer specifying maximum number of iterations allowed for EM algorithm. Default is 1000.
#' @param em.tol a numeric value setting tolerance for the convergence of EM algorithm. Default is 1e-6.
#' @param nlm.iterlim an integer specifying maximum number of iterations allowed for estimation with \code{nlm} function. Default is 500.
#' @param nlm.tol a numeric value setting tolerance for the convergence of \code{nlm} function. Default is 1e-8.
#' @param init.param initial parameters.
#' @param nrep number of trials.
#' @param test.iter an integer specifying maximum number of iterations allowed for testing parameters.
#' @param na.rm a logical value indicating whether to remove missing values (NA)
#' @param verbose a logical value indicating whether to display progress updates during the estimation process
#'
#' @return a \code{list} with control parameters for slcm estimation.
#'
#' @seealso [slcm]
#'
#' @export
slcmControl <- function(
   em.iterlim = 5000, em.tol = 1e-6,
   nlm.iterlim = 1000, nlm.tol = 1e-8,
   init.param = NULL,
   nrep = 1, test.iter = 100,
   na.rm = FALSE, verbose = TRUE
) {
   ctrl <- list(
      em.iterlim = em.iterlim, em.tol = em.tol,
      nlm.iterlim = nlm.iterlim, nlm.tol = nlm.tol,
      init.param = init.param,
      nrep = nrep, test.iter = test.iter,
      na.rm = na.rm, verbose = verbose
   )

   class(ctrl) <- "slcmControl"
   ctrl
}
