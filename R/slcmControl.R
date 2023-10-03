#' Control Parameters for slcm Estimation
#'
#' @param em.iterlim maximum number of iterations for EM algorithm. Default is 1000.
#' @param em.tol tolerance for the convergence of EM algorithm. Default is 1e-6.
#' @param nlm.iterlim maximum number of iterations for estimation with \code{nlm} function. Default is 500.
#' @param nlm.tol tolerance for the convergence of \code{nlm} function. Default is 1e-8.
#' @param init.param initial parameter
#' @param nrep number of trials
#' @param verbose a logical value whether to print estimation procedure.
#'
#' @return a \code{list} with control parameters for slcm estimation.
#'
#' @seealso [slcm]
#'
#' @export
slcmControl <- function(
   em.iterlim = 1000, em.tol = 1e-6,
   nlm.iterlim = 500, nlm.tol = 1e-8,
   init.param = NULL, nrep = 1,
   verbose = TRUE
) {
   ctrl <- list(
      em.iterlim = em.iterlim, em.tol = em.tol,
      nlm.iterlim = nlm.iterlim, nlm.tol = nlm.tol,
      init.param = init.param, nrep = nrep,
      verbose = verbose
   )

   class(ctrl) <- "slcmControl"
   ctrl
}
