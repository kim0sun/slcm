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
