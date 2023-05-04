#' @export
slcm.control <- function(
   em.iterlim = 1000, em.tol = 1e-6,
   nlm.iterlim = 500, nlm.tol = 1e-8,
   var.iterlim = 1000, var.tol = 1e-3, alpha = NULL,
   verbose = TRUE, new.iter = 1000,
   init.param = NULL, nrep = 1, contraints = NULL, ...
) {
   ctrl <- list(
      em.iterlim = em.iterlim, em.tol = em.tol,
      nlm.iterlim = nlm.iterlim, nlm.tol = nlm.tol,
      var.iterlim = var.iterlim, var.tol = var.tol,
      alpha = alpha,
      verbose = verbose, new.iter = new.iter,
      init.param = init.param, nrep = nrep,
      contraints = contraints
   )
   if (!is.null(init.param)) ctrl$nrep = 1

   class(ctrl) <- "slcm.control"
   ctrl
}
