#' @exportS3Method base::plot slcm
plot.slcm <- function(x, abbreviation = FALSE, dir = "TD",
                        equal_rank = NULL, font = "Helvetica", ...) {
   latent <- x$model$latent
   tree <- x$model$tree

   latent$child <- NA
   if (abbreviation)
      latent$children[latent$leaf] <-
         sapply(latent$children[latent$leaf], function(x)
            paste0("'", x[1], " ~ ", x[length(x)], "'"))

   parent <- row.names(latent)
   child <- sapply(latent$children, paste, collapse = ", ")

   node <- paste0(
      "node [shape = box]\n",
      paste(setdiff(unlist(latent$children), parent), collapse = ", "),
      "\n\n node [shape = oval]\n",
      paste(names(child), collapse = ", ")
   )
   path <- paste(paste(names(child), " -> {", child, "}"), collapse = "\n")
   if (!missing(equal_rank)) {
      equal_rank <- equal_rank[equal_rank %in% row.names(latent)]
      ranks <- paste0(
         "{rank = same; '", paste0(equal_rank, collapse = "';'"), "';}"
      )
   }

   text <- paste0(
      "digraph { \n  rankdir = '", dir, "';",
      "node[fontname = '", font, "']\n\n",
      node, "\n\n", path, "\n\n",
      if (!missing(equal_rank)) ranks, "\n}"
   )

   DiagrammeR::grViz(text, ...)
}

