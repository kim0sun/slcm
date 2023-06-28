#' @export
plot.slcm <- function(x, abbreviation = FALSE, dir = "TD",
                        equal_rank = NULL, font = "Helvetica", ...) {
   latent <- x$model$latent
   tree <- x$model$tree

   latent$child <- NA
   latent$child[latent$leaf] <- if (abbreviation)
      sapply(latent$children[latent$leaf], function(x)
         paste0("'", x[1], " ~ ", x[length(x)], "'")) else
            sapply(latent$children[latent$leaf], function(x)
               paste0("'", x, "'", collapse = ", "))
   latent$child[!latent$leaf] <-
      sapply(latent$children[!latent$leaf], function(x)
         paste0("'", x, "'", collapse = ", "))

   node <- paste0(
      "node [shape = box]\n",
      paste(tree$child[tree$leaf], collapse = ", "),
      "\n\n node [shape = oval]\n",
      paste(rownames(latent), collapse = ", ")
   )
   path <- paste(paste(rownames(latent), " -> {", latent$child, "}"), collapse = "\n")
   if (!missing(equal_rank)) {
      equal_rank <- equal_rank[equal_rank %in% row.names(latent)]
      ranks <- paste0("{rank = same; '", paste0(equal_rank, collapse = "';'"), "';}")
   }

   text <- paste0(
      "digraph { \n  rankdir = '", dir, "';",
      "node[fontname = '", font, "']\n\n",
      node, "\n\n", path, "\n\n",
      if (!missing(equal_rank)) ranks, "\n}"
   )

   DiagrammeR::grViz(text, ...)
}


