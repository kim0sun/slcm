get_lformula <- function(f) {
   lhs <- deparse(f[[2]])
   lhs_split <- strsplit(lhs, "[()]|\\[|\\]")[[1]]

   attr(f, "label") <- lhs_split[1]
   if (exists(lhs_split[2])) nc <- get(lhs_split[2])
   else nc <- lhs_split[2]
   attr(f, "nclass") <- as.numeric(nc)
   attr(f, "vars") <- labels(terms(f))

   f
}

proc_tree <- function(nc, iv) {
   tree <- data.frame(
      parent = iv[,1], "->", child = iv[,2],
      fix.empty.names = FALSE
   )
   for (i in seq_len(nrow(tree))) {
      par <- tree$parent[i]
      rank <- 1
      while (par %in% tree$child) {
         rank <- rank + 1
         par <- tree$parent[tree$child == par]
      }
      tree$root[i] <- par
      tree$rank[i] <- rank
      tree$leaf[i] <- !(tree$child[i] %in% tree$parent)
   }

   label <- unique(tree$parent)
   child <- unique(tree$child)
   leaf <- tapply(tree$leaf, tree$parent, any)
   rank <- tapply(tree$rank, tree$parent, max)
   root <- sapply(label, "%in%", tree$root)
   nvar <- tapply(tree$child, tree$parent, length)
   indv <- tapply(tree$child, tree$parent, c)
   latent <- data.frame(
      nclass = as.numeric(nc[label, drop = FALSE]),
      leaf = leaf[label, drop = FALSE],
      rank = rank[label, drop = FALSE],
      root = root[label, drop = FALSE],
      nvar = nvar[label, drop = FALSE],
      children = I(indv[label, drop = FALSE]),
      row.names = label
   )
   leaf <- label[latent$leaf]

   obs <- lapply(indv, function(x)
      x[!(x %in% label)])
   measure <- data.frame(
      nclass = as.numeric(nc[leaf, drop = FALSE]),
      nvar = sapply(obs, length)[leaf, drop = FALSE],
      indicator = I(obs[leaf, drop = FALSE]),
      row.names = leaf
   )

   struct <- tree[!tree$leaf, , drop = FALSE]
   struct <- struct[order(struct$rank, factor(struct$child, label)), 1:5]
   rownames(struct) <- NULL
   struct$parent <- factor(struct$parent, label)
   struct$child <- factor(struct$child, label)

   list(tree = tree,
        latent = latent,
        measure = measure,
        struct = struct)
}


identify_constr <- function(constr, model) {
   lt <- model$latent
   mr <- model$measure
   st <- model$struct
   label <- row.names(lt)
   leaf <- row.names(mr)

   const_leaf <- letters[seq_len(nrow(mr))]
   names(const_leaf) <- leaf
   const_edge <- LETTERS[seq_len(nrow(st))]

   if (!is.null(constr) && !is.list(constr))
      constr <- list(constr)


   for (i in seq_along(constr)) {
      constr[[i]] <- gsub(" ", "", constr[[i]])
      const <- deparse(constr[[i]])
      con <- strsplit(constr[[i]], "->|~")
      type <- unique(lengths(con))
      if (length(type) != 1) {
         warning("Wrong constraint: ", const, call. = FALSE)
         next
      }
      if (type == 1) {
         islf <- sapply(con, "%in%", leaf)
         if (sum(islf) > 1) {
            con <- con[islf]
            dims <- lapply(con, function(x)
               unlist(mr[as.character(x), 1:2]))
            if (length(unique(dims)) == 1) {
                const_leaf[unlist(con)] <- i
            } else {
               warning("Wrong constraint: ", const, call. = FALSE)
               next
            }
         } else {
            warning("Wrong constraint: ", const, call. = FALSE)
            next
         }
      } else if (type == 2) {
         nclass <- lapply(con, function(x) lt[x, "nclass"])
         if (length(unique(nclass)) == 1) {
            pl <- apply(st[, c("parent", "child")], 1, paste0, collapse = "~")
            pc <- sapply(con, paste0, collapse = "~")
            const_edge[match(pc, pl)] <- i
         } else {
            warning(i, "th constraint is wrong.")
            next
         }
      } else {
         warning("Wrong constraint: ", const, call. = FALSE)
         next
      }
   }

   reval <- function(x, c) c[factor(x, unique(x))]
   mr$constraint <- reval(const_leaf, letters)
   st$constraint <- reval(const_edge, LETTERS)

   list(tree = model$tree,
        latent = lt,
        measure = mr,
        struct = st,
        constraints = constr)
}

proc_formula <- function(formula, constr) {
   formula <- sapply(formula, get_lformula)
   label  <- sapply(formula, attr, "label")
   nclass <- lapply(formula, attr, "nclass")
   indvar <- lapply(formula, attr, "vars")

   nc <- cbind(rep(label, lengths(nclass)), unlist(nclass))
   iv <- cbind(rep(label, lengths(indvar)), unlist(indvar))

   nclass <- tapply(nc[,2], nc[,1], function(x) unique(x[!is.na(x)]))
   indvar <- tapply(iv[,2], iv[,1], function(x) unique(x[!is.na(x)]))
   err1 <- sapply(nclass, length)
   if (!all(err1 == 1)) stop("wrong number of classes for ", paste(
      paste0("`", names(err1)[err1 != 1], "`"), collapse = ", "))
   err2 <- table(unlist(indvar))
   if (!all(err2 == 1)) stop("wrong number of parents for ", paste(
      paste0("`", names(err2)[err2 != 1], "`"), collapse = ", "))

   model <- proc_tree(nclass, iv)
   identify_constr(constr, model)
}
