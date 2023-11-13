#' Construct structural latent class model
#'
#' Function for constructing structural latent class model.
#'
#' @param formula a formula for specifying latent structure. The details of model specification are given under 'Details'.
#' @param ... additional formulae of the model.
#' @param constraints a list containing constraints for the measurement invariance assumption. The details of posing constraints are given under 'Details'.
#'
#' @details
#' The \strong{`formula`} can be broadly categorized into three main types, each serving a unique purpose:
#'
#' 1. \strong{Defining Latent Class Variables with Manifest Indicators}: This is where you specify the relationship between a latent class variable and its manifest indicators. In these formulas, the latent class variable, denoted with square brackets or parentheses indicating the number of classes, is on the left-hand-side (lhs) and its manifest indicators are specified on right-hand-side (rhs). For example,
#' \preformatted{LC1[k] ~ x1 + x2 + x3
#' LC2[k] ~ y1 + y2 + y3
#' LC3(k) ~ z1 + z2 + z3}
#' Within these formulas, `k` denotes an integer reflecting the number of latent classes associated with the latent class variable.
#'
#' 2. \strong{Relating Latent Class Variables to Each Other}: Here, one latent class variable can be associated with or influenced by another latent class variable. The subsequent example implies that `LC2` is conditionally affected based on `LC1`.
#' \preformatted{LC1 ~ LC2}
#'
#' 3. \strong{Defining higher-level latent class variable}: Latent class variables can be identified by other latent class variables, instead of manifest indicators. Following example suggests that the `P` is measured by `LC1`, `LC2`, and `LC3` -- all of which are latent class variables.
#' \preformatted{P[k] ~ LC1 + LC2 + LC3}
#'
#' In all types of the formula, variables specified on the left-hand side (lhs) influence those on the right-hand side (rhs).
#'
#' The \strong{`constraints`} option offers a way to impose restrictions on certain probabilities in order to achieve more precise model inference. For instance, in Longitudinal Latent Class Analysis (LTA), it's imperative that latent class variables across various time-points convey identical meanings. With the `constraints` option, users can uphold measurement invariance in both the measurement and structural components of the model.
#'
#' 1. \strong{Measurement Invariance for Measurement Model}: By imposing measurement invariance in the measurement model, the software ensures that the probabilities associated with various latent class variables remain consistent. This consistency allows each latent class within these variables to retain an identical semantic meaning.
#' \preformatted{c("LC1", "LC2", "LC3")}
#' This command ensures that variables `LC1`, `LC2`, and `LC3` are semantically consistent, adhering to the measurement invariance assumption.
#'
#' 2. \strong{Measurement invariance for structural model}: In addition to the measurement model, the `constraints` option also supports the enforcement of constraints within the structural model, promoting consistent interpretations of transition probabilities.
#' \preformatted{c("P ~ LC1", "P -> LC2")}
#' This command implies that the transition probabilities from `P` to `LC1` and from `P` to `LC2` are consistent.
#'
#'
#' @returns
#' An object of class `slcm` with an element `model`. `model` includes following objects
#' \item{`tree`}{a `data.frame` describing whole parent-child relationships among latent class variables and manifest variables in the model}
#' \item{`latent`}{a `data.frame` containing descriptions for latent class variable in the model}
#' \item{`measure`}{a `data.frame` describing measurement part of the model}
#' \item{`struct`}{a `data.frame` describing structural part of the model}
#'
#' The object prints model description with four part.
#' 1. Latent variables: This delineates the latent class variables incorporated in the model, along with the number of classes for each variable. The root variable is marked by asterisk (`*`).
#' 2. Measurement model: Here, the manifest indicators for each latent class variable are presented. Any measurement constraints are indicated using lowercase alphabets. An identical alphabet signifies a consistent measurement structure.
#' 3. Structural model: This part describes the structural model by specifying the conditional dependency between latent class variables.
#' 4. Dependency constraints: This part suggests the constraints applied to the conditional dependencies between latent class variables. Each column marked with an uppercase alphabet symbolizes a consistent dependency structure.
#'
#'
#' @example examples/slcm.R
#'
#' @export
slcm = function(
   formula = NULL, ..., constraints = NULL
) {
   if (!is.list(formula)) formulae <- list(formula, ...)
   else formulae <- formula
   formulae <- formulae[sapply(formulae, inherits, "formula")]

   res <- list()
   res$model <- proc_formula(formulae, constraints)
   class(res) <- "slcm"
   res
}
