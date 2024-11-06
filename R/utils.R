#' Resample function
#'
#' @param x Vector
#' @param ... Additional arguments passed to \code{sample}
#'
#' @export

resample <- function(x, ...) x[sample.int(length(x), ...)]

#' Yield partial derivatives
#'
#' @param r Numeric. Intrinsic growth rate of species i.
#' @param a Numeric vector. Interaction coefficients.
#' @param i Integer. Species numeric ID.
#' @inheritParams stability
#'
#' @export

fun_partial <- function(r, a, i, x0, model) {

  # check input -------------------------------------------------------------
  if (length(r) != 1)
    stop(paste("Input 'r' must be a scalar:",
               "'r' has length", length(r)))

  if (length(i) != 1)
    stop(paste("'i' must be a scalar:",
               "'i' has length", length(i)))

  if (length(a) != length(x0))
    stop(paste("Invalid inputs in 'a' or 'x0':",
               "'a' has length =", length(a),
               "while 'x0' has length =", length(x0)))

  if (!(model %in% c("ricker", "bh", "glv")))
    stop(paste("Invalid model type: 'model' must be either of",
               "'ricker', 'bh', or 'glv'"))

  # model formula -----------------------------------------------------------
  ## declare vectorized parameters and variables
  v_a <- paste0("a[", seq_len(length(a)), "]")
  v_x <- paste0("x[", seq_len(length(x0)), "]")
  arg <- paste(c("r", "x", "a"), collapse = ", ")

  ## linear combination
  lcm <- paste(v_a, "*", v_x)

  ## function
  f <- NULL

  ## Generalized Lotka-Volterra model
  if (model == "glv") {
    ## get a model formula
    fm <- c("r", lcm)
    m <- paste0("x[", i, "]", " * ",
                "(",
                paste(fm, collapse = " + "),
                ")")

    ## function text for evaluation
    f_text <- parse(text = paste0("f <- function(", arg, ") {",
                                  m,
                                  "}"))

    eval(f_text)
  }

  ## Ricker model
  if (model == "ricker") {
    ## get a model formula
    fm <- c("r", lcm)
    m <- paste0("x[", i, "]", " * ",
                "exp(",
                paste(fm, collapse = " + "),
                ")")

    ## function text for evaluation
    f_text <- parse(text = paste0("f <- function(", arg, ") {",
                                  m,
                                  "}"))

    eval(f_text)
  }

  ## Beverton-Holt model
  if (model == "bh") {
    ## get a model formula
    m <- paste0("x[", i, "]", " * ", "exp(r)",
                " * ",
                "(1 + ",
                paste(lcm, collapse = " + "),
                ") ** -1")

    ## function text for evaluation
    f_text <- parse(text = paste0("f <- function(", arg, ") {",
                                  m,
                                  "}"))

    eval(f_text)
  }

  ## return partial derivatives evaluated at x0 (equilibrium)
  return(pracma::jacobian(f, x0 = x0, r = r, a = a))
}


#' Extra prey function
#'
#' @inheritParams ppm
#' @param alpha0 Numeric matrix. \code{n_species x n_species} interaction matrix
#' @param j Integer. Consumer's index
#' @param i0 Integer. Index for the first prey
#' @param tp Numeric. Initial trophic positions
#' @param kappa Integer. Number of extra prey items
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export

extra_prey <- function(alpha0,
                       j,
                       i0,
                       tp,
                       theta,
                       kappa,
                       cannibal = FALSE) {
  # alpha0: adjacency matrix defining trophic interactions
  # j:      consumer's index (j - 1 is the number of possible prey)
  # i0:     index of the first prey chosen (1 <= i0 < j)
  # tp:     vector of trophic positions for possible prey nodes (length = S)
  # kappa:  number of extra prey nodes in addition to the first prey i0
  # theta:  scale parameter for an exponential decay of prey preference

  # verify inputs
  if (length(tp) != ncol(alpha0)) stop(paste("'tp' length seems incorrect;",
                                             "the length must be",
                                             ncol(alpha0)))

  if (!(i0 < j && 1 <= i0)) stop("i0 must be a non-zero intger smaller than j")

  if (kappa >= j) stop("kappa must be an integer smaller than j")

  # probability of consumer j picks prey i
  lambda <- 1 / theta
  tp[j] <- tp[i0] + 1
  p_ij <- exp(- lambda * abs(tp[i0] - tp))

  # pick kappa prey species out of j - 1 nodes (without cannibalism) or j nodes
  if (cannibal) {
    i_index <- 1:j
    p_ij <- p_ij[1:j]
  } else {
    i_index <- 1:(j - 1)
    p_ij <- p_ij[1:(j - 1)]
  }

  # exclude i0 from resample() because i0 is the first pick
  # pick 'kappa' samples from index[-i0]
  # weighted by p_ij
  if (length(i_index) > 1) {
    i_pick <- resample(i_index[-i0],
                       size = kappa,
                       prob = p_ij[-i0])

    i <- c(i0, i_pick)
  } else {
    # when j = 2 with no cannibalism
    # no choice but i0 available as prey
    i <- i0
  }

  alpha0[i, j] <- 1

  return(alpha0)
}

#' Apply conversion efficiency and attack rate
#'
#' @inheritParams extra_prey
#' @param competition List for competition coefficients between producers.
#'  Specify minimum and maximum values for a uniform distribution.
#' @param attack List for attack rates.
#'  Specify minimum and maximum values for a uniform distribution.
#' @param convert List for conversion efficiency.
#'  Specify minimum and maximum values for a uniform distribution.
#' @param regulation List for self regulation (or intraspecific competition).
#'  Specify minimum and maximum values for a uniform distribution.
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export

to_alpha <- function(alpha0,
                     competition = list(min = 0,
                                        max = 0),
                     attack = list(min = 0,
                                   max = 1),
                     convert = list(min = 0,
                                    max = 1),
                     regulation = list(min = 0,
                                       max = 1)) {

  # identify basal species
  u_alpha0 <- abs(alpha0)
  u_alpha0[lower.tri(u_alpha0)] <- 0
  basal <- which(colSums(u_alpha0) == 0)

  # scale by # of resources
  su_alpha0 <- t(t(u_alpha0[, -basal]) / colSums(u_alpha0)[-basal])

  # generate random parameters
  ## matrix for attack rates
  a <- with(attack,
            matrix(stats::runif(ncol(alpha0)^2,
                                min = min,
                                max = max),
                   nrow = nrow(alpha0),
                   ncol = ncol(alpha0)))

  ## matrix for conversion efficiency
  b <- with(convert,
            matrix(stats::runif(ncol(alpha0)^2,
                                min = min,
                                max = max),
                   nrow = nrow(alpha0),
                   ncol = ncol(alpha0)))

  ## vector for intraspecific competition
  sr <- with(regulation,
             stats::runif(ncol(alpha0),
                          min = min,
                          max = max))

  # interaction matrix
  ## diag(alpha) represent net effects of cannibalism
  u_alpha0[, -basal] <- a[, -basal] * su_alpha0
  alpha <- -u_alpha0 + b * t(u_alpha0)

  ## competition between producers
  alpha[basal, basal] <- with(competition,
                              -stats::runif(n = length(basal)^2,
                                            min = min,
                                            max = max))

  if (length(basal) > 1) {
    diag(alpha[basal, basal]) <- 0
  } else {
    alpha[basal, basal] <- 0
  }

  ## intraspecific interaction
  ## added to net effects of cannibalism as "diag(alpha) - sr"
  diag(alpha) <- diag(alpha) - sr

  return(alpha)
}

