#' Resample function
#'
#' @param x Vector
#' @param ... Additional arguments passed to \code{sample}
#'
#' @export

resample <- function(x, ...) x[sample.int(length(x), ...)]

#' Compute partial derivatives of a single-species growth function
#'
#' This function evaluates the Jacobian (partial derivatives) of a species-specific
#' growth function at equilibrium densities, for GLV, Ricker, or Beverton-Holt models.
#'
#' @param r Numeric. Intrinsic growth rate of species i.
#' @param a Numeric vector. Interaction coefficients of species i with all species.
#' @param i Integer. Index of the focal species.
#' @inheritParams stability
#'
#' @return Numeric matrix. Jacobian of the focal species' growth function w.r.t all species.
#'
#' @export

fn_partial <- function(r, a, i, x0, model) {

  ## Ricker model
  if (model == "ricker") {
  # check inputs
  if (!is.numeric(r) || length(r) != 1) stop("r must be a scalar numeric")
  if (!is.numeric(i) || length(i) != 1) stop("i must be a scalar numeric")
  if (!is.numeric(a) || length(a) != length(x0)) stop("Length of 'a' and 'x0' must match")

  # define model function
  f <- switch(model,
              glv    = function(x) x[i] * (r + sum(a * x)),
              ricker = function(x) x[i] * exp(r + sum(a * x)),
              bh     = function(x) x[i] * exp(r) / (1 + sum(a * x)),
              stop("Unknown model"))

  # return Jacobian evaluated at x0
  pracma::jacobian(f, x0)
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

#' Equilibrium solver generator
#'
#' Returns a function that computes equilibrium abundances for a given
#' intrinsic growth vector \code{r} and interaction matrix \code{alpha},
#' using one of three population models: \code{"glv"}, \code{"ricker"}, or \code{"bh"}.
#'
#' @param model Character string specifying the model type.
#'   One of \code{"glv"}, \code{"ricker"}, or \code{"bh"}.
#'
#' @return A function with arguments \code{r} and \code{alpha} that returns
#'   the equilibrium vector \code{x0}.
#'
#' @export

fn_x0 <- function(model) {

  ## transformation of r differs by model
  r_transform <- switch(model,
                        # x0 = solve(alpha) %*% (-r)
                        glv    = function(r) -r,
                        ricker = function(r) -r,

                        # x0 = solve(alpha) %*% (exp(r) - 1)
                        bh     = function(r) exp(r) - 1
  )

  ## Main equilibrium solver shared by all models
  compute_x0 <- function(r, alpha) {

    if (det(alpha) == 0) return(NA)

    # base equilibrium
    rhs <- r_transform(r)
    x0  <- drop(solve(alpha) %*% rhs)

    # competitive exclusion adjustment
    if (any(x0 < 0)) {

      s_plus <- which(x0 > 0)
      s_neg  <- which(x0 <= 0)

      alpha_sub <- alpha[s_plus, s_plus]
      r_sub     <- r[s_plus]

      rhs_sub <- r_transform(r_sub)

      x0[s_plus] <- drop(solve(alpha_sub) %*% rhs_sub)
      x0[s_neg]  <- 0
    }

    x0
  }

  return(compute_x0)
}

#' Compute intrinsic growth rates for a given model
#'
#' @param model Model type: "glv", "ricker", or "bh".
#'
#' @return A function that computes r given alpha and x0.
#'
#' @export

fn_r <- function(model) {

  compute_r <- switch(model,
                      glv = function(alpha, x0) {
                        drop(-alpha %*% x0)
                      },
                      ricker = function(alpha, x0) {
                        drop(-alpha %*% x0)
                      },
                      bh = function(alpha, x0) {
                        log(1 + drop(alpha %*% x0))
                      }
  )

  return(compute_r)
}
