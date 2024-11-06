#' Find suitable intrinsic growth rates
#'
#' @param alpha Numeric matrix.
#'  `n_species x n_species` interaction matrix
#' @param k0 Numeric.
#'  Total carrying capacity for basal species combined
#' @param theta Numeric.
#'  Shape parameter of the Dirichlet distribution,
#'  which defines the relative equilibrium densities of basal species.
#' @param lambda0 Numeric.
#'  Initial value of a rate parameter for the exponential decay
#'  of equilibrium densities with trophic position.
#' @param interval Numeric.
#'  Increment of lambda value.
#' @param sigma Numeric.
#'  Degree of noise added to equilibrium densities
#' @param maxit Integer.
#'  Maximum number of iterations to find a suitable lambda
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export

findr <- function(alpha,
                  k0,
                  theta = 1,
                  lambda0 = 0,
                  interval = 0.01,
                  sigma = 0.01,
                  maxit = 1000) {

  # half interaction matrix
  alpha0 <- alpha
  alpha0[lower.tri(alpha0, diag = TRUE)] <- 0
  alpha0[alpha0 != 0] <- 1

  # basal species id
  id_basal <- which(colSums(alpha0) <= 0)
  n_basal <- length(id_basal)
  n_c <- ncol(alpha) - n_basal

  # basal equilibrium
  if (length(theta) != 1 && length(theta) != n_basal)
    stop("'theta' must be a scalar or have length of the number of basal species")

  if (any(theta <= 0))
    stop("'theta' must be positive")

  if (length(theta) == 1) {
    p_k <- drop(MCMCpack::rdirichlet(1, alpha = rep(theta, n_basal)))
  } else {
    p_k <- drop(MCMCpack::rdirichlet(1, alpha = theta))
  }

  f_k <- k0 * p_k

  # trophic position
  tp <- attr(alpha, "tp")

  ## k0 varies by basal species
  ## for basal species, use f_k
  ## for consumers, use mean k0 (`mean(f_k)`)
  w_k0 <- rep(mean(f_k), ncol(alpha))
  w_k0[1:n_basal] <- f_k

  # initialize lambda and r
  lambda <- lambda0
  eps <- exp(stats::rnorm(n_c, mean = 0, sd = sigma))

  x <- w_k0 * exp(-lambda * (tp - 1))
  x[-id_basal] <- x[-id_basal] * eps
  r <- drop(- alpha %*% x)

  # loop until all consumer's r < 0
  for (i in seq_len(maxit)) {
    lambda <- lambda + interval
    x <- w_k0 * exp(-lambda * (tp - 1))
    x[-id_basal] <- x[-id_basal] * eps
    r <- drop(- alpha %*% x)
    if (all(r[-id_basal] < 0)) break
  }

  if (any(r[-id_basal] >= 0))
    message("one or more consumer's r remain positive; increase 'maxit'?")

  cout <- cbind(r, x, tp)
  rownames(cout) <- seq_len(ncol(alpha))
  colnames(cout) <- c("r", "equilibrium", "tp")
  attr(cout, "lambda") <- lambda
  attr(cout, "iteration") <- i

  return(cout)
}

#' Preferential prey model
#'
#' @param n_species Integer.
#'  Number of species
#' @param n_basal Integer.
#'  Number of basal species
#' @param l Integer.
#'  Expected number of links in the upper triangle
#' @param theta Numeric.
#'  Scale parameter of an exponential distribution.
#'  Smaller values indicate greater trophic specialization.
#' @param cannibal Logical.
#'  If \code{TRUE}, cannibalism allowed
#' @param lower_tri Logical.
#'  If \code{TRUE}, lower triangle elements of the matrix will be returned
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export

ppm <- function(n_species,
                n_basal,
                l,
                theta,
                cannibal = FALSE,
                lower_tri = TRUE) {
  # n_species: number of species in the pool
  # n_basal: number of basal species
  # l: expected number of links
  # theta: scale parameter

  # verify inputs
  if (n_species < 3) stop("At least three species needed")
  if (n_basal < 1) stop("At least one basal species needed to construct a food web")
  if (n_basal >= n_species) stop("n_basal must be smaller than n_species")

  # number of consumers
  n_c <- n_species - n_basal

  # tp: trophic position
  # assign trophic position for basal species
  # assign -1 for consumers as initial values (to be updated)
  tp <- rep(-1, n_species)
  tp[seq_len(n_basal)] <- 1

  # beta parameter for beta distribution
  # determined so that E(L) = L
  if (l < n_c)
    stop("l must be at least equal to the number of consumers (n_species - n_basal)")

  if (cannibal) {
    max_l <- sum((n_basal + 1):n_species)

    if (l > max_l)
      stop(paste("maximum l is", max_l))

    if (n_species + n_basal < 1)
      stop("n_species + n_basal must be equal to or greater than 1")

    b <- ((n_species + n_basal - 1) * n_c) / (2 * (l - n_c)) - 1
  } else {
    max_l <- sum((n_basal + 1):n_species - 1)

    if (l > max_l)
      stop(paste("maximum l is", max_l))

    if (n_species + n_basal < 3)
      stop("n_species + n_basal must be equal to or greater than 3")

    b <- ((n_species + n_basal - 3) * n_c) / (2 * (l - n_c)) - 1
  }

  # vector for link proportions for all consumers
  v_xi <- stats::rbeta(n_c,
                       shape1 = 1,
                       shape2 = b)

  # vector for initial prey choice for all consumers
  v_i0 <- c(rep(-1, n_basal),
            sapply(X = (n_basal + 1):n_species,
                   FUN = function(j) resample(seq_len(j - 1), size = 1)))

  # realized number of prey nodes for all consumers
  # kappa follows a beta-binomial distribution
  # size parameter is the possible maximum number of prey nodes
  if (cannibal) {
    v_kappa <- c(rep(-1, n_basal),
                 stats::rbinom(n = n_c,
                               size = n_basal:(n_species - 1),
                               prob = v_xi))
  } else {
    v_kappa <- c(rep(-1, n_basal),
                 stats::rbinom(n = n_c,
                               size = (n_basal - 1):(n_species - 2),
                               prob = v_xi))
  }

  # alpha0: S x S interaction binary matrix
  # initialized with all zero
  # update the initial consumer's first prey
  alpha0 <- matrix(0, n_species, n_species)
  alpha0[v_i0[n_basal + 1], n_basal + 1] <- 1

  for (j in (n_basal + 1):n_species) {
    alpha0 <- extra_prey(alpha0 = alpha0,
                         j = j,
                         i0 = v_i0[j],
                         tp = tp,
                         kappa = v_kappa[j],
                         theta = theta,
                         cannibal = cannibal)

    v_n_prey <- colSums(alpha0)
    v_n_prey[v_n_prey == 0] <- 1

    tp_new <- (tp %*% alpha0) / v_n_prey + 1
    tp[j] <- tp_new[j]
  }

  if (lower_tri) {
    cout <- - alpha0 + t(alpha0)
  } else {
    cout <- - alpha0
  }

  attr(cout, "tp") <- tp

  return(cout)
}

#' Extra prey function
#'
#' @inheritParams ppm
#' @param alpha0 Numeric matrix.
#'  Adjacency matrix (\code{n_species x n_species})
#'  defining trophic interactions
#' @param j Integer.
#'  Consumer's index
#' @param i0 Integer.
#'  Index for the first prey
#' @param tp Numeric.
#'  Initial trophic positions
#' @param kappa Integer.
#'  Number of extra prey items
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
#'  Specify minimum (\code{min}) and maximum values (\code{max})
#'  for a uniform distribution.
#' @param attack List for attack rates.
#'  Specify minimum (\code{min}) and maximum values (\code{max})
#'  for a uniform distribution.
#' @param convert List for conversion efficiency.
#'  Specify minimum (\code{min}) and maximum values (\code{max})
#'  for a uniform distribution.
#' @param regulation List for self regulation (or intraspecific competition).
#'  Specify minimum (\code{min}) and maximum values (\code{max})
#'  for a uniform distribution.
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

#' Generate a food web based on the preferential prey model
#'
#' @inheritParams ppm
#' @inheritParams to_alpha
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export

foodweb <- function(n_species,
                    n_basal,
                    l,
                    theta,
                    cannibal = FALSE,
                    competition = list(min = 0,
                                       max = 0),
                    attack = list(min = 0,
                                  max = 1),
                    convert = list(min = 0,
                                   max = 1),
                    regulation = list(min = 0,
                                      max = 1)) {

  alpha0 <- ppm(n_species = n_species,
                n_basal = n_basal,
                theta = theta,
                l = l,
                cannibal = cannibal,
                lower_tri = FALSE)

  alpha <- to_alpha(alpha0 = alpha0,
                    competition = competition,
                    attack = attack,
                    convert = convert,
                    regulation = regulation)

  return(alpha)
}

#' Calculate a leading eigenvalue
#'
#' @param n_species Integer. Number of species.
#' @param r Numeric vector.
#'  Intrinsic growth rates for modeled species.
#' @param x0 Numeric vector.
#'  Equilibrium densities for modeled species.
#' @param alpha Numeric matrix.
#'  Interaction matrix for which linear stability is evaluated.
#' @param model Character string specifying a model type.
#'  Either \code{"ricker"} (Ricker),
#'  \code{"bh"} (Beverton-Holt), or
#'  \code{"glv"} (Generalized Lotka-Volterra).
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export

stability <- function(n_species,
                      r,
                      x0,
                      alpha,
                      model = "glv") {

  # check input -------------------------------------------------------------

  if (any(unique(dim(alpha)) != n_species))
    stop("dimension mismatch in alpha")

  if (!missing(r)) {
    if (length(r) != n_species)
      stop("dimension mismatch in alpha or r")
  }

  if (!missing(x0)) {
    if (length(x0) != n_species)
      stop("dimension mismatch in alpha or x0")
  }


  # get maximum absolute eigenvalue -----------------------------------------

  if (det(alpha) == 0) {
    ## if det(alpha) = 0, return NA
    max_lambda <- NA

  } else {

    if (model == "ricker" || model == "glv") {
      ## Ricker or GLV model

      if (missing(r) && !missing(x0)) {
        ## if equilibrium density provided,
        ## calculate r from alpha and x0
        r <- drop(-alpha %*% x0)
      }

      if (!missing(r) && missing(x0)) {
        ## if intrinsic growth provided,
        ## calculate x0 from alpha and r
        x0 <- drop(-solve(alpha) %*% r)
      }
    }

    if (model == "bh") {
      ## Beverton-Holt model

      if (missing(r) && !missing(x0)) {
        ## if equilibrium density provided,
        ## calculate r from alpha and x0
        r <- drop(log(1 + alpha %*% x0))
      }

      if (!missing(r) && missing(x0)) {
        ## if intrinsic growth provided,
        ## calculate x0 from alpha and r
        x0 <- drop(solve(alpha) %*% (exp(r) - 1))
      }
    }

    ## check negative equilibrium
    if (any(x0 < 0))
      stop("Negative equilibrium density is not allowed")

    ## get Jacobian matrix
    jm <- t(sapply(seq_len(n_species),
                   function(i) {
                     fun_partial(r = r[i],
                                 a = alpha[i, ],
                                 x0 = x0,
                                 i = i,
                                 model = model)
                   }))

    ## leading eigenvalue
    lambda <- eigen(jm)

    if (model == "ricker" || model == "bh") {
      ## discrete models
      max_lambda <- max(abs(Re(lambda$values)))
    } else {
      ## continuous models
      max_lambda <- max(Re(lambda$values))
    }

    attr(max_lambda, "jacobian") <- jm
    attr(max_lambda, "r") <- r
    attr(max_lambda, "x0") <- x0
  }

  return(max_lambda)
}
