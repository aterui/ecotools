
# two species case --------------------------------------------------------

x0 <- runif(2, 0, 1)
A <- rbind(c(runif(1, -1, 0), 0),
           c(runif(1, -1, 0), runif(1, -1, 0)))

test_that("stability, two-species system GLV", {

  ## glv;
  ## dx/dt = x (r[1] + A[1,1] x + A[1,2] y)
  ## dy/dt = y (r[2] + A[2,1] x + A[2,2] y)

  ## jacobian
  r <- drop(-(A %*% x0))
  j1_ <- c(r[1] + 2 * A[1, 1] * x0[1] + A[1, 2] * x0[2],
           A[1, 2] * x0[2])

  j2_ <- c(A[2, 1] * x0[1],
           r[2] + 2 * A[2, 2] * x0[2] + A[2, 1] * x0[1])

  J <- rbind(j1_, j2_)

  ## eigen values
  lm_max0 <- max(eigen(J)$values)
  lm_max1 <- stability(n_species = 2,
                       r = r,
                       x0 = x0,
                       alpha = A,
                       model = "glv")
  lm_max1 <- c(lm_max1) # remove attributes

  expect_equal(lm_max0, lm_max1, tolerance = 1e-6)
})

test_that("stability, two-species system Ricker", {

  ## ricker;
  ## x = x exp(r[1] + A[1,1] x + A[1,2] y)
  ## y = y exp(r[2] + A[2,1] x + A[2,2] y)

  ## jacobian
  r <- drop(-(A %*% x0))
  j1_ <- c((1 + x0[1] * A[1, 1]) * exp(r[1] + A[1, 1] * x0[1] + A[1, 2] * x0[2]),
           x0[1] * A[1, 2] * exp(r[1] + A[1, 1] * x0[1] + A[1, 2] * x0[2]))

  j2_ <- c(x0[2] * A[2, 1] * exp(r[2] + A[2, 2] * x0[2] + A[2, 1] * x0[1]),
           (1 + x0[2] * A[2, 2]) * exp(r[2] + A[2, 2] * x0[2] + A[2, 1] * x0[1]))

  J <- rbind(j1_, j2_)

  ## eigen values
  lm_max0 <- max(eigen(J)$values)
  lm_max1 <- stability(n_species = 2,
                       r = r,
                       x0 = x0,
                       alpha = A,
                       model = "ricker")
  lm_max1 <- c(lm_max1) # remove attributes

  expect_equal(lm_max0, lm_max1, tolerance = 1e-6)
})

test_that("stability, two-species system Beverton-Holt", {

  ## beverton;
  ## x = x exp(r[1]) / (1 + A[1, 1] * x + A[1, 2] * y)
  ## y = y exp(r[2]) / (1 + A[2, 1] * x + A[2, 2] * y)

  ## jacobian
  A0 <- abs(A)
  lambda <- drop(1 + A0 %*% x0)

  fx <- (1 + A0[1, 1] * x0[1] + A0[1, 2] * x0[2])
  fy <- (1 + A0[2, 1] * x0[1] + A0[2, 2] * x0[2])

  j1_ <- c(lambda[1] * fx^(-1) - (x0[1] * A0[1, 1] * lambda[1]) * fx^(-2),
           -(x0[1] * lambda[1] * A0[1, 2]) * fx^(-2))
  j2_ <- c(-(x0[2] * lambda[2] * A0[2, 1]) * fy^(-1),
           lambda[2] * fy^(-1) - (x0[2] * lambda[2] * A0[2, 2]) * fy^(-2))

  J <- rbind(j1_, j2_)

  ## eigen values
  lm_max0 <- max(eigen(J)$values)
  lm_max1 <- stability(n_species = 2,
                       r = log(lambda),
                       x0 = x0,
                       alpha = A0,
                       model = "bh")
  lm_max1 <- c(lm_max1) # remove attributes

  expect_equal(lm_max0, lm_max1, tolerance = 1e-6)
})


# more than 2 species -----------------------------------------------------

set.seed(123)
n_species_list <- c(3, 4)

for (n_species in n_species_list) {

  x0 <- runif(n_species, 0.1, 1)
  A <- matrix(runif(n_species^2, -1, 0), n_species, n_species)

  test_that(paste("stability,", n_species, "species system GLV"), {

    r <- drop(-(A %*% x0))

    # Jacobian
    J <- t(sapply(seq_len(n_species), function(i) {
      fun_partial(r = r[i],
                  a = A[i, ],
                  x0 = x0,
                  i = i,
                  model = "glv")
    }))

    lm_max0 <- max(eigen(J)$values)
    lm_max1 <- stability(n_species = n_species,
                         r = r,
                         x0 = x0,
                         alpha = A,
                         model = "glv")
    lm_max1 <- c(lm_max1) # remove attributes

    expect_equal(lm_max0, lm_max1, tolerance = 1e-6)
  })

  test_that(paste("stability,", n_species, "species system Ricker"), {

    r <- drop(-(A %*% x0))

    J <- t(sapply(seq_len(n_species), function(i) {
      fun_partial(r = r[i],
                  a = A[i, ],
                  x0 = x0,
                  i = i,
                  model = "ricker")
    }))

    lm_max0 <- max(eigen(J)$values)
    lm_max1 <- stability(n_species = n_species,
                         r = r,
                         x0 = x0,
                         alpha = A,
                         model = "ricker")
    lm_max1 <- c(lm_max1)

    expect_equal(lm_max0, lm_max1, tolerance = 1e-6)
  })

  test_that(paste("stability,", n_species, "species system Beverton-Holt"), {

    A0 <- abs(A)
    lambda_vec <- 1 + A0 %*% x0
    r <- log(lambda_vec)

    J <- t(sapply(seq_len(n_species), function(i) {
      fun_partial(r = r[i],
                  a = A0[i, ],
                  x0 = x0,
                  i = i,
                  model = "bh")
    }))

    lm_max0 <- max(eigen(J)$values)
    lm_max1 <- stability(n_species = n_species,
                         r = r,
                         x0 = x0,
                         alpha = A0,
                         model = "bh")
    lm_max1 <- c(lm_max1)

    expect_equal(lm_max0, lm_max1, tolerance = 1e-6)
  })
}
