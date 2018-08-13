make_temporal_weights <- function(N, TT) {
  if (TT > 1) {
    rows <- c()
    cols <- c()
    for (tt in 2:TT) {
      start_row <- (tt-1)*N +  1
      end_row <- start_row +  N - 1
      start_col <- start_row - N
      end_col <- end_row - N
      rows <- c(rows, start_row:end_row)
      cols <- c(cols, start_col:end_col)
    }
    TM <- sparseMatrix(i = rows, j = cols, x = TRUE, dims = c(N*TT, N*TT))
  } else {
    TM <- sparseMatrix(i = {}, j = {}, dims = c(N*TT, N*TT))
  }

  return(TM)
}

## Dummy class so we can use lmtest::coeftest
make_dummyclass <- function(x) {
  x <- list(coef = x)
  class(x) <- 'dummyclass'
  return(x)
}
coef.dummyclass <- function(x) {
  x$coef
}

## Main estimation class
RisSpatialProbit <- R6Class("RisSpatialProbit",
                  public = list(

                    spatial_dep = NULL,
                    temporal_dep = NULL,

                    N = NULL,
                    TT = NULL,
                    NT = NULL,

                    X = NULL,
                    y = NULL,

                    W_t = NULL,
                    W = NULL,
                    TM = NULL,

                    rho = NULL,
                    gamma = NULL,
                    beta = NULL,

                    rand_mat = NULL,

                    VC_mat = NULL,

                    initialize = function(X, y, W_t, N, TT,
                                          spatial_dep = TRUE, temporal_dep = TRUE, r = 100) {

                      ## Check inputs
                      if (TT == 1 & temporal_dep) {
                        warning("You requested temporal dependence estimation, but the data only cover a single period.
                                Temporal dependence will be set to FALSE.")
                        temporal_dep <- FALSE
                      }

                      ## Assign data
                      self$X <- X
                      self$W_t <- W_t
                      self$y <- y

                      ## Assign dependence switches
                      self$spatial_dep <- spatial_dep
                      self$temporal_dep <- temporal_dep

                      ## Assign data dimensions
                      self$N <- N
                      self$TT <- TT
                      self$NT <- N*TT

                      ## Make W, NTxNT spatial dep matrix
                      suppressWarnings((self$W <- kronecker(as(.sparseDiagonal(TT), 'dgCMatrix'), W_t)))
                      self$W <- as(self$W, "dgCMatrix")

                      ## Make TM, NTxNT temporal dep matrix
                      # NOTE: If TT = 1, returns NxN matrix of zeros
                      self$TM <- make_temporal_weights(N, TT)

                      ## Initialize interdependence parameters
                      self$rho <- 0
                      self$gamma <- 0

                      ## Initialize beta
                      self$beta <- rep(0, ncol(self$X))

                      ## Initialize importance matrix
                      rand_mat1 <- matrix(runif(N*TT*(r/2)),N*TT,r/2)
                      rand_mat2 <- 1-rand_mat1
                      rand_mat <- cbind(rand_mat1, rand_mat2)
                      self$rand_mat <- rand_mat

                    },

                    pack_theta = function() {
                      ##  Get vector containing all parameters to be estimated

                      theta <- self$beta
                      if (self$spatial_dep) {
                        # theta <- c(theta, atanh(self$rho))
                        theta <- c(theta, self$rho)
                      }
                      if (self$temporal_dep) {
                        # theta <- c(theta, atanh(self$gamma))
                        theta <- c(theta, self$gamma)
                      }

                      return(theta)
                    },

                    unpack_theta = function(theta) {
                      ## Takes the packed theta vector and returns a list with all parameters as attribues

                      K <- length(self$beta)
                      beta <- theta[1:K]
                      idx <- K+1
                      if (self$spatial_dep) {
                        # rho <- tanh(theta[idx])
                        rho <- theta[idx]
                        idx <- idx + 1
                      } else {
                        rho <- 0
                      }
                      if (self$temporal_dep) {
                        # gamma <- tanh(theta[idx])
                        gamma <- theta[idx]
                      } else {
                        gamma <- 0
                      }

                      out <- list(beta = beta, rho = rho, gamma = gamma)
                      return(out)
                    },

                    transform_theta = function(theta) {
                      ## Takes the packed theta vector and returns a theta vector where all
                      ## transformed params are returned to their original support
                      # i.e. from atanh(rho) to rho

                      transformed_theta <- rep(NA, length(theta))

                      K <- length(unlist(self$beta))
                      transformed_theta[1:K] <- theta[1:K]

                      idx <- K+1
                      if (self$spatial_dep) {
                        # transformed_theta[idx] <- tanh(theta[idx])
                        transformed_theta[idx] <- theta[idx]
                        idx <- idx + 1
                      }
                      if (self$temporal_dep) {
                        # transformed_theta[idx] <- tanh(theta[idx])
                        transformed_theta[idx] <- theta[idx]
                      }

                      return(transformed_theta)
                    },

                    llik = function(theta, verbose = FALSE) {

                      r <- ncol(self$rand_mat)
                      I <- diag(self$NT)

                      # Extract params
                      par <- self$unpack_theta(theta)
                      rho <- par$rho
                      gamma <- par$gamma
                      beta <- par$beta

                      # Safety checks
                      if (any(!is.finite(c(rho, gamma)))) {
                        return(1e9)
                      }

                      if (sum(abs(c(rho, gamma))) > 0.95) {
                        return(1e9)
                      }

                      # (Inverse) spatial multiplier
                      G = I - rho*self$W - gamma*self$TM
                      Z <- diag(as.vector(1-2*self$y))
                      vcov <- Z%*%solve(G, t(solve(G)))%*%t(Z)

                      # RIS matrices
                      ACH <- base::chol(solve(vcov), pivot = TRUE)
                      BCH <- solve(ACH)

                      GX <- base::solve(G, self$X%*%beta)
                      V <- -Z%*%GX

                      # Recursion
                      ln_prob <- ris_recursion(r, self$NT, BCH, self$rand_mat, as.vector(V))

                      # Log-likelihood
                      if (self$TT > 1) {
                        ll <- mean(colSums(ln_prob[(self$N+1):self$NT,]))
                      } else {
                        ll <- mean(colSums(ln_prob))
                      }

                      if (verbose) {
                        cat(ll, fill = TRUE)
                      }

                      return(-1*ll)
                    },

                    train = function(theta_init = NULL,
                                     method = c('Nelder-Mead', 'optimr', 'BFGS', 'newuoa'),
                                     vcov = TRUE,
                                     verbose = FALSE) {

                      method <- match.arg(method)

                      if (is.null(theta_init)) {
                        theta_init <- self$pack_theta()
                      }

                      if (method == 'Nelder-Mead') {
                        fit <- optim(par = theta_init, fn = self$llik, verbose = verbose, hessian = vcov)
                        if (vcov) {
                          self$compute_vcov(fit$hessian)
                        }
                      } else if (method == 'optimr') {
                        require(optimr)
                        fit <- optimr::optimr(par = theta_init, fn = self$llik, verbose = verbose, hessian = vcov)
                        if (vcov) {
                          self$compute_vcov(fit$hessian)
                        }
                      } else if (method == 'BFGS') {
                        fit <- optim(par = theta_init, fn = self$llik, verbose = verbose, method = 'BFGS', hessian = vcov)
                        if (vcov) {
                          self$compute_vcov(fit$hessian)
                        }
                      } else if (method == 'newuoa') {
                        require(nloptr)
                        fit <- nloptr::newuoa(x0 = theta_init, fn = self$llik, verbose = verbose)
                        if (vcov) {
                          self$compute_vcov(theta = fit$par) # Uses optimHessian
                        }
                      }

                      # Unpack vars
                      par <- self$unpack_theta(fit$par)
                      self$beta <- par$beta
                      self$rho <- par$rho
                      self$gamma <- par$gamma

                      # Return minimum negative likelihood
                      return(fit$value)
                    },

                    compute_vcov = function(hessian = NULL, theta = NULL) {

                      if (is.null(hessian)) {
                        if (is.null(theta)) {
                          theta <- self$pack_theta()
                        }
                        hessian <- optimHess(theta, fn = self$llik)
                      }

                      VC <- solve(Matrix(hessian))  # VC is inverse of observed Fisher Information
                      VC_pd <- as.matrix(nearPD(VC)$mat) # Force positive definite
                      self$VC_mat <- VC_pd

                    },

                    get_theta_names = function() {

                      theta_names <- c()

                      if (is.null(colnames(self$X))) {
                        beta_names <- paste0(paste0('beta_', 0:(ncol(self$X)-1)))
                      } else {
                        beta_names <- colnames(self$X)
                      }
                      theta_names <- c(theta_names, beta_names)

                      if (self$spatial_dep) {
                        theta_names <- c(theta_names, "rho")
                      }
                      if (self$temporal_dep) {
                        theta_names <- c(theta_names, "gamma")
                      }

                      return(theta_names)
                    },

                    coeftest = function(S = 1000, print_out = TRUE) {
                      ## Does two things:
                      # 1. Uses MC simulation to get VC matrix for transformed params
                      # 2. Prints pretty regression table (or returns list with coefs & SEs)

                      ## MC sim
                      theta <- self$pack_theta()
                      theta_trans <- self$transform_theta(theta)
                      theta_sim <- mvtnorm::rmvnorm(n = S, mean = theta, sigma = self$VC_mat)
                      theta_trans_sim <- t(apply(theta_sim, MARGIN = 1, FUN = self$transform_theta))
                      theta_trans_vc <- var(theta_trans_sim)

                      ## Pack in dummy class and print coef test
                      names(theta_trans) <- self$get_theta_names()
                      dc <- make_dummyclass(theta_trans)

                      if (print_out) {
                        lmtest::coeftest(x = dc, vcov. = theta_trans_vc)
                      } else {
                        out <- list(coef = theta_trans, se = sqrt(diag(theta_trans_vc)))
                        return(out)
                      }
                    }
                  ))

