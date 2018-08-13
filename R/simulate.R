## Make spatial weights matrix from lattice
# N must have natural root
make_W_t <- function(N) {

  ras <- raster::raster(matrix(1, sqrt(N), sqrt(N)))
  spdf <- raster::rasterToPolygons(ras)
  B <- rgeos::gIntersects(spdf, byid = TRUE)
  diag(B) <- FALSE
  B <- B*1
  W_t <- B / rowSums(B)
  W_t <- Matrix(W_t, sparse = TRUE)

  return(W_t)
}

## Monte Carlo simulation of DIM-STAR data
# N must have natural root
simulate_data <- function(N, TT,
                          rho, gamma,
                          beta,
                          X_params = c(0, 1)) {

  ## Make NxN spatial weights matrix
  W_t <- make_W_t(N)

  ## Make NT x NT spatial weights matrix list (one per outcome)
  suppressWarnings(W <- kronecker(as(.sparseDiagonal(TT), 'dgCMatrix'), W_t))

  ## Make NT x NT temporal weights matrix list (one per outcome)
  TM <- make_temporal_weights(N, TT)

  ## Compute A
  A <- .sparseDiagonal(N*TT) - rho*W - gamma*TM

  ## Make predictors
  K <- length(beta)
  X <- cbind(1, matrix(runif((K-1)*N*TT, X_params[1], X_params[2]), N*TT, K-1))
  Xbeta <- X%*%beta

  ## Sample ystar
  eps <- rnorm(N*TT, 0, 1)
  ystar <- solve(A, Xbeta + eps)

  ## Sample y
  y <- (ystar>0)*1

  ## Prepare return list
  out <- list(y = y, X = X, W_t = W_t, ystar = ystar, A = A, N = N, TT = TT)

  return(out)
}
