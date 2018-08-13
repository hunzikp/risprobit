#########################
# Set parameters

N <- 10*10
TT <- 1
G <- 1

temporal_dep <- FALSE
spatial_dep <- TRUE

rho <- ifelse(spatial_dep, 0.5, 0)
gamma <- ifelse(temporal_dep, 0.3, 0)

beta <- c(0, 3)


#########################
# Simulate data

set.seed(3)
out <- simulate_data(N, TT, rho = rho, gamma = gamma,
                     beta = beta, X_params = c(0, 1))
X <- out$X
y <- out$y
W_t <- out$W_t

#########################
# Train

set.seed(1)
ris <- RisSpatialProbit$new(X = X, y = y, W_t = W_t, N = N, TT = TT,
                            spatial_dep = spatial_dep, temporal_dep = temporal_dep,
                            r = 100)
ris$train(verbose = FALSE, method = 'newuoa', vcov = TRUE)
ris$coeftest()


library(spatialprobit)
fit <- sar_probit_mcmc(as.vector(ris$y), ris$X, Matrix(ris$W, sparse = TRUE))
summary(fit)
