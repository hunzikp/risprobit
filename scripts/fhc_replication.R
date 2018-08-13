library(risprobit)
library(Matrix)

## Read the data
data.df <- read.csv('data/hb_ksg_rep_Africa.csv', header = FALSE)
names(data.df)[1:2] <- c('cowid', 'year')
names(data.df)[4] <- 'incidence'
regressor_names <- c("neighpol", "neighpolsq", "neighlgdpl", "lnpop", "polity2l", "polity2sq", "postcoldw", "lgdp96l")
names(data.df)[8:15] <- regressor_names
X <- cbind(1, as.matrix(data.df[,8:15]))
y <- data.df$incidence

W <- as.matrix(read.csv('data/hb_ksg_contig_W.csv', header = FALSE))
denom <- ifelse(rowSums(W)==0, 1, rowSums(W))
W <- Matrix(W/denom, sparse = TRUE)

TM <- Matrix(as.matrix(read.csv('data/hb_ksg_TL.csv', header = FALSE)), sparse = TRUE)
TM <- as(TM, 'dgCMatrix')

## 'Naive' Probit
data.df$ytl <- as.vector(TM%*%y)
data.df$naive_ysl <- as.vector(W%*%y)

regr <- c(regressor_names, 'ytl', 'naive_ysl')
form <- paste("incidence ~",paste(regr, collapse = ' + '))
glm.fit <- glm(formula = as.formula(form), family = binomial(link = 'probit'), data = data.df[31:nrow(data.df),])
summary(glm.fit)


## Spatio-temporal Probit
N <- nrow(X)/2
TT <- 2
NT <- N*TT
W_dummy <- as(matrix(0, N, N), "dgCMatrix")
ris <- RisSpatialProbit$new(X = X, y = y, W_t = W_dummy, N = N, TT = TT,
                            spatial_dep = TRUE, temporal_dep = TRUE,
                            r = 100)
ris$W <- W
ris$TM <- TM
ris$train(method = 'BFGS', vcov = FALSE, verbose = TRUE)

