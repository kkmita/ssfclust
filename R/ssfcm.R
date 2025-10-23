#' Creates DHE (stands for "distances horizontally exploded") and DVE
#' (stands for "distances vertically exploded") matrices.
#'
#' @param A Matrix of size N x c.
#' @param vertical Boolean switch.
#' If `TRUE`, create DVE (vertical explosion).
#' If `FALSE`, create DHE (horizontal explosion).
#'
#' @return Matrix of size Nc x c
#'
#' @examples
#' A <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 2, byrow = TRUE)
#' AVE <- dheve(A, vertical = TRUE)
#' AHE <- dheve(A, vertical = FALSE)
#'
dheve <- function(A, vertical) {
  if (vertical == TRUE) {
    elements <- A[rep(1:nrow(A), each=ncol(A)), ]
  } else {
    elements <- matrix(c(t(A)))[, rep(1, ncol(A))]
  }
  return(elements)
}


#' Aggregates elements of DHE and DVE matrices in a step to build
#' evidence matrix E.
#'
#' @param dhe DHE matrix of size Nc x c.
#' @param dve DVE matrix of size Nc x c.
#'
#' @return Matrix of size Nc x 1.
#'
gamma_fcm <- function(dhe, dve) {
  # 1 / ((dhe/dve) %*% matrix(rep(1, ncol(dhe))))
  1 / rowSums(dhe/dve)
}


#' Rearranges elements of input matrix from a block matrix with vertical blocks
#' (column vectors) to a block matrix with horizontal blocks (row vectors).
#'
#' @param A Matrix of size Nc x 1.
#' @param c Number of columns in the wanted matrix.
#' Associated with the number of clusters.
#'
#' @return Matrix of size N x c.
#'
xi_fcm <- function(A, c) {
  matrix(A, ncol=c, byrow=TRUE)
}


#' Calculates data evidence matrix E from distances matrix D.
#'
#' @param D Distances matrix of size N x c.
#'
#' @return Matrix of size N x c.
#' @export
#'
calculate_evidence <- function(D) {
  dve <- dheve(D, vertical = TRUE)
  dhe <- dheve(D, vertical = FALSE)
  output <- xi_fcm(A = gamma_fcm(dhe, dve), c = ncol(dhe))
  return(output)
}


#' Estimated U matrix with memberships.
#'
#' @param X
#' a matrix *X* of dimension (N, p) containing predictor variables.
#'
#' @param V
#' a prototypes matrix of dimension (c, p)
#'
#' @param F_
#' the supervision  binary matrix of the same dimension as *U*.
#'
#' @param alpha
#' the scaling factor, a floating point > 0.
#'
#' @param function_dist
#' A function of two arguments: matrices X and V of the same
#' number of columns.
#' It should return a matrix of (nrow(X) x nrow(V)) of distances
#' between each row of X and all rows of V.
#' In case of Euclidean distance, the result should not be squared!
#'
estimate_U <-
  function(
    X,
    V,
    superF,
    alpha,
    function_dist
  ) {
    D <- function_dist(X, V)^2
    E <- calculate_evidence(D)

    U <- (E + alpha * superF) / (1 + alpha * rowSums(superF))

    return(U)

    # if (is.null(alpha)) {
    #   return(E)
    # } else {
    #   M <- matrix((1 / (1 + alpha)), ncol = 1)[, rep(1, ncol(F_))]
    #   F_alpha = F_*(alpha/(1+alpha))
    #
    #   return(M*E + F_alpha)
    # }
  }


#' Equation to calculate clusters' prototypes matrix $\hat{V}$.
#'
#' @param Phi Matrix with weights of size N x c.
#'
#' @param X Matrix with predictors of size N x p.
#'
#' @return Clusters' prototypes matrix of size c x p.
#' @export
#'
estimate_V <- function(Phi, X) {
  # Phi_tilde <- sweep(Phi, 2, colSums(Phi), "/")
  # return(t(t(X) %*% Phi_tilde))
  V <- (t(Phi) %*% X) / colSums(Phi)
  return(V)
}


#' Semi-Supervised Fuzzy C-Means model.
#'
#' @description
#' If *alpha* and *F_* are not supplied (their default values are `NULL`),
#' then a regular unsupervised Fuzzy C-Means algorithm is fitted.
#'
#' @param X
#' a matrix *X* with predictor variables.
#'
#' @param C
#' a number of clusters to find.
#'
#' @param U
#' optionally: a first memberships matrix to initialize the algorithm.
#' Used mainly for reproducibility to compare calculations with other packages
#' (e.g. in Python).
#'
#' @param function_dist
#' A function of two arguments: matrices X and V of the same
#' number of columns.
#' It should return a matrix of (nrow(X) x nrow(V)) of distances
#' between each row of X and all rows of V.
#' In case of Euclidean distance, the result should not be squared!
#'
#' @param alpha
#' the scaling factor, a floating point > 0.
#'
#' @param F_
#' the supervision  binary matrix of the same dimension as *U*.
#'
#' @export
#'
#' @examples
#' # simulate 50 obs from N2((5, 5)) and 50 obs from N2((7, 7))
#' library(MASS)
#' X <- rbind(
#'   MASS::mvrnorm(50, mu=c(5, 8), Sigma=matrix(c(3, 0, 0, 3), ncol=2)),
#'   MASS::mvrnorm(50, mu=c(7, 10), Sigma=matrix(c(3, 0, 0, 3), ncol=2))
#'   )
#'
#' # simulate supervision for 10% of each class
#' F_ <- matrix(0, nrow=100, ncol=2)
#' F_[sample(1:50, 10), 1] <- 1
#' F_[sample(51:100, 10), 2] <- 1
#'
#' model <- SSFCM(X=X, C=2)
#' model_ss <- SSFCM(X=X, C=2, alpha=1, F_=F_)
#'
#' acc_unsupervised_1 <- sum(apply(model$U, 1, which.max) == c(rep(1, 50), rep(2, 50)))
#' acc_unsupervised_2 <- sum(apply(model$U, 1, which.max) == c(rep(2, 50), rep(1, 50)))
#' acc_supervised <- sum(apply(model_ss$U, 1, which.max) == c(rep(1, 50), rep(2, 50)))
#'
SSFCM <- function(
    X,
    C,
    U = NULL,
    max_iter = 200,
    conv_criterion = 1e-4,
    function_dist = rdist::cdist,
    alpha = NULL,
    superF = NULL
) {

  if (!is.null(a) && length(a) != 1) {
    stop("'a' must be either NULL or a scalar (length 1 value).")
  }

  if (is.null(U)) {
    U <- matrix(runif(nrow(X)*C), ncol=C)
  }

  # Rows of U should sum up to 1
  U <- t(apply(U, 1, function(x) x / sum(x)))

  counter = 0
  U_history <- list()
  V_history <- list()
  Phi_history <- list()

  for (iter in 1:max_iter) {
    counter <- counter + 1
    U_previous_iter <- U

    Phi <- U_previous_iter^2 + (U_previous_iter-F)^2 * alpha * rowSums(superF)

    # Modify `Phi` if running semi-supervised FCM
    # if (!is.null(alpha) && prototypes) {
    #   U_alpha <- alpha * (U_previous_iter - F_)^2
    #   U_alpha[h_indices, ] <- 0
    #   Phi <- Phi + U_alpha
    # }

    Phi_history[[counter]] <- Phi

    V <- estimate_V(Phi, X)

    V_history[[counter]] <- V

    U <- estimate_U(
        X = X,
        V = V,
        superF = superF,
        alpha = alpha,
        function_dist = function_dist)

    U_history[[counter]] <- U

    conv_iter <- base::norm(U - U_previous_iter, type="F")

    if (conv_iter < conv_criterion) {
      break
    }
  }

  z <- list(
    U = U,
    V = V,
    function_dist = function_dist,
    counter = counter,
    V_history = V_history,
    U_history = U_history,
    Phi_history = Phi_history
  )

  class(z) <- "ssfcm"

  return(z)
}



#' Semi-Supervised Fuzzy C-Means model.
#'
#' @description
#' If *alpha* and *F_* are not supplied (their default values are `NULL`),
#' then a regular unsupervised Fuzzy C-Means algorithm is fitted.
#'
#' @param X
#' a matrix *X* with predictor variables.
#'
#' @param C
#' a number of clusters to find.
#'
#' @param U
#' optionally: a first memberships matrix to initialize the algorithm.
#' Used mainly for reproducibility to compare calculations with other packages
#' (e.g. in Python).
#'
#' @param function_dist
#' A function of two arguments: matrices X and V of the same
#' number of columns.
#' It should return a matrix of (nrow(X) x nrow(V)) of distances
#' between each row of X and all rows of V.
#' In case of Euclidean distance, the result should not be squared!
#'
#' @param alpha
#' the scaling factor, a floating point > 0.
#'
#' @param F_
#' the supervision  binary matrix of the same dimension as *U*.
#'
#' @export
#'
#' @examples
#' # simulate 50 obs from N2((5, 5)) and 50 obs from N2((7, 7))
#' library(MASS)
#' X <- rbind(
#'   MASS::mvrnorm(50, mu=c(5, 8), Sigma=matrix(c(3, 0, 0, 3), ncol=2)),
#'   MASS::mvrnorm(50, mu=c(7, 10), Sigma=matrix(c(3, 0, 0, 3), ncol=2))
#'   )
#'
#' # simulate supervision for 10% of each class
#' F_ <- matrix(0, nrow=100, ncol=2)
#' F_[sample(1:50, 10), 1] <- 1
#' F_[sample(51:100, 10), 2] <- 1
#'
#' model <- SSFCM(X=X, C=2)
#' model_ss <- SSFCM(X=X, C=2, alpha=1, F_=F_)
#'
#' acc_unsupervised_1 <- sum(apply(model$U, 1, which.max) == c(rep(1, 50), rep(2, 50)))
#' acc_unsupervised_2 <- sum(apply(model$U, 1, which.max) == c(rep(2, 50), rep(1, 50)))
#' acc_supervised <- sum(apply(model_ss$U, 1, which.max) == c(rep(1, 50), rep(2, 50)))
#'
FCM <- function(
    X,
    C,
    U = NULL,
    max_iter = 200,
    conv_criterion = 1e-4,
    function_dist = rdist::cdist
) {

  if (is.null(U)) {
    U <- matrix(runif(nrow(X)*C), ncol=C)
  }

  # Rows of U should sum up to 1
  U <- t(apply(U, 1, function(x) x / sum(x)))

  counter = 0
  U_history <- list()
  V_history <- list()
  Phi_history <- list()

  for (iter in 1:max_iter) {
    counter <- counter + 1
    U_previous_iter <- U

    Phi <- U_previous_iter^2

    # Modify `Phi` if running semi-supervised FCM
    # if (!is.null(alpha) && prototypes) {
    #   U_alpha <- alpha * (U_previous_iter - F_)^2
    #   U_alpha[h_indices, ] <- 0
    #   Phi <- Phi + U_alpha
    # }

    Phi_history[[counter]] <- Phi

    V <- estimate_V(Phi, X)

    V_history[[counter]] <- V

    D <- function_dist(X, V)^2
    U <- calculate_evidence(D)

    U_history[[counter]] <- U

    conv_iter <- base::norm(U - U_previous_iter, type="F")

    if (conv_iter < conv_criterion) {
      break
    }
  }

  z <- list(
    U = U,
    V = V,
    function_dist = function_dist,
    counter = counter,
    V_history = V_history,
    U_history = U_history,
    Phi_history = Phi_history
  )

  class(z) <- "fcm"

  return(z)
}

