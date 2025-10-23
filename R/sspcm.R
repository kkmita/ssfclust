#' #' Equation to calculate clusters' prototypes matrix $\hat{V}$.
#' #'
#' #' @param Phi Matrix with weights of size N x c.
#' #'
#' #' @param X Matrix with predictors of size N x p.
#' #'
#' #' @return Clusters' prototypes matrix of size c x p.
#' #' @export
#' #'
#' estimate_V <- function(Phi, X) {
#'   Phi_tilde <- sweep(Phi, 2, colSums(Phi), "/")
#'   return(t(t(X) %*% Phi_tilde))
#' }
#'
#'
#' #' Estimated T matrix with typicalities.
#' #'
#' #' @param X
#' #' a matrix *X* of dimension (N, p) containing predictor variables.
#' #'
#' #' @param V
#' #' a prototypes matrix of dimension (c, p)
#' #'
#' #' @param superF
#' #' the supervisionbinary matrix of the same dimension as *U*.
#' #'
#' #' @param alpha
#' #' a *N*-vector of observation-specific scaling factor values
#' #'
#' #' @param function_dist
#' #' A function of two arguments: matrices X and V of the same
#' #' number of columns.
#' #' It should return a matrix of (nrow(X) x nrow(V)) of distances
#' #' between each row of X and all rows of V.
#' #' In case of Euclidean distance, the result should not be squared!
#' #'
#' #' @param gammas
#' #' a *c*-vector of cluster-specific gamma parameter
#' #'
#' estimate_T <-
#'   function(
#'     X,
#'     V,
#'     superF,
#'     alpha,
#'     function_dist,
#'     gammas
#'   ) {
#'     D <- function_dist(X, V)^2
#'     G <- matrix(gammas, nrow = 1)[rep(1, nrow(superF)), ]
#'     M1 <- superF * alpha
#'     M2 <- matrix(rowSums(superF))[, rep(1, ncol(superF))] * alpha + 1
#'     Tm <- (G + (M1 * D) )/ (G + (M2 * D))
#'
#'     return(Tm)
#' }
#'
#'
#'
#' #' @export
#' SSPCM <- function(
#'     X,
#'     C,
#'     U = NULL,
#'     gammas = NULL,
#'     max_iter = 200,
#'     conv_criterion = 1e-4,
#'     function_dist = rdist::cdist,
#'     alpha = NULL,
#'     superF = NULL
#' ) {
#'   if (is.null(U)) {
#'     Tm <- matrix(runif(nrow(X)*C), ncol=C)
#'   } else{
#'     Tm <- U
#'   }
#'
#'   # Rows of U should sum up to 1
#'   Tm <- t(apply(Tm, 1, function(x) x / sum(x)))
#'
#'   if (is.null(gammas)) {
#'     gammas <- rep(1, C)
#'   }
#'
#'   counter = 0
#'   T_history <- list()
#'   V_history <- list()
#'   Phi_history <- list()
#'
#'   for (iter in 1:max_iter) {
#'     counter <- counter + 1
#'     Tm_previous_iter <- Tm
#'
#'     Phi <- Tm_previous_iter^2
#'
#'     # Modify `Phi` if running semi-supervised PCM
#'     if (!is.null(alpha)) {
#'       Tm_alpha <- alpha * (Tm_previous_iter - superF)^2
#'       Phi <- Phi + Tm_alpha
#'     }
#'
#'     Phi_history[[counter]] <- Phi
#'
#'     V <- estimate_V(Phi, X)
#'
#'     V_history[[counter]] <- V
#'
#'     Tm <- estimate_T(
#'       X = X,
#'       V = V,
#'       superF = superF,
#'       alpha = alpha,
#'       function_dist = function_dist,
#'       gammas = gammas
#'     )
#'
#'     T_history[[counter]] <- Tm
#'
#'     conv_iter <- base::norm(Tm - Tm_previous_iter, type="F")
#'
#'     if (conv_iter < conv_criterion) {
#'       break
#'     }
#'   }
#'
#'   z <- list(
#'     Tm = Tm,
#'     V = V,
#'     function_dist = function_dist,
#'     counter = counter,
#'     gammas = gammas,
#'     V_history = V_history,
#'     T_history = T_history,
#'     Phi_history = Phi_history
#'   )
#'
#'   class(z) <- "sspcm"
#'
#'   return(z)
#' }
#'
#'
#'
#' #' Soft assignment score function
#' #'
#' #' @param object
#' #' @param newdata
#' #'
#' #' @return
#' #'
#' #' @export
#' #' @examples
#' predict.sspcm <- function(object, newdata) {
#'   output <- estimate_Tm(
#'     X = newdata,
#'     V = object$V,
#'     superF = NULL,
#'     alpha = NULL,
#'     function_dist = object$function_dist,
#'     gammas = object$gammas
#'   )
#'   return(output)
#' }
