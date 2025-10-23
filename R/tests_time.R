library(ppclust)
library(fclust)
library(e1071)
library(ssfclust)
library(bench)


createX <- function(.N, .p) {
  X <- rbind(
    matrix(rnorm(.N), ncol = .p),
    matrix(rnorm(.N, mean = 1), ncol = .p),
    matrix(rnorm(.N, mean = 2), ncol = .p)
  )
  X
}


X <- createX(1e2, 2)

res <- bench::mark(
  ppclust = ppclust::fcm(x = X, centers = 3),
  fclust = fclust::FKM(X = X, k = 3),
  e1071 = e1071::cmeans(x = X, centers = 3),
  ssfclust = ssfclust::FCM(X = X, C = 3),
  check = FALSE
)


v_ppclust <- ppclust::fcm(x = X, centers = 3)$v |> round(2)
v_fclust <- fclust::FKM(X = X, k = 3)$H |> round(2)
v_e1071 <- e1071::cmeans(x = X, centers = 3)$centers |> round(2)
V_mine <- ssfclust::FCM(X = X, C = 3)$V |> round(2)



X <- rbind(
  matrix(rnorm(1e5), ncol = 5),
  matrix(rnorm(1e5, mean = 1), ncol = 5),
  matrix(rnorm(1e5, mean = 2), ncol = 5)
)


U0 <- matrix(runif(nrow(X)*3), ncol = 3)
U0 <- t(apply(U0, 1, function(x) x / sum(x)))


res <- bench::mark(
  ppclust = ppclust::fcm(x = X, centers = 3),
  # fclust = fclust::FKM(X = X, k = 3),
  e1071 = e1071::cmeans(x = X, centers = 3),
  ssfclust = ssfclust::FCM(X = X, C = 3),
  check = FALSE
)

t.mine1 <- Sys.time()
model_mine <- ssfclust::FCM(X = X, C = 3, U = U0)
t.mine2 <- Sys.time()

print(t.mine2 - t.mine1)


s1 <- Sys.time()
model_e = e1071::cmeans(x = X, centers = 3)
s2 <- Sys.time()

print(s2 - s1)

bench::mark(
  e1071 = e1071::cmeans(x = X, centers = 3),
  ssfclust = ssfclust::FCM(X = X, C = 3),
  fclust = fclust::FKM(X = X, k = 3),
  check = FALSE
)
