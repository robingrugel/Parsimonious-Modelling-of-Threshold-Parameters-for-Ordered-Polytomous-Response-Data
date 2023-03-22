library(purrr) # for cross()
library(haven) # for the convenient import of .sav
library(dplyr) # efficient data.frame operations
library(tidyr) # 
library(cmdstanr) # Stan interface through CmdStan
library(splines) # for easy B-spline computation
library(forcats) # for special operations with factors
library(loo) # für loo_compare() etc
library(purrr)
library(readr)
library(tables) # easy generation of latex tables
library(ggplot2)
library(RColorBrewer)
library(tidyr)
library(tikzDevice)

# seed for practical applications
seedygonzales <- 2024


grmcrclogit <- function(eta, tresh) {
  plogis(c(tresh, Inf) - eta) - plogis(c(-Inf, tresh) - eta)
}

con_thresh_matrix <- function(K = 100,
                              shape_mat = matrix(c(3, 2, 1.5, 1, 1,
                                                   0.5, 0.5, 0.5, 0.5, 1), ncol = 2),
                              scale_par = 10) {
  k <- 0:(K-2)
  apply(shape_mat, 1,
        function(x) {
          scale_par * (qbeta(k/(max(k)), shape1 = x[1], shape2 = x[2]) - 0.5)
        })
}

simulate_simple1 <- function(abil_par,
                             thresh_par,
                             disc_par) {
  
  I <- length(disc_par)
  P <- length(abil_par)
  R <- matrix(nrow = P, ncol = I)
  K <- nrow(thresh_par) + 1
  
  for (p in 1:P) {
    for (i in 1:I) {
      all_probs <- grmcrclogit(eta = disc_par[i] * abil_par[p],
                               tresh = thresh_par[, i])
      # print(paste0("i=",i,",p=",p))
      R[p, i] <- sample(1:K, 1, prob = all_probs)
    }
  }
  
  return(R)
}


p_list_raw_I6_K25 <- list(
  
  sim_K = 25,
  
  sim_I = 6,
  
  sim_P = c(25, 50),
  
  sim_disc_par = list(list(name = "excellent", 
                           value = rep(3, 6)),
                      list(name = "acceptable",
                           value = rep(1, 6))),
  
  sim_thresh_par = list(con_thresh_matrix(K = 25,
                                          scale_par = 20,
                                          shape_mat = matrix(c(4, 2, 0.25, 0.5, 7.5, 2.5,
                                                               0.5, 0.5, 0.25, 0.5, 7.5, 2.5), ncol = 2)))
)
p_list_I6_K25 <- p_list_raw_I6_K25 %>% cross()

p_list_raw_I6_K100 <- list(
  
  sim_K = 100,
  
  sim_I = 6,
  
  sim_P = c(25, 50),
  
  sim_disc_par = list(list(name = "excellent", 
                           value = rep(3, 6)),
                      list(name = "acceptable",
                           value = rep(1, 6))),
  
  sim_thresh_par = list(con_thresh_matrix(K = 100,
                                          scale_par = 20,
                                          shape_mat = matrix(c(4, 2, 0.25, 0.5, 7.5, 2.5,
                                                               0.5, 0.5, 0.25, 0.5, 7.5, 2.5), ncol = 2)))
)
p_list_I6_K100 <- p_list_raw_I6_K100 %>% cross()
p_list_I6 <- append(p_list_I6_K25, p_list_I6_K100)


p_list <- p_list_I6

