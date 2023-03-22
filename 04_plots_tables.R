source("00_basicfunctions.R")

options(tikzLatexPackages=c(getOption("tikzLatexPackages"),
"\\usepackage{amsfonts}",
"\\usepackage{bbm}"))

# example discrimination parameter for one item
# example ability parameters for one item
example_grm <- function(thresh_par = c(-2.1, -1.1, 0.5, 1.7),
                        abil_par = c(-2, -0.75, 2.5, 3),
                        disc_par = 1.3) {
    
    limitx <- 6
    x_all <- seq(-(limitx - 0.0001), (limitx - 0.0001), length.out = 4998)
    # y_all <- c(0, dnorm(x_all, 0, 1), 0) / dnorm(0, 0, 1)
    y_all <- c(0, dlogis(x_all, 0, 1), 0) / dlogis(0, 0, 1)
    x_all <- c(-limitx, x_all, limitx)
    
    basic_grm <- ggplot() +
      theme_bw() +
      xlim(range(abil_par) + c(0, 1)) + 
      ylim(c(min(disc_par * abil_par) - limitx, max(disc_par * abil_par) + limitx)) + 
      xlab("$\\theta_p$") + ylab("$\\Gamma_i$") +
      geom_abline(intercept = 0,
                  slope = disc_par,
                  color = "red") +
      geom_hline(yintercept = thresh_par,
                 alpha = 0.4,
                 linetype = 'dashed')
    
    for (a in abil_par) {
      dat_a <- data.frame(x = a + y_all, y = disc_par * a + x_all) %>%
        mutate(category = cut(y, 
                              breaks = c(-Inf, thresh_par, Inf), 
                              labels = 0:length(thresh_par)))
      
      dat_a <- dat_a %>% 
        add_row(x = a, y = thresh_par[1], category = "0") %>%
        add_row(x = a, y = thresh_par[length(thresh_par)], category = paste0(length(thresh_par)))
      
      for (t in 2:length(thresh_par)) {
        dat_a <- dat_a %>%
          add_row(x = a, y = thresh_par[t], category = paste0(t - 1)) %>%
          add_row(x = a, y = thresh_par[t - 1], category = paste0(t - 1))
      }
      
      # basic_grm <- basic_grm +
      #   geom_polygon_pattern(data = dat_a,
      #                        aes(x = x, y = y, pattern = category, fill = category), 
      #                        alpha = 0.2,
      #                        # color = "black", 
      #                        # pattern_fill = "black",
      #                        pattern_angle = 45,
      #                        pattern_density = 0.025,
      #                        pattern_spacing = 0.002)
      
      basic_grm <- basic_grm +
        geom_polygon(data = dat_a,
                     aes(x = x, y = y, fill = category),
                     alpha = 0.9, linetype = "solid")
    }
    
    basic_grm <- basic_grm +
      geom_segment(data = data.frame(x = abil_par,
                                     xend = abil_par + max(y_all),
                                     y = abil_par * disc_par,
                                     yend = abil_par * disc_par),
                   aes(x = x, xend = xend, y = y, yend = yend), 
                   linetype = "dotted")
    
    # basic_grm <- basic_grm + scale_pattern_fill_manual(values = 0:length(thresh_par))
    
    basic_grm <- basic_grm +
      scale_fill_manual(values = brewer.pal(length(thresh_par) + 1, name = "Blues"),
                        labels = 1:(length(thresh_par)+1)) +
      theme(legend.position="bottom", 
            legend.box.spacing = unit(3, "pt"))
    
    # basic_grm <- basic_grm + 
    #   annotation_custom(textGrob("xyz", x = -.05, gp = gpar(col = "red")), 
    #                     ymin=-3, ymax=-3, xmin=-Inf, xmax=Inf) +
    #   annotation_custom(linesGrob(x = c(-.02, .02),  gp = gpar(col = "red")) ,
    #                     ymin=-3, ymax=-3, xmin=-Inf, xmax=Inf)
    
    return(basic_grm)
}
# example_grm()

grm_crc <- function(thresh_par = c(-2.1, -1.1, 0.5, 1.7),
                    disc_par = 1.3,
                    lims = c(-4,4)) {
    grmcrc <- as.data.frame(t(sapply(seq(lims[1], lims[2], length.out = 1000),
                                     function(t)
                                       grmcrclogit(eta = disc_par * t, tresh = thresh_par))))
    
    grmcrc <- grmcrc %>% pivot_longer(cols = everything()) %>%
      mutate(theta = rep(seq(lims[1], lims[2], length.out = 1000), each = length(thresh_par) + 1))
    grmcrc$name <- as.factor(grmcrc$name)
    levels(grmcrc$name) <- paste0(1:(length(thresh_par)+1))
    
    ggplot(grmcrc, aes(x = theta, y = value, col = name)) + 
      geom_line(linewidth = 1.2) +
      xlab("$\\theta_p$") + 
      ylab("$p_{ik}(\\theta)$") + theme_bw() +
      geom_segment(data = data.frame(x = thresh_par,
                                     xend = thresh_par,
                                     y = 0,
                                     yend = 1),
                   aes(x = x, xend = xend, y = y, yend = yend), 
                   linetype = "dashed", color = "black", alpha = 0.4) +
      geom_point(data = data.frame(x = c(thresh_par[1] / disc_par,
                                         (thresh_par[-length(thresh_par)] + thresh_par[-1]) / (2 * disc_par),
                                         thresh_par[length(thresh_par)] / disc_par),
                                  y = c(0.5,
                                        diag(sapply((thresh_par[-length(thresh_par)] + thresh_par[-1]) / (2 * disc_par),
                                                function(t) grmcrclogit(eta = disc_par * t, tresh = thresh_par))[-c(1, length(thresh_par) + 1),]),
                                        0.5)),
                   aes(x = x, y = y), 
                   color = "black", alpha = 0.4, size = 3) +
      scale_color_manual(name = "category", values = brewer.pal(length(thresh_par) + 1, name = "Blues")) +
      theme(legend.position="bottom", 
            legend.box.spacing = unit(3, "pt"))
}
# grm_crc()


grm_crb <- function(thresh_par = c(-2.1, -1.1, 0.5, 1.7),
                    disc_par = 1.3,
                    abil_par = -0.75,
                    lims = c(-4,4)) {
  
  grmcrb <- as.data.frame(t(sapply(seq(lims[1], lims[2], length.out = 1000),
                                   function(t)
                                     plogis(c(thresh_par, Inf)- disc_par * t))))
  
  grmcrb <- grmcrb %>% pivot_longer(cols = everything(), names_to = "category") %>%
    mutate(theta = rep(seq(lims[1], lims[2], length.out = 1000), each = length(thresh_par) + 1))
  grmcrb$category <- as.factor(grmcrb$category)
  levels(grmcrb$category) <- paste0(1:(length(thresh_par) + 1))
  
  ggplot(grmcrb, aes(x = theta, y = value, col = category)) + geom_line(size = 1.2) +
    xlab("$\\theta_p$") + ylab("$P(Y_{i} \\leq k | \\theta_p)$") + theme_bw() +
    geom_segment(data = data.frame(x = abil_par,
                                   xend = abil_par,
                                   y = 0,
                                   yend = 1),
                 aes(x = x, xend = xend, y = y, yend = yend), 
                 linetype = "dashed", color = "black", alpha = 0.4) +
    geom_point(data = data.frame(x = rep(abil_par,(length(thresh_par) + 2)), 
                                y = plogis(c(-Inf, thresh_par, Inf) - disc_par * abil_par)),
                 aes(x = x, y = y), 
               color = "black", alpha = 0.4, size = 3) +
    scale_color_manual(values = brewer.pal(length(thresh_par) + 1, name = "Blues")) +
    theme(legend.position="bottom", 
          legend.box.spacing = unit(3, "pt"))
}
# grm_crb()


d1grmcrclogit <- function(disc_par,
                          abil_par,
                         thresh_par) {
  eta <- disc_par * abil_par
  - disc_par * (dlogis(c(thresh_par, Inf) - eta) - dlogis(c(-Inf, thresh_par) - eta))
}


d2grmcrclogit <- function(disc_par,
                          abil_par,
                          thresh_par) {
  eta <- disc_par * abil_par
  
  skal <- disc_par^2
  term1 <- (exp(eta + thresh_par[-1]) * (exp(eta) - exp(thresh_par[-1]))) / 
    (exp(eta) + exp(thresh_par[-1]))^3
  term2 <- (exp(thresh_par[-length(thresh_par)] + eta) * (exp(eta) + exp(thresh_par[-length(thresh_par)]))) / 
    (exp(eta) + exp(thresh_par[-length(thresh_par)]))^3
  
  c(skal * (exp(eta + thresh_par[1]) * (exp(eta) - exp(thresh_par[1]))) / 
      (exp(eta) + exp(thresh_par[1]))^3,
    skal * (term1 - term2),
    skal * -(exp(eta + thresh_par[length(thresh_par)]) * (exp(eta) - exp(thresh_par[length(thresh_par)]))) / 
      (exp(eta) + exp(thresh_par[length(thresh_par)]))^3)
}

catinfologit <- function(disc_par,
                         abil_par,
                         thresh_par) {
  pik <- grmcrclogit(disc_par * abil_par, thresh_par)
  pik1 <- d1grmcrclogit(disc_par, abil_par, thresh_par)
  pik2 <- d2grmcrclogit(disc_par, abil_par, thresh_par)
    
  exp(log(pik1^2) - log(pik)) - pik2
}


grm_info_logit <- function(thresh_par = c(-2.1, -1.1, 0.5, 1.7),
                           disc_par = 1.3,
                           lims = c(-4,4)) {
  
  grminfo <- t(sapply(seq(lims[1], lims[2], length.out = 1000),
                                   function(t) catinfologit(thresh_par = thresh_par,
                                                            disc_par = disc_par,
                                                            abil_par = t)))
  
  grminfo <- as.data.frame(cbind(grminfo, rowSums(grminfo)))
  grminfo <- grminfo %>% pivot_longer(cols = everything(), names_to = "category") %>%
    mutate(theta = rep(seq(lims[1], lims[2], length.out = 1000), each = length(thresh_par) + 2))
  
  grminfo$category <- as.factor(grminfo$category)
  levels(grminfo$category) <- c(paste0(1:(length(thresh_par) + 1)), "Item information")
  
  ggplot(grminfo, aes(x = theta, y = value, col = category)) + geom_line(linewidth = 1.2) +
    xlab("$\\theta_p$") + ylab("$I_{ik}(\\theta_p) / I_{i}(\\theta_p)$") + theme_bw() +
    scale_color_manual(values = c(brewer.pal(length(thresh_par) + 1, name = "Blues"), "red")) +
    theme(legend.position="bottom", 
          legend.box.spacing = unit(3, "pt"))
}
# grm_info_logit()


K <- 25
Bknots <- floor(sqrt(K))
Bdegree <- 3
knot_seq1 <- seq(1, K-1, length.out = Bknots)

B1 <- bs(
  seq(1,(K-1), length.out = 1000),
  knots = knot_seq1,
  degree = Bdegree
)
B1 <- B1[,-ncol(B1)]

set.seed(2341)
unifcoeff <- runif(ncol(B1), 0, 1)

bsplinebasis <- as.data.frame(cbind(B1, B1 %*% unifcoeff)) %>%
  pivot_longer(everything()) %>%
  mutate(k = rep(seq(1,(K-1), length.out = 1000),
                 each = Bknots + 2 + 1),
         name = as.factor(name))
levels(bsplinebasis$name) <- c(paste0("$B_{", 1:(Bknots + 2),",3,t}(k)$"), "$s(k)$")

B2 <- bs(
  1:(K-1),
    knots = knot_seq1,
    degree = Bdegree
  )
  B2 <- B2[,-ncol(B2)]
  
  bsplinebasis2 <- as.data.frame(cbind(B2, B2 %*% unifcoeff)) %>%
    pivot_longer(everything()) %>%
    mutate(k = rep(1:(K-1),
                   each = Bknots + 2 + 1),
           name = as.factor(name))
  levels(bsplinebasis2$name) <- c(paste0("$B_{", 1:(Bknots + 2),",3,t}(k)$"), "$s(k)$")

bspline <- ggplot(bsplinebasis,
       aes(x = k, y = value, col = name)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c(brewer.pal(length(Bknots) + 2 + Bdegree + 1, name = "Blues"), "red")) +
  xlab("$k$") +
  ylab("$B_{m,3,t}(k)$") +
  geom_point(data = bsplinebasis2, aes(x = k, y = value)) +
  theme_bw() +
  theme(legend.position="bottom", 
        legend.box.spacing = unit(0, "pt"),
        legend.title = element_blank())


unifcoeff2 <- sort(unifcoeff)
isplinebasis <- as.data.frame(cbind(B1, B1 %*% unifcoeff2)) %>%
  pivot_longer(everything()) %>%
  mutate(k = rep(seq(1,(K-1), length.out = 1000),
                 each = ncol(B1) +1),
         name = as.factor(name))
levels(isplinebasis$name) <- c(paste0("$B_{", 1:(Bknots + 2),",3,t}(k)$"), "$s(k)$")

isplinebasis2 <- as.data.frame(cbind(B2, B2 %*% unifcoeff2)) %>%
  pivot_longer(everything()) %>%
  mutate(k = rep(1:(K-1),
                 each = ncol(B1) +1),
         name = as.factor(name))
levels(isplinebasis2$name) <- c(paste0("$B_{", 1:(Bknots + 2),",3,t}(k)$"), "$s(k)$")

ispline <- ggplot(isplinebasis,
                   aes(x = k, y = value, col = name)) +
  geom_line(alpha = 0.8, linewidth = 1) +
  scale_color_manual(values = c(brewer.pal(length(Bknots) + 2 + Bdegree + 1, name = "Blues"), "red")) +
  xlab("$k$") +
  ylab("$B_{m,3,t}(k)$") +
  geom_point(data = isplinebasis2, aes(x = k, y = value)) +
  theme_bw() +
  theme(legend.position="bottom", 
        legend.box.spacing = unit(0, "pt"),
        legend.title = element_blank())


K_thresh <- 100
thresh_par <- con_thresh_matrix(K = K_thresh,
                                scale_par = 20,
                                shape_mat = matrix(c(4.5, 2, 0.25, 0.5, 7.5, 2.5,
                                                     0.5, 0.5, 0.25, 0.5, 7.5, 2.5), ncol = 2))
thresh_par <- as.data.frame(thresh_par)
colnames(thresh_par) <- paste0("Item ", 1:6)

thresh_plot <- ggplot(thresh_par %>% 
         mutate(k = 1:(K_thresh-1)) %>% 
         pivot_longer(cols = paste0("Item ", 1:6)),
       aes(x = k, y = value, col = name)) + geom_point(size= 1.6) +
  xlab("$k$") +
  ylab("$\\gamma_{ik}$") +
  geom_line(linewidth = 1, alpha = 0.3) +
  scale_color_manual(values = brewer.pal(6, name = "Blues")) + 
  theme_bw() +
  theme(legend.position="bottom",
        legend.title = element_blank(), 
        legend.box.spacing = unit(0, "pt"))



P <- 10000
abil_par <- qnorm((1:P)/((P + 1)))
thresh_par <- thresh_par

set.seed(1234)
sim_run_disc1 <- simulate_simple1(
  abil_par = abil_par,
  thresh_par = thresh_par,
  disc_par = rep(1, 6))
sim_run_disc1 <- as.data.frame(sim_run_disc1)
colnames(sim_run_disc1) <- paste0("Item ", 1:6)
sim_run_disc1  <- sim_run_disc1 %>% 
  mutate(p = 1:P) %>% 
  pivot_longer(cols = paste0("Item ", 1:6)) %>%
  mutate(disc = as.factor(1))

sim_run_disc3 <- simulate_simple1(
  abil_par = abil_par,
  thresh_par = thresh_par,
  disc_par = rep(3, 6))
sim_run_disc3 <- as.data.frame(sim_run_disc3)
colnames(sim_run_disc3) <- paste0("Item ", 1:6)
sim_run_disc3  <- sim_run_disc3 %>% 
  mutate(p = 1:P) %>% 
  pivot_longer(cols = paste0("Item ", 1:6)) %>%
  mutate(disc = as.factor(3))


ex_sim_hist <- ggplot() + 
  geom_bar(data = sim_run_disc1,
                 aes(x = value, fill = disc), alpha = 0.6) +
  geom_bar(data = sim_run_disc3,
                 aes(x = value, fill = disc), alpha = 0.6) +
  xlab("$k$") +
  # xlab("") +
  facet_wrap(. ~ name, scales = "free_y", ncol = 3) +
  scale_fill_manual(name = "$\\alpha_{i}$", 
                    values = brewer.pal(4, name = "Dark2")) +
  theme_bw() +
  theme(strip.background = element_rect(fill="white"),
        legend.position = "bottom", 
        legend.box.spacing = unit(0, "pt"))



all_results <- data.frame()
  for(p in c(1:length(p_list))) {
    path_p <- paste0("01_sim_final_march2023/res_sim_regular_gh_p", p, ".csv")
    new_results <- readr::read_csv2(path_p)
    new_results$prior <- "global"
    new_results$approach <- "regular"
    new_results$disc <- p_list[[p]]$sim_disc_par$name
    all_results <- rbind(all_results, new_results)
    
    path_p <- paste0("01_sim_final_march2023/res_sim_regular_ih_p", p, ".csv")
    new_results <- readr::read_csv2(path_p)
    new_results$prior <- "item"
    new_results$approach <- "regular"
    new_results$disc <- p_list[[p]]$sim_disc_par$name
    all_results <- rbind(all_results, new_results)
    
    path_p <- paste0("01_sim_final_march2023/res_sim_bspline_gh_p", p, ".csv")
    new_results <- readr::read_csv2(path_p)
    new_results$prior <- "global"
    new_results$approach <- "bsplineknots"
    new_results$disc <- p_list[[p]]$sim_disc_par$name
    all_results <- rbind(all_results, new_results)
    
    path_p <- paste0("01_sim_final_march2023/res_sim_bspline_ih_p", p, ".csv")
    new_results <- readr::read_csv2(path_p)
    new_results$prior <- "item"
    new_results$approach <- "bsplineknots"
    new_results$disc <- p_list[[p]]$sim_disc_par$name
    all_results <- rbind(all_results, new_results)
    
  }


which_models <- c("regular_gh", "regular_ih",
                  "bsplineknots_gh", "bsplineknots_ih")


all_results$fitted_model <- as.factor(all_results$fitted_model)
levels(all_results$fitted_model) <- c("B-Splineknots global", "B-Splineknots item",
                                      "regular GRM global", "regular GRM item")
all_results$K <- as.factor(all_results$K)
levels(all_results$K) <- paste0("$K=", levels(all_results$K),"$")
all_results$P <- as.factor(all_results$P)
levels(all_results$P) <- paste0("$P=", levels(all_results$P),"$")
all_results$disc <- as.factor(all_results$disc)
levels(all_results$disc) <- c("$\\alpha_i=1$", "$\\alpha_i=3$")
all_results$prior <- as.factor(all_results$prior)
levels(all_results$prior) <- c("$\\mu,\\sigma$", "$\\mu_i,\\sigma_i$")
all_results$approach <- as.factor(all_results$approach)
levels(all_results$approach) <- c("B-spline", "regular GRM")


summed_res <- all_results %>% mutate(resi = mean - true,
                                     ci_50_ind = true >= `25%` & true <= `75%`,
                                     ci_50_len = `75%` - `25%`,
                                     ci_90_ind = true >= `5%` & true <= `95%`,
                                     ci_90_len = `95%` - `5%`) %>%
  group_by(variable, approach, prior, K, I, P, disc,par_type,true) %>% 
  summarise(rmse = sqrt(sum(resi^2)),
            mean_cov50 = mean(ci_50_ind),
            mean_cov90 = mean(ci_90_ind),
            mean_sd = mean(sd),
            mean_len50 = mean(ci_50_len),
            mean_len90 = mean(ci_90_len),
            bias = mean(mean - true),
            var = var(mean))

abil_rmse <- ggplot(summed_res %>% filter(par_type == "ability"), 
                    aes(x = true, y = rmse, col = approach)) + 
  geom_point(aes(shape = approach), alpha = 0.7, size = 1.7) + 
  ggh4x::facet_grid2(P + disc ~  K + prior, scales = "free_y") +
  geom_hline(data = summed_res %>% filter(par_type == "ability") %>% 
               group_by(approach, K, I, P, disc,par_type, prior) %>% 
               summarise(armse = mean(rmse)),
             aes(yintercept = armse, col = approach)) +
  xlab("$\\theta_{p}$") +
  ylab("$\\widehat{\\textrm{RMSE}}(\\theta_{p})$") +
  scale_color_manual(values = brewer.pal(4, name = "Dark2")) + 
  theme_bw() +
  theme(strip.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.position = "bottom", 
        legend.box.spacing = unit(0, "pt"))

abil_bias <- ggplot(summed_res %>% filter(par_type == "ability"), 
       aes(x = true, y = bias, col = approach)) + 
  geom_point(aes(shape = approach), alpha = 0.7, size = 1.7) + 
  facet_grid(P + disc ~  K + prior, scales = "free_y") +
  geom_hline(yintercept = 0) +
  xlab("$\\theta_{p}$") +
  ylab("$\\widehat{\\textrm{Bias}}(\\theta_{p})$") +
  scale_color_manual(values = brewer.pal(4, name = "Dark2")) + 
  theme_bw() +
  theme(strip.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.position = "bottom", 
        legend.box.spacing = unit(0, "pt"))


abil_mean_cov50 <- ggplot(summed_res %>% filter(par_type == "ability"),
       aes(x = true, y = mean_cov50, col = approach)) +
  geom_point(aes(shape = approach), alpha = 0.7, size = 1.7) +
  facet_grid(P + disc ~  K + prior) +
  geom_hline(yintercept = 0.5) +
  xlab("$\\theta_{p}$") +
  ylab("$\\widehat{\\textrm{E}}(\\mathbbm{1}_{S_{\\textrm{cred},0.5}}(\\theta_p))$") +
  ylim(0:1) +
  scale_color_manual(values = brewer.pal(4, name = "Dark2")) +
  theme_bw() +
  theme(strip.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.position = "bottom", 
        legend.box.spacing = unit(0, "pt"))



disc_rmse <- ggplot(summed_res %>% filter(par_type == "discrimination") %>% 
                      ungroup() %>% 
                      mutate(w_item = as.factor(rep(1:6, each = length(p_list)*length(which_models)))), 
                    aes(x = w_item, y = rmse, col = approach)) + 
  geom_point(aes(shape = approach), alpha = 0.7, size = 1.7) + 
  ggh4x::facet_grid2(P + disc ~  K + prior, scales = "free_y") +
  geom_hline(data = summed_res %>% filter(par_type == "discrimination") %>% 
               group_by(approach, prior, K, I, P) %>% 
               summarise(armse = mean(rmse)),
             aes(yintercept = armse, col = approach)) +
  xlab("$i$") +
  ylab("$\\widehat{\\textrm{RMSE}}(\\alpha_{i})$") +
  scale_color_manual(values = brewer.pal(4, name = "Dark2")) + 
  theme_bw() +
  theme(strip.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.position = "bottom", 
        legend.box.spacing = unit(0, "pt"))

disc_mean_cov50 <- ggplot(summed_res %>% filter(par_type == "discrimination") %>% 
                            ungroup() %>% mutate(w_item = rep(1:6, each = length(p_list)*length(which_models))),
                          aes(x = w_item, y = mean_cov50, col = approach)) +
  geom_point(aes(shape = approach), alpha = 0.7, size = 1.7) +
  facet_grid(P + disc ~  K + prior, scales = "free") +
  geom_hline(yintercept = 0.5) +
  xlab("$\\alpha_{i}$") +
  ylab("$\\widehat{\\textrm{E}}(\\mathbbm{1}_{S_{\\textrm{cred},0.5}}(\\alpha_i))$") +
  ylim(0:1) +
  scale_color_manual(values = brewer.pal(4, name = "Dark2")) +
  theme_bw() +
  theme(strip.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.position = "bottom", 
        legend.box.spacing = unit(0, "pt"))

disc_bias <- ggplot(summed_res %>% filter(par_type == "discrimination") %>% 
         ungroup() %>% mutate(w_item = rep(1:6, each = length(p_list)*length(which_models))),
       aes(x = w_item, y = bias, col = approach)) +
  geom_point(aes(shape = approach), alpha = 0.7, size = 1.7) +
  facet_grid(P + disc ~  K + prior, scales = "free") +
  geom_hline(yintercept = 0) +
  xlab("$\\alpha_{i}$") +
  ylab("$\\widehat{\\textrm{Bias}}(\\alpha_{i})$") +
  scale_color_manual(values = brewer.pal(4, name = "Dark2")) +
  theme_bw() +
  theme(strip.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.position = "bottom", 
        legend.box.spacing = unit(0, "pt"))

thresh_rmse <- ggplot(summed_res %>% filter(par_type == "threshold") %>%
                        ungroup() %>% mutate(w_item = rep(1:6, 
                          each = length(p_list)/2 * (25-1 + 100-1) * length(which_models))) %>%
                        group_by(approach, K, I, P, disc, w_item, prior) %>%
                        summarise(armse = mean(rmse)), 
                    aes(x = w_item, y = armse, col = approach)) + 
  geom_point(aes(shape = approach), alpha = 0.7, size = 1.7) + 
  facet_grid(P + disc ~  K + prior, scales = "free_y") +
  geom_hline(data = summed_res %>% filter(par_type == "threshold") %>%
               group_by(prior, approach, K, P, disc,par_type) %>%
               summarise(armse = mean(rmse)),
             aes(yintercept = armse, col = approach)) +
  xlab("$i$") +
  ylab("$\\widehat{\\textrm{ARMSE}}(\\boldmath{\\gamma}_{i})$") +
  scale_color_manual(values = brewer.pal(4, name = "Dark2")) + 
  theme_bw() +
  theme(strip.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.position = "bottom", 
        legend.box.spacing = unit(0, "pt"))



abil_armse <- tabular(
  Heading("$P$") * P * Heading("$K$") * K * Heading("$\\alpha$") * disc ~
    Heading() *prior * Heading() * approach *  Heading()*  Heading() * identity * armse,
  data = summed_res %>% filter(par_type == "ability") %>%
    group_by(approach, prior, K, I, P, disc) %>%
    summarise(armse = mean(rmse)) %>%
    group_by(K, P, I, disc) %>%
    mutate(armse = round(armse, digits = 3)) %>%
    mutate(armse = ifelse(row_number() == which.min(armse), paste0(armse, " \\cellcolor{red!25}"), armse))
)

disc_armse <- tabular(
  Heading("$P$") * P * Heading("$K$") * K * Heading("$\\alpha$") * disc ~
    Heading() *prior * Heading() * approach *  Heading()*  Heading() * identity * armse,
  data = summed_res %>% filter(par_type == "discrimination") %>%
    group_by(approach, prior, K, I, P, disc) %>%
    summarise(armse = mean(rmse)) %>%
    group_by(K, P, I, disc) %>%
    mutate(armse = round(armse, digits = 3)) %>%
    mutate(armse = ifelse(row_number() == which.min(armse), paste0(armse, " \\cellcolor{red!25}"), armse))
)

thresh_armse <- tabular(
  Heading("$P$") * P * Heading("$K$") * K * Heading("$\\alpha$") * disc~
    Heading() *prior * Heading() * approach *  Heading()*  Heading() * identity * armse,
  data = summed_res %>% filter(par_type == "threshold") %>%
    group_by(approach, prior, K, I, P, disc) %>%
    summarise(armse = mean(rmse)) %>%
    group_by(K, P, disc) %>%
    mutate(armse = round(armse, digits = 2)) %>%
    mutate(armse = ifelse(row_number() == which.min(armse), paste0(armse, " \\cellcolor{red!25}"), armse))
)


# speeded C-test analysis
ctestdata <- read_sav("S-C-Test_crsm.sav")
names(ctestdata) <- c("group", "sex", paste0("Item ", 1:6))

ctestdata <- ctestdata %>% 
  mutate(group2 = fct_collapse(ctestdata %>% 
                                mutate(group = as.factor(group)) %>% 
                                pull(group), 'early learner' = c("1", "2"),
                               'late learner' = c("3", "4")))


complete_hist <- 
  ggplot(ctestdata %>% 
           pivot_longer(cols = paste0("Item ", 1:6), values_to = "k"), 
         aes(k)) + 
  geom_bar() + 
  facet_grid(group ~ name, scales = "free_y") + 
  theme_bw()+
  xlab("$k$")+
  theme(strip.background = element_rect(fill="white"),
        legend.position = "bottom", 
        legend.box.spacing = unit(0, "pt"))

grouped_hist <- ggplot(ctestdata %>% 
    pivot_longer(cols = paste0("Item ", 1:6)), aes(x = value)) + 
  geom_bar() + facet_grid(group2 ~ name)  + theme_bw()+
  xlab("$k$")+
  theme(strip.background = element_rect(fill="white"),
        legend.position = "bottom", 
        legend.box.spacing = unit(0, "pt"))


sampledmodel_regular2 <- as_cmdstan_fit(paste0("02_cmdstan_ctest_regular_2groups/", 
                                              list.files("02_cmdstan_ctest_regular_2groups")))
sampledmodel_regular4 <- as_cmdstan_fit(paste0("02_cmdstan_ctest_regular_4groups/", 
                                               list.files("02_cmdstan_ctest_regular_4groups")))

sampledmodel_bsplineknots2 <- as_cmdstan_fit(paste0("02_cmdstan_ctest_bspline_2groups/", 
                                                   list.files("02_cmdstan_ctest_bspline_2groups")))
sampledmodel_bsplineknots4 <- as_cmdstan_fit(paste0("02_cmdstan_ctest_bspline_4groups/", 
                                                    list.files("02_cmdstan_ctest_bspline_4groups")))

# sum(sampledmodel_regular2$summary(NULL, c("rhat"))$rhat > 1.1, na.rm = TRUE)
# sum(sampledmodel_regular4$summary(NULL, c("rhat"))$rhat > 1.1, na.rm = TRUE)
# sum(sampledmodel_bsplineknots2$summary(NULL, c("rhat"))$rhat > 1.1, na.rm = TRUE)
# sum(sampledmodel_bsplineknots4$summary(NULL, c("rhat"))$rhat > 1.1, na.rm = TRUE)
# sum(sampledmodel_regular2$summary(NULL, c("ess_bulk"))$ess_bulk < 40, na.rm = TRUE)
# sum(sampledmodel_regular4$summary(NULL, c("ess_bulk"))$ess_bulk < 40, na.rm = TRUE)
# sum(sampledmodel_bsplineknots2$summary(NULL, c("ess_bulk"))$ess_bulk < 40, na.rm = TRUE)
# sum(sampledmodel_bsplineknots4$summary(NULL, c("ess_bulk"))$ess_bulk < 40, na.rm = TRUE)


reg2draws <- sampledmodel_regular2$draws(variables = "Rpred")
reg2postpred_all <- data.frame()
set.seed(123)
for(d in 1:10) {
  reg2postpred_solo <- as.data.frame(matrix(sapply(1:(6 * nrow(ctestdata)), 
                                                  function(x) sample(as.vector(reg2draws[,,x]),1)), 
                                           byrow = FALSE, nrow = nrow(ctestdata)))
  colnames(reg2postpred_solo) <- colnames(as.data.frame(as.matrix(ctestdata[,-c(1, 2,9)] + 1)))
  reg2postpred_solo$type <- paste0("$\\mathbf{Y}^{\\textrm{rep}}_{", d, "}$")
  reg2postpred_solo$group <- ctestdata$group
  reg2postpred_all <- rbind(reg2postpred_all, reg2postpred_solo)
}
reg2postpred_all <- reg2postpred_all %>% pivot_longer(cols = paste0("Item ",1:6))

reg2true <- as.data.frame(as.matrix(ctestdata[,-c(1, 2, 9)] + 1))
reg2true$group <- ctestdata$group
reg2true <- reg2true %>% pivot_longer(cols = paste0("Item ",1:6))

reg2_ppplot_all <- ggplot(reg2postpred_all, aes(x = value)) + 
  geom_bar(alpha = 0.8) + 
  geom_bar(data = reg2true,
                 aes(x = value), fill = "red", alpha = 0.35) + 
  facet_grid(type ~ name) +
  theme_bw()  +
  xlab("$k$") +
  theme(strip.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.position = "bottom", 
        legend.box.spacing = unit(0, "pt"))

reg4draws <- sampledmodel_regular4$draws(variables = "Rpred")
reg4postpred_all <- data.frame()
set.seed(123)
for(d in 1:10) {
  reg4postpred_solo <- as.data.frame(matrix(sapply(1:(6 * nrow(ctestdata)), 
                                                   function(x) sample(as.vector(reg4draws[,,x]),1)), 
                                            byrow = FALSE, nrow = nrow(ctestdata)))
  colnames(reg4postpred_solo) <- colnames(as.data.frame(as.matrix(ctestdata[,-c(1, 2,9)] + 1)))
  reg4postpred_solo$type <- paste0("$\\mathbf{Y}^{\\textrm{rep}}_{", d, "}$")
  reg2postpred_solo$group <- ctestdata$group
  reg4postpred_all <- rbind(reg4postpred_all, reg4postpred_solo)
}
reg4postpred_all <- reg4postpred_all %>% pivot_longer(cols = paste0("Item ",1:6))

reg4true <- as.data.frame(as.matrix(ctestdata[,-c(1, 2, 9)] + 1))
reg4true$group <- ctestdata$group
reg4true <- reg4true %>% pivot_longer(cols = paste0("Item ",1:6))

reg4_ppplot_all <- ggplot(reg4postpred_all, aes(x = value)) + 
  geom_bar(alpha = 0.8) + 
  geom_bar(data = reg4true,
           aes(x = value), fill = "red", alpha = 0.35) + 
  facet_grid(type ~ name) +
  theme_bw()  +
  xlab("$k$") +
  theme(strip.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.position = "bottom", 
        legend.box.spacing = unit(0, "pt"))

bsp2draws <- sampledmodel_bsplineknots2$draws(variables = "Rpred")
bsp2postpred_all <- data.frame()
set.seed(123)
for(d in 1:10) {
  bsp2postpred_solo <- as.data.frame(matrix(sapply(1:(6 * nrow(ctestdata)), 
                                                  function(x) sample(as.vector(bsp2draws[,,x]),1)), 
                                           byrow = FALSE, nrow = nrow(ctestdata)))
  colnames(bsp2postpred_solo) <- colnames(as.data.frame(as.matrix(ctestdata[,-c(1, 2,9)] + 1)))
  bsp2postpred_solo$type <- paste0("$\\mathbf{Y}^{\\textrm{rep}}_{", d, "}$")
  bsp2postpred_solo$group <- ctestdata$group
  bsp2postpred_all <- rbind(bsp2postpred_all, bsp2postpred_solo)
}
bsp2postpred_all <- bsp2postpred_all %>% pivot_longer(cols = paste0("Item ",1:6))

bsp2true <- as.data.frame(as.matrix(ctestdata[,-c(1, 2,9)] + 1))
bsp2true$group <- ctestdata$group
bsp2true <- bsp2true %>% pivot_longer(cols = paste0("Item ",1:6))

bsp2_ppplot_all <- ggplot(bsp2postpred_all, aes(x = value)) + 
  geom_bar(alpha = 0.8) + 
  geom_bar(data = bsp2true,
           aes(x = value), fill = "red", alpha = 0.35) + 
  facet_grid(type ~ name) +
  theme_bw() +
  xlab("$k$") + 
  theme(strip.background = element_rect(fill="white"),
        legend.position = "bottom", 
        legend.box.spacing = unit(0, "pt"))

bsp4draws <- sampledmodel_bsplineknots4$draws(variables = "Rpred")
bsp4postpred_all <- data.frame()
set.seed(123)
for(d in 1:10) {
  bsp4postpred_solo <- as.data.frame(matrix(sapply(1:(6 * nrow(ctestdata)), 
                                                   function(x) sample(as.vector(bsp4draws[,,x]),1)), 
                                            byrow = FALSE, nrow = nrow(ctestdata)))
  colnames(bsp4postpred_solo) <- colnames(as.data.frame(as.matrix(ctestdata[,-c(1, 2,9)] + 1)))
  bsp4postpred_solo$type <- paste0("$\\mathbf{Y}^{\\textrm{rep}}_{", d, "}$")
  bsp4postpred_solo$group <- ctestdata$group
  bsp4postpred_all <- rbind(bsp4postpred_all, bsp4postpred_solo)
}
bsp4postpred_all <- bsp4postpred_all %>% pivot_longer(cols = paste0("Item ",1:6))

bsp4true <- as.data.frame(as.matrix(ctestdata[,-c(1, 2,9)] + 1))
bsp4true$group <- ctestdata$group
bsp4true <- bsp4true %>% pivot_longer(cols = paste0("Item ",1:6))

bsp4_ppplot_all <- ggplot(bsp4postpred_all, aes(x = value)) + 
  geom_bar(alpha = 0.8) + 
  geom_bar(data = bsp4true,
           aes(x = value), fill = "red", alpha = 0.35) + 
  facet_grid(type ~ name) +
  theme_bw() +
  xlab("$k$") + 
  theme(strip.background = element_rect(fill="white"),
        legend.position = "bottom", 
        legend.box.spacing = unit(0, "pt"))

pppvar <- c("meanRI","sdRI", "q25RI", "q75RI", "skewRI", "kurtRI")
pppnames <- c("mean", "sd", "$25\\%-\\textrm{quantile}$", "$75\\%-\\textrm{quantile}$", "skewness", "kurtosis")
pppval <- rbind(
  sampledmodel_regular2$summary(variables = pppvar) %>%
    select(variable, mean) %>% mutate(type = "regular",
                                      item = rep(paste0(1:6), 6),
                                      groups = 2,
                                      stat = rep(pppnames, each = 6)),
  sampledmodel_regular4$summary(variables = pppvar) %>%
    select(variable, mean) %>% mutate(type = "regular",
                                      item = rep(paste0(1:6), 6),
                                      groups = 4,
                                      stat = rep(pppnames, each = 6)),
  sampledmodel_bsplineknots2$summary(variables = pppvar) %>%
    select(variable, mean) %>% mutate(type = "B-Spline",
                                      item = rep(paste0(1:6), 6),
                                      groups = 2,
                                      stat = rep(pppnames, each = 6)),
  sampledmodel_bsplineknots4$summary(variables = pppvar) %>%
    select(variable, mean) %>% mutate(type = "B-Spline",
                                      item = rep(paste0(1:6), 6),
                                      groups = 4,
                                      stat = rep(pppnames, each = 6))
) %>% mutate(variable = as.factor(variable),
             type = as.factor(type),
             item = as.factor(item),
             groups = as.factor(groups),
             stat = as.factor(stat))

ppp_stats <- ggplot(pppval, aes(x = item, y = mean, col = type, shape = type)) + 
  geom_point() + 
  facet_grid(groups ~ stat) + 
  geom_hline(yintercept = c(0.05, 0.95)) +  theme_bw() +
  scale_color_manual(values = brewer.pal(3, name = "Dark2")) +
  ylim(c(0, 1)) +
  xlab("$i$") +
  ylab("$P(T(Y_{i}^{\\textrm{rep}}) \\geq T(Y_{i}))$") +
  theme(strip.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.position = "bottom", 
        legend.box.spacing = unit(0, "pt"))


emp_rel <- tabular(Heading() * identity * Heading() * emp_rel ~ Heading() * fitted_model *  Heading() * groups,
data = rbind(
  sampledmodel_regular2$summary(variables = c("theta"), "mean", "sd") %>% 
    mutate(se2 = sd^2,
           resi2 = (mean - mean(mean))^2) %>%
    summarise(emp_rel = 1 - mean(se2) / mean(resi2)) %>% 
    mutate(fitted_model = "regular GRM",
           groups = "2 groups"),
  sampledmodel_bsplineknots2$summary(variables = c("theta"), "mean", "sd") %>% 
    mutate(se2 = sd^2,
           resi2 = (mean - mean(mean))^2) %>%
    summarise(emp_rel = 1 - mean(se2) / mean(resi2)) %>% 
    mutate(fitted_model = "B-Spline",
           groups = "2 groups"),
  sampledmodel_regular4$summary(variables = c("theta"), "mean", "sd") %>% 
    mutate(se2 = sd^2,
           resi2 = (mean - mean(mean))^2) %>%
    summarise(emp_rel = 1 - mean(se2) / mean(resi2)) %>% 
    mutate(fitted_model = "regular GRM",
           groups = "4 groups"),
  sampledmodel_bsplineknots4$summary(variables = c("theta"), "mean", "sd") %>% 
    mutate(se2 = sd^2,
           resi2 = (mean - mean(mean))^2) %>%
    summarise(emp_rel = 1 - mean(se2) / mean(resi2)) %>% 
    mutate(fitted_model = "B-Spline",
           groups = "4 groups")
) %>% mutate(fitted_model = as.factor(fitted_model),
             groups = as.factor(groups)))


thetasum <- rbind(
  sampledmodel_regular2$summary(variables = "theta",
                                ~quantile(.x, probs = c(0.05, 0.25, 0.75, 0.95)), "mean", "sd") %>% 
    mutate(type = "regular GRM",
           examinee = 1:nrow(ctestdata),
           group = ctestdata$group,
           groups = "2 groups"),
  sampledmodel_bsplineknots2$summary(variables = "theta",
                                     ~quantile(.x, probs = c(0.05, 0.25, 0.75, 0.95)), "mean", "sd") %>% 
    mutate(type = "B-Spline",
           examinee = 1:nrow(ctestdata),
           group = ctestdata$group,
           groups = "2 groups"),
  sampledmodel_regular4$summary(variables = "theta",
                                ~quantile(.x, probs = c(0.05, 0.25, 0.75, 0.95)), "mean", "sd") %>% 
    mutate(type = "regular GRM",
           examinee = 1:nrow(ctestdata),
           group = ctestdata$group,
           groups = "4 groups"),
  sampledmodel_bsplineknots4$summary(variables = "theta",
                                     ~quantile(.x, probs = c(0.05, 0.25, 0.75, 0.95)), "mean", "sd") %>% 
    mutate(type = "B-Spline",
           examinee = 1:nrow(ctestdata),
           group = ctestdata$group,
           groups = "4 groups")
)


ctest_theta <- ggplot(thetasum, aes(x = examinee, y = mean, col = type)) + 
  geom_point(alpha = 0.7, size = 2) + 
  geom_errorbar(aes(ymin = `5%`, ymax = `95%`),
  # geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), 
                position = position_dodge(width = 1.5)) + 
  scale_color_manual(values = brewer.pal(3, name = "Dark2")) +
  xlab("$p$") +
  facet_grid(groups ~ group, scales = "free_x") +
  ylab("$\\theta_{p}$") +
  theme_bw() +
  theme(strip.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.position = "bottom", 
        legend.box.spacing = unit(0, "pt"))

theta2_sum <- rbind(
  sampledmodel_regular2$summary(variables = c("mu_theta", "sigma_theta"), c("mean", "sd")) %>% 
    mutate(type = "regular"),
  sampledmodel_bsplineknots2$summary(variables = c("mu_theta", "sigma_theta"), c("mean", "sd")) %>% 
    mutate(type = "bsplineknots") 
) %>% mutate(variable = as.factor(variable),
             type = as.factor(type)) %>% 
  pivot_longer(col = c("mean", "sd")) %>% mutate(name = as.factor(name))
levels(theta2_sum$variable) <- c("$\\mu_1$", "$\\sigma_1$")
levels(theta2_sum$type) <- c("B-spline", "regular GRM")

theta_hyper2 <- tabular(Heading() * variable * Heading() * identity * Heading() * value ~ 
                          Heading() * type * Heading() * name,
                        data = theta2_sum)

theta4_sum <- rbind(
  sampledmodel_regular4$summary(variables = c("mu_theta", "sigma_theta"), c("mean", "sd")) %>% 
    mutate(type = "regular"),
  sampledmodel_bsplineknots4$summary(variables = c("mu_theta", "sigma_theta"), c("mean", "sd")) %>% 
    mutate(type = "bsplineknots") 
) %>% mutate(variable = as.factor(variable),
             type = as.factor(type)) %>% 
  pivot_longer(col = c("mean", "sd")) %>% mutate(name = as.factor(name))
levels(theta4_sum$variable) <- c("$\\mu_1$", "$\\mu_2$", "$\\mu_3$",
                                 "$\\sigma_1$", "$\\sigma_2$", "$\\sigma_3$")
levels(theta4_sum$type) <- c("B-spline", "regular GRM")

theta_hyper4 <- tabular(Heading() *variable* Heading() * identity * Heading() * value ~ 
                          Heading() * type *Heading() * name,
                        data = theta4_sum)


alphasum <- rbind(
  sampledmodel_regular2$summary(variables = "alpha",
                               ~quantile(.x, probs = c(0.05, 0.25, 0.75, 0.95)), "mean", "sd") %>% 
    mutate(type = "regular", item = 1:6,
           groups = "2 groups"),
  sampledmodel_bsplineknots2$summary(variables = "alpha",
                                    ~quantile(.x, probs = c(0.05, 0.25, 0.75, 0.95)), "mean", "sd") %>% 
    mutate(type = "B-Spline", item = 1:6,
           groups = "2 groups"),
  sampledmodel_regular4$summary(variables = "alpha",
                               ~quantile(.x, probs = c(0.05, 0.25, 0.75, 0.95)), "mean", "sd") %>% 
    mutate(type = "regular", item = 1:6,
           groups = "4 groups"),
  sampledmodel_bsplineknots4$summary(variables = "alpha",
                                    ~quantile(.x, probs = c(0.05, 0.25, 0.75, 0.95)), "mean", "sd") %>% 
    mutate(type = "B-Spline", item = 1:6,
           groups = "4 groups")
)


ctest_alpha <- ggplot(alphasum, aes(x = item, y = mean, col = type)) + 
  geom_point(alpha = 0.7, size = 2) + 
  geom_errorbar(aes(ymin = `5%`, ymax = `95%`), 
                position = position_dodge(width = 0.05)) + 
  scale_color_manual(values = brewer.pal(3, name = "Dark2")) +
  xlab("$i$") +
  ylab("$\\alpha_{i}$") +
  facet_wrap(. ~ groups) +
  theme_bw() +
  theme(strip.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.position = "bottom", 
        legend.box.spacing = unit(0, "pt"))

gammasum <- rbind(
  sampledmodel_regular2$summary(variables = "gamma",
                                ~quantile(.x, probs = c(0.05, 0.25, 0.75, 0.95)), "mean", "sd") %>% 
    mutate(type = "regular", 
           item = paste0("Item ", rep(1:6, 25)),
           k = rep(1:25, each = 6),
           groups = "2 groups"),
  sampledmodel_bsplineknots2$summary(variables = "gamma",
                                     ~quantile(.x, probs = c(0.05, 0.25, 0.75, 0.95)), "mean", "sd") %>% 
    mutate(type = "B-Spline",
           item = paste0("Item ", rep(1:6, 25)),
           k = rep(1:25, each = 6),
           groups = "2 groups"),
  sampledmodel_regular4$summary(variables = "gamma",
                                ~quantile(.x, probs = c(0.05, 0.25, 0.75, 0.95)), "mean", "sd") %>% 
    mutate(type = "regular", 
           item = paste0("Item ", rep(1:6, 25)),
           k = rep(1:25, each = 6),
           groups = "4 groups"),
  sampledmodel_bsplineknots4$summary(variables = "gamma",
                                     ~quantile(.x, probs = c(0.05, 0.25, 0.75, 0.95)), "mean", "sd") %>% 
    mutate(type = "B-Spline",
           item = paste0("Item ", rep(1:6, 25)),
           k = rep(1:25, each = 6),
           groups = "4 groups")
)

ctest_gamma <- ggplot(gammasum, aes(x = k, y = mean, col = type)) + 
  geom_point(alpha = 0.7, size = 1.3) + 
  geom_errorbar(aes(ymin = `5%`, ymax = `95%`), 
                position = position_dodge(width = 0.2)) + 
  scale_color_manual(values = brewer.pal(3, name = "Dark2")) +
  xlab("$k$") +
  ylab("$\\gamma_{ik}$") +
  facet_grid(groups ~ item) +  theme_bw() +
  theme(strip.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.position = "bottom", 
        legend.box.spacing = unit(0, "pt"))


loo_regular2 <- sampledmodel_regular2$loo()
loo_bsplineknots2 <- sampledmodel_bsplineknots2$loo()
loo_regular4 <- sampledmodel_regular4$loo()
loo_bsplineknots4 <- sampledmodel_bsplineknots4$loo()
loo_comp <- loo::loo_compare(list("regular GRM (2 groups)" = loo_regular2,
                                  "B-spline (2 groups)" = loo_bsplineknots2,
                                  "regular GRM (4 groups)" = loo_regular4,
                                  "B-spline (4 groups)" = loo_bsplineknots4))

ctest_loo <- tabular(Heading() * groups ~ 
          Heading() * type *Heading() * identity * Heading() * elpd_loo_diff,
        data = data.frame(elpd_loo_diff = loo_comp[,1]) %>% 
          mutate(groups = as.factor(rep(c("4 groups", "2 groups"), 2)),
                 type = as.factor(rep(c("B-spline", "regular GRM"), each = 2))))




# reading motivation analysis
readdat <- readr::read_csv2("rohdaten_lesemotivation.csv")
readdat <- readdat %>% drop_na()

weird_grouplab <- (readdat[,2] == "treat") == (readdat[,13] == "treat")
readdat <- readdat[weird_grouplab,-c(1,2)]

readdat2 <- as.data.frame(rbind(as.matrix(readdat[,1:10] + 1),
                  as.matrix(readdat[,12:21] + 1)))
colnames(readdat2) <- paste0("Item ", 1:10)



reading_grm <- as_cmdstan_fit(paste0("03_cmdstan_readingmotivation_regular/", 
                                               list.files("03_cmdstan_readingmotivation_regular")))
reading_bsplineknots <- as_cmdstan_fit(paste0("03_cmdstan_readingmotivation_bspline/", 
                                               list.files("03_cmdstan_readingmotivation_bspline")))


# sum(reading_grm$summary(NULL, c("rhat"))$rhat > 1.1, na.rm = TRUE)
# sum(reading_bsplineknots$summary(NULL, c("rhat"))$rhat > 1.1, na.rm = TRUE)
# sum(reading_grm$summary(NULL, c("ess_bulk"))$ess_bulk < 40, na.rm = TRUE)
# sum(reading_bsplineknots$summary(NULL, c("ess_bulk"))$ess_bulk < 40, na.rm = TRUE)



pppvar <- c("meanRI","sdRI", "q25RI", "q75RI", "skewRI", "kurtRI")
pppnames <- c("mean", "sd", "$25\\%-\\textrm{quantile}$", "$75\\%-\\textrm{quantile}$", "skewness", "kurtosis")
pppval2 <- rbind(
  reading_grm$summary(variables = pppvar) %>%
    select(variable, mean) %>% mutate(type = "regular GRM",
                                      item = rep(paste0(1:10), 6),
                                      stat = rep(pppnames, each = 10)),
  reading_bsplineknots$summary(variables = pppvar) %>%
    select(variable, mean) %>% mutate(type = "B-spline",
                                      item = rep(paste0(1:10), 6),
                                      stat = rep(pppnames, each = 10))
) %>% mutate(variable = as.factor(variable),
             type = as.factor(type),
             item = factor(item, levels = paste0(1:10)),
             stat = as.factor(stat))


ppp_stats2 <- ggplot(pppval2, aes(x = item, y = mean, col = type, shape = type)) + 
  geom_point() + 
  facet_grid(. ~ stat) + 
  geom_hline(yintercept = c(0.05, 0.95)) +  theme_bw() +
  scale_color_manual(values = brewer.pal(3, name = "Dark2")) +
  ylim(c(0, 1)) +
  xlab("$i$") +
  ylab("$P(T(Y_{i}^{\\textrm{rep}}) \\geq T(Y_{i}))$") +
  theme(strip.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.position = "bottom", 
        legend.box.spacing = unit(0, "pt"))

abilsum <- rbind(
  reading_grm$summary(variables = "theta", ~quantile(.x, probs = c(0.05, 0.25, 0.75, 0.95)), "mean", "sd") %>% 
    mutate(type = "regular GRM",
           examinee = rep(1:(nrow(readdat)), 2),
           group = rep(readdat$actual_group, 2),
           pre_post = rep(c("pre", "post"), each = nrow(readdat))),
  reading_bsplineknots$summary(variables = "theta", ~quantile(.x, probs = c(0.05, 0.25, 0.75, 0.95)), "mean", "sd")%>% 
    mutate(type = "B-spline",
           examinee = rep(1:(nrow(readdat)), 2),
           group = rep(readdat$actual_group, 2),
           pre_post = rep(c("pre", "post"), each = nrow(readdat)))) %>%
  mutate(pre_post = factor(pre_post, levels = c("pre", "post")),
         group = as.factor(group))
levels(abilsum$group) <- c("control", "treatment")

reading_theta <- ggplot(abilsum, aes(x = examinee, y = mean, col = type)) + 
  geom_errorbar(aes(ymin = `5%`, ymax = `95%`), 
                position = position_dodge(width = 0.2)) + 
  geom_point(alpha = 0.7, size = 2) + theme_bw()+
  facet_grid(group ~ pre_post) +
  xlab("$p$") +
  ylab("$\\theta_p$") +
  scale_color_manual(values = brewer.pal(3, name = "Dark2")) +
  theme_bw() +
  theme(strip.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.position = "bottom", 
        legend.box.spacing = unit(0, "pt"))


theta_str <- rbind(reading_grm$summary(variables = c("mu_all","mu_theta_treat"), c("mean", "sd")),
reading_bsplineknots$summary(variables = c("mu_all","mu_theta_treat"), c("mean", "sd"))) %>% 
  mutate(type = as.factor(rep(c("regular GRM", "B-spline"), each = 2)),
         variable = as.factor(rep(c("$\\mu_{\\textrm{all}}$", "$\\mu_{t}$"), 2))) %>%
          pivot_longer(cols = c("mean", "sd")) %>%
         mutate(name = as.factor(name))
theta_str2 <- tabular(Heading() *variable* Heading() * identity * Heading() * value ~ 
          Heading() * type *Heading() * name,data = theta_str)


theta_strb <- rbind(reading_grm$summary(variables = c("Corr_treat[2,1]", "Cov_treat[2,2]",
                                  "Corr_control[2,1]", "Cov_control[2,2]"), 
                         c("mean", "sd")) %>% 
        mutate(type = as.factor("regular GRM"),
               variable = as.factor(c("$\\rho_{t}$", "$\\sigma_{t}^2$", "$\\rho_{c}$", "$\\sigma_{c}^2$"))),
reading_bsplineknots$summary(variables = c("Corr_treat[2,1]", "Cov_treat[2,2]",
                                  "Corr_control[2,1]", "Cov_control[2,2]"), 
                    c("mean", "sd")) %>%
       mutate(type = as.factor("B-spline"),
              variable = as.factor(c("$\\rho_{t}$", "$\\sigma_{t}^2$", "$\\rho_{c}$", "$\\sigma_{c}^2$")))) %>%
  pivot_longer(cols = c("mean", "sd")) %>%
  mutate(name = as.factor(name))

theta_str3 <- tabular(Heading() *variable* Heading() * identity * Heading() * value ~ 
                        Heading() * type *Heading() * name,data = theta_strb)


discsum <- rbind(
  reading_grm$summary(variables = "alpha", ~quantile(.x, probs = c(0.05, 0.25, 0.75, 0.95)), "mean", "sd") %>% 
    mutate(type = "regular GRM", item = 1:10),
  reading_bsplineknots$summary(variables = "alpha", ~quantile(.x, probs = c(0.05, 0.25, 0.75, 0.95)), "mean", "sd") %>% 
    mutate(type = "B-spline", item = 1:10)) %>%
  mutate(item = factor(item, levels = paste0(1:10)))


reading_alpha <- ggplot(discsum, aes(x = item, y = mean, col = type)) + 
  geom_errorbar(aes(ymin = `5%`, ymax = `95%`), 
                position = position_dodge(width = 0.1)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = brewer.pal(3, name = "Dark2")) +
  xlab("$i$") +
  ylab("$\\alpha_i$") +
  theme_bw() +
  theme(strip.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.position = "bottom", 
        legend.box.spacing = unit(0, "pt"))

gammasum <- rbind(
  reading_grm$summary(variables = "gamma", ~quantile(.x, probs = c(0.05, 0.95)), "mean", "sd") %>% 
    mutate(type = "regular GRM", item = rep(paste0("Item ", 1:10), 100), k = rep(1:100, each = 10)),
  reading_bsplineknots$summary(variables = "gamma", ~quantile(.x, probs = c(0.05, 0.95)), "mean", "sd") %>% 
    mutate(type = "B-spline", item = rep(paste0("Item ", 1:10), 100), k = rep(1:100, each = 10))
) %>% mutate(item = factor(item, levels = paste0("Item ", 1:10)))

reading_gamma <- ggplot(gammasum, aes(x = k, y = mean, col = type)) + 
  geom_point(alpha = 0.7, size = 1.3) + 
  facet_wrap(. ~ item, scales = "free_y")+
  geom_errorbar(aes(ymin = `5%`, ymax = `95%`), 
                position = position_dodge(width = 0.2)) + 
  xlab("$k$") +
  ylab("$\\gamma_{ik}$") + 
  scale_color_manual(values = brewer.pal(3, name = "Dark2")) +
  theme_bw() +
  theme(strip.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.position = "bottom", 
        legend.box.spacing = unit(0, "pt"))
