source("00_basicfunctions.R")

library(cmdstanr) # use sample() or sample_mpi() for simulations on cluster
library(splines)
library(rlecuyer)

model_spec <- "regular_gh"
p <- as.integer(Sys.getenv("PBS_ARRAYID"))
print(p)
num_cores <- 4
nsim <- 50


.lec.SetPackageSeed(c(6,13,73,4,52,1)) # rlecuyer equivalent for set.seed()
nstreams <- length(p_list) # number of rng streams
names <- paste0("stan_sim_stream",1:nstreams) # names for rng streams
.lec.CreateStream(names) # create rng streams
.lec.CurrentStream(names[p]) # select the p-th stream


modelpath <- paste0("01_sim_", model_spec, ".stan")

exe_path <- paste0("cmdstan_exe_p", p,"_",model_spec)
dir.create(exe_path)

cmdmodel <- cmdstan_model(stan_file = modelpath, compile = FALSE)
cmdmodel$compile(dir = exe_path)

compiled_model <- cmdstan_model(exe_file =  paste0(exe_path, "/01_sim_", model_spec))


all_iters <- data.frame()

RNGkind("L'Ecuyer-CMRG")
for(i in 1:nsim) {
  thresh_par <- p_list[[p]]$sim_thresh_par
  abil_par <- qnorm((1:p_list[[p]]$sim_P)/((p_list[[p]]$sim_P) + 1))
  
  sim_run <- simulate_simple1(abil_par = abil_par,
                              thresh_par = thresh_par,
                              disc_par = p_list[[p]]$sim_disc_par$value)
  
  
  data_list <- list(
    I = p_list[[p]]$sim_I,
    P = p_list[[p]]$sim_P,
    K = p_list[[p]]$sim_K,
    R = sim_run)
  
  sampledmodel <- compiled_model$sample(
    data = data_list,
    chains = 4,
    parallel_chains = 4,
    refresh = 100
  )
  
  iter_sum <- sampledmodel$summary(variables = c("theta", "alpha", "gamma"),
                       "mean", "median", "sd", "mad", 
                       ~quantile(.x, probs = c(0.05, 0.25, 0.75, 0.95)),
                       "rhat", "ess_bulk", "ess_tail")
  
  iter_sum$true <- c(
    abil_par,
    p_list[[p]]$sim_disc_par$value,
    as.vector(t(thresh_par)))
  
  iter_sum$par_type <- c(rep("ability", p_list[[p]]$sim_P),
                       rep("discrimination", p_list[[p]]$sim_I),
                       rep("threshold", p_list[[p]]$sim_I * (p_list[[p]]$sim_K - 1)))
  
  iter_sum$fitted_model <- model_spec
  iter_sum$iteration <- i
  iter_sum$I <- p_list[[p]]$sim_I
  iter_sum$P <- p_list[[p]]$sim_P
  iter_sum$K <- p_list[[p]]$sim_K

  iter_sum$no_div <- sum(sampledmodel$diagnostic_summary()$num_divergent)
  iter_sum$no_max_tree <- sum(sampledmodel$diagnostic_summary()$num_max_treedepth)
  
  all_iters <- rbind(all_iters, iter_sum)
}
write.csv2(all_iters, file = paste0("res_sim_", model_spec, "_p", p,".csv"))
unlink(exe_path, recursive = TRUE)

gc()