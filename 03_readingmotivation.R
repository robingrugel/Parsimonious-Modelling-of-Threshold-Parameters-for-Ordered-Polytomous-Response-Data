## run Stan models and save output for the reading motivation analysis

source("00_basicfunctions.R")

readdat <- readr::read_csv2("rohdaten_lesemotivation.csv")
readdat <- readdat %>% drop_na()

# adress labelling errors in the raw data files
weird_grouplab <- (readdat[,2] == "treat") == (readdat[,13] == "treat")
readdat <- readdat[weird_grouplab,-c(1,2)]

# run regular GRM model in Stan
reading_grm_stan <- cmdstan_model(stan_file = "03_readingmotivation_regular.stan")
reading_grm_samp <- reading_grm_stan$sample(
  data = list(
    I = 10,
    P = nrow(readdat),
    K = 101,
    R = rbind(as.matrix(readdat[,1:10] + 1),
              as.matrix(readdat[,12:21] + 1)),
    group = as.double(as.factor(readdat$actual_group))),
  chains = 4,
  parallel_chains = 4,
  seed = seedygonzales
)
dir.create("03_cmdstan_readingmotivation_regular")
reading_grm_samp$save_output_files(dir = "03_cmdstan_readingmotivation_regular")


# calculate the item-specific B-spline matrix based on smoothed sample 
# quantiles of the response pattern
Bknots <- floor(sqrt(101)) # number of inner knots
Bdegree <- 3 # degree of B-spline basis
# include extremes of the scale into the inner knot sequence
alt_knot <- unname(rbind(rep(1, 10),
                         apply(rbind(as.matrix(readdat[,1:10] + 1),
                                     as.matrix(readdat[,12:21] + 1)), 2,
                               function(x) {
                                 ydens <- density(x, adjust = 1, from = 0, to = 101)
                                 ydens$x[sapply((1:Bknots) / (Bknots + 2), function(q) which(cumsum(ydens$y) / sum(ydens$y) >= q)[1])]
                               }),
                         rep(101-1, 10)))

all_knots <- array(dim = c(10, 101-1,length(c(0, (1:Bknots) / (Bknots + 1), 1)) + Bdegree - 1))
for(a in 1:ncol(alt_knot)) {
  B_all <- bs(
    1:(101-1),
    knots = alt_knot[,a],
    degree = Bdegree
  )
  all_knots[a,,] <- B_all[,-ncol(B_all)]
}

# run B-spline version for the GRM model in Stan
reading_bspline_stan <- cmdstan_model(stan_file = "03_readingmotivation_bspline.stan")
reading_bspline_samp <- reading_bspline_stan$sample(
  data = list(
    I = 10,
    P = nrow(readdat),
    K = 101,
    R = rbind(as.matrix(readdat[,1:10] + 1),
              as.matrix(readdat[,12:21] + 1)),
    group = as.double(as.factor(readdat$actual_group)),
    M = dim(all_knots)[3],
    B = all_knots),
  chains = 4,
  parallel_chains = 4,
  seed = seedygonzales
)
dir.create("03_cmdstan_readingmotivation_bspline")
reading_bspline_samp$save_output_files(dir = "03_cmdstan_readingmotivation_bspline")
