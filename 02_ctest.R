## run Stan models and save output for the c-test analysis

source("00_basicfunctions.R")

ctestdata <- read_sav("S-C-Test_crsm.sav")
names(ctestdata) <- c("group", "sex", paste0("Item ", 1:6))

## create different aggregation of groups
ctestdata <- ctestdata %>% 
  mutate(group2 = fct_collapse(ctestdata %>% 
                                mutate(group = as.factor(group)) %>% 
                                pull(group), nativespeaker = c("1", "2"),
                              secondlanguage = c("3", "4")))


compiled_model_regular <- cmdstan_model(stan_file = "02_ctest_regular.stan")
compiled_model_bsplineknots <- cmdstan_model(stan_file = "02_ctest_bspline.stan")

sampledmodel_regular_4groups <- compiled_model_regular$sample(
  data = list(
    I = 6,
    P = nrow(ctestdata),
    K = 26,
    R = as.matrix(ctestdata[,-c(1, 2, 9)] + 1),
    no_group = max(as.double(ctestdata$group)),
    group = ctestdata$group),
  chains = 4,
  parallel_chains = 4,
  refresh = 100,
  seedygonzales
)
dir.create("02_cmdstan_ctest_regular_4groups")
sampledmodel_regular_4groups$save_output_files(dir = "02_cmdstan_ctest_regular_4groups")

sampledmodel_regular_2groups <- compiled_model_regular$sample(
  data = list(
    I = 6,
    P = nrow(ctestdata),
    K = 26,
    R = as.matrix(ctestdata[,-c(1, 2, 9)] + 1),
    no_group = max(as.double(ctestdata$group2)),
    group = as.double(ctestdata$group2)),
  chains = 4,
  parallel_chains = 4,
  refresh = 100,
  seedygonzales
)
dir.create("02_cmdstan_ctest_regular_2groups")
sampledmodel_regular_2groups$save_output_files(dir = "02_cmdstan_ctest_regular_2groups")

# calculate the item-specific B-spline matrix based on smoothed sample 
# quantiles of the response pattern
Bknots <- floor(sqrt(26)) # number of inner knots
Bdegree <- 3 # degree of B-spline basis
# include extremes of the scale into the inner knot sequence
alt_knot <- unname(rbind(rep(1, 6),
                         apply(as.matrix(ctestdata[,-c(1, 2, 9)] + 1), 2,
                               function(x) {
                                 ydens <- density(x, adjust = 1, from = 0, to = 26)
                                 ydens$x[sapply((1:Bknots) / (Bknots + 2), function(q) which(cumsum(ydens$y) / sum(ydens$y) >= q)[1])]
                               }),
                         rep(26-1, 6)))

all_knots <- array(dim = c(6, 26-1,length(c(0, (1:Bknots) / (Bknots + 1), 1)) + Bdegree - 1))
for(a in 1:ncol(alt_knot)) {
  B_all <- bs(
    1:(26-1),
    knots = alt_knot[,a],
    degree = Bdegree
  )
  all_knots[a,,] <- B_all[,-ncol(B_all)]
}

# run B-spline version for the GRM model in Stan for 4 groups approach
sampledmodel_bsplineknots_4groups <- compiled_model_bsplineknots$sample(
  data = list(
    I = 6,
    P = nrow(ctestdata),
    K = 26,
    R = as.matrix(ctestdata[,-c(1, 2, 9)] + 1),
    M = dim(all_knots)[3],
    B = all_knots,
    no_group = max(as.double(ctestdata$group)),
    group = ctestdata$group),
  chains = 4,
  parallel_chains = 4,
  refresh = 100,
  seedygonzales
)
dir.create("02_cmdstan_ctest_bspline_4groups")
sampledmodel_bsplineknots_4groups$save_output_files(dir = "02_cmdstan_ctest_bspline_4groups")

# run B-spline version for the GRM model in Stan for 2 groups approach
sampledmodel_bsplineknots_2groups <- compiled_model_bsplineknots$sample(
  data = list(
    I = 6,
    P = nrow(ctestdata),
    K = 26,
    R = as.matrix(ctestdata[,-c(1, 2, 9)] + 1),
    M = dim(all_knots)[3],
    B = all_knots,
    no_group = max(as.double(ctestdata$group2)),
    group = as.double(ctestdata$group2)),
  chains = 4,
  parallel_chains = 4,
  refresh = 100,
  seedygonzales
)
dir.create("02_cmdstan_ctest_bspline_2groups")
sampledmodel_bsplineknots_2groups$save_output_files(dir = "02_cmdstan_ctest_bspline_2groups")
