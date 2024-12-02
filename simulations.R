library(slendr)
init_env(quiet = TRUE)

library(parallel)
library(ggplot2)
library(dplyr)

NE_START <- c(6000, 2000, 500)
NE_END <- 200

GENERATION_TIME <- 32
SEQUENCE_LENGTH <- 10e6
RECOMBINATION_RATE <- 1e-8
MUTATION_RATE <- 1.65e-8

T_START <- 1
T_BOTTLE <- 5000
T_END <- T_BOTTLE + 50000

REPS <- 10

create_model <- function(Ne) {
  pop <- population("Beluga", time = T_START, N = Ne) %>%
    resize(time = T_BOTTLE, N = NE_END, how = "step")

  model <- compile_model(
    pop, generation_time = GENERATION_TIME,
    simulation_length = T_START + T_END,
    serialize = FALSE
  )

  times <- seq(500, T_END, by = 500)
  samples <- schedule_sampling(model, times = times, list(pop, 25))

  list(model = model, samples = samples)
}

simulate_ts <- function(model, samples, engine_fun = msprime) {
  ts <- engine_fun(model, sequence_length = SEQUENCE_LENGTH, recombination_rate = RECOMBINATION_RATE, samples = samples)

  if (identical(engine_fun, slim)) {
    Ne_start <- extract_parameters(model)$splits[, "N"]
    ts <- ts_recapitate(ts, Ne = Ne_start, recombination_rate = RECOMBINATION_RATE) %>%
      ts_simplify()
  }

  ts <- ts_mutate(ts, mutation_rate = MUTATION_RATE)

  ts
}

compute_pi <- function(ts) {
  result <- ts_samples(ts)
  result$pi <- ts_diversity(ts, result$name)$diversity
  result
}

t_start <- Sys.time()
df <- lapply(NE_START, function(Ne) {
  # config <- create_model(2000)
  config <- create_model(Ne)
  model <- config$model
  samples <- config$samples
  # plot_model(model)

  reps_df <- mclapply(1:REPS, function(rep_i) {
    ts <- simulate_ts(model, samples)
    pi <- compute_pi(ts)
    pi$rep <- rep_i
    pi$Ne_start <- Ne
    pi
  }, mc.cores = detectCores()) %>% do.call(rbind, .)

  reps_df
}) %>% do.call(rbind, .)
t_end <- Sys.time()
t_end - t_start

df %>%
mutate(time = time - T_BOTTLE) %>%
mutate(scenario = paste("Ne =", Ne_start, "→ Ne =", NE_END)) %>%
group_by(time, Ne_start, scenario, name) %>%
summarise(pi = mean(pi)) %>%
ggplot(aes(time, pi, group = interaction(as.factor(time), scenario), color = scenario)) +
  geom_boxplot(outlier.shape = NA) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "time after population crash [ky]",
       y = "nucleotide diversity") +
  scale_x_continuous(breaks = sort(c(-3000, seq(0, 50000, by = 5000)))) +
  guides(color = guide_legend("bottleneck scenario")) +
  coord_cartesian(ylim = c(0, 4.2e-4)) +
  # scale_x_continuous(breaks = unique(df$time)) +
  theme_bw()


