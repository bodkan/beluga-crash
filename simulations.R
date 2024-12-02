library(slendr)
init_env(quiet = TRUE)

library(parallel)
library(ggplot2)
library(dplyr)
library(tidyr)

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

grid_df <- expand_grid(Ne_start = NE_START, Ne_end = NE_END, rep = 1:10)

t_start <- Sys.time()

df <- mclapply(1:nrow(grid_df), function(rep_i) {
  params_df <- grid_df[rep_i, ]

  # config <- create_model(2000)
  config <- create_model(params_df$Ne_start)
  model <- config$model
  samples <- config$samples
  # plot_model(model)

  ts <- simulate_ts(model, samples)
  pi <- compute_pi(ts)

  params_df$result <- list(pi)

  params_df
}, mc.cores = detectCores()) %>%
  do.call(rbind, .)

t_end <- Sys.time()
t_end - t_start

final_df <- df %>%
  unnest(cols = result) %>%
  mutate(time = time - T_BOTTLE) %>%
  mutate(scenario = paste("Ne =", Ne_start, "→ Ne =", NE_END)) %>%
  group_by(time, Ne_start, scenario, name) %>%
  summarise(pi = mean(pi)) %>%
  ungroup()

final_df %>%
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


final_df %>%
group_by(time, scenario) %>%
summarise(
  mean_pi = mean(pi),
  se_pi = sd(pi) / sqrt(n()),
  lower_ci = mean_pi - qt(0.975, df = n() - 1) * se_pi,
  upper_ci = mean_pi + qt(0.975, df = n() - 1) * se_pi
) %>%
ggplot() +
  geom_ribbon(aes(time, ymin = lower_ci, ymax = upper_ci, fill = scenario), alpha = 0.5) +
  geom_line(aes(time, mean_pi, color = scenario)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "time after population crash [ky]",
       y = "nucleotide diversity") +
  scale_x_continuous(breaks = sort(c(-3000, seq(0, 50000, by = 5000)))) +
  guides(color = guide_legend("bottleneck scenario"),
         fill = guide_legend("bottleneck scenario")) +
  coord_cartesian(ylim = c(0, 4.2e-4)) +
  # scale_x_continuous(breaks = unique(df$time)) +
  theme_bw()
