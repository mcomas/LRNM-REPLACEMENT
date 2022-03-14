library(dplyr)

l_results = list.files("sim-01/data/", pattern = "evaluate.*RData", full.names = TRUE) |>
  lapply(function(x) try(readRDS(x)))
l_results = l_results[sapply(l_results, is.data.frame)]

ddata = list.files("sim-01/data/", pattern = "^data.*RData", full.names = TRUE) |>
  lapply(function(x){
    load(x)
    tibble(
      seed = sub("(.*)data_(.*)-seed_(.*).RData", "\\3", x)
    )
  }) |> bind_rows()

lns = readLines("sim-01.make")
DATA = scan(text = sub("^L_data = (.+)", "\\1", lns[grepl("^L_data = ", lns)]), what = 'character', quiet = TRUE)
SEED = as.integer(system(sub("^L_seed = \\$\\(shell (.+)\\)", "\\1", lns[grepl("^L_seed = ", lns)]), intern = TRUE))
COUNT = scan(text = sub("^L_count = (.+)", "\\1", lns[grepl("^L_count = ", lns)]), what = 'character', quiet = TRUE)
SIZE = as.integer(system(sub("^L_size = \\$\\(shell (.+)\\)", "\\1", lns[grepl("^L_size = ", lns)]), intern = TRUE))
REPLACEMENT = scan(text = sub("^L_replacement = (.+)", "\\1", lns[grepl("^L_replacement = ", lns)]), what = 'character', quiet = TRUE)
EVALUATE = scan(text = sub("^L_evaluate = (.+)", "\\1", lns[grepl("^L_evaluate = ", lns)]), what = 'character', quiet = TRUE)

replacement = function(GEN) sub("replacement_(.+)-count_(.+)-size_([0-9]+)-(.+)", "\\1", GEN)
size = function(GEN) as.numeric(sub("replacement_(.+)-count_(.+)-size_([0-9]+)-(.+)", "\\3", GEN))
dataset = function(GEN) sub("(.*)-data_(.*)-seed_(.*)", "\\2", GEN)
seed = function(GEN) sub("(.*)-data_(.*)-seed_(.*)", "\\3", GEN)

dresults = do.call(rbind, l_results) |>
  transform(
    dataset = dataset(gen),
    seed = seed(gen),
    replacement = factor(replacement(gen)),
    size = size(gen),
    metric = factor(metric)) |>
  subset(dataset %in% DATA) |>
  subset(replacement %in% REPLACEMENT) |>
  as_tibble() |>
  left_join(ddata, by = 'seed')

library(ggplot2)
set.seed(1)
p = ggplot(data=dresults) +
  geom_smooth(aes(x = size, y = value, col = replacement), se = TRUE) +
  geom_jitter(aes(x = size, y = value, col = replacement, group = paste(size,replacement)), alpha = 0.6, size = 0.5) +
  facet_grid(metric~dataset, scales = 'free_y') +
  scale_x_continuous(trans = "reverse") +
  theme(legend.position = 'top')

# library(mgcv)
# library(modelbased)
# m = gam(value~replacement*metric+s(size, k = 4, by = replacement)+s(size, k = 4, by = replacement:metric)+s(size, k = 4, by = metric), data = dresults)
# dplot1 = as.data.frame(estimate_means(m, at = c('replacement', 'size=c(30,50,70,90,110)', 'metric')))

library(dplyr)
dplot2 = dresults |>
  group_by(replacement, metric, size) |>
  summarise(m = mean(value), lo = m - 1.96*sd(value)/sqrt(n()), hi = m + 1.96*sd(value)/sqrt(n()))

ggplot(data=dplot2) +
  geom_errorbar(aes(x=size, y = m, ymin = lo, ymax=hi, col = replacement), 
                size = 1, width = 5, position = position_dodge()) +
  facet_wrap(~metric, ncol = 1, scales = 'free_y') +
  scale_x_continuous(trans = "reverse", breaks = unique(dplot2$size)) +
  theme(legend.position = 'top') + labs(y = '') 

p = ggplot(data=subset(dplot2, replacement %in% c('dm', 'lrnm', 'lrnm-cond-hermite'))) +
  geom_errorbar(aes(x=size, y = m, ymin = lo, ymax=hi, col = replacement), 
                size = 1, width = 5, position = position_dodge()) +
  facet_wrap(~metric, ncol = 1, scales = 'free_y') +
  scale_x_continuous(trans = "reverse", breaks = unique(dplot2$size)) +
  theme(legend.position = 'top') + labs(y = '') 
p
ggsave(p, width = 6, height = 5.4, filename = "overleaf/fig-sim01a.pdf")


ggplot(data=dresults) +
  geom_smooth(aes(x = cor*cor, y = value, col = replacement), se = TRUE) +
  geom_jitter(aes(x = cor*cor, y = value, col = replacement, group = paste(size,replacement)), alpha = 0.6, size = 0.5) +
  facet_grid(metric~dataset, scales = 'free_y') +
  scale_x_continuous(trans = "reverse") +
  theme(legend.position = 'top')
