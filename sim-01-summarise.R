l_results = list.files("sim-01/data/", pattern = "evaluate.*RData", full.names = TRUE) |>
  lapply(function(x) try(readRDS(x)))
l_results = l_results[sapply(l_results, is.data.frame)]

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

dresults = do.call(rbind, l_results) |>
  transform(
    dataset = dataset(gen),
    Replacement = replacement(gen),
    size = size(gen)) |>
  subset(dataset %in% DATA) |>
  subset(Replacement %in% REPLACEMENT)

library(ggplot2)
set.seed(1)
ggplot(data=dresults) +
  geom_smooth(aes(x = size, y = value, col = Replacement), se = TRUE) +
  geom_jitter(aes(x = size, y = value, col = Replacement, group = paste(size,Replacement)), alpha = 0.6, size = 0.5) +
  facet_grid(metric~dataset, scales = 'free_y') +
  scale_x_continuous(trans = "reverse") +
  theme(legend.position = 'top')

ggplot(data=subset(dresults, grepl("^lrn", Replacement))) +
  geom_smooth(aes(x = size, y = value, col = Replacement), se = TRUE) +
  geom_jitter(aes(x = size, y = value, col = Replacement, group = paste(size,Replacement)), alpha = 0.6, size = 0.5) +
  facet_grid(metric~dataset, scales = 'free_y') +
  scale_x_continuous(trans = "reverse") +
  theme(legend.position = 'top')
