l_results = list.files("sim-01/data/", pattern = "evaluate.*RData", full.names = TRUE) |>
  lapply(readRDS)
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
ggplot(data=dresults) +
  geom_smooth(aes(x = size, y = value, col = Replacement)) +
  # geom_boxplot(aes(x = size, y = value, fill = Replacement, group = paste(size,Replacement))) +
  facet_grid(metric~dataset, scales = 'free_y') +
  scale_x_continuous(trans = "reverse") +
  theme(legend.position = 'top')

