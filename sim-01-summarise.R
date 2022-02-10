l_results = list.files("sim-01/data/", pattern = "evaluate.*RData", full.names = TRUE) |>
  lapply(readRDS)

replacement = function(GEN) sub("replacement_(.+)-count_(.+)-size_([0-9]+)-(.+)", "\\1", GEN)
size = function(GEN) as.numeric(sub("replacement_(.+)-count_(.+)-size_([0-9]+)-(.+)", "\\3", GEN))
dataset = function(GEN) sub("(.*)-data_(.*)-seed_(.*)", "\\2", GEN)

dresults = do.call(rbind, l_results) |>
  transform(
    dataset = dataset(gen),
    Replacement = replacement(gen),
    size = size(gen))

library(ggplot2)
ggplot(data=dresults) +
  geom_smooth(aes(x = size, y = value, col = Replacement)) +
  # geom_boxplot(aes(x = size, y = value, fill = Replacement, group = paste(size,Replacement))) +
  facet_grid(metric~dataset, scales = 'free_y') +
  scale_x_continuous(trans = "reverse") +
  theme(legend.position = 'top')

