library(tidyverse)

l_results = list.files("sim-01b/data/", pattern = "evaluate.*rds", full.names = TRUE) %>%
  lapply(function(x) try(readRDS(x)))
l_results = l_results[sapply(l_results, is.data.frame)]

dresults = l_results %>%
  bind_rows() %>%
  extract(gen, c('replacement', 'count', 'size', 'data', 'dim', 'seed'), 
          'replacement_(.+)-count_(.+)-size_(.+)-data_(.+)-dim_(.+)-seed_(.+)', convert = TRUE) %>%
  mutate(rsize = if_else(count == 'uniform', size, dim * size))

dplot = dresults %>%
  group_by(replacement, metric, count, dim, size) %>%
  summarise(m = mean(value), lo = m - 1.96*sd(value)/sqrt(n()), hi = m + 1.96*sd(value)/sqrt(n()))

library(ggplot2)
p = ggplot(data=dplot) +
  geom_errorbar(aes(x=size, y = m, ymin = lo, ymax=hi, col = replacement), 
                size = 1, width = 5, position = position_dodge(width = 20)) +
  facet_grid(count~dim, scales = 'free_y') +
  # facet_wrap(~metric, ncol = 1, scales = 'free_y') +
  # scale_x_continuous(trans = "log", breaks = unique(dplot$size)) +
  theme(legend.position = 'top') + labs(y = '') 
p
