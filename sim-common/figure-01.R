library(tidyverse)

if(!exists("SIM")) SIM = 'sim-01a'

l_results = list.files(sprintf("%s/data/", SIM), pattern = "evaluate.*rds", full.names = TRUE) %>%
  lapply(function(x) try(readRDS(x)))
l_results = l_results[sapply(l_results, is.data.frame)]

dresults = l_results %>%
  bind_rows() %>%
  as_tibble() %>%
  filter(metric == 'STRESS') %>%
  extract(gen, c('replacement', 'count', 'size', 'data', 'dim', 'seed'), 
          'replacement_(.+)-count_(.+)-size_(.+)-data_(.+)-dim_(.+)-seed_(.+)', convert = TRUE) %>%
  mutate(rsize = if_else(count == 'uniform', size, dim * size),
         dim = fct_reorder(sprintf("d = %d", dim), dim))

dplot = dresults %>%
  group_by(replacement, metric, count, dim, size, rsize) %>%
  summarise(m = mean(value, na.rm=TRUE), lo = m - 1.96*sd(value, na.rm=TRUE)/sqrt(n()), hi = m + 1.96*sd(value, na.rm=TRUE)/sqrt(n()))

dplot = dresults %>%
  group_by(replacement, metric, count, dim, size, rsize) %>%
  summarise(m = median(value, na.rm=TRUE), 
            lo = quantile(value, 0.25, na.rm=TRUE), 
            hi = quantile(value, 0.75, na.rm=TRUE))


library(ggplot2)
p = ggplot(data=dplot) +
  geom_errorbar(aes(x=factor(rsize), y = m, ymin = lo, ymax=hi, col = replacement), 
                size = 1, width = 0.2, position = position_dodge(width =0.5)) +
  facet_grid(.~dim, scales = 'free_x', drop = TRUE) +
  # facet_wrap(~metric, ncol = 1, scales = 'free_y') +
  # scale_x_continuous(breaks = unique(dplot$rsize)) +
  theme_minimal() +
  theme(legend.position = 'top') + labs(col = '', x = 'n', y = 'STRESS')

ggsave(p, filename = sprintf("overleaf/%s-fig01.pdf", SIM), width = 7, height = 3)
