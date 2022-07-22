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
         replacement = fct_recode(replacement,
                                  "lrnm" = "lrnm-montecarlo",
                                  "lrnm-cond-1" ="lrnb-cond-1-hermite-new",
                                  "lrnm-cond" = "lrnm-cond-montecarlo"),
         d = dim,
         dim = fct_reorder(sprintf("d = %d", dim), dim))

dplot = dresults %>%
  mutate(
    pcor = pmap_dbl(list(data, d, seed), function(.data, .d, .seed){
      fname = sprintf("%s/data/data_%s-dim_%d-seed_%05d.RData", SIM, .data, .d, .seed)
      load(fname)
      pseudo_cor = function(S){
        l = eigen(S)$values
        v = cumsum(l) / sum(l)
        1-2*(1-quantile(v[-length(v)], 0.5))
      }
      pcor = pseudo_cor(SIGMA)
      unname(pcor)
    })
  )

library(ggplot2)
p = ggplot(data=dplot) +
  geom_smooth(aes(x=pcor, y = value, col = replacement)) +
  facet_grid(.~dim, scales = 'free_x', drop = TRUE) +
  # facet_wrap(~metric, ncol = 1, scales = 'free_y') +
  # scale_x_continuous(breaks = unique(dplot$rsize)) +
  theme_minimal() +
  theme(legend.position = 'top') + labs(col = '', x = 'pseudo-correlation', y = 'STRESS')

ggsave(p, filename = sprintf("overleaf/%s-fig01b.pdf", SIM), width = 7, height = 3)
