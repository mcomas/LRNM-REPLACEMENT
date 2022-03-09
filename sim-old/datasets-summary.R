lns = readLines("sim-01.make")
DATA = scan(text = sub("^L_data = (.+)", "\\1", lns[grepl("^L_data = ", lns)]), what = 'character', quiet = TRUE)
SEED = as.integer(system(sub("^L_seed = \\$\\(shell (.+)\\)", "\\1", lns[grepl("^L_seed = ", lns)]), intern = TRUE))
COUNT = scan(text = sub("^L_count = (.+)", "\\1", lns[grepl("^L_count = ", lns)]), what = 'character', quiet = TRUE)
SIZE = as.integer(system(sub("^L_size = \\$\\(shell (.+)\\)", "\\1", lns[grepl("^L_size = ", lns)]), intern = TRUE))

datasets = expand.grid(count = COUNT, size = SIZE, data = DATA, seed = SEED)
datasets = transform(datasets,
                     fname = sprintf("sim-01/data/count_%s-size_%05d-data_%s-seed_%05d.RData",
                                    count, size, data, seed))
dsummary = lapply(datasets$fname, function(fname_){
  load(fname_)
  is0 = X == 0
  is1 = X == 1
  data.frame(
    parts = ncol(X),
    prop.0 = mean(is0),
    prop.1 = mean(is1),
    complete.obs = mean(rowMeans(is0) == 0),
    complete.vars = mean(colMeans(is0) == 0)
  )
}) |> (function(l_res) do.call(rbind, l_res))()

datasets_summary = cbind(datasets, dsummary)

save(datasets_summary, file = 'sim-01/datasets-summary.RData')
