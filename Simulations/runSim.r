# run simmulated data through replicateFest
# use matrix created by Leslie
# matrix contains 90 clones expanded at 3 different levels (low, medium, high)
# from 3 different starting frequencies (low (0.5%), medium, high)

library(replicateFest)
library(dplyr)
library(ggplot2)
library(tidyr)


# read in data
simData = readRDS("countMatrixFromSim.rds")
dim(simData)
colnames(simData)

#====================
# run one experiment to test
#=====================
# replace with new function
res_exp = replicateFest:::runFromMatrix(simData, saveToFile = "simResults.csv", 
                                        refSamp="NoPeptide")

# get some statistics

tab = splitFileName(res_exp$clone)
# remove NA 
tab = tab[!is.na(tab$condition),]
dim(tab)
table(data.frame(frequency = tab[,2], expansion = tab[,3]))

#=============================
# run 10 simulated experiments
#=============================
# get list of RDS objects to run
files = list.files(path = "Simulated_matrices", pattern = "*.rds", 
                   full.names = TRUE)
# contigency table of results
sim_res = list()
# results for each matrix
sim_res_tables = list()
for(f in files)
{
  # read the object with simulated counts
  simData = readRDS(f)
   # run replicateFest
  res_exp = replicateFest:::runFromMatrix(simData, 
                                          refSamp="NoPeptide")
  
  # get some statistics
  tab = splitFileName(res_exp$clone)
  # remove NA 
  tab = tab[!is.na(tab$condition),]
  dim(tab)
  tab1 = table(data.frame(frequency = factor(tab[,2], levels = c("L1","L2","L3")), 
                          expansion = factor(tab[,3], levels = c("E1","E2","E3"))))
  sim_res[[basename(f)]] = tab1
  sim_res_tables [[basename(f)]] = res_exp
}

saveRDS(sim_res, file = "sim_res.rds")
saveRDS(sim_res_tables, file = "sim_res_tables.rds")


# Stack all 10 tables into an array and compute mean ± SD detection rate
rate_array <- simplify2array(lapply(sim_res, function(x) x / 10))

mean_mat <- apply(rate_array, 1:2, mean)
sd_mat   <- apply(rate_array, 1:2, sd)

# Format as "mean (SD)" table
result_table <- matrix(
  sprintf("%.2f (%.2f)", mean_mat, sd_mat),
  nrow = 3, ncol = 3,
  dimnames = dimnames(sim_res[[1]])
)

cat("Mean detection rate (SD) across 10 simulated replicates\n")
cat("Rows = starting frequency level; Columns = expansion factor\n\n")
print(as.data.frame(result_table))


# Build long-format data frame with mean and SD
rate_long <- expand.grid(
  frequency = factor(c("L1","L2","L3"), levels = c("L1","L2","L3")),
  expansion = factor(c("E1","E2","E3"), levels = c("E1","E2","E3"))
)
rate_long$mean_rate <- as.vector(mean_mat)
rate_long$sd_rate   <- as.vector(sd_mat)
rate_long$label     <- sprintf("%.2f\n(\u00b1%.2f)", rate_long$mean_rate, rate_long$sd_rate)

p = ggplot(rate_long, aes(x = expansion, y = rev(frequency), fill = mean_rate)) +
  geom_tile(color = "black", linewidth = 0.8) +
  geom_text(aes(label = label), size = 3.5) +
  scale_fill_gradient(low = "white", high = "#2166AC",
                      limits = c(0, 1), name = "Detection\nrate") +
  scale_x_discrete(labels = c("E1 (2\u00d7)", "E2 (4\u00d7)", "E3 (8\u00d7)")) +
  scale_y_discrete(labels = c("L3 (95th)", "L2 (75th)", "L1 (50th)")) +
  labs(
    x = "Expansion factor",
    y = "Starting frequency percentile",
    title = "Sensitivity: mean detection rate (\u00b1 SD) across 10 simulations"
  ) +
  theme_minimal(base_size = 13) +
  theme(panel.grid = element_blank())

print(p)

# --- Save table as CSV ---
table_df <- as.data.frame(result_table) |>
  tibble::rownames_to_column("frequency")

write.csv(table_df, "detection_rate_summary.csv", row.names = FALSE)


ggsave("detection_rate_heatmap.pdf", plot = p, width = 5, height = 4)
ggsave("detection_rate_heatmap.png", plot = p, width = 5, height = 4, dpi = 300)


