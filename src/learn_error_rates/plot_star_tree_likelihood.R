library(ggplot2)
setwd("~/Projects/mitochon/mt-scite-benchmark/src/learn_error_rates/")
n_cells = 100
n_mutations = 10
true_error_rate = 0.05
seed = 2
initial_mut_freq = 0.01
epsilon = 0.0001
# Do not use
# compute_normalized_likelihood <- function(dat, dat_star, epsilon= 0.0001) {
#   normalised_dat <- dat
#   for (row in 2:nrow(dat)) {
#     normalised_dat[row, -1] <- dat[row, -1] / (dat[row, -1] - dat_star[row, -1] + epsilon)
#   }
#   return(normalised_dat)
# }

compute_normalized_likelihood <- function(dat, dat_star) {
  normalised_dat <- dat
  for (row in 2:nrow(dat)) {
    normalised_dat[row, -1] <- dat[row, -1] /  dat_star[row, -1]
  }
  return(normalised_dat)
}

file_likelihoods = paste0("~/Projects/mitochon/mt-scite-benchmark/results/inference_output/",
                   n_mutations,"_120_", format(true_error_rate, scientific = F), "_500_", n_cells, "_", 
                   initial_mut_freq, "/seed_", seed, "/val_scores.txt")

file_star_likelihoods = paste0("~/Projects/mitochon/mt-scite-benchmark/results/inference_output/",
                     n_mutations,"_120_", format(true_error_rate, scientific = F), "_500_", n_cells, "_", 
                     initial_mut_freq, "/seed_", seed, "/val_scores_star_trees.txt")


dat = read.csv(header = F, file_likelihoods)
dat_star = read.csv(header = F, file_star_likelihoods)

new_dat = dat_star
#new_dat = compute_normalized_likelihood(dat, dat_star)

n_columns = ncol(new_dat)
error_rate = unlist(unname(new_dat[1, 2:n_columns]))
error_rate
mean_log_lik = apply(new_dat, MARGIN = 2, FUN = function(x){
  mean(x[2:31])
})[2:n_columns] # remove the first column because this is just the integers indicating the k-split
sd_log_lik = apply(new_dat, MARGIN = 2, FUN = function(x){
  sd(x[2:31])
})[2:n_columns]


# mean_log_lik_2 = apply(new_dat_2, MARGIN = 2, FUN = function(x){
#   mean(x[2:31])
# })[2:n_columns] # remove the first column because this is just the integers indicating the k-split
# sd_log_lik_2 = apply(new_dat_2, MARGIN = 2, FUN = function(x){
#   sd(x[2:31])
# })[2:n_columns]


dat_to_plot = data.frame(error_rate = error_rate, 
                         mean_log_lik = mean_log_lik,
                         sd_log_lik = sd_log_lik
                         #mean_2 = mean_log_lik_2,
                         #sd_2 = sd_log_lik_2
)
dat_to_plot



max_log_lik = min(mean_log_lik)
best_error_rate = error_rate[which(mean_log_lik == max_log_lik)]
best_error_rate

g= ggplot(dat_to_plot, aes(x=log(error_rate), y=mean_log_lik)) + 
  geom_point()+
  #geom_point(aes(y=mean_2), col="red")+
  #geom_errorbar(mapping = aes(ymin=mean_2-sd_2, ymax=mean_2+sd_2), col="red")+
  geom_errorbar(mapping = aes(ymin=mean_log_lik-sd_log_lik, ymax=mean_log_lik+sd_log_lik))+
  geom_vline(xintercept = log(true_error_rate), col = "darkgreen", linetype="dashed")+
  theme_classic()+
  geom_vline(xintercept = log(best_error_rate), col="red", linetype="dashed")+
  ylab("Star tree likelihood")+
  xlab("Log error rate")

g


ggsave(plot = g, width = 6.4, height = 4.8, units = "in",
       filename = paste0("plots/", n_mutations,"_120_", format(true_error_rate, scientific = F),
                         "_500_", n_cells, "_", initial_mut_freq, "_", seed, "_star_tree_log_lik.png"))

#### plot mutation matrices and check if they are different ###
#mut_dir = paste0("~/Projects/mitochon/mt-scite-benchmark/results/inference_output/",
#                 n_mutations,"_120_", format(true_error_rate, scientific = F), "_500_", n_cells, "_", initial_mut_freq, "/seed_", seed, "/")
#list.files(mut_dir)
#mut_file = paste0(mut_dir, "0.500000.csv")
#mut_5e5 = read.csv(header = F, sep = " ",
#                   file=mut_file)
#apply(mut_5e5, MARGIN = 2, max)