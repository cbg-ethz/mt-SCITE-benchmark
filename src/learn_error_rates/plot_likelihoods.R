library(ggplot2)
setwd("~/Projects/mitochon/mt-scite-benchmark/src/learn_error_rates/")
n_cells = 500
n_mutations = 50
true_error_rate = 0.005
seed = 10
initial_mut_freq = 0.1
epsilon = 0.0001

file_name = paste0("~/Projects/mitochon/mt-scite-benchmark/results/inference_output/",
                   n_mutations,"_120_", format(true_error_rate, scientific = F), "_500_", n_cells, "_", 
                   initial_mut_freq, "/seed_", seed, "/val_scores.txt")

file_name_2 = paste0("~/Projects/mitochon/mt-scite-benchmark/results/inference_output/",
                   n_mutations,"_120_", format(true_error_rate, scientific = F), "_500_", n_cells, "_", 
                   initial_mut_freq, "/seed_", seed, "/val_scores_star_trees.txt")


file_name
dat = read.csv(header = F, file_name)
dat_star = read.csv(header = F, file_name_2)

new_dat = dat
new_dat[2, ] = ((dat[2, ] - dat_star[2, ])/dat[2, ])
new_dat[3, ] = ((dat[3, ] - dat_star[3, ])/dat[3, ])
## Plot statistic

new_dat[2, ] = dat[2, ] / (dat[2, ] - dat_star[2, ] + epsilon)
new_dat[3, ] = dat[3, ] / (dat[3, ] - dat_star[3, ] + epsilon)
new_dat[4, ] = dat[4, ] / (dat[4, ] - dat_star[4, ] + epsilon)

# Plot signal in alignment only
#new_dat[2, ] =  (dat[2, ] - dat_star[2, ])
#new_dat[3, ] =  (dat[3, ] - dat_star[3, ])
#new_dat[4, ] =  (dat[4, ] - dat_star[4, ])

# Plot log lik only
#new_dat[2, ] =  (dat[2, ] )
#new_dat[3, ] =  (dat[3, ] )
#new_dat[4, ] =  (dat[4, ] )

n_columns = ncol(new_dat)
error_rate = unlist(unname(new_dat[1, 2:n_columns]))
error_rate
mean_log_lik = apply(new_dat, MARGIN = 2, FUN = function(x){
  mean(x[2:31])
})[2:n_columns] # remove the first column because this is just the integers indicating the k-split
sd_log_lik = apply(new_dat, MARGIN = 2, FUN = function(x){
  sd(x[2:31])
})[2:n_columns]


dat_to_plot = data.frame(error_rate = error_rate, 
                         mean_log_lik = mean_log_lik,
                         sd_log_lik = sd_log_lik
                         )
dat_to_plot



max_log_lik = max(mean_log_lik)
best_error_rate = error_rate[which(mean_log_lik == max_log_lik)]
best_error_rate

g= ggplot(dat_to_plot, aes(x=log(error_rate), y=mean_log_lik)) + 
  geom_point()+
  geom_errorbar(mapping = aes(ymin=mean_log_lik-sd_log_lik, ymax=mean_log_lik+sd_log_lik))+
  geom_vline(xintercept = log(true_error_rate), col = "darkgreen")+
  geom_hline(yintercept = 0, col = "pink") +
  geom_vline(xintercept = log(true_error_rate*10))+
  geom_vline(xintercept = log(true_error_rate/10)) +
  theme_bw()+
  geom_vline(xintercept = log(best_error_rate), col="blue")+
  ylab("Diff log lik")

g


ggsave(plot = g, width = 6.4, height = 4.8, units = "in",
       filename = paste0("plots/", n_mutations,"_120_", format(true_error_rate, scientific = F),
                         "_500_", n_cells, "_", initial_mut_freq, "_", seed, "_rel_diff_loglik.png"))

#### plot mutation matrices and check if they are different ###
mut_dir = paste0("~/Projects/mitochon/mt-scite-benchmark/results/inference_output/",
                   n_mutations,"_120_", format(true_error_rate, scientific = F), "_500_", n_cells, "_", initial_mut_freq, "/seed_", seed, "/")
list.files(mut_dir)
mut_file = paste0(mut_dir, "0.500000.csv")
mut_5e5 = read.csv(header = F, sep = " ",
               file=mut_file)
apply(mut_5e5, MARGIN = 2, max)

