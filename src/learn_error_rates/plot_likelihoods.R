library(ggplot2)
setwd("~/Projects/mitochon/mt-scite-benchmark/src/learn_error_rates/")
n_cells = 500
n_mutations = 10
true_error_rate = 0.005
seed = 1
initial_mut_freq = 0.1
epsilon = 0.0001


compute_normalized_likelihood <- function(dat, dat_star, epsilon) {
  normalised_dat <- dat
  for (row in 2:nrow(dat)) {
    normalised_dat[row, -1] <- dat[row, -1] / (dat[row, -1] - dat_star[row, -1] + epsilon)
  }
  return(normalised_dat)
}

# assuming x and y are log likelihoods
score_1 = function(lik_tree, lik_star){
  score_div = lik_tree / lik_star
  return(score_div)
}

score_2 = function(lik_tree, lik_star, epsilon){
  score = lik_tree / (lik_tree - lik_star + epsilon)
  return(score)
}


file_name = paste0("~/Projects/mitochon/mt-scite-benchmark/results/inference_output/",
                   n_mutations,"_120_", format(true_error_rate, scientific = F), "_500_", n_cells, "_", 
                   initial_mut_freq, "/seed_", seed, "/val_scores.txt")

file_name_2 = paste0("~/Projects/mitochon/mt-scite-benchmark/results/inference_output/",
                   n_mutations,"_120_", format(true_error_rate, scientific = F), "_500_", n_cells, "_", 
                   initial_mut_freq, "/seed_", seed, "/val_scores_star_trees.txt")


file_name
dat = read.csv(header = F, file_name)
dat_star = read.csv(header = F, file_name_2)

score_2_dat = compute_normalized_likelihood(dat, dat_star, epsilon)

score_1_dat = score_2_dat = dat

for (i in 2:nrow(dat)){
  score_1_dat[i, ] = score_1(dat[i, ], dat_star[i, ])
  score_2_dat[i, ] = score_2(lik_tree = dat[i, ], lik_star = dat_star[i, ], epsilon = epsilon)
}

# new_dat_2 = new_dat = dat
# new_dat[2, ] = ((dat[2, ] - dat_star[2, ])/dat[2, ])
# new_dat[3, ] = ((dat[3, ] - dat_star[3, ])/dat[3, ])
# ## Plot statistic
# 
# new_dat[2, ] = dat[2, ] / (dat[2, ] - dat_star[2, ] + epsilon)
# new_dat[3, ] = dat[3, ] / (dat[3, ] - dat_star[3, ] + epsilon)
# new_dat[4, ] = dat[4, ] / (dat[4, ] - dat_star[4, ] + epsilon)
# 
# 
# ## Plot rewritten statistic
# 
# new_dat_2[2, ] = 1 - dat_star[2, ] / dat[2, ]
# # Plot signal in alignment only
# #new_dat[2, ] =  (dat[2, ] - dat_star[2, ])
# #new_dat[3, ] =  (dat[3, ] - dat_star[3, ])
# #new_dat[4, ] =  (dat[4, ] - dat_star[4, ])
# 
# # Plot log lik only
# new_dat[2, ] =  (dat[2, ] /dat_star[2, ])
# new_dat[3, ] =  (dat[3, ] /dat_star[3, ])
# new_dat[4, ] =  (dat[4, ] /dat_star[4, ])

n_columns = ncol(dat)
n_rows = nrow(dat)
error_rate = unlist(unname(dat[1, 2:n_columns]))
error_rate

#mean_log_lik = apply(new_dat, MARGIN = 2, FUN = function(x){
#  mean(x[2:31])
#})[2:n_columns] # remove the first column because this is just the integers indicating the k-split
#sd_log_lik = apply(new_dat, MARGIN = 2, FUN = function(x){
#  sd(x[2:31])
#})[2:n_columns]

mean_score_1 = apply(score_1_dat, MARGIN = 2, FUN = function(x){
  mean(x[2:n_rows])
})[2:n_columns]
mean_score_2 = apply(score_2_dat, MARGIN = 2, FUN = function(x){
  mean(x[2:n_rows])
})[2:n_columns]
mean_lik_tree = apply(dat, MARGIN = 2, FUN = function(x){
  mean(x[2:n_rows])
})[2:n_columns]
mean_lik_star = apply(dat_star, MARGIN = 2, FUN = function(x){
  mean(x[2:n_rows])
})[2:n_columns]

# Plot this differently!
dat_to_plot = data.frame(error_rate = error_rate, 
                         lik_tree = mean_lik_tree,
                         #sd_log_lik = sd_log_lik,
                         lik_star = mean_lik_star,
                         score_1 = mean_score_1,
                         score_2 = mean_score_2
                         )
dat_to_plot



max_log_lik = max(mean_score_1)
best_error_rate = error_rate[which(mean_score_1 == max_log_lik)]
best_error_rate

max_log_lik = max(mean_score_2)
best_error_rate = error_rate[which(mean_score_2 == max_log_lik)]
best_error_rate

g= ggplot(dat_to_plot, aes(x=log(error_rate), y=lik_tree)) + 
  geom_point()+
  geom_jitter(aes(y=lik_star), width=0.05, col="red")+
  geom_jitter(aes(y=score_2), width=0.05, col="blue") + 

  geom_vline(xintercept = log(true_error_rate), col = "darkgreen")+
  theme_bw()+
  ylab("Lik div star")

  
g


ggsave(plot = g, width = 6.4, height = 4.8, units = "in",
       filename = paste0("plots/", n_mutations,"_120_", format(true_error_rate, scientific = F),
                         "_500_", n_cells, "_", initial_mut_freq, "_", seed, "_loglik_div_starlik.png"))


#geom_jitter(aes(y=score_1), width=0.05, col="lightgreen")+
#geom_jitter(aes(y=score_2), width=0.05, col="blue")+
#geom_errorbar(mapping = aes(ymin=mean_log_lik-sd_log_lik, ymax=mean_log_lik+sd_log_lik))+
#geom_hline(yintercept = 0, col = "pink") +
#geom_vline(xintercept = log(best_error_rate), col="blue")+
#ylim(-500, 2)
#scale_y_log10() # Set y-axis to log scale

#### plot mutation matrices and check if they are different ###
mut_dir = paste0("~/Projects/mitochon/mt-scite-benchmark/results/inference_output/",
                   n_mutations,"_120_", format(true_error_rate, scientific = F), "_500_", n_cells, "_", initial_mut_freq, "/seed_", seed, "/")
list.files(mut_dir)
mut_file = paste0(mut_dir, "0.500000.csv")
mut_5e5 = read.csv(header = F, sep = " ",
               file=mut_file)
apply(mut_5e5, MARGIN = 2, max)

####

x =0.1
y = 0.01

score_1 = function(x, y){
  score = log(x) / (log(x) - log(y))
  return(score)
}

score_2 = function(x, y){
  score = 1/( 1 - log(y) / log(x))
  return(score)
}

x = 0.001
y = 0.01

score_1(x,y)
score_2(x,y)