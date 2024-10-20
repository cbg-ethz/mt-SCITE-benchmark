
# determine the number of reads to use in benchmark based on the data 
# YFV2001_scRNAseq_sub1. Calculate the mean number of reads for all sites per cell, 
# that do not have 0 coverage. Then, average across all cells, take the mean of 
# means. 

# randomly pick 5 data files
data_dir = "../../data/data_johanna/YFV2001_scRNAseq_sub1/"
files = list.files(path = dat_dir)


mean = 0

for (random_file in files){
  dat = read.csv(file = paste0(data_dir, random_file), sep = "\t")
  dat_wo_0 = dat[ dat$Good_depth !=0, ]
  mean_random_file = mean(dat_wo_0$Good_depth)
  
  mean = mean + mean_random_file
}

mean / length(files) # 468.4557

# This is approximately 500, thus take num_reads as 500 in benchmark.