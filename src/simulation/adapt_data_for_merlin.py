import sys
import os
import numpy as np
import pandas as pd

from functions import transform_to_variant_matrix, transform_to_total_matrix

def adapt_data_for_merlin_inference(read_counts_file, output_dir, seed=1):

    # Ensure the output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
 
    # Write the read counts to a text file using numpy.savetxt for efficiency
    read_counts = np.load(read_counts_file)

    variant_matrix = transform_to_variant_matrix(read_counts)
    variant_file = os.path.join(output_dir, f"merlin_variant_matrix_{seed}.csv")
    variant_matrix.to_csv(variant_file, header=[f"cell{i+1}" for i in range(variant_matrix.shape[1])])
    
    
    total_matrix = transform_to_total_matrix(read_counts)
    total_file = os.path.join(output_dir, f"merlin_total_matrix_{seed}.csv")
    total_matrix.to_csv(total_file, header=[f"cell{i+1}" for i in range(total_matrix.shape[1])])


if __name__ == "__main__":

    # Get the parameters from the command line arguments
    read_counts_file = sys.argv[1]
    output_dir =  sys.argv[2]
    seed = sys.argv[3]


    adapt_data_for_merlin_inference(read_counts_file, output_dir, seed)
