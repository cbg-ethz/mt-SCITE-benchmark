import sys
import os
import numpy as np

from functions import transform_array_to_presence_absence_matrix

def adapt_data_for_scite_inference(read_counts_file, output_dir, seed=1):

    # Ensure the output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
 
    # Write the read counts to a text file using numpy.savetxt for efficiency
    read_counts = np.load(read_counts_file)

    input_dat_scite = transform_array_to_presence_absence_matrix(read_counts)

    output_file = os.path.join(output_dir, f"scite_input_{seed}.txt")
    np.savetxt(output_file, input_dat_scite, delimiter=" ", fmt="%d")


if __name__ == "__main__":

    # Get the parameters from the command line arguments
    read_counts_file = sys.argv[1]
    output_dir =  sys.argv[2]
    seed = sys.argv[3]


    adapt_data_for_scite_inference(read_counts_file, output_dir, seed)
