import sys
import os
import networkx as nx
import numpy as np
#sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src/simulation')))

from functions import reparameterize_beta, simulate_mutation_freq_tree, simulate_read_counts

def simulate_data(num_mutations, concentration, error_rate, num_reads, num_cells, output_dir, seed=1):

    # set seed
    np.random.seed(seed)
    
    # Ensure the output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Simulate the mutation frequency tree
    tree = simulate_mutation_freq_tree(num_mutations, concentration)
    
    # Write the tree to a file using networkx
    tree_file = os.path.join(output_dir, f"tree_{seed}.gml")
    nx.write_gml(tree, tree_file)
    
    # Simulate the read counts
    read_counts = simulate_read_counts(tree, error_rate, num_reads, num_cells)
    
    # Write the read counts to a text file using numpy.savetxt for efficiency
    read_counts_file = os.path.join(output_dir, f"read_counts_{seed}.npy")
    np.save(read_counts_file, read_counts)


if __name__ == "__main__":

    # Get the parameters from the command line arguments
    num_mutations = int(sys.argv[1])
    concentration = float(sys.argv[2])
    error_rate = float(sys.argv[3])
    num_reads = int(sys.argv[4])
    num_cells = int(sys.argv[5])
    output_dir = sys.argv[6]
    seed=int(sys.argv[7])

    
    # Run the simulation
    simulate_data(num_mutations, concentration, error_rate, num_reads, num_cells, output_dir, seed)
