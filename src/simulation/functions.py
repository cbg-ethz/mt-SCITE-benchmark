import networkx as nx
import numpy as np
import pandas as pd

def reparameterize_beta(mode, concentration):
    """
    Reparameterize the Beta distribution using mode and concentration.

    Parameters:
    mode (float): The mode of the Beta distribution (the central tendency of the distribution).
                  Must be between 0 and 1.
    concentration (float): The concentration parameter controlling the variance/spread of the distribution.
                           Must be greater than 2.

    Returns:
    (float, float): The alpha and beta shape parameters of the Beta distribution.
    """
    if not (0 < mode < 1):
        raise ValueError("Mode must be between 0 and 1.")
    if concentration <= 2:
        raise ValueError("Concentration must be greater than 2 to be valid for this parameterization.")

    alpha = mode * (concentration - 2) + 1
    beta = (1 - mode) * (concentration - 2) + 1

    return alpha, beta


def simulate_mutation_freq_tree(num_mutations, concentration, initial_mutation_freq=0.01):

    # input validation
    if not isinstance(num_mutations, int) or num_mutations < 1:
        raise ValueError("num_mutations must be a positive integer.")

    tree = nx.DiGraph()
    # Add root node with label 0
    mutation_freqs = [0] * num_mutations
    tree.add_node(0, mutation_freq=mutation_freqs) 

    ## Note: Here i from 1, ..., num_mutations +1, denotes the mutation nodes.
    ## i-1 from 0,  ..., num_mutations indexes the mutations in the node_freqs list
    for i in range(1, num_mutations + 1):
        # Randomly select an existing node as the parent
        parent = np.random.randint(i)
        #print("New loop round")
        #print(parent)

        # Set node mutation frequencies to parent's mutation frequencies
        node_freqs = tree.nodes[parent]['mutation_freq'].copy()
        #print(node_freqs)

        # Adjust mutation frequencies for descendant nodes, drawing from a Beta distribution 
        # centered on the parent frequency
        for mutation_nr, freq in enumerate(node_freqs):
            #print(mutation_nr)
            if freq == 0 or freq == 1:
                continue
            alpha, beta = reparameterize_beta(mode=freq, concentration=concentration)
            # Draw frequency for descendant node from the Beta distribution
            node_freqs[mutation_nr] = np.random.beta(alpha, beta)
            #print(node_freqs)

        #print(node_freqs)

        # Draw a new mutation frequency from Beta centered around intial mutation freq
        a, b = reparameterize_beta(mode=initial_mutation_freq, concentration=concentration)
        new_mutation_freq = np.random.beta(a, b)

        # Add the new mutation and its frequency
        node_freqs[i-1] = new_mutation_freq
        #print(node_freqs)

        # Add the new node with mutation frequencies 
        tree.add_node(i, mutation_freq=node_freqs)
        tree.add_edge(parent, i)

    #mutation_node_mapping = {i : i for i in range(num_mutations + 1)}
    return tree #, mutation_node_mapping


def simulate_read_counts(mutation_tree, error_rate, num_reads, num_cells):

    n_nodes =  mutation_tree.number_of_nodes()
    n_mutations = n_nodes - 1

    # init output tensor num_cells x num sites (or mutations) x 4
    read_counts = np.zeros((num_cells, n_mutations, 4), dtype=int)

    for cell in range(num_cells):

        # pick cell attachment node at random
        attachment_node = np.random.randint(n_nodes)
        mutation_freqs_node = mutation_tree.nodes[attachment_node]['mutation_freq']
        #print("attachment node")
        #print(attachment_node)

        # First, handle allele counts for mutations 
        for mutation_nr in range(n_mutations):

            # Draw frequency of mitochondrial mutation within cell (heteroplasmy)
            freq = mutation_freqs_node[mutation_nr] 

            # Handle allele counts for no mutations (sequencing errors only)
            if freq == 0:
                allele_counts = np.random.multinomial(num_reads, [1 - error_rate, error_rate / 3, error_rate / 3, error_rate / 3])

            # Simulate allele counts given mutation frequency, accounting for sequencing error
            else:
                # Ensure valid probabilities
                p0 =  1 - freq - error_rate  # Probability for allele 0
                #print(p0)
                p1 =  freq  # Probability for mutated allele

                # Check if p0 and p1 are valid (within the range [0, 1])
                if not (0 <= p0 <= 1):
                    raise RuntimeError(f"Invalid probability p0={p0}. Expected value between 0 and 1.")
                if not (0 <= p1 <= 1):
                    raise RuntimeError(f"Invalid probability p1={p1}. Expected value between 0 and 1.")

                p2 = p3 = error_rate / 2  # Probability for sequencing error

                allele_counts = np.random.multinomial(num_reads, [p0, p1, p2, p3])


            read_counts[cell, mutation_nr] = allele_counts

    return read_counts

def transform_array_to_presence_absence_matrix(array):
    """
    Transforms a 3D numpy array (cells x mutations x states) into a presence/absence matrix.
    Rows represent mutations, columns represent cells.
    Presence is marked as 1 if any of X, Y, Z is greater than 0 for a mutation in a given cell.
    
    Parameters:
    array (numpy array): 3D numpy array where shape is (cells, mutations, states), and states are [R, X, Y, Z].
    
    Returns:
    numpy array: A presence/absence matrix with mutations as rows and cells as columns. This will be used as the input to scite.
    """
    # Extract the X, Y, Z counts (columns 1, 2, 3 in the array)
    # Check if any of X, Y, Z > 0 for each mutation in each cell
    presence_absence_matrix = (array[:, :, 1:] > 0).any(axis=2).astype(int)
    
    return presence_absence_matrix.T  # Transpose to match mutations as rows and cells as columns


def transform_to_variant_matrix(read_count_matrix):
    """
    Transforms the original read count matrix into a variant matrix.
    
    Parameters:
    - read_count_matrix (np.array): A 3D numpy array where each cell contains 
                                    counts for each mutation state.
    
    Returns:
    - pd.DataFrame: A DataFrame representing the variant matrix, where each row 
                    is a mutation, each column is a cell, and each entry 
                    represents the number of variant reads.
    """
    num_cells = read_count_matrix.shape[0]
    num_mutations = read_count_matrix.shape[1]
    
    # Initialize an empty variant matrix
    variant_matrix = np.zeros((num_mutations, num_cells), dtype=int)
    
    # Populate the variant matrix
    for cell_idx in range(num_cells):
        for mut_idx in range(num_mutations):
            # Sum the variant reads (assuming these are the 2nd, 3rd, and 4th entries)
            variant_reads = read_count_matrix[cell_idx, mut_idx, 1:]
            variant_matrix[mut_idx, cell_idx] = np.sum(variant_reads)
    
    # Convert the matrix to a DataFrame with appropriate row and column names
    variant_df = pd.DataFrame(
        variant_matrix,
        index=[f"mut{i+1}" for i in range(num_mutations)],
        columns=[f"cell{j+1}" for j in range(num_cells)]
    )
    
    return variant_df


def transform_to_total_matrix(read_count_matrix):
    """
    Transforms the original read count matrix into a total matrix.
    
    Parameters:
    - read_count_matrix (np.array): A 3D numpy array where each cell contains 
                                    counts for each mutation state.
    
    Returns:
    - pd.DataFrame: A DataFrame representing the total matrix, where each row 
                    is a mutation, each column is a cell, and each entry 
                    represents the total number of reads for that mutation.
    """
    num_cells = read_count_matrix.shape[0]
    num_mutations = read_count_matrix.shape[1]
    
    # Initialize an empty total matrix
    total_matrix = np.zeros((num_mutations, num_cells), dtype=int)
    
    # Populate the total matrix
    for cell_idx in range(num_cells):
        for mut_idx in range(num_mutations):
            # Sum all reads for each mutation across all states
            total_reads = np.sum(read_count_matrix[cell_idx, mut_idx, :])
            total_matrix[mut_idx, cell_idx] = total_reads
    
    # Convert the matrix to a DataFrame with appropriate row and column names
    total_df = pd.DataFrame(
        total_matrix,
        index=[f"mut{i+1}" for i in range(num_mutations)],
        columns=[f"cell{j+1}" for j in range(num_cells)]
    )
    
    return total_df
