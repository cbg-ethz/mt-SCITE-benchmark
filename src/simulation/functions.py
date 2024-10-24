import networkx as nx
import numpy as np

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
