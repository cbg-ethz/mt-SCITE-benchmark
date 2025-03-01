import marimo

__generated_with = "0.8.3"
app = marimo.App(width="medium")


@app.cell
def __():
    import marimo as mo
    import networkx as nx
    import numpy as np
    import matplotlib.pyplot as plt
    import pandas as pd
    return mo, np, nx, pd, plt


@app.cell
def __():
    # old code
    ''' 
    def simulate_mutation_tree(num_mutations):
        """
    Simulates a mutation tree with m + 1 nodes, where each additional node
    represents a new mutation attached to a randomly chosen parent node.

    Parameters:
        num_mutations (int): The number of mutation nodes (excluding the root) to generate.

    Returns:
        tree (networkx.DiGraph): A directed graph representing the mutation tree, 
            where node 0 is the root, and each other node corresponds to a mutation.
        mutation_node_mapping (dict): A dictionary mapping each node to its mutation.

    Raises:
        ValueError: If 'm' is not a positive integer.
    """

        # input validation
        if not isinstance(num_mutations, int) or num_mutations < 1:
            raise ValueError("num_mutations must be a positive integer.")

        tree = nx.DiGraph()
        # Add root node with label 0
        tree.add_node(0, mutations=[]) 

        for i in range(1, num_mutations + 1):
            # Randomly select an existing node as the parent
            parent = np.random.randint(i)

            # Add the new mutation node and create an edge between parent and child
            tree.add_node(i, mutations=[i])
            tree.add_edge(parent, i)

        mutation_node_mapping = {i : i for i in range(num_mutations + 1)}
        return tree, mutation_node_mapping


    def map_node_to_accumulated_mutations(tree, mutation_node_mapping):
        """
        Maps each node in the tree to all mutations accumulated along its path from the root.

        Parameters:
            tree (networkx.DiGraph): The mutation tree generated by simulate_mutation_tree.
            mutation_node_mapping (dict): A dictionary mapping each node to its mutation.

        Returns:
            accumulated_mutations (dict): A dictionary mapping each node to a list of mutations
                accumulated along its path from the root.
        """
        accumulated_mutations = {}

        # Traverse each node in the tree
        for node in tree.nodes:
            mutations = []
            current_node = node

            # Traverse upwards to the root
            while current_node != 0:  # root node is 0
                mutations.append(mutation_node_mapping[current_node])
                current_node = list(tree.predecessors(current_node))[0]  # move to parent node

            accumulated_mutations[node] = list(reversed(mutations))  # reverse to get root-to-leaf order

        return accumulated_mutations


    def simulate_read_counts(mutation_tree, error_rate, map_node_to_accumulated_mutations, num_reads, num_cells):

        """
        Simulates read counts for a population of cells given a mutation tree.

        Parameters:
            mutation_tree (networkx.DiGraph): The mutation tree with parent-child relationships.
            error_rate (float): The sequencing error rate.
            map_node_to_accumulated_mutations (dict): A dictionary mapping nodes to a list of mutations on the path from the root.
            num_reads (int): The total number of reads per site (with/o mutation).
            num_cells (int): The total number of cells to simulate.

        Returns:
            read_counts (np.array): A tensor of shape (num_cells, num_mutations, 4), representing allele counts for the 4 nucleotides reference, mutation and 2 alternative nucleotides.
        """

        n_nodes =  mutation_tree.number_of_nodes()

        # init output tensor num_cells x num sites (or mutations) x 4
        read_counts = np.zeros((num_cells, n_nodes-1, 4), dtype=int)

        for cell in range(num_cells):

            # pick cell attachment node at random
            attachment_node = np.random.randint(n_nodes)
            print("attachment node")
            print(attachment_node)

            # First, handle allele counts for mutations 
            for mutation in map_node_to_accumulated_mutations[attachment_node]:

                # Draw frequency of mitochondrial mutation within cell (heteroplasmy)
                freq = np.random.uniform(0.1, 0.9) 
                 # np.random.beta(2, 2) = 0.98 then problems with error rate = 0.05!

                # Simulate allele counts given mutation frequency, accounting for sequencing error
                allele_counts = np.random.multinomial(
                    num_reads, [1 - freq - error_rate, freq, error_rate / 2, error_rate / 2]
                )

                read_counts[cell, mutation-1] = allele_counts

            # Second, handle allele counts for no mutations (sequencing errors only)
            for mutation in range(1, n_nodes):
                print(mutation)
                if mutation not in map_node_to_accumulated_mutations[attachment_node]:
                    allele_counts = np.random.multinomial(
                        num_reads, [1 - error_rate, error_rate / 3, error_rate / 3, error_rate / 3]
                    )
                    read_counts[cell, mutation-1] = allele_counts
                    print(read_counts)

        return read_counts

    '''
    return


@app.cell
def __(np, nx):
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

        return tree #, mutation_node_mapping
    return reparameterize_beta, simulate_mutation_freq_tree


@app.cell
def __(np):
    err_cells = np.random.randint(10, size=2)
    print(err_cells)
    if 0 in err_cells:
        print ("tre")
    return err_cells,


@app.cell
def __():
    n_err_mutations = 2
    n_err_samples = 2
    num_cells = 100
    return n_err_mutations, n_err_samples, num_cells


@app.cell
def __(n_err_mutations, n_err_samples, np, num_cells):
    err_mutations_with_cells = np.zeros((n_err_mutations, n_err_samples), dtype=int)
    err_mutations_with_cells
    for err_mutation in range(n_err_mutations):
        cells_with_err_mutations = np.random.randint(num_cells, size=n_err_samples)
        err_mutations_with_cells[err_mutation, ] = cells_with_err_mutations
        
    return cells_with_err_mutations, err_mutation, err_mutations_with_cells


@app.cell
def __(err_mutations_with_cells):
    err_mutations_with_cells
    return


@app.cell
def __(np):
    def sample_allele_counts(freq, error_rate, num_reads):
        p0 =  1 - freq - error_rate  # Probability for allele 0
        p1 =  freq  # Probability for mutated allele

        # Check if p0 and p1 are valid (within the range [0, 1])
        if not (0 <= p0 <= 1):
            raise RuntimeError(f"Invalid probability p0={p0}. Expected value between 0 and 1.")
        if not (0 <= p1 <= 1):
            raise RuntimeError(f"Invalid probability p1={p1}. Expected value between 0 and 1.")

        p2 = p3 = error_rate / 2  # Probability for sequencing error

        allele_counts = np.random.multinomial(num_reads, [p0, p1, p2, p3])

        return allele_counts
        
    def simulate_read_counts(mutation_tree, error_rate, num_reads, num_cells, n_err_mutations, n_err_samples):

        n_nodes =  mutation_tree.number_of_nodes()
        n_mutations = n_nodes - 1

        # init output tensor num_cells x num sites (or mutations) x 4
        read_counts = np.zeros((num_cells, n_mutations+n_err_mutations, 4), dtype=int)
        mutation_probabilities = np.zeros((num_cells, n_mutations), 
                                          dtype=float)

        # Determine, which non-clonal/erroneous mutations are present in which cells
        err_mutations_with_cells = np.zeros((n_err_mutations, n_err_samples), dtype=int)
        err_mutations_with_cells
        for err_mutation in range(n_err_mutations):
            cells_with_err_mutations = np.random.randint(num_cells, size=n_err_samples)
            err_mutations_with_cells[err_mutation, ] = cells_with_err_mutations
        

        for cell in range(num_cells):

            # pick cell attachment node at random
            attachment_node = np.random.randint(n_nodes)
            
            # Create read counts for clonal mutations
            mutation_freqs_node = mutation_tree.nodes[attachment_node]['mutation_freq']

            mutation_probabilities[cell, ] = mutation_freqs_node
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
                    
                    allele_counts = sample_allele_counts(freq, error_rate, num_reads)


                read_counts[cell, mutation_nr] = allele_counts

            # Create read counts for non-clonal mutations
            for n_err_mutation in range(n_err_mutations):
                cells_with_err_mutations = err_mutations_with_cells[n_err_mutation, ]

                # if this cell has an erroneous mutation, sample its read count
                if cell in cells_with_err_mutations:
                    # use fixed initial mutation frequency
                    mutation_freq = 0.01
                    allele_counts = sample_allele_counts(mutation_freq, error_rate, num_reads)

                read_counts[cell, n_mutations+n_err_mutation] = allele_counts
                

        return [read_counts, mutation_probabilities]
    return sample_allele_counts, simulate_read_counts


@app.cell
def __(mf_tree, simulate_read_counts):
    dat, mut_probs = simulate_read_counts(mutation_tree=mf_tree, error_rate=0.0009, num_reads=500, num_cells=10, n_err_mutations=2, n_err_samples=10)
    return dat, mut_probs


@app.cell
def __(dat):
    dat[:, -3:, ]
    return


@app.cell
def __(dat):
    dat
    return


@app.cell
def __(dat, np):
    np.shape(dat)
    return


@app.cell
def __(dat):
    dat[0]
    return


@app.cell
def __(dat, np):
    allele_frequencies = np.sum(dat, axis=2)  # Sum over the four alleles to get total reads per mutation
    normalized_frequencies = dat / np.expand_dims(allele_frequencies, axis=2)  # Normalize allele counts
    return allele_frequencies, normalized_frequencies


@app.cell
def __(allele_frequencies):
    allele_frequencies
    return


@app.cell
def __():
    return


@app.cell
def __(mo):
    mo.md(r"""# A little example""")
    return


@app.cell
def __(mo):
    mo.md(r"""## For showing the frequency simulation""")
    return


@app.cell
def __(np, nx, simulate_mutation_freq_tree):
    np.random.seed(1)
    mf_tree = simulate_mutation_freq_tree(num_mutations=50, concentration=120, initial_mutation_freq=0.01)
    mut_freq = nx.get_node_attributes(mf_tree, 'mutation_freq')
    mut_freq
    return mf_tree, mut_freq


@app.cell
def __(mf_tree, nx, plt):
    nx.draw_networkx(mf_tree, with_labels=True, arrows=True)
    plt.show()
    return


@app.cell
def __(mo):
    mo.md(r"""## For checking the read count simulation""")
    return


@app.cell
def __(mf_tree, simulate_read_counts):
    read_counts = simulate_read_counts(mutation_tree=mf_tree, error_rate=0.001, num_reads=500, num_cells=100)
    return read_counts,


@app.cell
def __():
    return


@app.cell
def __(mo):
    mo.md(r"""# Prep data for different inference methods""")
    return


@app.cell
def __():
    return


@app.cell
def __(np):
    # Example data with the provided 2x5x4 array
    array_two_cells = np.array([
        # Cell 1
        [[458, 41, 1, 0],  # Mutation 1
         [479, 20, 0, 1],  # Mutation 2
         [500, 0, 0, 0],   # Mutation 3
         [500, 0, 0, 0],   # Mutation 4
         [500, 0, 0, 0]],  # Mutation 5
        # Cell 2
        [[458, 38, 1, 3],   # Mutation 1
         [470, 2, 18, 10],   # Mutation 2
         [500, 0, 0, 0],   # Mutation 3
         [500, 0, 0, 0],   # Mutation 4
         [500, 0, 0, 0]]   # Mutation 5
    ])

    np.save("test_data.npy", array_two_cells)
    return array_two_cells,


@app.cell
def __(mo):
    mo.md(r"""## For scite""")
    return


@app.cell
def __(array_two_cells):
    def transform_array_to_presence_absence_matrix(array):
        """
        x  into a presence/absence matrix.
        Rows represent mutations, columns represent cells.
        Presence is marked as 1 if any of X, Y, Z is greater than 0 for a mutation in a given cell.

        Parameters:
        array (numpy array): 3D numpy array where shape is (cells, mutations, states), and states are [R, X, Y, Z].

        Returns:
        numpy array: A presence/absence matrix with mutations as rows and cells as columns.
        """
        # Extract the X, Y, Z counts (columns 1, 2, 3 in the array)
        # Check if any of X, Y, Z > 0 for each mutation in each cell
        presence_absence_matrix = (array[:, :, 1:] > 0).any(axis=2).astype(int)

        return presence_absence_matrix.T  # Transpose to match mutations as rows and cells as columns


    # Apply the function to the example array
    presence_absence_matrix = transform_array_to_presence_absence_matrix(array_two_cells)
    print(presence_absence_matrix)
    return (
        presence_absence_matrix,
        transform_array_to_presence_absence_matrix,
    )


@app.cell
def __(mo):
    mo.md(r"""For MERLIN""")
    return


@app.cell
def __(np, pd):
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
    return transform_to_total_matrix, transform_to_variant_matrix


@app.cell
def __(array_two_cells, transform_to_variant_matrix):
    transform_to_variant_matrix(array_two_cells)
    return


@app.cell
def __(array_two_cells, transform_to_total_matrix):
    total_matrix_df = transform_to_total_matrix(array_two_cells)
    total_matrix_df
    return total_matrix_df,


@app.cell
def __():
    #total_matrix_df.to_csv("total_matrix.csv", header=[f"cell{i+1}" for i in range(total_matrix_df.shape[1])])
    return


@app.cell
def __():
    return


@app.cell
def __():
    return


@app.cell
def __():
    return


if __name__ == "__main__":
    app.run()
