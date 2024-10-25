import marimo

__generated_with = "0.8.14"
app = marimo.App(width="medium")


@app.cell
def __(__file__):
    import marimo as mo
    import pytest
    import os
    import sys
    import pytest
    import networkx as nx
    import numpy as np
    import matplotlib.pyplot as plt

    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src/simulation')))

    from functions import reparameterize_beta, simulate_mutation_freq_tree, simulate_read_counts
    return (
        mo,
        np,
        nx,
        os,
        plt,
        pytest,
        reparameterize_beta,
        simulate_mutation_freq_tree,
        simulate_read_counts,
        sys,
    )


@app.cell
def __(mo):
    mo.md(r"""This marimo python file exists for interactive checking of the test functions""")
    return


@app.cell
def __(mo):
    mo.md(r"""Test beta distribution""")
    return


@app.cell
def __(pytest, reparameterize_beta):
    def test_reparameterize_beta_valid_input():
        # Test with a valid mode and concentration
        mode = 0.5
        concentration = 5
        alpha, beta = reparameterize_beta(mode, concentration)
        assert alpha == pytest.approx(2.5, rel=1e-2)
        assert beta == pytest.approx(2.5, rel=1e-2)

    def test_reparameterize_beta_invalid_mode():
        # Test with invalid mode (outside 0 to 1 range)
        with pytest.raises(ValueError, match="Mode must be between 0 and 1."):
            reparameterize_beta(-0.1, 5)

        with pytest.raises(ValueError, match="Mode must be between 0 and 1."):
            reparameterize_beta(1.1, 5)

    def test_reparameterize_beta_invalid_concentration():
        # Test with invalid concentration (less than or equal to 2)
        with pytest.raises(ValueError, match="Concentration must be greater than 2 to be valid for this parameterization."):
            reparameterize_beta(0.5, 2)

        with pytest.raises(ValueError, match="Concentration must be greater than 2 to be valid for this parameterization."):
            reparameterize_beta(0.5, 1.5)

    def test_reparameterize_beta_edge_case():
        # Test edge cases with mode near 0 and 1
        mode = 0.01
        concentration = 10
        alpha, beta = reparameterize_beta(mode, concentration)
        assert alpha == pytest.approx(1.08, rel=1e-2)
        assert beta == pytest.approx(8.92, rel=1e-2)

        mode = 0.99
        alpha, beta = reparameterize_beta(mode, concentration)
        assert alpha == pytest.approx(8.92, rel=1e-2)
        assert beta == pytest.approx(1.08, rel=1e-2)
    return (
        test_reparameterize_beta_edge_case,
        test_reparameterize_beta_invalid_concentration,
        test_reparameterize_beta_invalid_mode,
        test_reparameterize_beta_valid_input,
    )


@app.cell
def __(mo):
    mo.md(r"""Test mutation tree simulation""")
    return


@app.cell
def __(np, nx, pytest, simulate_mutation_freq_tree):
    def test_simulate_mutation_freq_tree_valid_input():
        np.random.seed(42)  # For reproducibility
        num_mutations = 5
        concentration = 10
        initial_mutation_freq = 0.01

        tree = simulate_mutation_freq_tree(num_mutations, concentration, initial_mutation_freq)

        # Test that the tree has the correct number of nodes and edges
        assert len(tree.nodes) == num_mutations + 1  # Root node + num_mutations
        assert len(tree.edges) == num_mutations      # One edge per mutation

        # Test that each node has the correct attributes
        for node in tree.nodes:
            mutation_freqs = tree.nodes[node]['mutation_freq']
            assert isinstance(mutation_freqs, list)
            assert len(mutation_freqs) == num_mutations
            # Frequencies should be between 0 and 1
            for freq in mutation_freqs:
                assert 0 <= freq <= 1

    def test_simulate_mutation_freq_tree_invalid_num_mutations():
        concentration = 10
        with pytest.raises(ValueError, match="num_mutations must be a positive integer."):
            simulate_mutation_freq_tree(0, concentration)
        with pytest.raises(ValueError, match="num_mutations must be a positive integer."):
            simulate_mutation_freq_tree(-5, concentration)
        with pytest.raises(ValueError, match="num_mutations must be a positive integer."):
            simulate_mutation_freq_tree(3.5, concentration)

    def test_simulate_mutation_freq_tree_mutation_freqs_in_range():
        np.random.seed(42)
        num_mutations = 3
        concentration = 5
        tree = simulate_mutation_freq_tree(num_mutations, concentration)
        for node in tree.nodes:
            mutation_freqs = tree.nodes[node]['mutation_freq']
            for freq in mutation_freqs:
                assert 0 <= freq <= 1

    def test_simulate_mutation_freq_tree_tree_structure():
        num_mutations = 10
        concentration = 5
        tree = simulate_mutation_freq_tree(num_mutations, concentration)
        # Check if the graph is connected and acyclic (i.e., a tree)
        assert nx.is_tree(tree)
        # Check that there are no isolated nodes
        assert nx.number_connected_components(tree.to_undirected()) == 1

    # HEERE 
    def test_simulate_mutation_freq_tree_initial_mutation_freq_effect():
        np.random.seed(42)
        num_mutations = 5
        concentration = 100
        initial_mutation_freq = 0.5

        tree = simulate_mutation_freq_tree(num_mutations, concentration, initial_mutation_freq)

        # Check that the new mutations are centered around initial_mutation_freq
        for node in tree.nodes:
            mutation_freqs = tree.nodes[node]['mutation_freq']
            # The frequencies of new mutations should be around initial_mutation_freq
            new_mutation_freq = mutation_freqs[node - 1] if node != 0 else None
            if new_mutation_freq is not None and new_mutation_freq != 0:
                assert pytest.approx(new_mutation_freq, 0.1) == initial_mutation_freq

    ## THis could be used for choosing a value for the concentration parameter!!
    def test_simulate_mutation_freq_tree_similar_parent_child_freqs():
        np.random.seed(42)
        num_mutations = 10
        concentration = 500
        initial_mutation_freq = 0.5

        tree = simulate_mutation_freq_tree(num_mutations, concentration, initial_mutation_freq)

        # Compare the mutation frequencies between parent and child nodes
        for node in tree.nodes:
            if node == 0:
                continue  # Skip the root node

            parent = list(tree.predecessors(node))[0]  # Get the parent node
            parent_freqs = tree.nodes[parent]['mutation_freq']
            node_freqs = tree.nodes[node]['mutation_freq']

            # Check if the mutation frequencies are similar for existing mutations between parent and child
            for i in range(num_mutations):
                if parent_freqs[i] == 0 or parent_freqs[i] == 1:
                    # Skip mutations that haven't been initialized or are fixed
                    continue
                assert pytest.approx(node_freqs[i], 0.1) == parent_freqs[i], \
                    f"Mutation frequency of node {node} for mutation {i} does not match parent node {parent}"


    def test_simulate_mutation_freq_tree_random_seed_consistency():
        num_mutations = 5
        concentration = 5
        initial_mutation_freq = 0.01

        np.random.seed(42)
        tree1 = simulate_mutation_freq_tree(num_mutations, concentration, initial_mutation_freq)

        np.random.seed(42)
        tree2 = simulate_mutation_freq_tree(num_mutations, concentration, initial_mutation_freq)

        # Compare the mutation frequencies of the two trees
        for node in tree1.nodes:
            freqs1 = tree1.nodes[node]['mutation_freq']
            freqs2 = tree2.nodes[node]['mutation_freq']
            assert freqs1 == freqs2

    def test_simulate_mutation_freq_tree_large_num_mutations():
        num_mutations = 100
        concentration = 10
        tree = simulate_mutation_freq_tree(num_mutations, concentration)
        assert len(tree.nodes) == num_mutations + 1
        assert len(tree.edges) == num_mutations
    return (
        test_simulate_mutation_freq_tree_initial_mutation_freq_effect,
        test_simulate_mutation_freq_tree_invalid_num_mutations,
        test_simulate_mutation_freq_tree_large_num_mutations,
        test_simulate_mutation_freq_tree_mutation_freqs_in_range,
        test_simulate_mutation_freq_tree_random_seed_consistency,
        test_simulate_mutation_freq_tree_similar_parent_child_freqs,
        test_simulate_mutation_freq_tree_tree_structure,
        test_simulate_mutation_freq_tree_valid_input,
    )


@app.cell
def __():
    return


@app.cell
def __(np, plt):
    from scipy.stats import beta

    def beta_quantiles(mode, concentration):
        # Compute alpha and beta from mode and concentration
        alpha = mode * (concentration - 2) + 1
        beta_param = (1 - mode) * (concentration - 2) + 1

        # Calculate the 5% and 95% quantiles
        quantile_5 = beta.ppf(0.05, alpha, beta_param)
        quantile_95 = beta.ppf(0.95, alpha, beta_param)

        return quantile_5, quantile_95

    def plot_beta_pdf(mode, concentration):
        """
        Plot the PDF of a Beta distribution for given alpha and beta parameters.

        Parameters:
        alpha (float): The alpha shape parameter of the Beta distribution.
        beta_param (float): The beta shape parameter of the Beta distribution.
        """
        alpha = mode * (concentration - 2) + 1
        beta_param = (1 - mode) * (concentration - 2) + 1

        # Generate values between 0 and 1
        x = np.linspace(0, 1, 1000)
        # Compute the Beta distribution PDF
        y = beta.pdf(x, alpha, beta_param)

        # Plot the PDF
        plt.figure(figsize=(8, 6))
        plt.plot(x, y, label=f'Beta PDF (α={alpha}, β={beta_param})', color='blue')
        plt.title('Beta Distribution PDF')
        plt.xlabel('x')
        plt.ylabel('Density')
        plt.legend()
        plt.grid(True)
        plt.show()
    return beta, beta_quantiles, plot_beta_pdf


@app.cell
def __(mo):
    mo.md(r"""Example for centered beta""")
    return


@app.cell
def __(beta_quantiles):
    # Example for centered concentration
    mode_centered = 0.01
    concentration_centered = 120
    q5_centered, q95_centered = beta_quantiles(mode_centered, concentration_centered)
    print(f"5% quantile: {q5_centered}")
    print(f"95% quantile: {q95_centered}")
    return concentration_centered, mode_centered, q5_centered, q95_centered


@app.cell
def __(concentration_centered, mode_centered, plot_beta_pdf):
    plot_beta_pdf(mode_centered, concentration_centered)
    return


@app.cell
def __(mo):
    mo.md(r"""Example for uncetered beta""")
    return


@app.cell
def __(beta_quantiles):
    # Example for uncentered concentration
    mode_uncentered = 0.01
    concentration_uncentered = 30
    q5_uncentered, q95_uncentered = beta_quantiles(mode_uncentered, concentration_uncentered)
    print(f"5% quantile: {q5_uncentered}")
    print(f"95% quantile: {q95_uncentered}")
    return (
        concentration_uncentered,
        mode_uncentered,
        q5_uncentered,
        q95_uncentered,
    )


@app.cell
def __(concentration_centered, mode_centered, plot_beta_pdf):
    plot_beta_pdf(mode_centered, concentration_centered)
    return


@app.cell
def __(mo):
    mo.md(r"""Test read count simulation""")
    return


@app.cell
def __(np, nx, pytest, simulate_read_counts):
    # Example mutation tree setup
    def create_dummy_mutation_tree():
        # Create a simple mutation tree
        tree = nx.DiGraph()

        # Add nodes with mutation frequencies for testing
        tree.add_node(0, mutation_freq=[0, 0])  # Root node
        tree.add_node(1, mutation_freq=[0.5, 0])  # Child node
        tree.add_node(2, mutation_freq=[0.5, 0.3])  # Child node

        tree.add_edge(0, 1)
        tree.add_edge(1, 2)

        return tree

    # Test the simulate_read_counts function
    def test_simulate_read_counts():
        np.random.seed(42)  # Set seed for reproducibility

        # Create a dummy mutation tree with 2 mutations
        mutation_tree = create_dummy_mutation_tree()

        # Define parameters
        error_rate = 0.01
        num_reads = 100
        num_cells = 5

        # Call the function
        read_counts = simulate_read_counts(mutation_tree, error_rate, num_reads, num_cells)

        # Check that the shape of the read_counts tensor is correct
        assert read_counts.shape == (num_cells, 2, 4), "Unexpected shape of read_counts tensor."

        # Check that the allele counts sum up to the number of reads for each site
        for cell in range(num_cells):
            for mutation_nr in range(2):  # Assuming 2 mutations in the tree
                assert read_counts[cell, mutation_nr].sum() == num_reads, f"Allele counts for cell {cell}, mutation {mutation_nr} do not sum to {num_reads}."

        # Check that the counts are integers
        assert np.issubdtype(read_counts.dtype, np.integer), "Allele counts should be integers."

    # Additional test to check invalid probability handling
    def test_invalid_probability_handling():
        np.random.seed(42)

        # Create a dummy mutation tree
        mutation_tree = create_dummy_mutation_tree()

        # Introduce an invalid mutation frequency
        mutation_tree.nodes[1]['mutation_freq'][0] = 1.5  # Invalid mutation frequency

        error_rate = 0.01
        num_reads = 100
        num_cells = 5

        # Expecting a RuntimeError for invalid probabilities
        with pytest.raises(RuntimeError, match="Invalid probability p0="):
            simulate_read_counts(mutation_tree, error_rate, num_reads, num_cells)
    return (
        create_dummy_mutation_tree,
        test_invalid_probability_handling,
        test_simulate_read_counts,
    )


if __name__ == "__main__":
    app.run()
