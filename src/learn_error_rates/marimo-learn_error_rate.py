import marimo

__generated_with = "0.9.14"
app = marimo.App(width="medium")


@app.cell
def __():
    import marimo as mo
    import numpy as np
    import os
    import subprocess
    return mo, np, os, subprocess


@app.cell
def __(mo):
    mo.md(r"""# Exploratory coding""")
    return


@app.cell
def __(mo):
    mo.md(r"""Aim is to compare the learned trees to the true tree""")
    return


@app.cell
def __():
    input_dir = "../../results/simulated_data/10_120_0.005_500_100/"
    output_dir = "../../results/inference_output/10_120_0.005_500_100/"
    #to be defined:
    # input files for seed 1 10
    #for seed_nr in range(1, 11):
    #    generate_mutation_matrices(input_dir, output_dir)
    return input_dir, output_dir


@app.cell
def __(input_dir, os, output_dir):
    seed_nr = 1 
    input_file = os.path.join(input_dir, f"read_counts_{seed_nr}.npy")
    input_file
    output_dir_seed = os.path.join(output_dir, f"seed_{seed_nr}")
    output_dir_seed
    #os.mkdir(output_dir_seed)
    return input_file, output_dir_seed, seed_nr


@app.cell
def __():
    # run shell command
    # Note that this does not run, as the final setup in the snakemake pipeline is different
    #for error_rate in error_rates:
    #    print(error_rate)
    #    output_file = os.path.join(output_dir_seed, f"{error_rate:.6f}.csv")
    #    print(output_file)
    #    command = [
    #        "python", "../../src/simulation/adapt_data_for_mtscite.py",
    #        input_file, output_file, str(error_rate)
    #    ]
    #    print(command)

    #    try:
            # Run the command
    #        subprocess.run(command, check=True)
    #        print(f"Successfully computed mutation probability matrix given error rate {error_rate}")
    #    except subprocess.CalledProcessError as e:
            # Handle errors during script execution
    #        print(f"Error running script for error rate {error_rate}: {e}")
    return


@app.cell
def __(create_scite_mutation_dict):
    import pydot
    import networkx as nx
    import matplotlib.pyplot as plt


    def get_tree_from_gv(dot_file):

        graphs = pydot.graph_from_dot_file(dot_file)
        # pydot.graph_from_dot_file() can return multiple graphs, pick the first one
        pydot_graph = graphs[0]

        # Convert pydot graph to a NetworkX graph
        tree = nx.nx_pydot.from_pydot(pydot_graph)
        n_mutation = tree.number_of_nodes() - 1

        mutation_dict = create_scite_mutation_dict(n_mutation)
        tree_relabelled = nx.relabel_nodes(tree, mutation_dict)

        return tree_relabelled

    f = "../../results/inference_output/10_120_0.005_500_500/seed_1/learned_0.050000_0_map0.gv"
    f_true = "../../results/simulated_data/10_120_0.005_500_500/tree_1.gml"

    tree_1 = get_tree_from_gv(f)
    tree_true = nx.read_gml(f_true)
    return f, f_true, get_tree_from_gv, nx, plt, pydot, tree_1, tree_true


@app.cell
def __(tree_1):
    tree_1.number_of_nodes()
    return


@app.cell
def __(nx, plt, tree_1):
    # Visualize the graph (optional)
    nx.draw(tree_1, with_labels=True, node_size=700, node_color="lightblue")
    plt.show()
    return


@app.cell
def __(nx, plt, tree_true):
    nx.draw(tree_true, with_labels=True, node_size=700, node_color="lightblue")
    plt.show()
    return


@app.cell
def __(np):
    # Aim: recreating the plot 10_120_0.005_500_500_2_loglik.png
    # For the tree accuracy

    def compute_parent_child_distance(E1, E2):
            # Convert the lists to sets
            set_E1 = set(E1)
            set_E2 = set(E2)

            # Calculate the symmetric difference
            symmetric_difference = set_E1.symmetric_difference(set_E2)

            # Return the count of edges in the symmetric difference
            distance = len(symmetric_difference)
            return distance

            ## Example
            #e1 = [('0', '1'), ('0', '2'), ('1', '3')]
            #e2 = [('0', '1'), ('0', '2'), ('2', '3')]
            #print(compute_parent_child_distance(e1, e2))
            #print(set(e1).symmetric_difference(set(e2)))
            #return compute_parent_child_distance
    def compute_normalised_parent_child_distance(E1, E2):

        symmetric_difference = compute_parent_child_distance(E1, E2)

        maximum_possible_distance = len(E1) * 2

        normalised_distance = symmetric_difference / maximum_possible_distance

        return normalised_distance


    def compute_average_distances(true_trees, inferred_trees_list):
            # Each element in inferred_trees_list corresponds to one method's inferred trees
            avg_distances = []
            sd_distances = []

            # For each inference method
            for inferred_trees in inferred_trees_list:
                distances = []
                #Estimate the distances to the true trees.
                for i in range(0, 10):
                    true_tree = true_trees[i]
                    inferred_tree = inferred_trees[i]
                    
                    distance = compute_normalised_parent_child_distance(true_tree.edges, 
                                                                        inferred_tree.edges) 
                    #print(distance)
                    distances.append(distance)
        
                avg_distances.append(np.mean(distances))
                sd_distances.append(np.std(distances))
            #for inferred_trees in inferred_trees_list:
            #distances = []
            #    for true_tree in true_trees:
                    # Compute distance to each inferred tree for each true tree
            #        for inferred_tree in inferred_trees:
            #            distance = compute_parent_child_distance(true_tree.edges, inferred_tree.edges)  # Assuming this function is defined
            #            distances.append(distance)
            #    # Calculate average distance for this method
            #    avg_distances.append(np.mean(distances))
            #    sd_distances.append(np.std(distances))
            return [avg_distances, sd_distances]


    def create_scite_mutation_dict(n):
        # Create dictionary with node n+1 mapped to 0
        mutation_dict = {str(n+1): '0'}

        return mutation_dict

    def generate_error_rates(true_error_rate):
        """Generate 20 error rates spaced logarithmically between min and max."""
        error_rate_min = true_error_rate / 10
        error_rate_max = true_error_rate * 10
        error_rates = np.round(
            np.logspace(np.log10(error_rate_min), np.log10(error_rate_max), num=20),
            decimals=6
        )
        return error_rates
    return (
        compute_average_distances,
        compute_normalised_parent_child_distance,
        compute_parent_child_distance,
        create_scite_mutation_dict,
        generate_error_rates,
    )


@app.cell
def __(mo):
    mo.md(
        r"""
        # Systematic Approach

        For different conditions, I plot a topological distance (parent child distance) between the true tree and the tree learned during cross validation.
        """
    )
    return


@app.cell
def __(nx, os, plt):
    path_to_truth = "../../results/simulated_data/"
    true_error_rate = 0.005
    condition = f"10_120_{true_error_rate}_500_500_0.1"
    true_tree_file = os.path.join(path_to_truth, condition, "tree_1.gml")
    true_tree = nx.read_gml(true_tree_file)
    print(true_tree)
    nx.draw_networkx(true_tree, pos=nx.shell_layout(true_tree))
    plt.show()
    return (
        condition,
        path_to_truth,
        true_error_rate,
        true_tree,
        true_tree_file,
    )


@app.cell
def __(condition, generate_error_rates, os, true_error_rate):
    inference_output_dir = "../../results/inference_output/"
    seed_dir = "seed_1"
    inf_error_rates = generate_error_rates(true_error_rate)

    inferred_trees_dir = os.path.join(inference_output_dir, condition, seed_dir )

    template = "learned_[inferred_error_rate]_[rep]_map0.gv"
    inferred_trees = []
    inferred_trees_dir
    return (
        inf_error_rates,
        inference_output_dir,
        inferred_trees,
        inferred_trees_dir,
        seed_dir,
        template,
    )


@app.cell
def __(nx, plt, tree):
    nx.draw_networkx(tree, pos=nx.shell_layout(tree))
    plt.show()
    return


@app.cell
def __(
    compute_normalised_parent_child_distance,
    get_tree_from_gv,
    inf_error_rates,
    inferred_trees_dir,
    np,
    os,
    template,
    true_tree,
):
    avg_distances = []
    sd_distances = []

    for inf_error_rate in inf_error_rates:

        formatted_error_rate = f"{inf_error_rate:8f}"

        template_inf_rate = template.replace("[inferred_error_rate]", str(formatted_error_rate))
        distances = []

        repetitions = range(0,3)

        for rep in repetitions:
            template_fully_specified = template_inf_rate.replace("[rep]", str(rep))
            tree_file = os.path.join(inferred_trees_dir, template_fully_specified)

            tree = get_tree_from_gv(tree_file)
            distance = compute_normalised_parent_child_distance(tree.edges(data=False), true_tree.edges)
            distances.append(distance)

        avg_distances.append(np.mean(distances))
        sd_distances.append(np.std(distances))
    return (
        avg_distances,
        distance,
        distances,
        formatted_error_rate,
        inf_error_rate,
        rep,
        repetitions,
        sd_distances,
        template_fully_specified,
        template_inf_rate,
        tree,
        tree_file,
    )


@app.cell
def __(avg_distances, sd_distances):
    avg_distances
    sd_distances
    return


@app.cell
def __(condition, os):
    figure_file = os.path.join("plots", f"{condition}_tree_topol_dist.png")
    figure_file
    return (figure_file,)


@app.cell
def __(
    avg_distances,
    figure_file,
    inf_error_rates,
    plt,
    sd_distances,
    true_error_rate,
):
    plt.errorbar(inf_error_rates, avg_distances, yerr=sd_distances, fmt='o', capsize=5, label='Mean distances')

    # Customize plot
    plt.xscale('log')  # Use logarithmic scale for error rates
    plt.xlabel('Log Error Rates')
    plt.ylabel('Normalised Parent-Child Distance')
    plt.axvline(true_error_rate, color='green', linestyle='--', label='True error rate')

    #plt.title('Mean Distances with Error Bars')
    plt.legend()
    plt.grid(True)

    plt.savefig(figure_file, dpi=300, bbox_inches='tight')  # Save with high resolution


    # Show the plot
    plt.show()
    return


@app.cell
def __(mo):
    mo.md(
        r"""
        One more point of validation is the following. For the higher error rates, we would expect the mutation probabiliy matrix to have really low probability for any mutation.

        Maybe it would help to disentangle tree reconstruction and mutation matrix computation to compare the distance between the true and estimated mutation matrix.
        """
    )
    return


@app.cell
def __(mo):
    mo.md(r"""Normalise by score on a star tree""")
    return


@app.cell
def __():
    from networkx.drawing.nx_pydot import write_dot
    return (write_dot,)


@app.cell
def __(nx, write_dot):
    num_mutations = 10
    star_tree = nx.DiGraph()

    center_node = num_mutations + 1
    star_tree.add_edges_from((center_node, i) for i in range(1, num_mutations + 1))
    write_dot(star_tree, "star_tree.gv")
    return center_node, num_mutations, star_tree


@app.cell
def __(nx, plt, star_tree):
    nx.draw_networkx(star_tree)
    plt.show()
    return


@app.cell
def __():
    # test that scite can read the tree DONE
    # test that scite reads the tree correctly, i.e. that it interprets the node labels as correct mutations numbers
    # run mt scite with new normalisatoin
    return


if __name__ == "__main__":
    app.run()
