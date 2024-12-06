import marimo

__generated_with = "0.9.14"
app = marimo.App(width="medium")


@app.cell
def __():
    import marimo as mo
    import numpy as np
    import networkx as nx
    import os
    import matplotlib.pyplot as plt
    from Bio import Phylo
    import pandas as pd
    import plotnine
    return Phylo, mo, np, nx, os, pd, plotnine, plt


@app.cell
def __():
    # import data
    #dat_dir = "../../results/inference_output/"
    #conditions = os.listdir(dat_dir)
    #conditions
    return


@app.cell
def __():
    #condition = conditions[7]
    #print(condition)
    #path_to_condition = os.path.join(dat_dir, condition)
    return


@app.cell
def __():
    return


@app.cell
def __(mo):
    mo.md(r"""# Get all trees to the same format""")
    return


@app.cell
def __(mo):
    mo.md(r"""## Ground truth tree""")
    return


@app.cell
def __(condition, nx, os, plt):
    path_to_truth = "../../results/simulated_data/"
    true_tree_file = os.path.join(path_to_truth, condition, "tree_1.gml")
    true_tree = nx.read_gml(true_tree_file)
    print(true_tree)
    nx.draw_networkx(true_tree, pos=nx.shell_layout(true_tree))
    plt.show()
    return path_to_truth, true_tree, true_tree_file


@app.cell
def __(true_tree):
    print(true_tree.edges)
    return


@app.cell
def __(mo):
    mo.md(r"""Merlin input""")
    return


@app.cell
def __(create_merlin_mutation_dict, nx, os, path_to_condition):
    # import merlin tree
    merlin_graph = nx.read_edgelist(os.path.join(path_to_condition, "merlin_1_clone_tree_edge_list.txt"), 
                                   create_using = nx.DiGraph, delimiter=", ")
    mutation_dict_merlin = create_merlin_mutation_dict(50)
    merlin_graph = nx.relabel_nodes(merlin_graph, mapping=mutation_dict_merlin)
    #print(mutation_dict_merlin)
    print(merlin_graph.edges)
    return merlin_graph, mutation_dict_merlin


@app.cell
def __(merlin_graph, nx, plt):
    nx.draw_networkx(merlin_graph, pos=nx.shell_layout(merlin_graph))
    plt.show()
    return


@app.cell
def __(mo):
    mo.md(r"""## mt scite tree""")
    return


@app.cell
def __():
    return


@app.cell
def __():
    return


@app.cell
def __():
    #from Bio import Phylo
    #import networkx as nx
    from typing import Optional, Union
    return Optional, Union


@app.cell
def __(Optional, Phylo, nx):
    def create_scite_mutation_dict(n):
        # Create dictionary with node n+1 mapped to 0
        mutation_dict = {str(n+1): '0'}

        return mutation_dict

    def get_node_label(clade: Phylo.BaseTree.Clade) -> str:
            """
            Get the label for a clade, assuming 'name' for tips and 'confidence' for internal nodes. 
            This is 

            Args:
                clade (Phylo.BaseTree.Clade): A clade in the phylogenetic tree.

            Returns:
                str: The label for the clade.
            """
            return str(clade.name or clade.confidence or f"Node_{id(clade)}")

    def add_edge_with_attributes(graph: nx.DiGraph, parent_label: str, child_label: str, clade: Phylo.BaseTree.Clade):
            """
            Add an edge with branch attributes between parent and child nodes in the graph.

            Args:
                graph (nx.DiGraph): The NetworkX graph.
                parent_label (str): The label for the parent node.
                child_label (str): The label for the child node.
                clade (Phylo.BaseTree.Clade): The clade corresponding to the child node.
            """
            edge_data = {'weight': clade.branch_length or 1.0}
            if hasattr(clade, "color") and clade.color:
                edge_data['color'] = clade.color.to_hex()
            if hasattr(clade, "width") and clade.width is not None:
                edge_data['width'] = clade.width
            graph.add_edge(parent_label, child_label, **edge_data)

    def build_graph(clade: Phylo.BaseTree.Clade, graph: nx.DiGraph, 
                    parent_label: Optional[str] = None):
            """
            Recursively build the graph from the tree structure.

            Args:
                clade (Phylo.BaseTree.Clade): The current clade being processed.
                parent_label (Optional[str]): The label of the parent node. None for the root.
            """
            clade_label = get_node_label(clade)
            graph.add_node(clade_label)

            if parent_label:
                add_edge_with_attributes(graph, parent_label, clade_label, clade)

            for child in clade.clades:
                build_graph(child, graph, clade_label)


    def tree_to_networkx(tree: Phylo.BaseTree.Tree) -> nx.DiGraph:
        """
        Convert a Biopython Tree object to a NetworkX graph, including branch attributes.

        This function uses 'name' for tip labels and 'confidence' for internal node labels,
        if 'name' is not available. Additionally, it includes branch attributes like 'weight',
        'color', and 'width' where available.

        Args:
            tree (Phylo.BaseTree.Tree): A Biopython Tree object.

        Returns:
            Union[nx.DiGraph, nx.Graph]: A directed or undirected NetworkX graph, depending
                on whether the input tree is rooted.

        Raises:
            ValueError: If the input tree is not a valid Biopython Tree object.
        """
        if not isinstance(tree, Phylo.BaseTree.Tree):
            raise ValueError("Input must be a Biopython Phylo Tree object")

        # Create a directed graph 
        graph = nx.DiGraph() 

        # Start building the graph from the root
        root_label = get_node_label(tree.root)
        build_graph(clade=tree.root, graph=graph)

        return graph
    return (
        add_edge_with_attributes,
        build_graph,
        create_scite_mutation_dict,
        get_node_label,
        tree_to_networkx,
    )


@app.cell
def __(Phylo, os, path_to_condition):
    # import mt scite tree
    mt_scite_tree = Phylo.read(os.path.join(path_to_condition, "mt_scite_mutation_prob_0.0005_1_map0.newick"), format="newick", rooted=True)
    Phylo.draw(mt_scite_tree)
    print(mt_scite_tree)
    return (mt_scite_tree,)


@app.cell
def __():
    return


@app.cell
def __(create_scite_mutation_dict, mt_scite_tree, nx, tree_to_networkx):
    mt_scite_graph = tree_to_networkx(mt_scite_tree)
    mt_scite_dict = create_scite_mutation_dict(50)
    mt_scite_graph = nx.relabel_nodes(mt_scite_graph, mt_scite_dict)
    return mt_scite_dict, mt_scite_graph


@app.cell
def __(mt_scite_graph, nx, plt):
    nx.draw(mt_scite_graph, with_labels=True)
    plt.show()
    return


@app.cell
def __(mo):
    mo.md(r"""## Scite Tree""")
    return


@app.cell
def __(Phylo, os, path_to_condition):
    scite_tree = Phylo.read(os.path.join(path_to_condition, "scite_0.005_1_ml0.newick"), format="newick", rooted=True)
    Phylo.draw(scite_tree)
    print(scite_tree)
    return (scite_tree,)


@app.cell
def __(create_scite_mutation_dict, nx, scite_tree, tree_to_networkx):
    scite_graph = tree_to_networkx(scite_tree)
    scite_dict = create_scite_mutation_dict(10)
    scite_graph = nx.relabel_nodes(scite_graph, scite_dict)
    return scite_dict, scite_graph


@app.cell
def __(nx, plt, scite_graph):
    nx.draw_networkx(scite_graph)
    plt.show()
    return


@app.cell
def __(mo):
    mo.md(r"""#Compare trees using parent child distance""")
    return


@app.cell
def __(scite_graph, true_tree):
    print(true_tree.edges)
    print(scite_graph.edges)
    return


@app.cell
def __():
    def compute_parent_child_distance(E1, E2):
        # Convert the lists to sets
        set_E1 = set(E1)
        set_E2 = set(E2)

        # Calculate the symmetric difference
        symmetric_difference = set_E1.symmetric_difference(set_E2)

        # Return the count of edges in the symmetric difference
        distance = len(symmetric_difference)
        return distance

    def compute_normalised_parent_child_distance(E1, E2):

        symmetric_difference = compute_parent_child_distance(E1, E2)

        maximum_possible_distance = len(E1) * 2

        normalised_distance = symmetric_difference / maximum_possible_distance

        return normalised_distance


    ## Example
    e1 = [('0', '1'), ('0', '2'), ('1', '3')]
    e2 = [('0', '1'), ('0', '2'), ('2', '3')]
    print(compute_parent_child_distance(e1, e2))
    print(compute_normalised_parent_child_distance(e1, e2))
    print(set(e1).symmetric_difference(set(e2)))
    return (
        compute_normalised_parent_child_distance,
        compute_parent_child_distance,
        e1,
        e2,
    )


@app.cell
def __(compute_parent_child_distance, merlin_graph, true_tree):
    compute_parent_child_distance(true_tree.edges, merlin_graph.edges)
    return


@app.cell
def __(compute_parent_child_distance, mt_scite_graph, true_tree):
    compute_parent_child_distance(true_tree.edges, mt_scite_graph.edges)
    return


@app.cell
def __(compute_parent_child_distance, scite_graph, true_tree):
    compute_parent_child_distance(true_tree.edges, scite_graph.edges)
    return


@app.cell
def __(mo):
    mo.md(
        r"""
        However, the trees are not yet matched!!

        """
    )
    return


@app.cell
def __(mo):
    mo.md(r"""#Plot results""")
    return


@app.cell
def __(path_to_condition):
    path_to_condition
    return


@app.cell
def __():
    def create_merlin_mutation_dict(n):
        # Create dictionary with 'root' mapped to 0
        mutation_dict = {"root": 0}

        # Add mappings for each mutation
        for i in range(1, n + 1):
            mutation_dict[f"mut{i}"] = f"{i}"

        return mutation_dict
    return (create_merlin_mutation_dict,)


@app.cell
def __(
    Phylo,
    compute_normalised_parent_child_distance,
    create_merlin_mutation_dict,
    create_scite_mutation_dict,
    mutation_dict_merlin,
    np,
    nx,
    tree_to_networkx,
):
    # Define parsing functions for each method
    def parse_merlin(file_path, num_mutations):

        # Read graph
        merlin_graph = nx.read_edgelist(file_path, create_using = nx.DiGraph,
                                        delimiter=", ")

        # Remap node labels such that they are the same as in the true tree
        mutation_dict = create_merlin_mutation_dict(num_mutations)
        merlin_graph = nx.relabel_nodes(merlin_graph, mapping=mutation_dict_merlin, copy=False)

        return merlin_graph  

    def parse_scite_mtscite(file_path, num_mutations):
        # read tree from newick file
        scite_tree = Phylo.read(file_path, format="newick", rooted=True)

        # transform to network to enable comparison to ground truth tree
        scite_graph = tree_to_networkx(scite_tree)

        # Relabel nodes such that they are the same as in the true tree
        scite_dict = create_scite_mutation_dict(num_mutations)
        scite_graph = nx.relabel_nodes(scite_graph, scite_dict)

        return scite_graph


    def parse_true(file_path, num_mutations):

        true_tree = nx.read_gml(file_path)

        return true_tree

    # Dictionary mapping method names to parsing functions
    parsing_functions = {
        "merlin": parse_merlin,
        "mtscite": parse_scite_mtscite,
        "scite": parse_scite_mtscite,
        "true": parse_true
    }

    # Generalized function to load trees for a given method
    def load_trees(path_template, seed_range, method, num_mutations=10):
        trees = []
        # Get the parsing function based on the method
        parsing_function = parsing_functions.get(method)

        if parsing_function is None:
            raise ValueError(f"Unknown method: {method}")

        for seed in seed_range:
            # Replace [seed] in the path template with the actual seed number
            file_path = path_template.replace("[seed]", str(seed))

            # Parse the tree file and convert to NetworkX graph using the chosen parsing function
            tree_graph = parsing_function(file_path, num_mutations)
            trees.append(tree_graph)

        return trees

    def compute_average_distances(true_trees, inferred_trees_list):
        # Each element in inferred_trees_list corresponds to one method's inferred trees
        avg_distances = []
        sd_distances = []

        for inferred_trees in inferred_trees_list:
            distances = []
            for true_tree in true_trees:
                # Compute distance to each inferred tree for each true tree
                for inferred_tree in inferred_trees:
                    distance = compute_normalised_parent_child_distance(true_tree.edges, inferred_tree.edges)  
                    distances.append(distance)
            # Calculate average distance for this method
            avg_distances.append(np.mean(distances))
            sd_distances.append(np.std(distances))
        return [avg_distances, sd_distances]
    return (
        compute_average_distances,
        load_trees,
        parse_merlin,
        parse_scite_mtscite,
        parse_true,
        parsing_functions,
    )


@app.cell
def __():
    return


@app.cell
def __():
    return


@app.cell
def __(os):
    true_error_rate = 0.005
    inferrered_error_rates = [0.05, 0.005, 0.0005]
    condition = f"10_120_{true_error_rate}_500_500_0.1"

    # Define file path templates for ground truth trees
    simulation_output_dir = "../../results/simulated_data/"
    path_to_true_trees = os.path.join(simulation_output_dir, condition)
    true_path_template = os.path.join(path_to_true_trees, "tree_[seed].gml")

    # Define file path templates for each method
    dat_dir = "../../results/inference_output/"
    path_to_condition = os.path.join(dat_dir, condition)

    scite_path_template = os.path.join(path_to_condition, f"scite_{true_error_rate}_[seed]_ml0.newick")

    mtscite_path_template_0 = os.path.join(path_to_condition, f"mt_scite_mutation_prob_{inferrered_error_rates[0]}_[seed]_map0.newick")

    mtscite_path_template_1 = os.path.join(path_to_condition, f"mt_scite_mutation_prob_{inferrered_error_rates[1]}_[seed]_map0.newick")

    mtscite_path_template_2 = os.path.join(path_to_condition, f"mt_scite_mutation_prob_{inferrered_error_rates[2]}_[seed]_map0.newick")

    merlin_path_template = os.path.join(path_to_condition, "merlin_[seed]_clone_tree_edge_list.txt")
    return (
        condition,
        dat_dir,
        inferrered_error_rates,
        merlin_path_template,
        mtscite_path_template_0,
        mtscite_path_template_1,
        mtscite_path_template_2,
        path_to_condition,
        path_to_true_trees,
        scite_path_template,
        simulation_output_dir,
        true_error_rate,
        true_path_template,
    )


@app.cell
def __(true_path_template, true_trees):
    true_path_template
    true_trees
    return


@app.cell
def __(
    load_trees,
    merlin_path_template,
    mtscite_path_template_0,
    mtscite_path_template_1,
    mtscite_path_template_2,
    scite_path_template,
    true_path_template,
):
    # Load trees for each method using the generalized function
    seed_range = range(1, 11)
    scite_trees = load_trees(scite_path_template, seed_range, "scite")
    mtscite_trees_0 = load_trees(mtscite_path_template_0, seed_range, "mtscite")
    mtscite_trees_1 = load_trees(mtscite_path_template_1, seed_range, "mtscite")
    mtscite_trees_2 = load_trees(mtscite_path_template_2, seed_range, "mtscite")
    true_trees = load_trees(true_path_template, seed_range, "true")
    merlin_trees = load_trees(merlin_path_template, seed_range, "merlin")
    return (
        merlin_trees,
        mtscite_trees_0,
        mtscite_trees_1,
        mtscite_trees_2,
        scite_trees,
        seed_range,
        true_trees,
    )


@app.cell
def __(merlin_trees):
    merlin_trees[0].edges

    return


@app.cell
def __(true_trees):
    true_trees[0].edges
    return


@app.cell
def __(mtscite_trees_0):
    mtscite_trees_0[0].edges
    return


@app.cell
def __(scite_trees):
    scite_trees[0].edges
    return


@app.cell
def __():
    #list(set(merlin_trees[0].edges) & set(true_trees[0].edges))
    return


@app.cell
def __(compute_parent_child_distance, merlin_trees, true_trees):
    compute_parent_child_distance(merlin_trees[0].edges, true_trees[0].edges)
    return


@app.cell
def __(merlin_trees, nx, plt):
    nx.draw_networkx(merlin_trees[0])
    plt.show()
    return


@app.cell
def __(nx, plt, true_trees):
    nx.draw_networkx(true_trees[0])
    plt.show()
    return


@app.cell
def __(
    compute_average_distances,
    inferrered_error_rates,
    merlin_trees,
    mtscite_trees_0,
    mtscite_trees_1,
    mtscite_trees_2,
    np,
    pd,
    scite_trees,
    true_error_rate,
    true_trees,
):
    # Assume true_trees, merlin_trees, scite_trees, and mtscite_trees are available
    # inferred_trees_list contains trees from each method
    #merlin_trees = [merlin_graph]
    #scite_trees = [scite_graph]
    #mtscite_trees = [mt_scite_graph]
    #true_trees = [true_tree]

    inferred_trees_list = [merlin_trees, scite_trees, mtscite_trees_0, mtscite_trees_1, mtscite_trees_2]
    method_names = ['Merlin', 'SCITE', 'mtSCITE', 'mtSCITE', 'mtSCITE']

    error_rates = [np.nan, true_error_rate] +  inferrered_error_rates

    # Step 2: Calculate average distances for each method
    avg_distances, sd_distances = compute_average_distances(true_trees, inferred_trees_list)

    data = {
        "Method":method_names,
        "Error_Rate": error_rates,
        "Avg_Distance": avg_distances,
        "Sd_Distance": sd_distances
    }
    df = pd.DataFrame(data)
    return (
        avg_distances,
        data,
        df,
        error_rates,
        inferred_trees_list,
        method_names,
        sd_distances,
    )


@app.cell
def __():
    return


@app.cell
def __():
    from plotnine import ggplot, aes, geom_bar, geom_errorbar, labs, theme, position_dodge, scale_fill_manual
    return (
        aes,
        geom_bar,
        geom_errorbar,
        ggplot,
        labs,
        position_dodge,
        scale_fill_manual,
        theme,
    )


@app.cell
def __(
    aes,
    condition,
    df,
    geom_bar,
    geom_errorbar,
    ggplot,
    labs,
    position_dodge,
    scale_fill_manual,
):
    plot_2 = (
        ggplot(df, aes(x="Method", y="Avg_Distance", 
                       fill= "factor(Error_Rate)",
                       group= "factor(Error_Rate)"))+
        geom_bar(stat = "identity", width=0.7, position = position_dodge(0.7))+

        geom_errorbar(aes(ymin="Avg_Distance - 0.5 * Sd_Distance",
                          ymax="Avg_Distance + 0.5 * Sd_Distance"), 
                      width=0.2, position=position_dodge(0.7))+

        labs(
            x="Method",
            y="Average Parent-Child Distance",
            fill="Fixed Error rates"
        )
        + scale_fill_manual(values=["#999999", "#56B4E9", "#E69F00", "#009E73"])


    )
    plot_2.show()
    plot_2.save(f"../../results/plots/{condition}_fixedErrors_tree_topology_comparison.png", width=10, height=6, dpi=300)
    return (plot_2,)


@app.cell
def __(condition):
    condition
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
