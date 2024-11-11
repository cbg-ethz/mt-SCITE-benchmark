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
    return Phylo, mo, np, nx, os, plt


@app.cell
def __(os):
    # import data
    dat_dir = "../../results/inference_output/"
    conditions = os.listdir(dat_dir)
    conditions
    return conditions, dat_dir


@app.cell
def __(conditions, dat_dir, os):
    condition = conditions[7]
    condition
    path_to_condition = os.path.join(dat_dir, condition)

    return condition, path_to_condition


@app.cell
def __(mo):
    mo.md(
        r"""
        Merlin input

        """
    )
    return


@app.cell
def __(nx, os, path_to_condition):
    # import merlin tree
    merlin_tree = nx.read_edgelist(os.path.join(path_to_condition, "merlin_1_clone_tree_edge_list.txt"), create_using = nx.DiGraph, delimiter=", ")

    print(merlin_tree.edges)

    return (merlin_tree,)


@app.cell
def __(merlin_tree, nx, plt):
    nx.draw_networkx(merlin_tree, pos=nx.shell_layout(merlin_tree))
    plt.show()
    return


@app.cell
def __(mo):
    mo.md(r"""mt scite tree""")
    return


@app.cell
def __():
    #from Bio import Phylo
    #import networkx as nx
    from typing import Optional, Union
    return Optional, Union


@app.cell
def __(Optional, Phylo, nx):



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

        def build_graph(clade: Phylo.BaseTree.Clade, parent_label: Optional[str] = None):
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
                build_graph(child, clade_label)

        # Start building the graph from the root
        root_label = get_node_label(tree.root)
        graph.add_node(root_label)
        build_graph(tree.root, root_label)

        return graph

    return (tree_to_networkx,)


@app.cell
def __(Phylo, os, path_to_condition):
    # import mt scite tree
    mt_scite_tree = Phylo.read(os.path.join(path_to_condition, "mt_scite_mutation_prob_0.0001_1_map0.newick"), format="newick", rooted=True)
    Phylo.draw(mt_scite_tree)
    print(mt_scite_tree)
    return (mt_scite_tree,)


@app.cell
def __(mt_scite_tree, tree_to_networkx):
    mt_scite_graph = tree_to_networkx(mt_scite_tree)
    return (mt_scite_graph,)


@app.cell
def __(mt_scite_graph, nx, plt):
    nx.draw(mt_scite_graph, with_labels=True)
    plt.show()
    return


if __name__ == "__main__":
    app.run()
