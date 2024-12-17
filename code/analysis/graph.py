import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
from itertools import product
from typing import Iterable, Callable
import os


def make_graph(vertices: Iterable, is_edge: Callable) -> nx.Graph:
    """given a set a verticies and a decriptor function, construct a graph
    - vertices:Iterable, the set of vertices of the graph,
    -is_edge:Callable, take 2 vertices.
        Return a bool if there is an edge and a dictionnary of metadata
    """
    graph = nx.Graph()
    graph.add_nodes_from(vertices)
    for a, b in product(vertices, vertices):
        validation, metadata = is_edge(a, b)
        if validation:
            graph.add_edge(a, b, **metadata)
    return graph


def dumy_dataset(source: str) -> nx.Graph:
    graph = nx.Graph()
    source = pd.read_csv(source, delimiter='\t\t')
    graph.add_edges_from(
        ((index.Sp_sending, index.Sp_receiving)
            for index in graph.intertuples())
    )
    return graph


def analyse_graph(
    G: nx.Graph,
    plot_file: str = "graph_plot.png",
    graphml_file: str = "graph.graphml",
    degree_dist_file: str = "degree_distribution.png"
) -> None:
    """
    Analyze and visualize the graph `G`. This function performs these tasks:

    1. Plots the graph and saves it to a PNG file.
    2. Saves the graph in GraphML format.
    3. Plots the degree distribution of the graph and saves it to a PNG file.

    Parameters:
    G (networkx.Graph): The graph to analyze.
    plot_file (str): The filename to save the graph plot.
    graphml_file (str): The filename to save the graph.
    degree_dist_file (str): The filename to save the degree distribution plot.

    Returns:
    None
    """

    # Ensure the result directory exists
    os.makedirs("result", exist_ok=True)
    plot_file = os.path.join("result", plot_file)
    graphml_file = os.path.join("result", graphml_file)
    degree_dist_file = os.path.join("result", degree_dist_file)

    # Plot the graph and save it to a file using the object-oriented API
    fig, ax = plt.subplots(figsize=(8, 6))  # Create a figure and axis
    nx.draw(
        G,
        with_labels=True,
        node_color="skyblue",
        font_weight="bold",
        node_size=500,
        ax=ax,
    )
    ax.set_title("Graph Visualization")
    fig.savefig(plot_file)
    plt.close(fig)  # Close the figure to release memory

    # Save the graph in GraphML format
    nx.write_graphml(G, graphml_file)

    # Plot the degree distribution using the object-oriented API
    degrees = [G.degree(n) for n in G.nodes()]  # Get the degree of each node
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.hist(
        degrees, bins=range(min(degrees), max(degrees) + 1), edgecolor="black"
    )
    ax.set_title("Degree Distribution")
    ax.set_xlabel("Degree")
    ax.set_ylabel("Frequency")
    fig.savefig(degree_dist_file)
    plt.close(fig)  # Close the figure to release memory
