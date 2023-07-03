import networkx as nx
import numpy as np
import string
G = nx.from_numpy_matrix(X.values)
G = nx.relabel_nodes(G, dict(zip(range(len(G.nodes())),string.ascii_uppercase)))
G = nx.drawing.nx_agraph.to_agraph(G)
!conda install -c conda-forge pygraphviz
G = nx.drawing.nx_agraph.to_agraph(G)
G.node_attr.update(color="red", style="filled")
G.edge_attr.update(color="blue", width="2.0")
