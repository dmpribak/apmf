import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from graph_functions import *
from ff import *
import time

G = nx.DiGraph()
G.add_nodes_from([0,1,2,3,4,5])
G.add_weighted_edges_from([(0,1,3),(1,2,4),(1,3,3),(2,4,3),(3,4,4),(4,5,2)])

G2 = nx.Graph()
G2.add_nodes_from([0,1,2,3,4])
G2.add_weighted_edges_from([(0,1,1),(0,4,1),(1,2,1),(1,4,1),(2,3,1),(3,4,1)])

G3 = nx.complete_graph(20)
for u,v in G3.edges():
    set_weight(u,v,G3,1)

#show_graph(G2)
#gomory_hu(G2,plot=True)
time_graph(G3)
