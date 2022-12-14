import networkx as nx
import numpy as np
from graph_functions import * 
import matplotlib.pyplot as plt
import time

def ford_fulkerson(s,t,G):
    Gf = init_residual(G)
    path = dfs(s, t, Gf)
    f = 0
    while(path):
        min_weight = weight(*path[0],Gf)
        for u,v in path:
            if weight(u,v,Gf) < min_weight:
                min_weight = weight(u,v,Gf)
        for u,v in path:
            w = weight(u,v,Gf)
            w2 = weight(v,u,Gf)
            set_weight(u,v,Gf,w-min_weight)
            set_weight(v,u,Gf,w2+min_weight)
        f += min_weight
        path = dfs(s,t,Gf)
    return f,Gf

def all_pairs_max_flow(G, algorithm,output=True):
    for u in G.nodes:
        for v in G.nodes:
            if u <= v:
                continue
            flow,_ = algorithm(u,v,G)
            if flow and output:
                print("(%d,%d): %d" % (u,v,flow))

def time_graph(G):
    print("All-Pairs Ford-Fulkerson:")
    now = time.perf_counter()
    all_pairs_max_flow(G,ford_fulkerson,output=False)
    print(time.perf_counter()-now)

    print("All-Pairs Dinitz:")
    now = time.perf_counter()
    all_pairs_max_flow(G,dinitz,output=False)
    print(time.perf_counter()-now)

    now = time.perf_counter()
    print("Gomory-Hu")
    T = gomory_hu(G)
    print(time.perf_counter()-now)

def show_graph(G):
    print("All-Pairs Ford-Fulkerson:")
    all_pairs_max_flow(G,ford_fulkerson,output=True)

    print("All-Pairs Dinitz:")
    now = time.perf_counter()
    all_pairs_max_flow(G,dinitz,output=True)

    now = time.perf_counter()
    print("Gomory-Hu")
    T = gomory_hu(G)
    plot_edge_weights(G)
    plt.show()
    plot_edge_weights(T)
    plt.show()

def dinitz(s,t,G):
    Gf = init_residual(G)
    f = 0
    edges = bfs(s,t,Gf)
    while edges:
        Gp = nx.create_empty_copy(Gf)
        Gp.add_edges_from(edges)
        f += dinitz_blocking_flow(s,t,Gp,Gf)
        edges = bfs(s,t,Gf)
    return f,Gf

def gomory_hu(G,plot=False):
    vertex_mapping = {0:G.nodes()}
    T = nx.Graph()
    T.add_nodes_from([0])
    blob = find_blob(vertex_mapping)

    while blob != -1:
        neighbors = [n for n in T.neighbors(blob)]
        neighbor_weights = {blob2:weight(blob,blob2,T) for blob2 in neighbors}
        Tp,vertex_mapping_p,offset = expand_blob(G,T,vertex_mapping,blob) 
        s = 0
        t = 1
        f,Tpf = ford_fulkerson(s,t,Tp)
        A = min_cut(s,t,Tpf)
        
        new_A = []
        new_B = []
        for key,value in vertex_mapping_p.items():
            if key >= len(vertex_mapping[blob]):
                break
            if key in A:
                new_A.append(value[0])
            else:
                new_B.append(value[0])
        
        T.remove_node(blob)
        vertex_mapping[blob] = new_A
        T.add_node(blob)

        b_idx = T.number_of_nodes()
        vertex_mapping[b_idx] = new_B
        T.add_node(b_idx)

        for n in neighbors:
            if n-offset in A:
                if blob-offset in A:
                    T.add_edge(blob, n, weight=neighbor_weights[n])
                else:
                    T.add_edge(b_idx, n, weight=neighbor_weights[n])
            else:
                if blob-offset not in A:
                    T.add_edge(blob, n, weight=neighbor_weights[n])
                else: 
                    T.add_edge(b_idx, n, weight=neighbor_weights[n])

        T.add_edge(blob,b_idx,weight=f)
        if plot:
            plot_edge_weights(T)
            plt.show()
        
        blob = find_blob(vertex_mapping)

    T2 = nx.Graph()
    for u in T.nodes():
        T2.add_node(vertex_mapping[u][0])
    for u,v in T.edges():
        T2.add_edge(vertex_mapping[u][0],vertex_mapping[v][0],weight=weight(u,v,T))


    return T2
