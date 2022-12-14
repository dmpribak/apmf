import networkx as nx
import numpy as np

import copy

def weight(u, v, G):
    return G[u][v]["weight"]

def set_weight(u, v, G, w):
    G[u][v]["weight"] = w

def plot_edge_weights(G): 
    pos = nx.spring_layout(G)
    nx.draw(G,pos,with_labels=True)
    nx.draw_networkx_edge_labels(G,pos,edge_labels=nx.get_edge_attributes(G,"weight"))

def init_residual(G):
    Gf = nx.DiGraph()
    Gf.add_nodes_from(G.nodes())

    if G.is_directed():
        for u,v,d in G.edges(data=True):
            Gf.add_edge(u,v,weight=weight(u,v,G))
            if not G.has_edge(v,u):
                Gf.add_edge(v, u, weight=0)
    else:
        for u,v,d in G.edges(data=True):
            Gf.add_edge(u,v,weight=weight(u,v,G))
            Gf.add_edge(v,u,weight=weight(u,v,G))
        
    return Gf

"""
Returns 0 if no path between nodes
"""
def dfs(s, t, G):
    stack = [s]
    found = np.zeros(G.number_of_nodes(),dtype=np.int32)
    found[s] = 1
    path = []
    parent = np.zeros(G.number_of_nodes(),dtype=np.int32)
    while stack:
        loc = stack.pop()
        neighbors = G.neighbors(loc)
        for v in neighbors:
            if found[v] or (weight(loc, v, G) == 0):
                continue
            if v == t:
                backtrack = t
                parent[t] = loc
                while backtrack != s:
                    path.append((parent[backtrack],backtrack))
                    backtrack = parent[backtrack]
                path.reverse()
                return path
            stack.append(v)
            parent[v] = loc
            found[v] = 1
    return 0

def bfs(s, t, G):
    frontier = []
    current = [s]
    parents = [[] for i in range(G.number_of_nodes())]
    found = np.zeros(G.number_of_nodes(),dtype=np.int32)
    found[s]=1
    found2 = np.zeros(G.number_of_nodes(),dtype=np.int32)
    edges = []
    target_reached = False
    while current:
        for u in current:
            neighbors = G.neighbors(u)
            for v in neighbors:
                if ((found[v] and (found[v] < found[u]+1)) or (weight(u,v,G) == 0)):
                    continue
                if v == t:
                    target_reached = True
                parents[v].append(u)
                if not found[v]:
                    found[v] = found[u]+1
                    frontier.append(v)
        if target_reached:
            vertices = [t]
            while vertices:
                v = vertices.pop()
                for u in parents[v]:
                    edges.append((u,v))
                    if not found2[u]:
                        vertices.append(u)
                        found2[u] = 1
            return edges
        current = frontier
        frontier = []
    return 0

def dinitz_blocking_flow(s,t,Gp,Gf):
    stack = [s]
    found = np.zeros(Gf.number_of_nodes(),dtype=np.int32)
    found[s] = 1
    path = []
    parent = np.zeros(Gf.number_of_nodes(),dtype=np.int32)
    f = 0
    while stack:
        loc = stack.pop()
        neighbors = Gp.neighbors(loc)
        for v in neighbors:
            if found[v] or (weight(loc, v, Gf) == 0):
                continue
            if v == t:
                backtrack = t
                parent[t] = loc
                while backtrack != s:
                    path.append((parent[backtrack],backtrack))
                    backtrack = parent[backtrack]
                path.reverse()
                min_weight = path_min_weight(path,Gf)
                block_path(path,Gf,min_weight)
                f += min_weight
                continue
            stack.append(v)
            parent[v] = loc
            found[v] = 1
    return f

def path_min_weight(path, G):
    min_weight = weight(*path[0],G)
    for u,v in path:
        if weight(u,v,G) < min_weight:
            min_weight = weight(u,v,G)
    return min_weight

def block_path(path, G, min_weight):
    for u,v in path:
        w = weight(u,v,G)
        w2 = weight(v,u,G)
        set_weight(u,v,G,w-min_weight)
        set_weight(v,u,G,w2+min_weight)

def find_blob(vertex_mapping):
    for key,val in vertex_mapping.items():
        if len(val) > 1:
            return key
    return -1

def expand_blob(G,T,vertex_mapping,blob):
    s = vertex_mapping[blob][0]
    t = vertex_mapping[blob][1]
    offset = T.number_of_nodes()
    Tc = T.copy()
    Tc.remove_node(blob) # blob is still present in vertex_mapping!!!!
    Tp = nx.Graph()

    vertex_mapping_p = {}
    for i,u in enumerate(vertex_mapping[blob]): # break apart blob and add to Tp
        Tp.add_node(i)
        vertex_mapping_p[i] = [u]
    
    ccs = nx.connected_components(Tc)
    for i,cc in enumerate(ccs): # Contract nodes
        current_blob = i+len(vertex_mapping[blob])
        Tp.add_node(current_blob) 
        vertex_mapping_p[current_blob] = [] 
        for key in cc: # each cc is a collection of blobs from T
            for u in vertex_mapping[key]:
                vertex_mapping_p[current_blob].append(u)


    for i,u in enumerate(vertex_mapping[blob]): # add proper edges to Tp
        for key,vertices in vertex_mapping_p.items():
            if key+offset == blob:
                continue
            
            total_flow = 0
            for v in vertices:
                if G.has_edge(u,v):
                    total_flow += weight(u,v,G)
            if total_flow:
                 Tp.add_edge(i,key,weight=total_flow)
    
    return Tp,vertex_mapping_p,offset

def min_cut(s, t, Gf):
    stack = [s]
    found = np.zeros(Gf.number_of_nodes(),dtype=np.int32)
    found[s] = 1
    A = [s]
    while stack:
        loc = stack.pop()
        neighbors = Gf.neighbors(loc)
        for v in neighbors:
            if found[v] or (weight(loc, v, Gf) == 0):
                continue
            stack.append(v)
            found[v] = 1
            A.append(v)
    return A
