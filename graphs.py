def graph_from_mesh(nodes, triangles, node_coords=False, edge_length=False):
    '''
    Creates a networkx graph from a mesh
    '''
    
    import numpy as np
    import networkx as nx
    
    # initiate empty graph
    G=nx.Graph()
    
    # add node indices
    G.add_nodes_from(range(len(nodes)))
    
    # add edges
    G.add_edges_from(triangles[:,[0,1]])
    G.add_edges_from(triangles[:,[0,2]])
    G.add_edges_from(triangles[:,[1,2]])
    
    # caution! this adds a key 'coords' to the node
    # which will also be picked up by .neighbors methods
    
    if edge_length:
        for e in G.edges_iter():
            #G[e[0]][e[1]]['length']=np.linalg.norm(nodes[e[0]]-nodes[e[1]])
            G[e[0]][e[1]]=np.linalg.norm(nodes[e[0]]-nodes[e[1]])
            G[e[1]][e[0]]=np.linalg.norm(nodes[e[0]]-nodes[e[1]])
            
    if node_coords:
        for n in G.nodes_iter():
            G[n]['coords']=nodes[n]

    return G


def subcluster_graph(nodes, triangles, clustering):
    
    '''
    Finds the non-connected components of each cluster
    Currently works with networkx 1.6, not with 1.10
    '''
    
    import numpy as np
    import networkx as nx
    
    # initiate empty graph
    G=nx.Graph()
    
    # add node indices
    G.add_nodes_from(range(len(nodes)))
    
    # add edges
    G.add_edges_from(triangles[:,[0,1]])
    G.add_edges_from(triangles[:,[0,2]])
    G.add_edges_from(triangles[:,[1,2]])
        
        
    # make an array with one column per cluster
    # entries will be 0 if this node is not in the respective subcluster
    # otherwise they refer to the number of the sublcuster
    # therefore, the subclusters have to be one-indexed
    subclust_arr = np.zeros((clustering.shape[0], int(clustering.max()+1)))
    #for each cluster
    for nclust in range(int(clustering.max()+1)):
        # make a subgraph only containing nodes from this cluster
        clustgraph = nx.subgraph(G, np.where(clustering==nclust)[0])
        # find the seperate connected components of this subgraph
        subclust = nx.connected_components(clustgraph)
        
        for i in range(len(subclust)):
            for j in subclust[i]:
                subclust_arr[j][nclust] = i+1

    return subclust_arr

def node_cluster_id(node, subclustering):
    import numpy as np
    if np.where(subclustering[node]!=0)[0].shape[0] == 0:
        return None
    else:
        main_clust = int(np.where(subclustering[node]!=0)[0][0])
        sub_clust = int(subclustering[node][main_clust])
        clust_id = str(main_clust)+'_'+str(sub_clust)
        return clust_id


def adjacent_subcluster(nodes, triangles, subclustering):
    import numpy as np
    G = graph_from_mesh(nodes, triangles)
    
    unique_pairs=[]
    for e in G.edges():
        ca = node_cluster_id(e[0], subclustering)
        cb = node_cluster_id(e[1], subclustering)
        if (ca is not None and cb is not None and ca != cb):
            if (not (ca, cb) in unique_pairs and not (cb, ca) in unique_pairs):
                unique_pairs.append((ca, cb))
    
    return unique_pairs


def dijkstra(G, start, end=None):
    '''
    adapted to work with networkx graphs from:
    https://github.com/joyrexus/dijkstra/blob/master/dijkstra.py
    '''
    '''
    dijkstra's algorithm determines the length from `start` to every other 
    vertex in the graph.
    The graph argument `G` should be a dict indexed by nodes.  The value 
    of each item `G[v]` should also a dict indexed by successor nodes.
    In other words, for any node `v`, `G[v]` is itself a dict, indexed 
    by the successors of `v`.  For any directed edge `v -> w`, `G[v][w]` 
    is the length of the edge from `v` to `w`.
        graph = {'a': {'b': 1}, 
                 'b': {'c': 2, 'b': 5}, 
                 'c': {'d': 1},
                 'd': {}}
    Returns two dicts, `dist` and `pred`:
        dist, pred = dijkstra(graph, start='a') 
    
    `dist` is a dict mapping each node to its shortest distance from the
    specified starting node:
        assert dist == {'a': 0, 'c': 3, 'b': 1, 'd': 4}
    `pred` is a dict mapping each node to its predecessor node on the
    shortest path from the specified starting node:
        assert pred == {'b': 'a', 'c': 'b', 'd': 'c'}
    
    '''
    
    from pqdict import PQDict
    import networkx as nx
    
    inf = float('inf')
    D = {start: 0}          # mapping of nodes to their dist from start
    Q = PQDict(D)           # priority queue for tracking min shortest path
    P = {}                  # mapping of nodes to their direct predecessors
    #U = set(G.keys())       # unexplored nodes
    U = set(G.nodes())

    while U:                                    # nodes yet to explore
        (v, d) = Q.popitem()                    # node w/ min dist d on frontier
        D[v] = d                                # est dijkstra greedy score
        U.remove(v)                             # remove from unexplored
        if v == end: break

        # now consider the edges from v with an unexplored head -
        # we may need to update the dist of unexplored successors 
        for w in G[v]:                          # successors to v
            if w in U:                          # then w is a frontier node
                d = D[v] + G[v][w]              # dgs: dist of start -> v -> w
                if d < Q.get(w, inf):
                    Q[w] = d                    # set/update dgs
                    P[w] = v                    # set/update predecessor

    return D, P



def shortest_path(G, start, end):
    
    from pqdict import PQDict
    import networkx as nx
    
    dist, pred = dijkstra(G, start, end)
    v = end
    path = [v]
    while v != start:
        v = pred[v]
        path.append(v)        
    path.reverse()
    return path


def sorted_path(G, path):
    import networkx as nx
    subG = nx.subgraph(G, path)
    ends = []
    for node in path:
        if nx.degree(subG)[node] == 1:
            ends.append(node)
    sorted_path = shortest_path(subG, ends[0], ends[1])
    return sorted_path