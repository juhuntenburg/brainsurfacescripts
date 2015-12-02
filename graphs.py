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
            G[e[0]][e[1]]['length']=np.linalg.norm(nodes[e[0]]-nodes[e[1]])
            
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
    main_clust = int(np.where(subclust[node]!=0)[0][0])
    sub_clust = int(subclust[node][main_clust])
    clust_id = str(main_clust)+'_'+str(sub_clust)
    return clust_id


def adjacent_subcluster(nodes, triangles, subclustering):
    
    G = graph_from_mesh(nodes, triangles)
    
    unique_pairs=[]
    for e in G.edges():
        ca = node_cluster_id(e[0], subclustering)
        cb = node_cluster_id(e[1], subclustering)
        if ca != cb:
            if (not (ca, cb) in unique_pairs and not (cb, ca) in unique_pairs):
                unique_pairs.append((ca, cb))
    
    return unique_pairs