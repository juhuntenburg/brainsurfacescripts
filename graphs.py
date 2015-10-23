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
    if node_coords:
        for n in G.nodes_iter():
            G[n]['coords']=nodes[n]
    
    if edge_length:
        for e in G.edges_iter():
            G[e[0]][e[1]]['length']=np.linalg.norm(nodes[e[0]]-nodes[e[1]])

    return G
