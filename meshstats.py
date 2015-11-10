#calculate edge lengths and face areas from vertex and face arrays
def edges_areas_angles(vertices, faces):
    
    edges=np.empty((faces.shape[0], 3))
    angles=np.empty((faces.shape[0], 3))
    areas=np.empty((faces.shape[0], 1))
    
    for f in range(faces.shape[0]):
        # find all vertices belonging to this face
        vertex_coords=vertices[faces[f]]
        
        # make edge vectors
        edge_vec=np.array([vertex_coords[1]-vertex_coords[0],  # vector from vertex 0 to 1
                            vertex_coords[2]-vertex_coords[0],  # vector from vertex 0 to 2
                            vertex_coords[2]-vertex_coords[1]]) # vector from vertex 1 to 2
        
        # calculate edge lengths
        edge_len=np.linalg.norm(edge_vec, axis=1)
        
        # calculate unit vectors
        edge_len=edge_len[:,np.newaxis]
        edge_unit=edge_vec/edge_len
        
        # calculate angles
        angle = np.array([np.arccos(np.dot(edge_unit[0], edge_unit[1])),
                           np.arccos(np.dot(edge_unit[0], edge_unit[2])),
                           np.arccos(np.dot(edge_unit[1], edge_unit[2]))])
        
        # calculate areas
        area = (.5*edge_len[0]*edge_len[1]*np.sin(angle[0]))[0]
        
        # write into large arrays
        edges[f]=np.squeeze(edge_len)
        angles[f]=np.squeeze(angle)
        areas[f]=area
        
    return edges.flatten(), areas.flatten(), angles.flatten()


# calculate some basic stats on mesh
def meshstats(in_array):
    
    array_mean = np.nanmean(in_array)
    array_sd = np.nanstd(in_array)
    array_z=sp.stats.zscore(np.nan_to_num(in_array))

    return array_mean, array_sd, array_z