import numpy as np

def calculate_normals(vertices, faces):
    '''
    Calculate the vertex normals of a given mesh, as average of face normals
    weighted by face area. Normals are normalized to length 1. adapted from:
    https://sites.google.com/site/dlampetest/python/calculating-normals-of-a
    -triangle-mesh-using-numpy
    '''
    
    triangles = vertices[faces]
    face_normals = np.cross(triangles[::, 1 ] - triangles[::, 0]  , 
                            triangles[::, 2 ] - triangles[::, 0])
    # face_normals /= np.linalg.norm(face_normals, axis=1)[:,np.newaxis]
    # weighting by surface area of the triangle =  half-length of the normal
    
    face_normals /= 2  
    vertex_normals = np.zeros(vertices.shape, dtype=vertices.dtype)
    vertex_normals[faces[:, 0]] += face_normals
    vertex_normals[faces[:, 1]] += face_normals
    vertex_normals[faces[:, 2]] += face_normals
    vertex_normals /= np.linalg.norm(vertex_normals, axis=1)[:, np.newaxis]
    
    return vertex_normals


def compare_normals(normals_a, normals_b):
    '''
    Calculate the angles between two sets of normals
    '''
    diff_rad = np.zeros((normals_a.shape[0],))
    for i in range(normals_a.shape[0]):
        diff_rad[i] = np.arccos(np.dot(normals_a[i], normals_b[i]))
    
    diff_deg = diff_rad  * (180/np.pi)
    return diff_rad, diff_deg


def add_neighbours(node, length, graph, labels, tree):
    '''
    Find unlabelled neighbours of a node in the graph and add them to 
    AVLtree
    '''
    import numpy as np
    # find direct neighbours of the node
    neighbours = np.array(graph.neighbors(node))
    # check that they don't already have a label
    unlabelled = neighbours[np.where(labels[neighbours][:, 1] == -1)[0]]
    # insert source neighbour pair with edge length to tree
    for u in unlabelled:
        # new_length = length + graph[node][u]['length']
        new_length = length + graph[node][u]
        tree.insert(new_length, (node, u))
    
    return tree

def find_voronoi_seeds(simple_vertices, simple_faces, 
                       complex_vertices, complex_faces, 
                       cutoff_angle=(np.pi/2)):
    '''
    Finds those points on the complex mesh that correspond best to the
    simple mesh (taking into accound euclidian distance and direction of normals)
    while forcing a one-to-one mapping
    '''
    from bintrees import FastAVLTree
    import scipy.spatial as spatial
    
    
    # calculate normals for simple and complex vertices
    simple_normals = calculate_normals(simple_vertices, simple_faces)
    complex_normals = calculate_normals(complex_vertices, complex_faces)
    
    # prepare array to store seeds
    voronoi_seed_idx = np.zeros((simple_vertices.shape[0],), dtype='int64')-1
    missing = np.where(voronoi_seed_idx==-1)[0].shape[0]
    
    # initialize with all vertices and small number of neighbours
    remaining_idxs = range(simple_vertices.shape[0])
    neighbours = 100
    
    while missing > 0:
    
        print 'producing nearest neighbours', neighbours
        # find nearest neighbours of simple vertices on complex mesh using kdtree
        inaccuracy, mapping  = spatial.KDTree(complex_vertices).query(simple_vertices[remaining_idxs], k=neighbours)
        
        # create tidy long-format lists
        simple_idxs = np.asarray([neighbours*[simple_idx] for simple_idx in remaining_idxs]).flatten()
        candidate_idxs = mapping.flatten()
        diff_euclid = inaccuracy.flatten()
        
        # for each vertex pair calculate the angle between their normals
        diff_normals, _ = compare_normals(simple_normals[simple_idxs], complex_normals[candidate_idxs])
        print 'candidates', diff_normals.shape[0]
        # remove those pairs that have an angle / distance above cutoff
        #mask = np.unique(np.concatenate((np.where(diff_euclid>cutoff_euclid)[0], np.where(diff_normals>cutoff_rad)[0])))
        mask = np.unique(np.where(diff_normals>cutoff_rad)[0])
        diff_normals = np.delete(diff_normals, mask)
        diff_euclid = np.delete(diff_euclid, mask)
        simple_idxs = np.delete(simple_idxs, mask)
        candidate_idxs = np.delete(candidate_idxs, mask)
        
        print 'remaining candidates', diff_normals.shape[0]
        # calculate scores for each vertex pair
        scores = (diff_normals-np.mean(diff_normals)) + (diff_euclid-np.mean(diff_euclid))
        
        print 'producing tree'
        # make a binary search tree from the scores and vertex pairs, 
        # organisation is key: score, values: tuple(simple_vertex, candiate_complex_vertex)
        tree = FastAVLTree(zip(scores, zip(simple_idxs, candidate_idxs)))
        
        while tree.count > 0:
        
            min_item =  tree.pop_min()
            simple_idx = min_item[1][0]
            candidate_idx = min_item[1][1]
        
            if (voronoi_seed_idx[simple_idx] == -1):
                if candidate_idx not in voronoi_seed_idx:
                    voronoi_seed_idx[simple_idx] = candidate_idx
                else:
                    pass
            else:
                pass
        
            missing = np.where(voronoi_seed_idx==-1)[0].shape[0]
            if missing == 0:
                break
        
        # if the tree is empty, but there are still seeds missing, increase the number of nearest neighbours
        # and repeat procedure, but only for those simple vertices that have not been matched yet
        print 'missing ,', missing
        remaining_idxs = np.where(voronoi_seed_idx==-1)[0]
        neighbours *= 5

    return voronoi_seed_idx, inaccuracy


def competetive_fast_marching(vertices, graph, seeds):
    '''
    Label all vertices on highres mesh to the closest seed vertex
    using a balanced binary search tree
    '''
    import numpy as np
    import sys
    from bintrees import FastAVLTree
    # make a labelling container to be filled with the search tree
    # first column are the vertex indices of the complex mesh
    # second column are the labels from the simple mesh
    # (-1 for all but the corresponding points for now)
    labels = np.zeros((vertices.shape[0], 2), dtype='int64') - 1
    labels[:, 0] = range(vertices.shape[0])
    for i in range(seeds.shape[0]):
        labels[seeds[i]][1] = i
    # initiate AVLTree for binary search
    tree = FastAVLTree()
    # organisation of the tree will be
    # key: edge length; value: tuple of vertices (source, target)
    # add all neighbours of the voronoi seeds
    for v in seeds:
        add_neighbours(v, 0, graph, labels, tree)
    # Competetive fast marching starting from voronoi seeds
    printcount = 0
    while tree.count > 0:
        
        printcount += 1
        
        # pdb.set_trace()
        # pop the item with minimum edge length
        min_item = tree.pop_min()
        length = min_item[0]
        source = min_item[1][0]
        target = min_item[1][1]
        # if target no label yet (but source does!), assign label of source
        if labels[target][1] == -1:
            if labels[source][1] == -1:
                sys.exit('Source has no label, something went wrong!')
            else:
                # assign label of source to target
                labels[target][1] = labels[source][1]
        
            # test if labelling is complete
            if any(labels[:, 1] == -1):
                # if not, add neighbours of target to tree
                add_neighbours(target, length, graph, labels, tree)
            else:
                break
        
        # if the target already has a label the item is just popped out of the 
        # tree and nothing else happens
        else:
            pass
        
        # for monitoring the progress
        if np.mod(printcount, 100) == 0.0:
            print 'tree ' + str(tree.count)
            print 'labels ' + str(np.where(labels[:, 1] == -1)[0].shape[0])
    
    return labels


def sample_simple(highres_data, labels):
    
    import numpy as np
    '''
    Computes the mean of data from highres mesh that is assigned to the same
    label (typical simple mesh vertices).
    '''
    # create new empty lowres data array
    lowres_data = np.empty((int(labels.max() + 1), highres_data.shape[1]))
    # find all vertices on highres and mean
    for l in range(int(labels.max() + 1)):
        patch = np.where(labels == l)[0]
        patch_data = highres_data[patch]
        patch_mean = np.mean(patch_data, axis=0)
        lowres_data[l] = patch_mean
    
    return lowres_data


def sample_volume(nii_file, vertices):
    '''
    Samples volumetric data on surface mesh vertices
    '''
    import numpy as np
    import nibabel as nb
    
    img = nb.load(nii_file)
    affine = img.get_affine()
    data = img.get_data()
    
    # for each vertex in the highres mesh find voxel it maps to
    dim = -(np.round([affine[0, 0], affine[1, 1], affine[2, 2]], 1))
    idx = np.asarray(np.round(vertices / dim), dtype='int64')
    data_mesh = data[idx[:, 0], idx[:, 1], idx[:, 2]]
    
    return data_mesh
    
