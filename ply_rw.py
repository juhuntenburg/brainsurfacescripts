# -*- coding: utf-8 -*-

'''
function to read ply files
---------------------------
returns arrays of vertices and faces
'''
# TO DO: add reading comments and texture

def read_ply(file):
    import numpy as np
    import pandas as pd
    print "reading ply format"
    # read full file and drop empty lines
    ply_df = pd.read_csv(file, header=None, engine='python')
    ply_df = ply_df.dropna()
    # extract number of vertices and faces, and row that marks the end of header
    number_vertices = int(ply_df[ply_df[0].str.contains('element vertex')][0].iloc[0].split()[2])
    number_faces = int(ply_df[ply_df[0].str.contains('element face')][0].iloc[0].split()[2])
    end_header = ply_df[ply_df[0].str.contains('end_header')].index.tolist()[0]
    # read vertex coordinates into dict
    vertex_df = pd.read_csv(file, skiprows=range(end_header + 1),
                            nrows=number_vertices, sep='\s*', header=None,
                            engine='python')
    vertex_array = np.array(vertex_df)
    # vertices = {'values' : vertex_array, 'number' : number_vertices}
    # read face indices into dict
    face_df = pd.read_csv(file, skiprows=range(end_header + number_vertices + 1),
                          nrows=number_faces, sep='\s*', header=None,
                          engine='python')
    face_array = np.array(face_df.iloc[:, 1:4])
    # faces = {'values' : face_array, 'number' : number_faces}
    
    return vertex_array, face_array


'''
function to write ply files
---------------------------
uses arrays of vertices and faces
to construct ply file
can optionally take a comment
'''
# TO DO: add write texture

def write_ply(filename, vertices, faces, comment=None):
    import numpy as np
    import pandas as pd
    print "writing ply format"
    # infer number of vertices and faces
    number_vertices = vertices.shape[0]
    number_faces = faces.shape[0]
    # make header dataframe
    header = ['ply',
            'format ascii 1.0',
            'comment %s' % comment,
            'element vertex %i' % number_vertices,
            'property float x',
            'property float y',
            'property float z',
            'element face %i' % number_faces,
            'property list uchar int vertex_indices',
            'end_header'
             ]
    header_df = pd.DataFrame(header)
    # make dataframe from vertices
    vertex_df = pd.DataFrame(vertices)
    # make dataframe from faces, adding first row of 3s (indicating triangles)
    triangles = np.reshape(3 * (np.ones(number_faces)), (number_faces, 1))
    triangles = triangles.astype(int)
    faces = faces.astype(int)
    faces_df = pd.DataFrame(np.concatenate((triangles, faces), axis=1))
    # write dfs to csv
    header_df.to_csv(filename, header=None, index=False)
    with open(filename, 'a') as f:
        vertex_df.to_csv(f, header=False, index=False,
                         float_format='%.3f', sep=' ')
    with open(filename, 'a') as f:
        faces_df.to_csv(f, header=False, index=False,
                        float_format='%.0f', sep=' ')
