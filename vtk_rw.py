# -*- coding: utf-8 -*-

'''
function to read vtk files
---------------------------
returns arrays of vertices and faces
'''
# TO DO: add reading comments and texture

def read_vtk(file):
    import numpy as np
    import pandas as pd
    print "reading vtk format"
    # read full file while dropping empty lines 
    vtk_df = pd.read_csv(file, header=None, engine='python')
    vtk_df = vtk_df.dropna()
    # extract number of vertices and faces
    number_vertices = int(vtk_df[vtk_df[0].str.contains('POINTS')][0].iloc[0].split()[1])
    number_faces = int(vtk_df[vtk_df[0].str.contains('POLYGONS')][0].iloc[0].split()[1])
    # read vertices into df and array
    start_vertices = (vtk_df[vtk_df[0].str.contains('POINTS')].index.tolist()[0]) + 1
    vertex_df = pd.read_csv(file, skiprows=range(start_vertices),
                            nrows=number_vertices, sep='\s*', header=None,
                            engine='python')
    if np.array(vertex_df).shape[1] == 3:
        vertex_array = np.array(vertex_df)
    # when the vtk format has 9 indices per line, it has to be reshaped
    elif np.array(vertex_df).shape[1] == 9:
        vertex_df = pd.read_csv(file, skiprows=range(start_vertices),
                                nrows=number_vertices / 3 + 1,
                                sep='\s*', header=None,
                                engine='python')
        vertex_array = np.array(vertex_df.iloc[0:1, 0:3])
        vertex_array = np.append(vertex_array, vertex_df.iloc[0:1, 3:6], axis=0)
        vertex_array = np.append(vertex_array, vertex_df.iloc[0:1, 6:9], axis=0)
        for row in range(1, (number_vertices / 3 + 1)):
            for col in [0, 3, 6]:
                vertex_array = np.append(vertex_array,
                                         np.array(vertex_df.iloc[row:(row + 1),
                                                              col:(col + 3)]),
                                         axis=0) 
        # strip rows containing nans
        vertex_array = vertex_array[ ~np.isnan(vertex_array) ].reshape(number_vertices, 3)
    else:
        print "vertex indices out of shape"
    # vertices = {'values' : vertex_array, 'number' : number_vertices}
    # read faces into df and array
    start_faces = (vtk_df[vtk_df[0].str.contains('POLYGONS')].index.tolist()[0]) + 1
    face_df = pd.read_csv(file, skiprows=range(start_faces),
                          nrows=number_faces, sep='\s*', header=None,
                          engine='python')
    face_array = np.array(face_df.iloc[:, 1:4])
    # faces = {'values' : face_array, 'number' : number_faces}
    
    return vertex_array, face_array
        

'''
function to write vtk files
---------------------------
uses arrays of vertices and faces
to construct ply file
can optionally take a comment
'''
# TO DO: add write texture

def write_vtk(filename, vertices, faces, comment=None):
    import numpy as np
    import pandas as pd
    print "writing vtk format"
    # infer number of vertices and faces
    number_vertices = vertices.shape[0]
    number_faces = faces.shape[0]
    # make header and subheader dataframe
    header = ['# vtk DataFile Version 3.0',
            '%s' % comment,
            'ASCII',
            'DATASET POLYDATA',
            'POINTS %i float' % number_vertices
             ]
    header_df = pd.DataFrame(header)
    sub_header = ['POLYGONS %i %i' % (number_faces, 4 * number_faces)]
    sub_header_df = pd.DataFrame(sub_header)
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
        sub_header_df.to_csv(f, header=False, index=False)
    with open(filename, 'a') as f:
        faces_df.to_csv(f, header=False, index=False,
                        float_format='%.0f', sep=' ')
