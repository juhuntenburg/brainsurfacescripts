# -*- coding: utf-8 -*-
'''
Functions to read and write vtk files
------------------------------------
* read takes vtk file and returns vertex and face array
* write takes vertex and faces array and optional comment and returns vtk file
* reading and writing of texture not supported yet

TO DO: add reading comments
'''

def read_vtk(file):
    # read full file while dropping empty lines 
    vtk_df=pd.read_csv(file, header=None)
    vtk_df=vtk_df.dropna()
    # extract number of vertices and faces
    number_vertices=int(vtk_df[vtk_df[0].str.contains('POINTS')][0].iloc[0].split()[1])
    number_faces=int(vtk_df[vtk_df[0].str.contains('POLYGONS')][0].iloc[0].split()[1])
    # read vertices into df and array
    start_vertices= (vtk_df[vtk_df[0].str.contains('POINTS')].index.tolist()[0])+1
    vertex_df=pd.read_csv(file, skiprows=range(start_vertices), nrows=number_vertices, sep='\s*', header=None)
    if np.array(vertex_df).shape[1]==3:
        vertex_array=np.array(vertex_df)
    # sometimes the vtk format is weird with 9 indices per line, then it has to be reshaped
    elif np.array(vertex_df).shape[1]==9:
        vertex_df=pd.read_csv(file, skiprows=range(start_vertices), nrows=number_vertices/3+1, sep='\s*', header=None)
        vertex_array=np.array(vertex_df.iloc[0:1,0:3])
        vertex_array=np.append(vertex_array, vertex_df.iloc[0:1,3:6], axis=0)
        vertex_array=np.append(vertex_array, vertex_df.iloc[0:1,6:9], axis=0)
        for row in range(1,(number_vertices/3+1)):
            for col in [0,3,6]:
                vertex_array=np.append(vertex_array, array(vertex_df.iloc[row:(row+1),col:(col+3)]),axis=0) 
        # strip rows containing nans
        vertex_array=vertex_array[ ~np.isnan(vertex_array) ].reshape(number_vertices,3)
    else:
        print "vertex indices out of shape"
    # read faces into df and array
    start_faces= (vtk_df[vtk_df[0].str.contains('POLYGONS')].index.tolist()[0])+1
    face_df=pd.read_csv(file, skiprows=range(start_faces), nrows=number_faces, sep='\s*', header=None)
    face_array=np.array(face_df.iloc[:,1:4])
    # read data into df and array if exists
    if vtk_df[vtk_df[0].str.contains('POINT_DATA')].index.tolist()!=[]:
        start_data=(vtk_df[vtk_df[0].str.contains('POINT_DATA')].index.tolist()[0])+3
        number_data = number_vertices
        data_df=pd.read_csv(file, skiprows=range(start_data), nrows=number_data, sep='\s*', header=None)
        data_array=np.array(data_df)
    else:
        data = np.empty(0)
    
    return vertex_array, face_array, data_array



def write_vtk(filename, vertices, faces, data=None, comment=None):
    # infer number of vertices and faces
    number_vertices=vertices.shape[0]
    number_faces=faces.shape[0]
    number_data=data.shape[0]
    # make header and subheader dataframe
    header=['# vtk DataFile Version 3.0',
            '%s'%comment,
            'ASCII',
            'DATASET POLYDATA',
            'POINTS %i float'%number_vertices
             ]
    header_df=pd.DataFrame(header)
    sub_header=['POLYGONS %i %i'%(number_faces, 4*number_faces)]
    sub_header_df=pd.DataFrame(sub_header)    
    # make dataframe from vertices
    vertex_df=pd.DataFrame(vertices)
    # make dataframe from faces, appending first row of 3's (indicating the polygons are triangles)
    triangles=np.reshape(3*(np.ones(number_faces)), (number_faces,1))
    triangles=triangles.astype(int)
    faces=faces.astype(int)
    faces_df=pd.DataFrame(np.concatenate((triangles,faces),axis=1))
    # write dfs to csv
    header_df.to_csv(filename, header=None, index=False)
    with open(filename, 'a') as f:
        vertex_df.to_csv(f, header=False, index=False, float_format='%.3f', sep=' ')
    with open(filename, 'a') as f:
        sub_header_df.to_csv(f, header=False, index=False)
    with open(filename, 'a') as f:
        faces_df.to_csv(f, header=False, index=False, float_format='%.0f', sep=' ')
    # if there is data append second subheader and data
    if data!=None:
        datapoints=data.shape[1]
        sub_header2=['POINT_DATA %i'%(number_data),
                     'SCALARS EmbedVertex float %i'%(datapoints),
                     'LOOKUP_TABLE default']
        sub_header_df2=pd.DataFrame(sub_header2)
        data_df=pd.DataFrame(data)
        with open(filename, 'a') as f:
            sub_header_df2.to_csv(f, header=False, index=False)
        with open(filename, 'a') as f:
            data_df.to_csv(f, header=False, index=False, float_format='%.16f', sep=' ')
