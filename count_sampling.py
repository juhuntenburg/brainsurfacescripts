# -*- coding: utf-8 -*-
import nibabel as nb
import numpy as np
from math import sqrt
from vtk_rw import read_vtk
from ply_rw import read_ply

'''
Counts how often each voxel in input volume is sampled by input mesh.
Mesh and volume need to be in the same space.
'''
# TO DO: allow transformation input, find a better way to get pixdims?

vol_file = raw_input('volume file: ')  # /scr/ilz3/myelinconnect/old/KSYT_sub010/preprocessed/rest1_1/coregister/epi2highres_nonlin.nii.gz'
surf = raw_input('surface file: ')  # /scr/ilz3/myelinconnect/surfaces/simplification_sub010/rh_avg_simple_a.ply'
out = raw_input('output file: ')  # /scr/ilz3/myelinconnect/surfaces/simplification_sub010/rh_avg_simple_a_ply_qa.nii.gz'


# load information from volume
vol = nb.load(vol_file)
vol_data = vol.get_data()
vol_hdr = vol.get_header()
vol_affine = vol.get_affine()

# read vertices from surface
if surf[-3:] == 'vtk':
    v, f = read_vtk(surf)
elif surf[-3:] == 'ply':
    v, f = read_ply(surf)
else:
    print "only vtk and ply format supported yet"


# scale and translate surface into volume space
pixdims = list(vol_hdr['pixdim'][1:4])

v_scaled = np.empty_like(v)
v_scaled[:, 0] = ((1. / pixdims[0]) * v[:, 0])
v_scaled[:, 1] = (-(1. / pixdims[1]) * v[:, 1])
v_scaled[:, 2] = (-(1. / pixdims[2]) * v[:, 2])

trans = [0, (vol.shape[1] - 1), (vol.shape[2] - 1)]

v_trans = np.empty_like(v_scaled)
v_trans[:, 0] = v_scaled[:, 0] + trans[0]
v_trans[:, 1] = v_scaled[:, 1] + trans[1]
v_trans[:, 2] = v_scaled[:, 2] + trans[2]

# for each vertex, find the closest voxel, count how often each voxel is sampled
count = np.empty_like(vol_data)
print "counting vertex sampling"

for i in range(v_trans.shape[0]):
    
    # get coordinates of vertex
    x = v_trans[i, 0]
    y = v_trans[i, 1]
    z = v_trans[i, 2]

    # define 8 neighbouring voxels
    n1 = (np.floor(x), np.floor(y), np.floor(z))
    n2 = (n1[0] + 1, n1[1], n1[2])
    n3 = (n1[0], n1[1] + 1, n1[2])
    n4 = (n1[0], n1[1], n1[2] + 1)
    n5 = (n1[0] + 1, n1[1] + 1, n1[2])
    n6 = (n1[0], n1[1] + 1, n1[2] + 1)
    n7 = (n1[0] + 1, n1[1], n1[2] + 1)
    n8 = (n1[0] + 1, n1[1] + 1, n1[2] + 1)
    
    neighbours = [n1, n2, n3, n4, n5, n6, n7, n8]
    
    # calculate distances
    distances = []
    for n in range(len(neighbours)):
        d = sqrt((x - neighbours[n][0]) ** 2 + (y - neighbours[n][1]) ** 2 + 
                 (z - neighbours[n][2]) ** 2)
        distances.append(d)
        
    # find nearest neighbour
    m = min(distances)
    nn = neighbours[[i for i, j in enumerate(distances) if j == m][0]]
    
    # add one to voxel in that position
    count[nn[0], nn[1], nn[2]] += 1


# save count to nifti
count_img = nb.Nifti1Image(count, vol_affine, vol_hdr)
print "saving image"
nb.save(count_img, out)
