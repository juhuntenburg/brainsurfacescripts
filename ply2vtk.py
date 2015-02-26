from ply_rw import read_ply
from vtk_rw import write_vtk

'''
Meta script to convert ply to vtk file
'''

in_ply = raw_input("input ply file: ")
out_vtk = raw_input("output vtk file: ")
comment = raw_input("comment (optional): ")

if comment == "":
    comment = 'vtk format converted from ply'

v, f = read_ply(in_ply)
write_vtk(out_vtk, v, f, comment)
