from vtk_rw import read_vtk
from ply_rw import write_ply

'''
Meta script to convert vtk to ply file
'''

in_vtk = raw_input("input vtk file: ")
out_ply = raw_input("output ply file: ")
comment = raw_input("comment (optional): ")

if comment == "":
    comment = 'ply format converted from vtk'

v, f, d= read_vtk(in_vtk)
write_ply(out_ply, v, f, comment)
