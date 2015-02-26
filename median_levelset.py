import nibabel as nb
import os
import numpy as np
from math import sqrt
import pandas as pd

'''
Finding median levelset surface in a group of subjects
----------------------------------------------------------
* Computes the mean distance of each surface to the average in a narrow band
* Original idea Christine Tardif
'''

# desired path for output csv file
out = '/scr/ilz3/myelinconnect/surfaces/median.csv'

# get list of right and left hemisphere surfaces
rh_surf_dir = '/scr/ilz3/myelinconnect/surfaces/surf_rh/orig/mid_surface/'
rh_surf_files = os.listdir(rh_surf_dir)

lh_surf_dir = '/scr/ilz3/myelinconnect/surfaces/surf_lh/orig/mid_surface/'
lh_surf_files = os.listdir(lh_surf_dir)

# create average right and left hemisphere surfaces
os.chdir(rh_surf_dir)
rh_avg = np.empty_like(nb.load(rh_surf_files[0]).get_data())
for surf in rh_surf_files:
    surf_data = nb.load(surf).get_data()
    rh_avg += surf_data
rh_avg = rh_avg / len(rh_surf_files)

os.chdir(lh_surf_dir)
lh_avg = np.empty_like(nb.load(lh_surf_files[0]).get_data())
for surf in lh_surf_files:
    surf_data = nb.load(surf).get_data()
    lh_avg += surf_data
lh_avg = lh_avg / len(lh_surf_files)


# compute distances in narrowband of sqrt(3)
os.chdir(rh_surf_dir)
rh_distances = dict()
for surf in rh_surf_files:
    surf_data = nb.load(surf).get_data()
    mask = (surf_data <= sqrt(3)) & (surf_data >= -sqrt(3))
    narrowband = surf_data[mask]
    avg_narrowband = rh_avg[mask]
    distance = abs(avg_narrowband - narrowband)
    mean_distance = np.mean(distance)
    rh_distances[surf[0:4]] = mean_distance

os.chdir(lh_surf_dir)    
lh_distances = dict()
for surf in lh_surf_files:
    surf_data = nb.load(surf).get_data()
    mask = (surf_data <= sqrt(3)) & (surf_data >= -sqrt(3))
    narrowband = surf_data[mask]
    avg_narrowband = lh_avg[mask]
    distance = abs(avg_narrowband - narrowband)
    mean_distance = np.mean(distance)
    lh_distances[surf[0:4]] = mean_distance

# compute mean distances across left and right hemisphere for each subject
lh_dist = []
rh_dist = []
dist = []
subjects = []
for sub in rh_distances.keys():
    rh_val = rh_distances[sub]
    rh_dist.append(rh_val)
    lh_val = lh_distances[sub]
    lh_dist.append(lh_val)
    mean_val = (rh_val + lh_val) / 2
    dist.append(mean_val)
    subjects.append(sub)

# write data to csv file
df = pd.DataFrame(zip(subjects, rh_dist, lh_dist, dist),
                  columns=["subject", "rh distance",
                             "lh distance", "mean_distance"])
df.to_csv(out)
