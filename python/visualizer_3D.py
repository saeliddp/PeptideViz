# -*- coding: utf-8 -*-
"""
Created on Sun Nov 24 22:32:15 2019
@author: saeli
TO VISUALIZE STANDARD DEVIATIONS, SEE LINE 27
"""
import build_3D_delta_lists
from mpl_toolkits.mplot3d import Axes3D
import pickle
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

aas = "ACDEFGHIKLMNPQRSTVWY"

with open('meansxyzm.pickle', 'rb') as fr:
    meansxyzm = pickle.load(fr)
with open('stdevsxyzm.pickle', 'rb') as fr:
    stdevsxyzm = pickle.load(fr)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

plasma = cm.get_cmap('plasma', 256)

# To visualize means, leave as is. To visualize standard deviations,
# comment out lines 29-31 and uncomment lines 32-34
s3d = ax.scatter3D(meansxyzm[0], meansxyzm[1], meansxyzm[2], c=meansxyzm[3], cmap=plasma)
title = 'Average Change in Binding Affinity for Mutations'
cbar_label = 'Average Change in Binding Affinity'
#s3d = ax.scatter3D(stdevsxyzm[0], stdevsxyzm[1], stdevsxyzm[2], c=stdevsxyzm[3], cmap=plasma)
#title = 'Standard Deviation of Change in Binding Affinity for Mutations'
#cbar_label = 'Standard Deviation of Change in Binding Affinity'

ax.set_xlabel('Initial AA')
ax.set_ylabel('AA After Mutation')
ax.set_zlabel('Position of Mutation')
fig.colorbar(s3d, ax=ax, label=cbar_label)
fig.suptitle(title)

# define ticks and labels
xts = []
yts = []
zts = []
xls = []
yls = []
zls = []

for i in range(20):
    xts.append(i)
    yts.append(i)
    xls.append(aas[i])
    yls.append(aas[i])

for i in range(12):
    zts.append(i)
    zls.append(i + 1)

ax.set_zticks(zts)
ax.set_yticks(yts)
ax.set_xticks(xts)
ax.set_xticklabels(xls)
ax.set_yticklabels(yls)
ax.set_zticklabels(zls)

plt.show()