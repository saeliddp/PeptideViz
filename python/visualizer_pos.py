# -*- coding: utf-8 -*-
"""
Created on Sun Nov 24 19:39:01 2019
@author: saeli
TO VISUALIZE STANDARD DEVIATIONS, SEE LINE 54
"""
import build_pos_delta_lists
import means_stdevs_pos
import pickle
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

labels = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 
              'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

with open('means_pos.pickle', 'rb') as fr:
    means = pickle.load(fr)
with open('stdevs_pos.pickle', 'rb') as fr:
    stdevs = pickle.load(fr)
    
def buildValidArray(data):
    max_val = 0
    min_val = 0
    output_arr = []
    ind = 0
    numzeroes = 0
    for label1 in range(12):
        output_arr.append([])
        for label2 in labels:
            val = data[label1 + 1][label2]
            if val is not None:
                output_arr[ind].append(val)
                max_val = max(max_val, val)
                min_val = min(min_val, val)
            else:
                output_arr[ind].append(0)
        ind += 1
    return [max_val, min_val, output_arr]

ticks_x = []
for i in range(20):
    ticks_x.append(i + 0.5)
ticks_y = []
for i in range(12):
    ticks_y.append(i + 0.5)
meansv = buildValidArray(means)
stdevsv = buildValidArray(stdevs)
viridis = cm.get_cmap('viridis', 256)
data = meansv[2]
fig, axs = plt.subplots(1, 1, figsize=(5, 4), constrained_layout=True)

# To visualize means, leave as is. To visualize standard deviations,
# comment out lines 57, 67, 68 and uncomment lines 58, 69, 70

psm = axs.pcolormesh(data, cmap=viridis, rasterized=True, vmin=meansv[1], vmax=meansv[0])
#psm = axs.pcolormesh(data, cmap=viridis, rasterized=True, vmin=stdevsv[1], vmax=stdevsv[0])

plt.xticks(ticks=ticks_x)
plt.yticks(ticks=ticks_y)
axs.set_xticklabels(labels)
axs.set_xlabel('Amino Acid after Mutation')
axs.set_yticklabels([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])
axs.set_ylabel('Position of Mutation')

fig.suptitle('Average Change in Binding Affinity after One Mutation')
fig.colorbar(psm, ax=axs, label='Average Change')
#fig.suptitle('Standard Deviation of Change in Binding Affinity')
#fig.colorbar(psm, ax=axs, label='Standard Deviation')
plt.show()