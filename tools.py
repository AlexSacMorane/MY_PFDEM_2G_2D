# -*- encoding=utf-8 -*-

from pathlib import Path
import shutil, os
import numpy as np

#------------------------------------------------------------------------------------------------------------------------------------------ #

def create_folder(name):
    '''
    Create a new folder. If it already exists, it is erased.
    '''
    if Path(name).exists():
        shutil.rmtree(name)
    os.mkdir(name)

#------------------------------------------------------------------------------------------------------------------------------------------ #

def tuplet_to_list(tuplet):
    '''
    Convert a tuplet into lists.
    '''
    L_x = []
    L_y = []
    p_center = np.array([0,0])
    n_mean = 0
    for v in tuplet:
        L_x.append(v[0])
        L_y.append(v[1])
        p_center = p_center + np.array([v[0], v[1]])
        n_mean = n_mean + 1
    L_x.append(L_x[0])
    L_y.append(L_y[0])
    p_center = p_center/n_mean
    # translate center to the point (0,0)
    for i in range(len(L_x)):
        L_x[i] = L_x[i] - p_center[0]
        L_y[i] = L_y[i] - p_center[1]
    return L_x, L_y
