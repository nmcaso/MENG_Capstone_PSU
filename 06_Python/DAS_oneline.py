"""
Paul Klippel
April 2022
SIMBA lab
This is a no loop algorithm after initial delay matrix setup for Delay and Sum reconstructions
"""

from scipy.io import loadmat
import numpy as np
import plotly.express as plt
import time as time

'''
Function: setup
inputs:
    filepath: exact file location for the .mat datafile
    xlims: min and max values for the x-axis in the reconstruction region
    ylims: min and max values for the y-axis in the reconstruction region
    res: pixel resolution in the reconstruction region
    
output:
    return_obj: Object containing all necessary variables
        .sens_count: number of transducer elements
        .rfdata: time vector with pressure functions for each element
        .sens_xlocs: x-position for each element
        .sens_ylocs: y-position for each element
        .c0: medium speed of sound
        .dt: inverse of the sampling frequency
        .xMat: all x-values for each reconstructed pixel
        .yMat: all y-values for each reconstructed pixel
        .allIndex: indeces corresponding to each pixel in the layered delay matrix for all elements
            l x m x n, l - corresponding to each element layer, m x n - corresponding to the reconstruction region
'''


class VariableSetup:
    def __init__(self, filepath, xlims, ylims, res):
        dat = loadmat(filepath)
        self.sens_count = np.squeeze(dat['ele'].astype(int))
        self.rfdata = np.array(dat['rfdata'])
        self.sens_xlocs = np.array(dat['x0']) * 1e-3
        self.sens_ylocs = np.array(dat['z0']) * 1e-3
        self.c0 = dat['sos'] * 1e3
        self.dt = 1 / (dat['fs'] * 1e6)
        xlength = np.arange(xlims[0], xlims[1], res).size
        ylength = np.arange(ylims[0], ylims[1], res).size
        self.xMat = np.transpose(np.tile(np.arange(xlims[0], xlims[1], res), (ylength, 1)))
        self.yMat = np.tile(np.arange(ylims[0], ylims[1], res), (xlength, 1))
        self.allIndex = np.zeros((self.sens_count, xlength, ylength), dtype='uint32')

    def initialize_index(self):
        shift = 0
        for i in np.arange(0, self.sens_count):
            self.allIndex[i, :, :] = np.round(np.sqrt(np.square(self.xMat - self.sens_xlocs[0, i]) +
                                         np.square(self.yMat - self.sens_ylocs[0, i])) / self.c0 / self.dt)+shift
            shift += len(self.rfdata)


'''
Main
'''
xlimits = [-0.015, 0.015]
ylimits = [-0.015, 0.015]
dx = 31e-6

setup_vars = VariableSetup("leaf.mat", xlimits, ylimits, dx)
setup_vars.initialize_index()

t_start = time.time()
flat_rfdata = setup_vars.rfdata.flatten('F')  # flattened rfdata for linear indexing
img = np.mean(flat_rfdata[setup_vars.allIndex.astype(int)], 0)
print('total time: ', time.time() - t_start)

fig = plt.imshow(img, labels=dict(x='m', y='m'), x=setup_vars.xMat[:, 0].astype(str), y=setup_vars.yMat[0, :].astype(str))
fig.show()
