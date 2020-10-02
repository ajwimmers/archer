import matplotlib
matplotlib.use('agg') # AJW 16Aug2018
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable

plt.close('all')

matplotlib.rcParams['axes.titlesize'] = 10
matplotlib.rcParams['axes.labelsize'] = 10
matplotlib.rcParams['legend.fontsize'] = 10
matplotlib.rcParams['font.size'] = 10

import numpy as np


from satpy import available_readers, Scene
from glob import glob
available_readers()

filenames = glob('./GIMGO-SVI01*')
scn = Scene(reader='viirs_sdr', filenames=filenames)

scn.load(['I01'])
scn.save_datasets(writer='cf', datasets=['I01'], filename='viirs_sdr_i01_damien_test.nc', 
    exclude_attrs=['raw_metadata'])

scn.load(['I04'])
scn.save_datasets(writer='cf', datasets=['I04'], filename='viirs_sdr_i04_damien_test.nc', 
    exclude_attrs=['raw_metadata'])

scn.load(['I05'])
scn.save_datasets(writer='cf', datasets=['I05'], filename='viirs_sdr_i05_damien_test.nc', 
    exclude_attrs=['raw_metadata'])

#Save test figure
plt.figure()
plt.imshow(scn['I05'][500:2200, 3000:4500], vmin=180, vmax=300, cmap='jet_r')
plt.colorbar()
plt.savefig('test_pysat.png', dpi=100)


filenames = glob('./GDNBO-SVDNB*')
scn = Scene(reader='viirs_sdr', filenames=filenames)

scn.load(['DNB'])
scn.save_datasets(writer='cf', datasets=['DNB'], filename='viirs_sdr_dnb_damien_test.nc', 
    exclude_attrs=['raw_metadata'])

