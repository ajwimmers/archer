# -*- coding: utf-8 -*-
"""
Created on Thu Feb 25 17:10:01 2016

@author: wimmers
"""

from __future__ import print_function
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
#from importlib import reload

import archer.utilities.ScoreFuncs as sftbx


def pcolorCenterShift(xGrid, yGrid, zGrid):
    
    xGrid = np.concatenate((1.5*xGrid[0:1,:] - 0.5*xGrid[1:2,:], 0.5*xGrid[0:-1,:] + 0.5*xGrid[1:,:], 1.5*xGrid[-1:,:] - 0.5*xGrid[-2:-1,:]), axis=0)
    yGrid = np.concatenate((1.5*yGrid[0:1,:] - 0.5*yGrid[1:2,:], 0.5*yGrid[0:-1,:] + 0.5*yGrid[1:,:], 1.5*yGrid[-1:,:] - 0.5*yGrid[-2:-1,:]), axis=0)
    
    xGrid = np.concatenate((1.5*xGrid[:,0:1] - 0.5*xGrid[:,1:2], 0.5*xGrid[:,0:-1] + 0.5*xGrid[:,1:], 1.5*xGrid[:,-1:] - 0.5*xGrid[:,-2:-1]), axis=1)
    yGrid = np.concatenate((1.5*yGrid[:,0:1] - 0.5*yGrid[:,1:2], 0.5*yGrid[:,0:-1] + 0.5*yGrid[:,1:], 1.5*yGrid[:,-1:] - 0.5*yGrid[:,-2:-1]), axis=1)
    
    zGrid = np.concatenate((zGrid, np.nan * zGrid[-1:,:]), axis=0)
    zGrid = np.concatenate((zGrid, np.nan * zGrid[:,-1:]), axis=1)
    
    return xGrid, yGrid, zGrid
    

def pcolorCenterShiftDataOnly(zGrid):
    
    zGrid = np.concatenate((zGrid, np.nan * zGrid[-1:,:]), axis=0)
    zGrid = np.concatenate((zGrid, np.nan * zGrid[:,-1:]), axis=1)
    
    return zGrid
    

def discrete_cmap(N, base_cmap=None):
    
    """Create an N-bin discrete colormap from the specified input map"""
    # Credit to: Jake VanderPlas
    # https://gist.github.com/jakevdp/91077b0cae40f8f8244a
    
    # Usage example: plt.scatter(x, y, c=c, s=50, cmap=discrete_cmap(N, 'cubehelix'))

    base = cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, N))
    
    # Force the first element to be gray (for satellite source plots):
    color_list[0] = (8./15., 8./15., 8./15., 1.0)
    
    cmap_name = base.name + str(N)
    
    return base.from_list(cmap_name, color_list, N)


def plot_diag_3panel(image, attrib, in_dict, out_dict, score_dict, display_filename=None):


    # Channel-specific settings
    if attrib['archer_channel_type'] == '89GHz':

        disp_grid0 = image['bt_grid']
        color_lo = [160, 160, 160]
        color_hi = [280, 280, 280]
        color_label = 'Brightness temp, K'

    elif attrib['archer_channel_type'] == '37GHz':

        disp_grid0 = image['bt_grid']
        color_lo = [160, 160, 160]
        color_hi = [280, 280, 280]
        color_label = 'Brightness temp, K'

    elif attrib['archer_channel_type'] == '183GHz':

        disp_grid0 = image['bt_grid']
        color_lo = [160, 160, 160]
        color_hi = [280, 280, 280]
        color_label = 'Brightness temp, K'

    elif attrib['archer_channel_type'] == 'IR':

        disp_grid0 = image['bt_grid']
        color_lo = [180, 180, 180]
        color_hi = [300, 300, 300]
        color_label = 'Brightness temp, K'

    elif attrib['archer_channel_type'] == 'SWIR':

        disp_grid0 = image['bt_grid']
        color_lo = [180, 180, 180]
        color_hi = [300, 300, 300]
        color_label = 'Brightness temp, K'

    elif attrib['archer_channel_type'] == 'Vis' or attrib['archer_channel_type'] == 'DNB':

        disp_grid0 = image['data_grid']
        color_lo = [0,   150, 150]
        color_hi = [255, 320, 320]
        color_label = 'Pseudo brightness temp, K'

    else:
        print('No real channel type to display!')


    # Shorthand terms
    op_lon = in_dict['op_lon']
    op_lat = in_dict['op_lat']
    op_vmax = in_dict['op_vmax']
    alpha = out_dict['alpha_parameter']
    ring_weight = 1 #in_dict['ring_weight'] # Helps you see something

    # Turn interactive plotting off
    plt.ioff()

    fig = plt.figure()

    # Fix for pcolor offsets
    lon_mx_disp, lat_mx_disp, bt_mx_disp = pcolorCenterShift(
        in_dict['lon_mx'], in_dict['lat_mx'], disp_grid0)

    gs = gridspec.GridSpec(1, 3, width_ratios=(1, 1, 1))

    ax0 = plt.subplot(gs[0,0])
    # This explains the nan-related runtime error for the following line: https://github.com/matplotlib/matplotlib/issues/9892
    img0 = plt.pcolormesh(lon_mx_disp, lat_mx_disp, bt_mx_disp, vmin=color_lo[0], vmax=color_hi[0], cmap='jet_r')
    #img0 = plt.scatter(lon_mx, lat_mx, c=bt_noco_mx, vmin=160, vmax=280, cmap='jet_r', edgecolors='none', marker=',')
    spiral_score_grid_disp = score_dict['spiral_score_grid'] + 0
    spiral_score_grid_disp[spiral_score_grid_disp < -1e8] = np.nan
    plt.contour(score_dict['lon_grid1'], score_dict['lat_grid1'], spiral_score_grid_disp, 
        np.arange(0,50), colors='gray')
    plt.plot(op_lon, op_lat, 'w+')
    plt.title('Spiral score contours, alpha = %4.1f' % (alpha))

    aspect_lon = np.cos(np.pi/180*op_lat) # Aspect ratio
    plt.xlim(op_lon-2.5/aspect_lon, op_lon+2.5/aspect_lon)
    plt.ylim(op_lat-2.5, op_lat+2.5)

    # Colorbar, h/t https://joseph-long.com/writing/colorbars/
    divider = make_axes_locatable(ax0)
    cax0 = divider.append_axes("bottom", size="5%", pad=0.30)
    fig.colorbar(img0, cax=cax0, orientation='horizontal')

    # The next two are plotted on the product grid:
    # Fix for pcolor offsets
    lon_grid_disp, lat_grid_disp, bt_grid_disp = pcolorCenterShift(
        score_dict['lon_grid1'], score_dict['lat_grid1'], score_dict['data_grid1'] )

    ax1 = plt.subplot(gs[0,1])
    # This explains the nan-related runtime error for the following line: https://github.com/matplotlib/matplotlib/issues/9892
    img1 = plt.pcolormesh(lon_grid_disp, lat_grid_disp, bt_grid_disp, vmin=color_lo[1], vmax=color_hi[1], cmap='jet_r')
    plt.xlim(np.min(lon_grid_disp[0,:]), np.max(lon_grid_disp[0,:]))
    plt.ylim(np.min(lat_grid_disp[:,0]), np.max(lat_grid_disp[:,0]))

    ring_score_grid_disp = score_dict['ring_score_grid'] + 0
    ring_score_grid_disp[ring_score_grid_disp == 0] = np.nan
    plt.contour(score_dict['lon_grid1'], score_dict['lat_grid1'], ring_weight * ring_score_grid_disp, 
        np.arange(0,50), colors='gray')
    #plt.plot(op_lon, op_lat, 'w+')
    plt.title(('%s      Vmax = %d kt\nRing score contours') % 
        (attrib['archer_channel_type'], op_vmax) )

    ring_score_max_idx = np.argmax(score_dict['ring_score_grid'])
    i_rs_max, j_rs_max = sftbx.ind2sub(np.shape(score_dict['ring_score_grid']), ring_score_max_idx)
    rs_max_lon = score_dict['lon_grid1'][i_rs_max, j_rs_max]
    rs_max_lat = score_dict['lat_grid1'][i_rs_max, j_rs_max]
    plt.plot(rs_max_lon, rs_max_lat, 'ws', markerfacecolor='none')


    # Colorbar, h/t https://joseph-long.com/writing/colorbars/
    divider = make_axes_locatable(ax1)
    cax1 = divider.append_axes("bottom", size="5%", pad=0.30)
    cbar = fig.colorbar(img1, cax=cax1, orientation='horizontal')
    cbar.set_label(color_label)


    ax2 = plt.subplot(gs[0,2])
    # This explains the nan-related runtime error for the following line: https://github.com/matplotlib/matplotlib/issues/9892
    img2 = plt.pcolormesh(lon_grid_disp, lat_grid_disp, bt_grid_disp, vmin=color_lo[2], vmax=color_hi[2], cmap='jet_r')
    plt.xlim(np.min(lon_grid_disp[0,:]), np.max(lon_grid_disp[0,:]))
    plt.ylim(np.min(lat_grid_disp[:,0]), np.max(lat_grid_disp[:,0]))

    if out_dict['uses_target']:
        combo_score_grid_disp = score_dict['combo_score_grid']
        with np.errstate(invalid='ignore'):
            combo_score_grid_disp[combo_score_grid_disp < -1e8] = np.nan
        plt.contour(score_dict['lon_grid1'], score_dict['lat_grid1'], combo_score_grid_disp, 
            np.arange(0,50), colors='gray')
        plt.plot(op_lon, op_lat, 'w+')
        plt.plot(out_dict['center_lon'], out_dict['center_lat'], 'ws', markerfacecolor='none')

        if not np.isnan(out_dict['center_lon']):
            ring_ang_arr = np.arange(0, 361, 5)
            lon_scale = np.cos(out_dict['center_lat'] * np.pi/180)
            ring_radius_deg = out_dict['ring_radius_deg']
            ring_lon_arr = out_dict['center_lon'] + np.cos(ring_ang_arr * np.pi/180) * ring_radius_deg / lon_scale
            ring_lat_arr = out_dict['center_lat'] + np.sin(ring_ang_arr * np.pi/180) * ring_radius_deg 
            plt.plot(ring_lon_arr, ring_lat_arr, 'm-')

    plt.title('Combined score contours')

    # Colorbar, h/t https://joseph-long.com/writing/colorbars/
    divider = make_axes_locatable(ax2)
    cax2 = divider.append_axes("bottom", size="5%", pad=0.30)
    fig.colorbar(img2, cax=cax2, orientation='horizontal')

    #print(np.round(ring_weight * score_dict['ring_score_grid'] * 100)[::4, ::4].astype('int'))

    # Colorbar, h/t https://joseph-long.com/writing/colorbars/
    divider = make_axes_locatable(ax2)
    cax2 = divider.append_axes("bottom", size="5%", pad=0.30)
    fig.colorbar(img2, cax=cax2, orientation='horizontal')


    fig.set_size_inches(15, 5)
    #plt.show()
    plt.savefig(display_filename, dpi=100)

    return 1