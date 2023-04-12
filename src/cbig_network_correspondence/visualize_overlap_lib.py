from os import path
import numpy as np
import scipy.io
import seaborn as sns
import matplotlib.pyplot as plt
import copy
from pathlib import Path


# define path
PROJECT_PATH = path.dirname(path.abspath(__file__))
NET_NAME_PATH = path.join(PROJECT_PATH, "data", "network_names")
RESULTS_PATH = path.join(PROJECT_PATH, "data", "overlap_results")
ATLAS_LIST_PATH = path.join(PROJECT_PATH, "data", "atlas_list")

def pair_match(overlap_mat):
    """
    Match each reference network to networks from another atlas, find the best
    matched pair and reorder the network order to shape the overlap_mat to be
    a diagnol shape.
    """
    if overlap_mat.shape[0] > overlap_mat.shape[1]:
        overlap_mat_new = overlap_mat[0:overlap_mat.shape[1],:]
    else:
        overlap_mat_new = overlap_mat[:]   
    order_idx = [0] * overlap_mat_new.shape[0]
    for i in range(overlap_mat_new.shape[0]):
        order_idx[i] = np.argmax(overlap_mat_new[i,:])
        overlap_mat_new[:,order_idx[i]] = -1
    if overlap_mat.shape[0] < overlap_mat.shape[1]:
        remain_idx = [i for i in np.arange(overlap_mat.shape[1]) if i not in order_idx]
        order_idx = list(order_idx) + list(remain_idx)
    return order_idx

def construct_network_name(atlas_input):
    """
    Finds the network names for a given atlas.
    1) atlas_input is a atlas name of available atlases: 
    reads in the atlas network names from the network name directory.
    2) atlas_input is a list:
    reads in the list and treats it as the list of network names.
    3) atlas_input is a path to a file containing network name list:
    reads in the file. Each row of the file is treated as a network name. The ordering
    of network names should be consistent with the input data. Specifically, if the
    input data is a hard parcellation, the network names should follow ascending order
    of the parcellation labels. If the input data is a set of soft parcellation spatial
    maps or a set of metric data files, the network names should follow the natural sort
    order of all file names.
    4) atlas_input is a single network name:
    reads in the single network name directly. The input data should be one single network
    or one single metric data.
    """
    atlas_list_file = open(ATLAS_LIST_PATH,"r", encoding="utf8")
    atlas_names = atlas_list_file.read().splitlines()
    if atlas_input in atlas_names:
        atlas_name_file = path.join(NET_NAME_PATH, atlas_input)
        atlas_name_file = open(atlas_name_file,"r", encoding="utf8")
        network_names = atlas_name_file.read().splitlines()
    elif isinstance(atlas_input, list) and all(isinstance(elem, str) for elem in atlas_input):
        network_names = atlas_input
    elif path.isfile(atlas_input):
        atlas_name_file = open(atlas_input,"r", encoding="utf8")
        network_names = atlas_name_file.read().splitlines()
    elif isinstance(atlas_input, str) and len(atlas_input) > 0:
        network_names = []
        network_names.append(atlas_input)
    else:
        network_names = None
    return network_names

def draw_overlap_mat(overlap_data, ref_atlas_name, other_atlas_name, minv, maxv, figfile):
    """
    Draw a given overlap matrix. 
    """
    overlap_data_reorder = copy.deepcopy(overlap_data)
    order_idx = pair_match(overlap_data_reorder)
    overlap_data = overlap_data[:, order_idx]
    reference_name = construct_network_name(ref_atlas_name)
    other_name = construct_network_name(other_atlas_name)

    other_name_new = [other_name[i] for i in order_idx]

    _, axes = plt.subplots(figsize=(8,8))
    sns.set(font_scale=1)
    sns.set_style("ticks")
    im = sns.heatmap(
            ax=axes,
            data=overlap_data,
            cmap="rocket",
            square=True,
            vmin=minv,
            vmax=maxv,
            linewidth=0.05,
            xticklabels=other_name_new,
            yticklabels=reference_name,
            cbar=False
            )

    frame_len=3.5
    axes.axhline(y=0, color='k',linewidth=frame_len)
    axes.axhline(y=overlap_data.shape[0], color='k',linewidth=frame_len)
    axes.axvline(x=0, color='k',linewidth=frame_len)
    axes.axvline(x=overlap_data.shape[1], color='k',linewidth=frame_len)

    axes.xaxis.tick_top()
    axes.tick_params(axis='x', which='major',labelsize=14,labelrotation=90,pad=10)
    axes.tick_params(axis='y', which='major',labelsize=14,pad=10)
        
    mappable = im.get_children()[0]
    cb=plt.colorbar(mappable,
                ax = axes,
                orientation = 'horizontal',
                location = 'bottom',
                label='Dice overlap',
                aspect=35.5,
                fraction=0.025,
                pad=0.08
                )
    cb.outline.set_linewidth(1.5)
    cb.outline.set_color('black')
    cb.ax.locator_params(nbins=10)
    cb.ax.tick_params(labelsize=10)
    cb.ax.tick_params(width=1.5)
    
    plt.tight_layout()
    plt.ioff()
    plt.savefig(figfile, dpi=200, bbox_inches='tight')

def draw_overlap_atlases(ref_atlas_name, other_atlas_name, minv, maxv, figfile):
    """
    Draw overlap matrix for any pair of existing atlases
    """
    overlap_file = path.join(RESULTS_PATH, ref_atlas_name, other_atlas_name + '.mat')
    data = scipy.io.loadmat(overlap_file)
    overlap_data = np.array(data.get('overlap'))
    draw_overlap_mat(overlap_data, ref_atlas_name, other_atlas_name, minv, maxv, figfile)
        
            
