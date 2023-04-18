from os import path
import numpy as np
import scipy.io
import seaborn as sns
import matplotlib.pyplot as plt
import copy
from nilearn import plotting
from nilearn import image
import matplotlib.colors as colors


# define path
PROJECT_PATH = path.dirname(path.abspath(__file__))
NET_NAME_PATH = path.join(PROJECT_PATH, "data", "network_names")
RESULTS_PATH = path.join(PROJECT_PATH, "data", "overlap_results")
ATLAS_LIST_PATH = path.join(PROJECT_PATH, "data", "atlas_list")
TEMPLATE_PATH = path.join(PROJECT_PATH, "data", "templates")

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
        

def draw_overlap_map_fs_LR_32k(overlap_data, outfig):
    """
    Draw surface maps for overlappping networks on fs_LR_32k surface.
    """
    colors4_list = [(0.7, 0.7, 0.7), (0.89, 0.39, 0.45),
                    (0.21, 0.18, 0.75), (0.55, 0.29, 0.6)]
    cmap4 = colors.ListedColormap(colors4_list)
    lh_surf_mesh = path.join(TEMPLATE_PATH, 'fs_LR_32k', 'lh.very_inflated')
    rh_surf_mesh = path.join(TEMPLATE_PATH, 'fs_LR_32k', 'rh.very_inflated')
    
    lh_fig_l = plotting.plot_surf(surf_mesh=lh_surf_mesh,
                                  surf_map=overlap_data[:32492],
                                  hemi='left', view='lateral', cmap=cmap4)
    lh_fig_m = plotting.plot_surf(surf_mesh=lh_surf_mesh,
                                  surf_map=overlap_data[:32492],
                                  hemi='left', view='medial', cmap=cmap4)
    rh_fig_l = plotting.plot_surf(surf_mesh=rh_surf_mesh,
                                  surf_map=overlap_data[32492:],
                                  hemi='right', view='lateral', cmap=cmap4)
    rh_fig_m = plotting.plot_surf(surf_mesh=rh_surf_mesh,
                                  surf_map=overlap_data[32492:],
                                  hemi='right', view='medial', cmap=cmap4)
   
    lh_fig_l.savefig(path.join(outfig, "_lh_lateral"), dpi=300)
    lh_fig_m.savefig(path.join(outfig, "_lh_medial"), dpi=300)
    rh_fig_l.savefig(path.join(outfig, "_rh_lateral"), dpi=300)
    rh_fig_m.savefig(path.join(outfig, "_rh_medial"), dpi=300)


def draw_overlap_map_fsaverage6(overlap_data, outfig):
    """
    Draw surface maps for overlappping networks on fsaverage6 surface.
    """
    colors4_list = [(0.7, 0.7, 0.7), (0.89, 0.39, 0.45),
                    (0.21, 0.18, 0.75), (0.55, 0.29, 0.6)]
    cmap4 = colors.ListedColormap(colors4_list)
    lh_surf_mesh = path.join(TEMPLATE_PATH, 'fsaverage', 'lh.inflated')
    rh_surf_mesh = path.join(TEMPLATE_PATH, 'fsaverage', 'rh.inflated')

    lh_fig_l = plotting.plot_surf(surf_mesh=lh_surf_mesh,
                                  surf_map=overlap_data[:40962],
                                  hemi='left', view='lateral', cmap=cmap4)
    lh_fig_m = plotting.plot_surf(surf_mesh=lh_surf_mesh,
                                  surf_map=overlap_data[:40962],
                                  hemi='left', view='medial', cmap=cmap4)
    rh_fig_l = plotting.plot_surf(surf_mesh=rh_surf_mesh,
                                  surf_map=overlap_data[40962:],
                                  hemi='right', view='lateral', cmap=cmap4)
    rh_fig_m = plotting.plot_surf(surf_mesh=rh_surf_mesh,
                                  surf_map=overlap_data[40962:],
                                  hemi='right', view='medial', cmap=cmap4)

    lh_fig_l.savefig(path.join(outfig, "_lh_lateral"), dpi=300)
    lh_fig_m.savefig(path.join(outfig, "_lh_medial"), dpi=300)
    rh_fig_l.savefig(path.join(outfig, "_rh_lateral"), dpi=300)
    rh_fig_m.savefig(path.join(outfig, "_rh_medial"), dpi=300)


def draw_overlap_map_FSLMNI2mm(overlap_data, outfig, coords=[-4, -31, 18]):
    """
    Draw volumetric maps for overlappping networks on FSLMNI2mm.
    coords is the coordinates of the point where the cut is performed.
    """
    colors3_list = [(0.89, 0.39, 0.45), (0.21, 0.18, 0.75), (0.55, 0.29, 0.6)]
    cmap3 = colors.ListedColormap(colors3_list)
    bg_image = image.load_img(path.join(TEMPLATE_PATH, 'FSLMNI2mm.nii.gz'))
    roi_image_header = image.load_img(path.join(TEMPLATE_PATH,
                                                'FSLMNI2mm_header.nii.gz'))
    new_image = image.new_img_like(roi_image_header, overlap_data)
    figv = plotting.plot_roi(roi_img=new_image, bg_img=bg_image,
                             cut_coords=coords, cmap=cmap3)
    figv.savefig(path.join(outfig, "_FSLMNI2mm"), dpi=300)


def draw_overlap_map_LairdColin2mm(overlap_data, outfig, coords=[-4, -31, 18]):
    """
    Draw volumetric maps for overlappping networks on LairdColin2mm.
    coords is the coordinates of the point where the cut is performed.
    Note that LairdColin2mm is the space of Laird atlas, we use 
    LairdColin1mm as the background image. 
    """
    colors3_list = [(0.89, 0.39, 0.45), (0.21, 0.18, 0.75), (0.55, 0.29, 0.6)]
    cmap3 = colors.ListedColormap(colors3_list)
    bg_image = image.load_img(path.join(TEMPLATE_PATH, 'LairdColin1mm.nii.gz'))
    roi_image_header = image.load_img(path.join(TEMPLATE_PATH,
                                                'FSLMNI2mm_header.nii.gz'))
    new_image = image.new_img_like(roi_image_header, overlap_data)
    figv = plotting.plot_roi(roi_img=new_image, bg_img=bg_image,
                             cut_coords=coords, cmap=cmap3)
    figv.savefig(path.join(outfig, "_LairdColin2mm"), dpi=300)


def draw_overlap_map_ShenColin1mm(overlap_data, outfig, coords=[-4, -31, 18]):
    """
    Draw volumetric maps for overlappping networks on ShenColin1mm.
    coords is the coordinates of the point where the cut is performed.
    """
    colors3_list = [(0.89, 0.39, 0.45), (0.21, 0.18, 0.75), (0.55, 0.29, 0.6)]
    cmap3 = colors.ListedColormap(colors3_list)
    bg_image = image.load_img(path.join(TEMPLATE_PATH, 'ShenColin1mm.nii.gz'))
    roi_image_header = image.load_img(path.join(TEMPLATE_PATH,
                                                'ShenColin1mm_header.nii.gz'))
    new_image = image.new_img_like(roi_image_header, overlap_data)
    figv = plotting.plot_roi(roi_img=new_image, bg_img=bg_image,
                             cut_coords=coords, cmap=cmap3)
    figv.savefig(path.join(outfig, "_ShenColin1mm"), dpi=300)


def draw_overlap_map(space_name, overlap_data, outfig, coords=[-4, -31, 18]):
    """
    Draw brain maps for overlappping networks. The brain map will be in the
    same space as the reference data/atlas space.
    """
    if space_name == "fs_LR_32k":
        draw_overlap_map_fs_LR_32k(overlap_data, outfig)
    elif space_name == "fsaverage6":
        draw_overlap_map_fsaverage6(overlap_data, outfig)
    elif space_name == "FSLMNI2mm":
        draw_overlap_map_FSLMNI2mm(overlap_data, outfig, coords)
    elif space_name == "LairdColin2mm":
        draw_overlap_map_LairdColin2mm(overlap_data, outfig, coords)
    elif space_name == "ShenColin1mm":
        draw_overlap_map_ShenColin1mm(overlap_data, outfig, coords)
    else:
        print("We current don't provide the function to draw brain \
              maps for" + space_name + "space. Please \
              check nilearn and use the function in nilearn to draw \
              the brain maps.")