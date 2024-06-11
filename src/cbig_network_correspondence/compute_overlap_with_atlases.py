""""Modules used for calculating overlap between data and atlases"""
import copy
from os import path
from os import listdir
from os import makedirs
from natsort import natsorted
import numpy as np
import nibabel as nib
from brainspace.null_models import SpinPermutations
import pickle
from . import grab_data_info as grab_info
from . import visualize_overlap_lib as vis
from . import visualize_report_lib as vis_report
from scipy.io import loadmat
from regfusion import vol_to_fsaverage
import pandas as pd

# define path
PROJECT_PATH = path.dirname(path.abspath(__file__))
ATLASES_PATH = path.join(PROJECT_PATH, "data", "atlases")
NETWORK_ASSIGN_PATH = path.join(PROJECT_PATH, "data", "network_assignment")
RESULTS_PATH = path.join(PROJECT_PATH, "data", "overlap_results")
NET_NAME_PATH = path.join(PROJECT_PATH, "data", "network_names")
ATLAS_LIST_PATH = path.join(PROJECT_PATH, "data", "atlas_list")
ROTAION_PATH = path.join(PROJECT_PATH, "data", "spin_rotations")


class DataParams:
    """
    This defines the parameter structures for a given data
    """
    def __init__(self, config, data_path):
        self.data_path = data_path
        self.project_path = PROJECT_PATH
        self.atlases_path = ATLASES_PATH
        self.network_assign_path = NETWORK_ASSIGN_PATH

        #if config is a file path
        if path.isfile(config):
            self.config = grab_info.read_config(config)
        #if config is a config object
        else:
            self.config = config

class RefAtlasParams:
    """
    This defines the parameter structures for an existing atlas.
    This atlas is served as the reference atlas.
    """
    def __init__(self, atlas_name):
        atlas_config = path.join(PROJECT_PATH,
                                 "data",
                                 "atlas_config",
                                 atlas_name)
        self.config = grab_info.read_config(atlas_config)
        self.data_path = path.join(ATLASES_PATH,
                                   self.config.space,
                                   self.config.category)
        self.project_path = PROJECT_PATH
        self.atlases_path = ATLASES_PATH
        self.network_assign_path = NETWORK_ASSIGN_PATH


class AtlasParams:
    """
    This defines the parameter structures for an existing atlas.
    This atlas is not served as the reference atlas.
    """
    def __init__(self, atlas_name, atlas_space):
        atlas_config = path.join(PROJECT_PATH,
                                 "data",
                                 "atlas_config",
                                 atlas_name)
        self.config = grab_info.read_config(atlas_config)
        self.data_path = path.join(ATLASES_PATH,
                                   atlas_space,
                                   self.config.category)
        self.project_path = PROJECT_PATH
        self.atlases_path = ATLASES_PATH
        self.network_assign_path = NETWORK_ASSIGN_PATH


def sort_files(data_path, filename):
    """
    Reads in files under the path and sort the files based on a natural order.
    """
    if isinstance(data_path, list):
        files_sorted_res = data_path
        return files_sorted_res 
    elif path.isdir(data_path):
        allfiles = listdir(data_path)
        files = []
        for eachfile in allfiles:
            if filename in eachfile:
                files.append(eachfile)
        files_sorted = natsorted(files)
        files_sorted_res = [path.join(data_path, str(i)) for i in files_sorted]
        return files_sorted_res
    elif path.isfile(data_path):
        files_sorted_res = [data_path]
        return files_sorted_res
    else:
        #print error message
        print("ERROR! Please provide a valid data path!")
        sys.exit(1)


def read_file(data_path):
    """
    Reads in data with given path

    data_path:
        the location to the data
    """

    if path.isfile(data_path):
        if ".mat" in data_path:
            mat_data = loadmat(data_path)
            for key in mat_data:
                # find variable in the mat file with name starting with "lh" or "rh"
                if key.startswith('lh'):
                    mat_data_lh = mat_data.get(key)
                elif key.startswith('rh'):
                    mat_data_rh = mat_data.get(key)
            # Transpose lh_labels and rh_labels if needed
            if mat_data_lh.shape[0] < mat_data_lh.shape[1]:
                mat_data_lh = mat_data_lh.T
            if mat_data_rh.shape[0] < mat_data_rh.shape[1]:
                mat_data_rh = mat_data_rh.T
            data = np.concatenate((mat_data_lh, mat_data_rh))
        elif ".npy" in data_path:
            data = np.load(data_path)
        elif ".nii" in data_path:
            nii_img = nib.load(data_path)
            data = nii_img.get_fdata()
            if data.shape[-1] == 1:
                # Remove the singleton dimension
                data = np.squeeze(data, axis=-1)

    return data


def read_atlas(params):
    """
    Reads in atlas based on config information

    self.project_path = project_path
    self.atlases_path = atlases_path
    self.network_assign_path = network_assign_path
    params.config:
    - params.config.category:
        the categopry of atlas or data
    - params.config.name:
        the name of atlas or data
    - params.config.space:
        the space of atlas or data
    - params.config.type:
        the data type of atlas or data
    - params.config.netassign: [optional]
        the network assignment flag of atlas or data
    - params.config.threshold: [optional]
        the threshold of soft parcellation or metric data
    """

    if params.config.type == 'Hard':
        print('This is a hard parcellation.')
        params.data_sort_path = sort_files(params.data_path, params.config.name)
        data = read_file(params.data_sort_path[0])
        # The hard parcellation is an areal-level parcellation
        # We need to do network assignment based on the povided mapping
        if params.config.netassign:
            datatmp = copy.deepcopy(data)
            mapfile = path.join(params.network_assign_path, params.config.name
                                + ".mat")
            mapfile_data = loadmat(mapfile)
            mapping = mapfile_data.get('mapping')
            #transpose mapping if needed
            if mapping.shape[0] > mapping.shape[1]:
                mapping = mapping.T
            for roi in range(1, np.max(mapping.shape)+1):
                data[datatmp == roi] = mapping[0][roi-1]
            print("Finish network assignmeent!")

        # For data in the volumetric space
        # Apply cortical mask
        if "mm" in params.config.space:
            corticalmask = path.join(params.project_path, "data",
                                     "cortical_masks",
                                     params.config.space + ".nii.gz")
            mask = read_file(corticalmask)
            data[mask == 0] = 0

    elif params.config.type in ('Soft', 'Metric'):
        print('This is a soft parcellation or metric data.')
        data = []
        # Read in cortical mask for data in volumetric space
        if "mm" in params.config.space:
            corticalmask = path.join(params.project_path, "data",
                                     "cortical_masks",
                                     params.config.space + ".nii.gz")
            mask = read_file(corticalmask)
        params.data_sort_path = sort_files(params.data_path, params.config.name)
        for curr_data_file in params.data_sort_path:
            curr_data = read_file(curr_data_file)
            # Apply cortical mask
            if "mm" in params.config.space:
                curr_data[mask == 0] = 0
            # Apply threshold
            if params.config.threshold is not None:
                th_l = params.config.threshold[0]
                th_h = params.config.threshold[1]
                curr_data = np.where(
                    np.logical_and(curr_data > th_l, curr_data < th_h), 1, 0)
            data.append(curr_data)

    return data

def proj_atlas(params):
    """
    Project FSLMNI152 atlas to the surface space

    self.project_path = project_path
    self.atlases_path = atlases_path
    self.network_assign_path = network_assign_path
    params.config:
    - params.config.category:
        the categopry of atlas or data
    - params.config.name:
        the name of atlas or data
    - params.config.space:
        the space of atlas or data
    - params.config.type:
        the data type of atlas or data
    - params.config.netassign: [optional]
        the network assignment flag of atlas or data
    - params.config.threshold: [optional]
        the threshold of soft parcellation or metric data
    """

    if params.config.type == 'Hard':
        print('This is a hard parcellation.')
        params.data_sort_path = sort_files(params.data_path, params.config.name)
        #get the parent directory of the data
        out_path = path.dirname(params.data_sort_path[0])
        temp_path = path.join(out_path, "temp")
        #create temp directory if it does not exist
        if not path.exists(temp_path):
            makedirs(temp_path)
        lh, rh = vol_to_fsaverage(params.data_sort_path[0], temp_path)
        lh_data = nib.load(lh)
        rh_data = nib.load(rh)
        lh_data_fs = lh_data.get_fdata()
        rh_data_fs = rh_data.get_fdata()
        #extract first 40962 elements from lh_data_fs and rh_data_fs
        lh_data_fs = lh_data_fs[:40962]
        rh_data_fs = rh_data_fs[:40962]
        data = np.concatenate((lh_data_fs, rh_data_fs), axis=0)
        data = np.squeeze(data)
        data = np.reshape(data, (data.shape[0], 1))
        out_path = path.join(temp_path, params.config.name + '.npy')
        np.save(out_path, data)

    elif params.config.type in ('Soft', 'Metric'):
        print('This is a soft parcellation or metric data.')
        data = []        
        params.data_sort_path = sort_files(params.data_path, params.config.name)
        out_path = path.dirname(params.data_sort_path[0])
        temp_path = path.join(out_path, "temp")
        out_path = path.join(out_path, "temp", "surf_data")
        #create temp directory if it does not exist
        if not path.exists(out_path):
            makedirs(out_path)
        for curr_data_file in params.data_sort_path:
            lh, rh = vol_to_fsaverage(curr_data_file, temp_path)
            lh_data = nib.load(lh)
            rh_data = nib.load(rh)
            lh_data_fs = lh_data.get_fdata()
            rh_data_fs = rh_data.get_fdata()
            #extract first 40962 elements from lh_data_fs and rh_data_fs
            lh_data_fs = lh_data_fs[:40962]
            rh_data_fs = rh_data_fs[:40962]
            curr_data = np.concatenate((lh_data_fs, rh_data_fs), axis=0)
            curr_data = np.squeeze(curr_data)
            curr_data = np.reshape(curr_data, (curr_data.shape[0], 1))
            np.save(path.join(out_path, path.basename(curr_data_file) + params.config.name + '.npy'), curr_data)

    return out_path

def if_atlas_exist(atlas_name):
    atlas_list_file = open(ATLAS_LIST_PATH,"r", encoding="utf8")
    atlas_names = atlas_list_file.read().splitlines()
    if atlas_name in atlas_names:
        return True
    else:
        return False

def dice_coeff(net1, net2):
    """
    Computes Dice coefficient between two arrays
    """
    intersection = np.sum(net1[net2 == 1])
    coefficient = (2.0 * intersection) / (np.sum(net1) + np.sum(net2))
    return coefficient

def permutation_pvalue(overlap_val, overlap_val_spin):
    """
    Computes spin permuation p-value for overlap matrix
    """
    diff_val = [overlap_val - overlap_val_spin[:,:,i] for i in range(overlap_val_spin.shape[2])]
    p_val = (np.sum(np.array(diff_val) < 0, axis=0) + 1)/(overlap_val_spin.shape[2] + 1)
    
    return p_val

def compute_overlap(ref_params, other_params):
    """
    Computes network overlap between reference data or atlas with
    another atlas.
    """
    other_params.config.space = copy.deepcopy(ref_params.config.space)
    ref_data = read_atlas(ref_params)
    other_data = read_atlas(other_params)
    if ref_params.config.type == 'Hard':
        # Reference atlas/data is a hard parcellation
        ref_net_uni = np.unique(ref_data.flatten())
        # Exclude regions labeled as 0
        ref_net_uni = np.setdiff1d(ref_net_uni, np.array([0]))
        if other_params.config.type == 'Hard':
            # Other atlas is also a hard parcellation

            other_net_uni = np.unique(other_data.flatten())
            # Exclude regions labeled as 0
            other_net_uni = np.setdiff1d(other_net_uni, np.array([0]))
            idx_i = 0
            overlap_val = np.zeros((len(ref_net_uni), len(other_net_uni)))
            for each_ref_net in ref_net_uni:
                curr_ref_net = ref_data == each_ref_net
                idx_j = 0
                for each_other_net in other_net_uni:
                    curr_other_net = other_data == each_other_net
                    overlap_val[idx_i][idx_j] = dice_coeff(curr_ref_net,
                                                           curr_other_net)
                    idx_j = idx_j + 1
                idx_i = idx_i + 1
        elif other_params.config.type in ('Soft', 'Metric'):
            other_net_uni = range(0, len(other_data))
            idx_i = 0
            overlap_val = np.zeros((len(ref_net_uni), len(other_net_uni)))
            for each_ref_net in ref_net_uni:
                curr_ref_net = ref_data == each_ref_net
                idx_j = 0
                for each_other_net in other_net_uni:
                    curr_other_net = other_data[each_other_net]
                    overlap_val[idx_i][idx_j] = dice_coeff(curr_ref_net,
                                                           curr_other_net)
                    idx_j = idx_j + 1
                idx_i = idx_i + 1
    elif ref_params.config.type in ('Soft', 'Metric'):
        # Reference atlas/data is a soft parcellation or any metric data
        ref_net_uni = range(0, len(ref_data))
        if other_params.config.type == 'Hard':
            # Other atlas is also a hard parcellation
            other_net_uni = np.unique(other_data.flatten())
            # Exclude regions labeled as 0
            other_net_uni = np.setdiff1d(other_net_uni, np.array([0]))
            idx_i = 0
            overlap_val = np.zeros((len(ref_net_uni), len(other_net_uni)))
            for each_ref_net in ref_net_uni:
                curr_ref_net = ref_data[each_ref_net]
                idx_j = 0
                for each_other_net in other_net_uni:
                    curr_other_net = other_data == each_other_net
                    overlap_val[idx_i][idx_j] = dice_coeff(curr_ref_net,
                                                           curr_other_net)
                    idx_j = idx_j + 1
                idx_i = idx_i + 1
        elif other_params.config.type in ('Soft', 'Metric'):
            other_net_uni = range(0, len(other_data))
            idx_i = 0
            overlap_val = np.zeros((len(ref_net_uni), len(other_net_uni)))
            for each_ref_net in ref_net_uni:
                curr_ref_net = ref_data[each_ref_net]
                idx_j = 0
                for each_other_net in other_net_uni:
                    curr_other_net = other_data[each_other_net]
                    overlap_val[idx_i][idx_j] = dice_coeff(curr_ref_net,
                                                           curr_other_net)
                    idx_j = idx_j + 1
                idx_i = idx_i + 1

    return overlap_val

def compute_overlap_spin(ref_params, other_params):
    """
    Computes network overlap between reference data or atlas with
    another atlas.
    """
    # prepare spin rotations
    rotation_file = path.join(ROTAION_PATH, ref_params.config.space, "1000_spin_permutations_state0.pkl")
    spin_rots = pickle.load(open(rotation_file, "rb"))
    if ref_params.config.space == "fsaverage6":
        hemi_size = 40962
    elif ref_params.config.space == "fs_LR_32k":
        hemi_size = 32492
    
    
    other_params.config.space = copy.deepcopy(ref_params.config.space)
    orig_ref_data = read_atlas(ref_params)
    other_data = read_atlas(other_params)

    # roate ref_data
    # if ref_data is a list, then it is a soft parcellation or metric data, rotate each element
    if isinstance(orig_ref_data, list):
        ref_data_rotated = []
        for each_ref_data in orig_ref_data:
            ref_data_rotated.append(np.hstack(spin_rots.randomize(each_ref_data[:hemi_size],each_ref_data[hemi_size:])))        
    else:
        ref_data_rotated = np.hstack(spin_rots.randomize(orig_ref_data[:hemi_size],orig_ref_data[hemi_size:]))
        
    if ref_params.config.type == 'Hard':
        # Reference atlas/data is a hard parcellation
        ref_net_uni = np.unique(orig_ref_data.flatten())
        # Exclude regions labeled as 0
        ref_net_uni = np.setdiff1d(ref_net_uni, np.array([0]))
        if other_params.config.type == 'Hard':
            # Other atlas is also a hard parcellation

            other_net_uni = np.unique(other_data.flatten())
            # Exclude regions labeled as 0
            other_net_uni = np.setdiff1d(other_net_uni, np.array([0]))
            idx_i = 0
            overlap_val = np.zeros((len(ref_net_uni), len(other_net_uni),1000))
            for each_ref_net in ref_net_uni:
                curr_ref_net = ref_data_rotated == each_ref_net
                idx_j = 0
                for each_other_net in other_net_uni:
                    curr_other_net = other_data == each_other_net
                    overlap_val[idx_i][idx_j] = [dice_coeff(curr_ref_net[i_rot,:],curr_other_net) for i_rot in range(curr_ref_net.shape[0])]
                    idx_j = idx_j + 1
                idx_i = idx_i + 1
        elif other_params.config.type in ('Soft', 'Metric'):
            other_net_uni = range(0, len(other_data))
            idx_i = 0
            overlap_val = np.zeros((len(ref_net_uni), len(other_net_uni)))
            for each_ref_net in ref_net_uni:
                curr_ref_net = ref_data_rotated == each_ref_net
                idx_j = 0
                for each_other_net in other_net_uni:
                    curr_other_net = other_data[each_other_net]
                    overlap_val[idx_i][idx_j] = [dice_coeff(curr_ref_net[i_rot,:],curr_other_net) for i_rot in range(curr_ref_net.shape[0])]

                    idx_j = idx_j + 1
                idx_i = idx_i + 1
    elif ref_params.config.type in ('Soft', 'Metric'):
        # Reference atlas/data is a soft parcellation or any metric data
        ref_net_uni = range(0, len(orig_ref_data))
        if other_params.config.type == 'Hard':
            # Other atlas is also a hard parcellation
            other_net_uni = np.unique(other_data.flatten())
            # Exclude regions labeled as 0
            other_net_uni = np.setdiff1d(other_net_uni, np.array([0]))
            idx_i = 0
            overlap_val = np.zeros((len(ref_net_uni), len(other_net_uni), 1000))
            for each_ref_net in ref_net_uni:
                curr_ref_net = ref_data_rotated[each_ref_net]
                idx_j = 0
                for each_other_net in other_net_uni:
                    curr_other_net = other_data == each_other_net
                    overlap_val[idx_i][idx_j] = [dice_coeff(curr_ref_net[i_rot,:],curr_other_net) for i_rot in range(curr_ref_net.shape[0])]
                    idx_j = idx_j + 1
                idx_i = idx_i + 1
        elif other_params.config.type in ('Soft', 'Metric'):
            other_net_uni = range(0, len(other_data))
            idx_i = 0
            overlap_val = np.zeros((len(ref_net_uni), len(other_net_uni), 1000))
            for each_ref_net in ref_net_uni:
                curr_ref_net = ref_data_rotated[each_ref_net]
                idx_j = 0
                for each_other_net in other_net_uni:
                    curr_other_net = other_data[each_other_net]
                    overlap_val[idx_i][idx_j] = [dice_coeff(curr_ref_net[i_rot,:],curr_other_net) for i_rot in range(curr_ref_net.shape[0])]
                    idx_j = idx_j + 1
                idx_i = idx_i + 1

    return overlap_val

def compute_overlap_data(data_params, atlas_name):
    """
    Computes overlap matrix for a given data and an existing atlas
    """
    # Reads in the atlas in the same space as the data
    atlas_params = AtlasParams(atlas_name, data_params.config.space)

    # Compute the overlap matrix between data and the atlas
    overlap_val = compute_overlap(data_params, atlas_params)

    return overlap_val

def compute_overlap_data_all(data_params, atlas_names=None):
    """
    Computes overlap matrix for a given data and all existing atlases
    """
    # Reads in the atlas in the same space as the data
    if atlas_names is None:
        atlas_list_file = open(ATLAS_LIST_PATH,"r", encoding="utf8")
        atlas_names = atlas_list_file.read().splitlines()
    overlap_val = {}
    for atlas_name in atlas_names:
        print("Computing overlap with " + atlas_name)
        if not if_atlas_exist(atlas_name):
            print("Atlas " + atlas_name + " does not exist.")
            sys.exit(1)
        atlas_params = AtlasParams(atlas_name, data_params.config.space)
        # Compute the overlap matrix between data and the atlas
        overlap_val[atlas_name] = compute_overlap(data_params, atlas_params)
    return overlap_val

def permute_overlap_data_all(data_params, atlas_names):
    """
    Computes overlap matrix for a given data and all existing atlases
    """
    # if data is in volume project it to surface
    if data_params.config.space == "FSLMNI2mm":
        data_params.config.space = "fsaverage6"
        #data_params.config.name = "input_data"
        data_path = proj_atlas(data_params)
        data_params.data_path = data_path

    if atlas_names is None:
        # Reads in the atlas in the same space as the data
        atlas_list_file = open(ATLAS_LIST_PATH,"r", encoding="utf8")
        atlas_names = atlas_list_file.read().splitlines()
    elif isinstance(atlas_names, str):
        atlas_names = [atlas_names]
    p_val = {}
    for atlas_name in atlas_names:
        print("Performing permutation with " + atlas_name)
        if not if_atlas_exist(atlas_name):
            print("Atlas " + atlas_name + " does not exist.")
            sys.exit(1)

        atlas_params = AtlasParams(atlas_name, data_params.config.space)
        
        # Compute the overlap matrix between data and the atlas
        overlap_val = compute_overlap(data_params, atlas_params)
        overlap_val_spin = compute_overlap_spin(data_params, atlas_params)
        p_val[atlas_name] = permutation_pvalue(overlap_val, overlap_val_spin)
    return p_val
    
def load_overlap_atlases(ref_atlas_name, other_atlas_name):
    """
    Load overlap matrix for any pair of existing atlases
    """
    overlap_file = path.join(RESULTS_PATH,
                             ref_atlas_name,
                             other_atlas_name + '.mat')
    data = loadmat(overlap_file)
    overlap_val = np.array(data.get('overlap'))

    return overlap_val


def locate_label_name(params, specify_net_name):
    """
    Returns the binarized data for a specific network of an atlas
    """
    data = read_atlas(params)
    network_names = vis.construct_network_name(params.config.name)
    if specify_net_name in network_names:
        index = network_names.index(specify_net_name)
        if params.config.type == 'Hard':
            net_uni = np.unique(data.flatten())
            # Exclude regions labeled as 0
            net_uni = np.setdiff1d(net_uni, np.array([0]))
            return (data == net_uni[index]).astype(int)
        else:
            return np.squeeze(data[index])
    elif specify_net_name is None and len(data) == 1:
        return np.squeeze(data)
    else:
        print('Please specify a valid network name.')
        return None


def obtain_overlap(ref_params, other_params, ref_net_name, other_net_name):
    """
    Computes network overlap between reference network or reference data
    and a network from another atlas. Reference network will be labeled as 1,
    other network will be labeled as 2, the overlap area will be labeled as 3.
    """
    other_params.config.space = copy.deepcopy(ref_params.config.space)
    net_data = locate_label_name(ref_params, ref_net_name)
    other_data = locate_label_name(other_params, other_net_name)
    other_data[other_data != 0] = 2
    overlap_data = net_data + other_data
    return overlap_data


def obtain_overlap_data(data_params, atlas_name, atlas_net_name):
    """
    Computes network overlap between binarized data and a network from another
    atlas. The binarized data will be labeled as 1, other network will be
    labeled as 2, the overlap area will be labeled as 3.
    """
    # Reads in the atlas in the same space as the data
    atlas_params = AtlasParams(atlas_name, data_params.config.space)
    overlap_data = obtain_overlap(data_params,
                                  atlas_params,
                                  data_params.config.name,
                                  atlas_net_name)
    return overlap_data


def obtain_overlap_atlases(ref_atlas_name, other_atlas_name,
                           ref_net_name, other_net_name):
    """
    Computes network overlap between a network from reference atlas and
    a network from another atlas. Reference network will be labeled as 1,
    other network will be labeled as 2, the overlap area will be labeled as 3.
    """
    ref_params = RefAtlasParams(ref_atlas_name)
    other_params = AtlasParams(other_atlas_name, ref_params.config.space)
    overlap_data = obtain_overlap(ref_params,
                                  other_params,
                                  ref_net_name,
                                  other_net_name)
    return overlap_data


def data_network_correspondence(ref_params, atlas_names_list):
    overlap_mat_all = compute_overlap_data_all(ref_params, atlas_names_list)
    # check the dictionary overlap_mat_all items
    # for item inside overlap_mat_all.items(), key is the atlas name, value is the overlap matrix
    # the overlap matrix should be a 1xN matrix, where N is the number of networks in the atlas
    # throw error message if the matrix is not 1xN
    for key, value in overlap_mat_all.items():
        if value.shape[0] != 1:
            print("ERROR! The overlap matrix for atlas " + key + " is not a 1xN matrix.")
            sys.exit(1)
    p_val = permute_overlap_data_all(ref_params, atlas_names_list)
    network_names = vis.construct_network_name_all(atlas_names_list)

    df1 = pd.DataFrame([(key, value) for key, values in overlap_mat_all.items() for value in values], columns=['group', 'dice'])
    #flatten df1
    df1 = df1.explode('dice')
    df1 = df1.reset_index(drop=True)
    df2 = pd.DataFrame([(key, value) for key, values in network_names.items() for value in values], columns=['group', 'name'])
    df3 = pd.DataFrame([(key, value) for key, values in p_val.items() for value in values], columns=['group', 'p_value'])
    df3 = df3.explode('p_value')
    df3 = df3.reset_index(drop=True)

    df = pd.concat([df1,df2['name'],df3['p_value']],axis=1)
    df = df[['group','name','dice','p_value']]

    return df

def atlas_network_correspondence(ref_params, atlas_names_list):
    # if ref_params is a string of atlas name
    if isinstance(ref_params, str):
        ref_params = RefAtlasParams(ref_params)
    overlap_mat = compute_overlap_data_all(ref_params, atlas_names_list)
    p_val = permute_overlap_data_all(ref_params, atlas_names_list)
    combined_data = {key: {'Dice': overlap_mat[key], 'P value': p_val[key]} for key in overlap_mat}
    
    return combined_data

def network_correspondence(ref_input, atlas_names_list,out_dir,ref_networks=None):
    # if ref_params is a string of atlas name
    if isinstance(ref_input, str):
        ref_params = RefAtlasParams(ref_input)
    elif isinstance(ref_input, DataParams):
        ref_params = ref_input

    overlap_mat_all = compute_overlap_data_all(ref_params, atlas_names_list)
    # check the dictionary overlap_mat_all items
    # for item inside overlap_mat_all.items(), key is the atlas name, value is the overlap matrix
    # the overlap matrix should be a 1xN matrix, where N is the number of networks in the atlas
    # throw error message if the matrix is not 1xN
    single_flag = True
    for key, value in overlap_mat_all.items():
        if value.shape[0] != 1:
            single_flag = False
    p_val = permute_overlap_data_all(ref_params, atlas_names_list)

    if single_flag:
        print("Single dimension data is provided. Visualize network correspondence as circular chart.")
        network_names = vis.construct_network_name_all(atlas_names_list)
        df1 = pd.DataFrame([(key, value) for key, values in overlap_mat_all.items() for value in values], columns=['group', 'dice'])
        #flatten df1
        df1 = df1.explode('dice')
        df1 = df1.reset_index(drop=True)
        df2 = pd.DataFrame([(key, value) for key, values in network_names.items() for value in values], columns=['group', 'name'])
        df3 = pd.DataFrame([(key, value) for key, values in p_val.items() for value in values], columns=['group', 'p_value'])
        df3 = df3.explode('p_value')
        df3 = df3.reset_index(drop=True)

        df = pd.concat([df1,df2['name'],df3['p_value']],axis=1)
        df = df[['group','name','dice','p_value']]
        vis_report.report_data_network_correspondence(df,out_dir)
    else:
        print("Multiple dimension data are provided. Visualize network correspondence as heatmap.")
        combined_data = {key: {'Dice': overlap_mat_all[key], 'P value': p_val[key]} for key in overlap_mat_all}
        if ref_networks is not None:
            ref_networks_data = ref_networks
        else:
            print("ref_networks is not provided. Use ref_input to get the reference network data.")
            if isinstance(ref_input, str):
                ref_networks_data = ref_input
            elif isinstance(ref_input, DataParams):
                ref_networks_data = ref_input.config.name
            else:
                # print error
                print("ERROR! Cannot get reference networks based on ref_input. Please provide a valid reference network data.")
                sys.exit(1)
        vis_report.report_atlas_network_correspondence(combined_data,ref_networks_data,out_dir)
        
