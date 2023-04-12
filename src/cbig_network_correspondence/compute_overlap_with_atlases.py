""""Modules used for calculating overlap between data and atlases"""
import copy
from os import path
from os import listdir
from natsort import natsorted
import numpy as np
import nibabel as nib
from . import grab_data_info as grab_info
from scipy.io import loadmat
import importlib
importlib.reload(grab_info)

# define path
PROJECT_PATH = path.dirname(path.abspath(__file__))
ATLASES_PATH = path.join(PROJECT_PATH, "data", "atlases")
NETWORK_ASSIGN_PATH = path.join(PROJECT_PATH, "data", "network_assignment")
RESULTS_PATH = path.join(PROJECT_PATH, "data", "overlap_results")

class DataParams:
    """
    This defines the parameter structures for a given data
    """
    def __init__(self, config, data_path):
        self.data_path = data_path
        self.project_path = PROJECT_PATH
        self.atlases_path = ATLASES_PATH
        self.network_assign_path = NETWORK_ASSIGN_PATH
        self.config = grab_info.read_config(config)

class RefAtlasParams:
    """
    This defines the parameter structures for an existing atlas.
    This atlas is served as the reference atlas.
    """
    def __init__(self, atlas_name):
        atlas_config = path.join(PROJECT_PATH, "data", "atlas_config", atlas_name)
        self.config = grab_info.read_config(atlas_config)
        self.data_path = path.join(ATLASES_PATH, self.config.space, self.config.category)
        self.project_path = PROJECT_PATH
        self.atlases_path = ATLASES_PATH
        self.network_assign_path = NETWORK_ASSIGN_PATH

class AtlasParams:
    """
    This defines the parameter structures for an existing atlas.
    This atlas is not served as the reference atlas.
    """
    def __init__(self, atlas_name, atlas_space):
        atlas_config = path.join(PROJECT_PATH, "data", "atlas_config", atlas_name)
        self.config = grab_info.read_config(atlas_config)
        self.data_path = path.join(ATLASES_PATH, atlas_space, self.config.category)
        self.project_path = PROJECT_PATH
        self.atlases_path = ATLASES_PATH
        self.network_assign_path = NETWORK_ASSIGN_PATH

def sort_files(data_path, filename):
    """
    Reads in files under the path and sort the files based on a natural order.
    """
    if path.isdir(data_path):
        allfiles = listdir(data_path)
        files = []
        for eachfile in allfiles:
            if filename in eachfile:
                files.append(eachfile)
        files_sorted = natsorted(files)
        files_sorted_full = [path.join(data_path, str(i)) for i in files_sorted]
        return files_sorted_full

    if path.isfile(data_path):
        return [data_path]

def read_file(data_path):
    """
    Reads in data with given path

    data_path:
        the location to the data
    """

    if path.isfile(data_path):
        if ".mat" in data_path:
            mat_data = loadmat(data_path)
            mat_data_lh = mat_data.get('lh_labels')
            mat_data_rh = mat_data.get('rh_labels')
            data = np.concatenate((mat_data_lh, mat_data_rh))

        elif ".nii" in data_path:
            nii_img = nib.load(data_path)
            data = nii_img.get_fdata()

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
        params.data_path = sort_files(params.data_path, params.config.name)
        data = read_file(params.data_path[0])
        # The hard parcellation is an areal-level parcellation
        # We need to do network assignment based on the povided mapping
        if params.config.netassign:
            datatmp = copy.deepcopy(data)
            mapfile = path.join(params.network_assign_path, params.config.name + ".mat")
            mapfile_data = loadmat(mapfile)
            mapping= mapfile_data.get('mapping')
            for roi in range(1,mapping.shape[1]+1):
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
        print('This is a soft parcellation.')
        data = []
        # Read in cortical mask for data in volumetric space
        if "mm" in params.config.space:
            corticalmask = path.join(params.project_path, "data",
                "cortical_masks",
                params.config.space + ".nii.gz")
            mask = read_file(corticalmask)
        params.data_path = sort_files(params.data_path, params.config.name)
        for curr_data_file in params.data_path:
            curr_data = read_file(curr_data_file)
            # Apply cortical mask
            if "mm" in params.config.space:
                curr_data[mask == 0] = 0
            # Apply threshold
            if params.config.threshold is not None:
                th_l = params.config.threshold[0]
                th_h = params.config.threshold[1]            
                curr_data = np.where(np.logical_and(curr_data > th_l, curr_data < th_h), 1, 0)
            data.append(curr_data)

    return data

def dice_coeff(net1, net2):
    """
    Computes Dice coefficient between two arrays
    """
    intersection = np.sum(net1[net2 == 1])
    coefficient = (2.0 * intersection) / (np.sum(net1) + np.sum(net2))
    return coefficient


def compute_overlap(ref_params,other_params):
    """
    Computes network overlap between reference data or atlas with another atlas.
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
                    overlap_val[idx_i][idx_j] = dice_coeff(curr_ref_net, curr_other_net)
                    idx_j = idx_j + 1
                idx_i = idx_i + 1
        elif other_params.config.type in ('Soft', 'Metric'):
            other_net_uni = range(0,len(other_data))
            idx_i = 0
            overlap_val = np.zeros((len(ref_net_uni), len(other_net_uni)))
            for each_ref_net in ref_net_uni:
                curr_ref_net = ref_data == each_ref_net
                idx_j = 0
                for each_other_net in other_net_uni:
                    curr_other_net = other_data[each_other_net]
                    overlap_val[idx_i][idx_j] = dice_coeff(curr_ref_net, curr_other_net)
                    idx_j = idx_j + 1
                idx_i = idx_i + 1
    elif ref_params.config.type in ('Soft', 'Metric'):
        # Reference atlas/data is a soft parcellation or any metric data
        ref_net_uni = range(0,len(ref_data))
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
                    overlap_val[idx_i][idx_j] = dice_coeff(curr_ref_net, curr_other_net)
                    idx_j = idx_j + 1
                idx_i = idx_i + 1
        elif other_params.config.type in ('Soft', 'Metric'):
            other_net_uni = range(0,len(other_data))
            idx_i = 0
            overlap_val = np.zeros((len(ref_net_uni), len(other_net_uni)))
            for each_ref_net in ref_net_uni:
                curr_ref_net = ref_data[each_ref_net]
                idx_j = 0
                for each_other_net in other_net_uni:
                    curr_other_net = other_data[each_other_net]
                    overlap_val[idx_i][idx_j] = dice_coeff(curr_ref_net, curr_other_net)
                    idx_j = idx_j + 1
                idx_i = idx_i + 1

    return overlap_val

def compute_overlap_data(data_params, atlas_name):
    """
    Computes overlap matrix for a given data and an existing atlas
    """
    ## Reads in the atlas in the same space as the data
    atlas_params = AtlasParams(atlas_name,data_params.config.space)

    ## Compute the overlap matrix between data and the atlas
    overlap_val = compute_overlap(data_params,atlas_params)

    return overlap_val

def load_overlap_atlases(ref_atlas_name, other_atlas_name):
    """
    Load overlap matrix for any pair of existing atlases
    """
    overlap_file = path.join(RESULTS_PATH, ref_atlas_name, other_atlas_name + '.mat')
    data = loadmat(overlap_file)
    overlap_val = np.array(data.get('overlap'))

    return overlap_val