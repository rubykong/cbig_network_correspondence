# cbig_network_correspondence

[![PyPI](https://img.shields.io/pypi/v/cbig_network_correspondence?style=flat-square)](https://pypi.python.org/pypi/cbig_network_correspondence/)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/cbig_network_correspondence?style=flat-square)](https://pypi.python.org/pypi/cbig_network_correspondence/)
[![PyPI - License](https://img.shields.io/pypi/l/cbig_network_correspondence?style=flat-square)](https://pypi.python.org/pypi/cbig_network_correspondence/)
[![Coookiecutter - Wolt](https://img.shields.io/badge/cookiecutter-Wolt-00c2e8?style=flat-square&logo=cookiecutter&logoColor=D4AA00&link=https://github.com/woltapp/wolt-python-package-cookiecutter)](https://github.com/woltapp/wolt-python-package-cookiecutter)


---

**Documentation**: [https://rubykong.github.io/cbig_network_correspondence](https://rubykong.github.io/cbig_network_correspondence)

**Source Code**: [https://github.com/rubykong/cbig_network_correspondence](https://github.com/rubykong/cbig_network_correspondence)

**PyPI**: [https://pypi.org/project/cbig_network_correspondence/](https://pypi.org/project/cbig_network_correspondence/)

---

This toolbox was used to explore the network correspondence between networks across different atlases. 

## Installation

## Python installation

Check here for how to install python: https://realpython.com/installing-python/

### Conda environment

We provide NCT python environment here: [NCT_env](https://github.com/rubykong/cbig_network_correspondence/blob/master/NCT_env.yml). You can create a conda environment using the following command:

```
conda env create -f NCT_env.yml
```

To use this environment, activate the environment using the following command:

```
conda activate NCT_env
```

## NCT installation
To install:

```sh
pip install cbig_network_correspondence
```

To upgrade:

```sh
pip install --upgrade cbig_network_correspondence
```

If you don't have pip installed, this [pip installation](https://pip.pypa.io/en/stable/installation/) can guide you through the process.

## Usage


### Quick Example!

Test the toolbox quickly using the example data! Check the last section in [Usage.ipynb](https://github.com/rubykong/cbig_network_correspondence/blob/master/Usage.ipynb) for the results.

```
import cbig_network_correspondence as cnc

atlas_names_list = ["EG17","MG360J12"]
example = cnc.load_example
# the example config file
print(example.example_config)
# the path to the example data
print(example.example_nii)
ref_params = cnc.compute_overlap_with_atlases.DataParams(example.example_config, example.example_nii)
cnc.compute_overlap_with_atlases.network_correspondence(ref_params, atlas_names_list,"~/NCTexample/example_results")

```


### Tutorial

We provide a tutorial for how to use this toolbox. The tutorial is available in the Usage.ipynb notebook:
+ [Usage.ipynb](https://github.com/rubykong/cbig_network_correspondence/blob/master/Usage.ipynb)

### Atlases we included

Check this list for the atlases we have included. We use abbreviations for the atlases. Here we provide the decription for each atlas. If you want us to include your atlas, please contact me (roo.cone@gmail.com).

| Abbreviation | Source                         | Description                                            |
|--------------|--------------------------------|--------------------------------------------------------|
| MG360J12     | Glasser2016; Ji2019            | Mattew Glasser2016 360-ROI with Ji2019 12 Cole-Anticevic networks |
| HCPICA       | Smith2009; Smith2013           | HCP 25-node ICA maps                                   |
| AL20         | Laird2011                      | Angela Laird2011 20-node ICA maps                      |
| AS200K17     | Schaefer2018; Kong2021         | Alex Schaefer2018 200-ROI with Kong2021 17 networks    |
| AS200Y17     | Schaefer2018; Yeo2011          | Alex Schaefer2018 200-ROI with Yeo2011 17 networks     |
| AS400K17     | Schaefer2018; Kong2021         | Alex Schaefer2018 400-ROI with Kong2021 17 networks    |
| AS400Y17     | Schaefer2018; Yeo2011          | Alex Schaefer2018 400-ROI with Yeo2011 17 networks     |
| XS268_8      | Shen2013                       | Xilin Shen2013 268-ROI with 8 networks                 |
| XS368_8      | Shen2013                       | Xilin Shen2013 368-ROI with 8 networks                 |
| WS90_14      | Shirer2012                     | William Shirer2012 90-ROI 14 networks                  |
| UKBICA       | Smith2009; Alfaro-Almagro2018  | UKBiobank 25-node ICA maps                             |
| EG286_12     | Power2011; Gordon2016          | Evan Gordon2016 286-ROI with 12 networks               |
| TL12         | Power2011; Laumann2015         | Tim Laumann2015 12 networks (Power2011)                |
| EG17         | Power2011; Gordon2017          | Evan Gordon2017 17 networks                            |
| TY7          | Yeo2011                        | Thomas Yeo 7 networks                                  |
| TY17         | Yeo2011                        | Thomas Yeo 17 networks                                 |
| XY200K17     | Yan2023; Kong2021              | Xiaoxuan Yan2023 200-ROI with Kong2021 17 networks     |
| XY200Y17     | Yan2023; Yeo2011               | Xiaoxuan Yan2023 200-ROI with Yeo2011 17 networks      |
| XY400K17     | Yan2023; Kong2021              | Xiaoxuan Yan2023 400-ROI with Kong2021 17 networks     |
| XY400Y17     | Yan2023; Yeo2011               | Xiaoxuan Yan2023 400-ROI with Yeo2011 17 networks      |


These atlases should be automatically downloaded when you install the toolbox. If you want to download the atlases manually, you can find the atlases here: https://github.com/rubykong/cbig_network_correspondence_data. The atlas data will be automatically downloaded in your python package directory: `cbig_network_correspondence/src/cbig_network_correspondence/data`.


### Input data

The toolbox can also explore the network correspondence of your own input data. The input data can be different format. We support input data in fs_LR_32k, fsaverage6, and FSLMNI2mm.

#### Format

Here are the formats we support:

For data in fs_LR_32k or fsaverage6 space:

+ `.mat` file: The file should contain two variables starting with `lh` and `rh`. For example, the variable names could be `lh_data`, `rh_data` or `lh_labels`, `rh_labels`. Each variable should be a 32492x1 vector for `fs_LR_32k` space or a 40962x1 vector for `fsaverage6` space. 

+ `.npy` file: The file should contain a numpy array with shape (64984, 1) for `fs_LR_32k` space or (81924, 1) for `fsaverage6` space.

For data in FSLMNI2mm space:

+ `.nii` or `.nii.gz` file: The file should be a nifti file with the same dimension as the FSL MNI152 2mm template.


#### Pass in the data as ...

+ A single path to a single file
  - `/data_dir/my_task_contrast_map.mat`
  - `/data_dir/my_task_contrast_map.npy`
  - `/data_dir/my_task_contrast_map.nii.gz`

+ A list of paths to multiple files
  - `['/data_dir/my_task_contrast_map1.mat', '/data_dir/my_task_contrast_map2.mat']`
  - `['/data_dir/my_task_contrast_map1.npy', '/data_dir/my_task_contrast_map2.npy']`
  - `['/data_dir/my_task_contrast_map1.nii.gz', '/data_dir/my_task_contrast_map2.nii.gz']`

+ A directory contains the input files
  - `/data_dir`
  under this directory, the toolbox will search for the files containing the `Data_Name` in the config file (see the Config file section for more details). For example, if the `Data_Name` is `my_task_contrast_map`, the toolbox will search for the files containing `my_task_contrast_map` in the file name

+ A string indicating the name of a existing atlas
  - `AS400Y17`
  - `EG17`
  - `TY17`

#### Reference data space

To explore the network correspondence between a given input data and existing atlases, the toolbox will use the input data as the reference data. The toolbox will use existing atlases in the same space as the input data. For example, if the input data is in `fs_LR_32k` space, the toolbox will use the existing atlases in `fs_LR_32k` space. If the input data is in `FSLMNI2mm` space, the toolbox will use the existing atlases in `FSLMNI2mm` space. 

#### Reference data names

User can specify the name for the given input data, if not, the toolbox will use the `Data_Name` in the config file as the reference data name. 

Reference data names could be important for a multi-dimension input data. For example, if the input data is a set of task contrast maps, the user might want to specify the reference data names for each task contrast map. Otherwise, the toolbox will name each task contrast map as `<Data_Name>1`, `<Data_Name>2`, etc. 

If the input data is the user's own atlas with K networks, the user should specify the reference data names for each network. Otherwise the toolbox will name each network as `<Data_Name>1`, `<Data_Name>2`, etc.

#### Config file

In this toolbox, the information for each atlas or input data is constructed as a config file. We have included the config files for existing atlases. The config files will be automatically downloaded when you install the toolbox. If you want to check the config files manually, you can find the config files in your python package directory: `cbig_network_correspondence/src/cbig_network_correspondence/data/atlas_config`. If you want to explore network correspondence of your own input data/atlas, you can create your own config file, check next sections for how to create a config file. 

The config file should be a text file. 

Here is an example of the config file for the Schaefer 400-ROI atlas with Yeo 17-network network assignemnt `AS400Y17`:
```
[data_info]
Data_Category: YeoLab
Data_Name: AS400Y17
Data_Space: fsaverage6
Data_Type: Hard
Data_NetworkAssignment: 1
```

Here is an example of the config file for a example metric data `Example`:
```
[data_info]
Data_Name: Example
Data_Space: FSLMNI2mm
Data_Type: Metric
Data_Threshold: [5,Inf]
```

Here is the decription for the config file:

+ `Data_Category` [optional]
  
  We included it for existing atlases because some atlases were from the same lab. In this case, we put the some atlases into the same category. This is not necessary for your own atlas or your own metric data.
   
+ `Data_Name`

  The name of the atlas or metric data. If the input data is a set of files instead of a single file, user can pass in the directory to the files. The toolbox will search for the files containing the `Data_Name` in the file name. For example, if the `Data_Name` is `my_task_contrast_map`, the toolbox will search for the files containing `my_task_contrast_map` in the file name.

+ `Data_Space`

  The space name of the space where the atlas/data is in. It can be a surface name or a volumetric space name. We currenly have atlases in the following spaces: `fs_LR_32k` surface space, `fsaverage6` surface space, `FSLMNI2mm` volumetric space, `LairdColin2mm` volumetric space, `ShenColin1mm` volumetric space. If the user wants to explore the network correspondence of their own data, the data should be in `fs_LR_32k`, `fsaverage6`, or `FSLMNI2mm` space.
  
+ `Data_Type`

  The type of the atlas/data. We have two types of atlas: `Hard` and `Soft`. The `Hard` atlas is the atlas that each ROI is assigned to a specific network. The `Soft` atlas is the atlas that each ROI is assigned to a probability of belonging to each network. We also allow for a metric data such as task contrast maps or probability maps, in this case `Data_Type` should be `Metric`. 

+ `Data_NetworkAssignment` [optional]

  Sometimes the atlas has fine-grained ROIs and these ROIs were assigned to large-scale networks. For example, the Schaefer2018 400 ROIs were assigned to the Yeo 17-network atlas. In this case the user can provide a mapping between ROIs to networks, here you can find the mapping for some of the existing atlases: https://github.com/rubykong/cbig_network_correspondence_data/tree/master/network_assignment

+ `Data_Threshold` [optional]

  For soft parcellation and metric data, we will need to apply a threshold to binarize the data. The threshold should be set as `[lower bound, upper bound]`. If the upper bound is `Inf`, then the threshold will be applied as `>= lower bound`. If the lower bound is `-Inf`, then the threshold will be applied as `<= upper bound`. If the `Data_Type` is `Soft` or `Metric`, but the `Data_Threshold` is not provided, the toolbox will use a default threshold `[0, Inf]`. `Data_Threshold` is not necessary for `Hard` atlas.

After importing the toolbox, here are the relevant commands. After decribing the commands, we provide examples for how to use this toolbox. We also provided a jupyter notebook for the usage examples: Usage.ipynb. 


### Usage scenarios

This toolbox has two main usage scenarios:

1. Use this toolbox to explore network correspondence between existing atlases.
2. Provide users' own input data, and explore the network correspondence between the input data and existing atlases.

**We have provided detailed examples for each usage scenario in the Usage.ipynb notebook. Here we provide brief descriptions for inputs.**

#### 1. Explore network correspondence between existing atlases

![NCT_workflow1](https://github.com/rubykong/cbig_network_correspondence/assets/20438248/5046486c-b192-4aa4-90e4-825d962fc94c)

```
import cbig_network_correspondence as cnc

cnc.compute_overlap_with_atlases.network_correspondence(reference_atlas, atlas_names_list, output_dir)
```

+ `reference_atlas`: str
  - The name of the reference atlas. The reference atlas should be one of the existing atlases. The toolbox will use the reference atlas to explore the network correspondence between the reference atlas and other atlases in the `atlas_names_list`.
+ `atlas_names_list`: list
  - A list of the names of the atlases. The toolbox will explore the network correspondence between the reference atlas and the atlases in the `atlas_names_list`.
+ `output_dir`: str
  - The directory to save the output figures and csv files.

Note: network name should be the abbreviation of the atlas. For example, `AS400Y17` for the Schaefer 400-ROI atlas with Yeo 17-network network assignemnt.

#### 2. Explore network correspondence between input data and existing atlases

![NCT_workflow2](https://github.com/rubykong/cbig_network_correspondence/assets/20438248/55f9f944-246d-4146-b530-978b7e951ddd)

```
import cbig_network_correspondence as cnc

# construct DataParams object based on the data file path and config
ref_params = cnc.compute_overlap_with_atlases.DataParams(config, file_path)

# compute the overlap with atlases and save the results
cnc.compute_overlap_with_atlases.network_correspondence(ref_params, atlas_names_list,output_dir,ref_data_names)

```

+ `config`: str
  - The path to the config file for the input data.
+ `file_path`: str/list
  - The path to a single input data file
  - A list of paths to multiple input data files
  - A directory contains the input data files
+ `atlas_names_list`: list
  - A list of the names of the atlases. The toolbox will explore the network correspondence between the input data and the atlases in the `atlas_names_list`.
+ `output_dir`: str
  - The directory to save the output figures and csv files.
+ `ref_data_names`: str/list
  - The name of the input data. If `ref_data_names` is not provided, the toolbox will use the `Data_Name` in the config file as the reference data name. If the input data is multi-dimension, each dimension will be named as `<Data_Name>1`, `<Data_Name>2`,... if the `ref_data_names` is not provided.


### Significance test

The toolbox provides a significance test to test the significance of the network correspondence between the given input data and an existing atlas. The toolbox will use the spin-test (Gordon2017;Alexander-Bloch2018)to test the significance of the network correspondence. Specically, SpinPermutations from BrainSpace python package is used here. For data in FSLMNI2mm space, we project the data to the fsaverage6 space for the spin-test. 


### Results

NCT provides different types of plots to help the user explore/visualize the network correspondences. 

+ Overlap Heatmap
  - If input data has multiple dimensions, e.g., a atlas with K networks, a set of task contrast maps. For a specified list of existing atlases, the toolbox will provide heatmaps to show the overlap between the input data and the existing atlases. Each atlas will have a heatmap. If input data has K dimension, and the other atlas has N networks, the heatmap will be a KxN heatmap. Each row represents the overlap between the input data and an existing atlas. Each column represents the overlap between the input data and a network in the existing atlas. The color represents the overlap value. The asterisk indicates the significance of the overlap. 

+ Network Clock
  - If input data is single dimension, e.g., a task contrast map, a probability map. For a specified list of existing atlases, the toolbox will provide a network clock to show the overlap between the input data and the existing atlases. The network clock will show the overlap between the input data and all networks. Networks from the same atlas are grouped together and colored by the same color. Only the networks with significant overlap will be indicated by network name in the network clock. The clock ring represents the overlap value. The network font size represents the overlap value. 

+ Network Radar
  - If input data is single dimension, e.g., a task contrast map, a probability map. For a specified list of existing atlases, the toolbox will provide a network radar to show the overlap between the input data and the existing atlases. Each atlas will have a radar map. The asterisk indicates the significance of the overlap.

+ Summary Table
  - The toolbox will provide a summary table to show the exact overlap values and p-values of the overlap between the input data and the existing atlases. Summary tables will be saved as figures and csv files.

## Development

* Clone this repository
* Requirements:
  * [Poetry](https://python-poetry.org/)
  * Python 3.9+
* Create a virtual environment and install the dependencies

```sh
poetry install
```

* Activate the virtual environment

```sh
poetry shell
```

---

This project was generated using the [wolt-python-package-cookiecutter](https://github.com/woltapp/wolt-python-package-cookiecutter) template.
