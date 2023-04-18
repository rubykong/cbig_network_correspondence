# cbig_network_correspondence

[![PyPI](https://img.shields.io/pypi/v/cbig_network_correspondence?style=flat-square)](https://pypi.python.org/pypi/cbig_network_correspondence/)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/cbig_network_correspondence?style=flat-square)](https://pypi.python.org/pypi/cbig_network_correspondence/)
[![PyPI - License](https://img.shields.io/pypi/l/cbig_network_correspondence?style=flat-square)](https://pypi.python.org/pypi/cbig_network_correspondence/)
[![Coookiecutter - Wolt](https://img.shields.io/badge/cookiecutter-Wolt-00c2e8?style=flat-square&logo=cookiecutter&logoColor=D4AA00&link=https://github.com/woltapp/wolt-python-package-cookiecutter)](https://github.com/woltapp/wolt-python-package-cookiecutter)


---

**Source Code**: [https://github.com/rubykong/cbig_network_correspondence](https://github.com/rubykong/cbig_network_correspondence)

**PyPI**: [https://pypi.org/project/cbig_network_correspondence/](https://pypi.org/project/cbig_network_correspondence/)

---

This toolbox was used to explore the network correspondence between networks across different atlases. 

## Installation

```sh
pip install cbig_network_correspondence
```

## Usage

This toolbox has the following functions:

+ Compute and visualize network overlap between two existing alases
+ Compute and visualize overlap between any metric data (e.g. task contrast map, probability map, other atlas) and an existing atlas
+ Visualze the overlapping area between two atlases or a metric data and an atlas

We will provide detailed usage of these function in the following sections.

### Atlases we included

Check this list for the atlases we have included. We use abbreviations for the atlases. Here we provide the decription for each atlas. If you want us to include your atlas, please contact me (roo.cone@gmail.com).

![networkcorrespondence](https://user-images.githubusercontent.com/20438248/232702654-3f334164-bb13-41af-a7aa-8f29bae9fca7.jpg)

These atlases should be automatically downloaded when you install the toolbox. If you want to download the atlases manually, you can find the atlases here: https://github.com/rubykong/cbig_network_correspondence_data. The atlas data will be automatically downloaded in your python package directory: `cbig_network_correspondence/src/cbig_network_correspondence/data`.

### Config file

In this toolbox, the information for each atlas was constructed as config file. We have included the config files for each atlas. The config files will be automatically downloaded when you install the toolbox. If you want to check the config files manually, you can find the config files in your python package directory: `cbig_network_correspondence/src/cbig_network_correspondence/data/atlas_config`. If you want to use your own atlas, you can create your own config file, check next sections for details. 

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

Here is an example of the config file for a example metric data `Exemplar1(pos)`:
```
[data_info]
Data_Name: Exemplar1(pos)
Data_Space: FSLMNI2mm
Data_Type: Metric
Data_Threshold: [5,Inf]
```

Here is the decription for the config file:

+ `Data_Category` [optional]
  
  We included it for existing atlases because some atlases were from the same lab. In this case, we put the some atlases into the same category. This is not necessary for your own atlas or your own metric data.
   
+ `Data_Name`

  The name of the atlas or metric data.

+ `Data_Space`

  The space name of the space where the atlas/data is in. It can be a surface name or a volumetric space name. We currenly have atlases in the following spaces: `fs_LR_32k` surface space, `fsaverage6` surface space, `FSLMNI2mm` volumetric space, `LairdColin2mm` volumetric space, `ShenColin1mm` volumetric space.
  
+ `Data_Type`

  The type of the atlas/data. We have two types of atlas: `Hard` and `Soft`. The `Hard` atlas is the atlas that each ROI is assigned to a specific network. The `Soft` atlas is the atlas that each ROI is assigned to a probability of belonging to each network. We also allow for a metric data such as task contrast maps or probability maps, in this case `Data_Type` should be `Metric`. 

+ `Data_NetworkAssignment` [optional]

  Sometimes the atlas has fine-grained ROIs and these ROIs were assigned to large-scale networks. For example, the Schaefer2018 400 ROIs were assigned to the Yeo 17-network atlas. In this case the user can provide a mapping between ROIs to networks, here you can find the mapping for some of the existing atlases: https://github.com/rubykong/cbig_network_correspondence_data/tree/master/network_assignment

+ `Data_Threshold` [optional]

  For soft parcellation and metric data, we will need to apply a threshold to binarize the data. The threshold should be set as `[lower bound, upper bound]`. If the upper bound is `Inf`, then the threshold will be applied as `>= lower bound`. If the lower bound is `-Inf`, then the threshold will be applied as `<= upper bound`. If the `Data_Type` is `Soft` or `Metric`, but the `Data_Threshold` is not provided, the toolbox will use a default threshold `[0, Inf]`. `Data_Threshold` is not necessary for `Hard` atlas.

### Compute and visualize network overlap between two existing alases

In this toolbox, we computes the network overlap between any pair of existing atlases. We treat each atlas as the reference atlas, and projected other atlases to the same space as the reference atlas. Then we compute the network overlap between the reference atlas and the projected atlas.

Here I will show you how to use this toolbox to compute and visualize the network overlap between two existing atlases.

If you want to check the network overlap matrix between the Glasser 2016 atlas with Ji2019 networks and the Yeo 17 networks, you can use the following code:

```python
import  cbig_network_correspondence as cnc

# When Yeo 17-network atlas is the reference atlas, the Glasser2016 atlas will be projected to fsaverage6 surface space
overlap_mat1 = cnc.compute_overlap_with_atlases.load_overlap_atlases("TY17", "MG360J12")

# When Glasser2016 atlas atlas is the reference atlas, theYeo 17-network atlas will be projected to fs_LR_32k surface space
overlap_mat2 = cnc.compute_overlap_with_atlases.load_overlap_atlases("MG360J12", "TY17")
```

If you want to visualize the network overlap matrix, you can use the following code:

```python
cnc.visualize_overlap_lib.draw_overlap_mat(overlap_mat1, "TY17", "MG360J12", 0, 1, "<path_to_figure>/TY17_MG360J12")
cnc.visualize_overlap_lib.draw_overlap_mat(overlap_mat2, "MG360J12", "TY17", 0, 1, "<path_to_figure>/MG360J12_TY17")
```

However, if you have no interest to check the overlap matrix and just want to see the visualized network overlap between two atlases, you can use the following code:

```python
import  cbig_network_correspondence as cnc

cnc.visualize_overlap_lib.draw_overlap_atlases("MG360J12", "TY17", 0, 1, "<path_to_figure>/MG360J12_TY17")
```

If you want to check the brain maps of the overlappping regions between two atlases, you can use the following code:

+ To check the overlapping regions between the `Default` network from `MG360J12` atlas and the `DefaultA` network from `TY17` atlas:

```python
import  cbig_network_correspondence as cnc

# When Yeo 17-network atlas is the reference atlas
overlap_data = cnc.compute_overlap_with_atlases.obtain_overlap_atlases("TY17", "MG360J12", "Default", "DefaultA")
cnc.visualize_overlap_lib.draw_overlap_map("fsaverage6", overlap_data, "<path_to_figures>/TY17_MG360J12_DefaultA_Default")

# When Glasser2016 atlas atlas is the reference atlas
overlap_data = cnc.compute_overlap_with_atlases.obtain_overlap_atlases("MG360J12", "TY17", "Default", "DefaultA")
# By default we use coords=[-4, -31, 18] to visualize volumetric data, you can change it to any other coordinates
cnc.visualize_overlap_lib.draw_overlap_map("fs_LR_32k", overlap_data, "<path_to_figures>/MG360J12_TY17_Default_DefaultA", [-4, -31, 18])

```

### Compute and visualize overlap between any metric data and an existing atlas

In this toolbox, we also allow users to upload their own data and compute the network overlap between the data and an existing atlas. We treat the data as the reference, and projected the atlas to the same space as the reference data. For a metric data, we binarize the data use `Data_Threshold`,then we compute the network overlap between the binarized data and the projected atlas. Here we provide an example data to show you how to use this toolbox.

You can load the example data by using the following code:

```python
import  cbig_network_correspondence as cnc

example = cnc.load_example
# the example config file
print(example.example_config)
# the path to the example data
print(example.example_nii)
# construct parameters for example data
data_params = cnc.compute_overlap_with_atlases.DataParams(example.example_config, example.example_nii)
```

In this example, we will compute the overlap between thie example data and the Glasser2016 atlas. Since the example data is in FSLMNI2mm space, we will project the Glasser2016 atlas to FSLMNI2mm space and compute the overlap matrix.

```python
overlap_mat = cnc.compute_overlap_with_atlases.compute_overlap_data(data_params,"MG360J12")
```

After that, we can visualize the overlap matrix by using the following code:

```python
cnc.visualize_overlap_lib.draw_overlap_mat(overlap_mat, "Exampler1(pos)", "MG360J12", 0, 1, "<path_to_figures>/exampler1_pos_MG360J12")
```

Based on the visualized overlap matrix, we can notice that the example data overlapped the most with the `Default` network from Glasser2016. If we want to check the overlapping regions between the example data and the `Default` network from Glasser2016, we can use the following code:

```python
overlap_data = cnc.compute_overlap_with_atlases.obtain_overlap_data(data_params, 'MG360J12', 'Default')
cnc.visualize_overlap_lib.draw_overlap_map("FSLMNI2mm", overlap_data, "<path_to_figures>/example_MG360J12_Default")
```



## Development

* Clone this repository
* Requirements:
  * [Poetry](https://python-poetry.org/)
  * Python 3.7+
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
