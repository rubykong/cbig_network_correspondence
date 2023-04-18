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

After importing the toolbox, here are the relevant commands. After decribing the commands, we provide examples for how to use this toolbox. We also provided a jupyter notebook for the examples: ExampleUsage.ipynb. Please note that the output figure path were hard-coded in the notebook, so you might need to change the output figure path to your own path.

```python
import cbig_network_correspondence as cnc
```

```
overlap_mat = cnc.compute_overlap_with_atlases.load_overlap_atlases(ref_atlas_name, other_atlas_name)

This function will load the pre-computed network overlap matrix between the reference atlas and the other atlas.

Input:
- ref_atlas_name: the name of the reference atlas
- other_atlas_name: the name of the other atlas

Output:
- overlap_mat: the network overlap matrix between the reference atlas and the other atlas
```

```
cnc.visualize_overlap_lib.draw_overlap_mat(overlap_mat, ref_atlas_name, other_atlas_name, minv, maxv, outputfigfile)

This function will draw the network overlap matrix between the reference atlas and the other atlas. In this function, the user can pass in the loaded overlap matrix from cnc.compute_overlap_with_atlases.load_overlap_atlases.

Input:
- overlap_mat: the network overlap matrix between the reference atlas and the other atlas
- ref_atlas_name: the name of the reference atlas
- other_atlas_name: the name of the other atlas
- minv: the minimum value of the colorbar
- maxv: the maximum value of the colorbar
- outputfigfile: the output figure file name

```
```
cnc.visualize_overlap_lib.draw_overlap_atlases(ref_atlas_name, other_atlas_name, minv, maxv, outputfigfile)

This function will draw the network overlap matrix between the reference atlas and the other atlas. This function will automatically load the pre-computed network overlap matrix between the reference atlas and the other atlas.

Input:
- ref_atlas_name: the name of the reference atlas
- other_atlas_name: the name of the other atlas
- minv: the minimum value of the colorbar
- maxv: the maximum value of the colorbar
- outputfigfile: the output figure file name
```

```
overlap_data = cnc.compute_overlap_with_atlases.obtain_overlap_atlases(ref_atlas_name, other_atlas_name, ref_network_name, other_network_name)

This function generate the overlapping area between networks from two atlases. This is used for visualization purpose.

Input:
- ref_atlas_name: the name of the reference atlas
- other_atlas_name: the name of the other atlas
- ref_network_name: the name of the reference network
- other_network_name: the name of the other network

Note that you can find the network names for each atlas here:
https://github.com/rubykong/cbig_network_correspondence_data/tree/master/network_names

Output:
- overlap_data: a numpy array in the same space as the reference atlas. The array contains three labels: 1, 2, 3.
1 indicates the reference network, 2 indicates the other network, 3 indicates the overlapping area between the two networks.
```

```
cnc.visualize_overlap_lib.draw_overlap_map(data_space, overlap_data, outputfigfile, coords)

This function visualize the overlapping area between networks from two atlases. The reference network will be colored in pink, the other network will be colored in blue, and the overlapping area will be colored in purple.

Input:
- data_space: the space of the data. It can be `fsaverage6`, `fs_LR_32k`, `FSLMNI2mm`, `LairdColin2mm`, `ShenColin1mm`
- overlap_data: a numpy array in the same space as the reference atlas. The array contains three labels: 1, 2, 3. This is the output of cnc.compute_overlap_with_atlases.obtain_overlap_atlases
- outputfigfile: the output figure file name
- coords: [for volumetric data only] the coordinates of the point where the cut is performed. BY default, we use [-4, -31, 18].
```


```
example = cnc.load_example

This function find the path to the example config and example nifti image

Output:
- example_config: the path to the example config file
- example_nii: the path to the example nifti image
```

```
data_params = cnc.compute_overlap_with_atlases.DataParams(config_file, data_path)

This function reads in the config information and the data path to construct the parameters for the data

Input:
- config_file: the path to the config file
- data_path: the path to the data
```

```
overlap_data = cnc.compute_overlap_with_atlases.obtain_overlap_data(data_params, atlas_name, network_name)

This function generate the overlapping area between the binarized data and an existing atlas. This is used for visualization purpose. The calculation is done in the same space as the data.

Input:
- data_params: the parameters for the data. This is the output of the function cnc.compute_overlap_with_atlases.DataParams
- atlas_name: the name of the atlas. For example, "MG360J12".
- network_name: the name of the network from the atlas. For example, "Default".

```

### Compute and visualize network overlap between two existing alases

In this toolbox, we computes the network overlap between any pair of existing atlases. We treat each atlas as the reference atlas, and projected other atlases to the same space as the reference atlas. Then we compute the network overlap between the reference atlas and the projected atlas.

Here we provide examples if you want to check the network overlap matrix between the Glasser 2016 atlas with Ji2019 networks and the Yeo 17 networks, you can use the following code:

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

The output figures will be save under the path `<path_to_figure>/TY17_MG360J12.png` and `<path_to_figure>/MG360J12_TY17.png`. Here are the corresponding figures. The left side are reference network names, the upper part are the other atlas network names

+ `<path_to_figure>/TY17_MG360J12.png` 

![image](https://user-images.githubusercontent.com/20438248/232752322-0fbe141f-3879-4c6d-9a04-1dd372753ea0.png)

+ `<path_to_figure>/MG360J12_TY17.png`

![image](https://user-images.githubusercontent.com/20438248/232752400-a162a152-043f-48e9-b45d-171bfea6b952.png)


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
overlap_data = cnc.compute_overlap_with_atlases.obtain_overlap_atlases("TY17", "MG360J12", "DefaultA", "Default")
# Yeo 17-network atlas is in fsaverage6 surface space
cnc.visualize_overlap_lib.draw_overlap_map("fsaverage6", overlap_data, "<path_to_figures>/TY17_MG360J12_DefaultA_Default")

# When Glasser2016 atlas atlas is the reference atlas
overlap_data = cnc.compute_overlap_with_atlases.obtain_overlap_atlases("MG360J12", "TY17", "Default", "DefaultA")
# Glasser2016 atlas is in fs_LR_32k surface space
cnc.visualize_overlap_lib.draw_overlap_map("fs_LR_32k", overlap_data, "<path_to_figures>/MG360J12_TY17_Default_DefaultA")

```

Here are the results:
![Overlap_surface](https://user-images.githubusercontent.com/20438248/232758349-ad5157cd-d6c2-4d4c-880c-7dda3a36fe0b.jpg)


### Compute and visualize overlap between any metric data and an existing atlas

In this toolbox, we also allow users to upload their own data and compute the network overlap between the data and an existing atlas. We treat the data as the reference, and projected the atlas to the same space as the reference data. For a metric data, we binarize the data use `Data_Threshold`, then we compute the network overlap between the binarized data and the projected atlas. Here we provide an example data to show you how to use this toolbox.

After importing the toolbox, here are the relevant commands. After decribing the commands, we provide examples for how to use this toolbox to compute and visualize the overlap between any metric data and an existing atlas.

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
overlap_mat = cnc.compute_overlap_with_atlases.compute_overlap_data(data_params, "MG360J12")
```

After that, we can visualize the overlap matrix by using the following code:

```python
cnc.visualize_overlap_lib.draw_overlap_mat(overlap_mat, "Exampler1(pos)", "MG360J12", 0, 1, "<path_to_figures>/exampler1_pos_MG360J12")
```

Here is how the overlap matrix looks like:
![image](https://user-images.githubusercontent.com/20438248/232759475-2eae4dad-0f7c-40d0-afe3-bd4731466542.png)

Based on the visualized overlap matrix, we can notice that the example data overlapped the most with the `Default` network from Glasser2016. If we want to check the overlapping regions between the example data and the `Default` network from Glasser2016, we can use the following code:

```python
overlap_data = cnc.compute_overlap_with_atlases.obtain_overlap_data(data_params, 'MG360J12', 'Default')
# By default we use coords=[-4, -31, 18] to visualize volumetric data, you can change it to any other coordinates
cnc.visualize_overlap_lib.draw_overlap_map("FSLMNI2mm", overlap_data, "<path_to_figures>/example_MG360J12_Default", [-4, -31, 18])
```

Here is the result:
![image](https://user-images.githubusercontent.com/20438248/232760522-f1327b70-2bc7-4b74-bdf2-30274146bfa6.png)



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
