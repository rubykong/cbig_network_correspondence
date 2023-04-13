from . import compute_overlap_with_atlases
from . import visualize_overlap_lib
from . import load_example
from git import Repo
from os import path

data_dir = path.join(path.dirname(path.abspath(__file__)), 'data')
if not path.exists(data_dir):
    print("Download data")
    Repo.clone_from('https://github.com/rubykong/cbig_network_correspondence_data', data_dir)
else:
    repo_dir = Repo(data_dir)
    repo_dir.remotes.origin.pull('master')
