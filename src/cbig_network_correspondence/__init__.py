from . import compute_overlap_with_atlases
from . import visualize_overlap_lib
from . import load_example
import os
os.environ["GIT_PYTHON_REFRESH"] = "quiet"
from git import Repo

data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
if not path.exists(data_dir):
    print("Download data")
    Repo.clone_from('https://github.com/rubykong/cbig_network_correspondence_data', data_dir)
else:
    repo_dir = Repo(data_dir)
    repo_dir.remotes.origin.pull('master')
