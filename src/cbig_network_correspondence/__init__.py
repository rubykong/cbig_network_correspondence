#from . import compute_overlap_with_atlases
#from . import visualize_overlap_lib
#from . import load_example
import asyncio
import os

output_dir = os.path.dirname(os.path.abspath(__file__))
os.system("gh-folder-download --url https://github.com/rubykong/cbig_network_correspondence_data/tree/master/data --output " + output_dir)