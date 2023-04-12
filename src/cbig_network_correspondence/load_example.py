from os import path

# define path
PROJECT_PATH = path.dirname(path.abspath(__file__))
EXAMPLE_PATH = path.join(PROJECT_PATH, "data", "examples")
example_config = path.join(EXAMPLE_PATH, "Exemplar1_pos")
example_nii = path.join(EXAMPLE_PATH, "Exemplar1.nii.gz")


def __init__(self):
    """
    This defines the parameter structures for a given data
    """
    self.data_path = example_nii
    self.config_file = example_config