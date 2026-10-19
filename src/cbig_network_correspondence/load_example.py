from . import datasets

_FILES = {
    "example_config": "examples/Exemplar1_pos",
    "example_nii": "examples/Exemplar1.nii.gz",
}


def __getattr__(name):
    if name in _FILES:
        return datasets.fetch_file(_FILES[name])
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
