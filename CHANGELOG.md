# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Fixed
- improve the readme and usage.ipynb to explain more about the downloading atlases section

## [0.4.0] - 2026-10-19
### Added
- release atlases from GINNA (Gillig2025), Jackson2021, Lipkin2022

### Changed
- Data is no longer fetched by `git clone`/`git pull` as a side effect of `import cbig_network_correspondence`. Instead, individual data files are downloaded lazily (on first actual use) and cached locally via [Pooch](https://www.fatiando.org/pooch/), pinned to a tagged release (`v1.1.0`) of [cbig_network_correspondence_data](https://github.com/rubykong/cbig_network_correspondence_data) and verified by checksum. `spin_rotations` (large, static) is fetched separately from the rest of the data, only when a permutation test is actually run.
- The cache location can be overridden with the `CBIG_NETWORK_CORRESPONDENCE_DATA_DIR` environment variable. A new `cnc.datasets.prefetch_all()` helper can warm the whole cache ahead of time (e.g. on a cluster login node before submitting offline jobs).
- Replaced the `gitpython` dependency with `pooch`.

### Fixed
- `import cbig_network_correspondence` still ran the pre-migration code path (`git clone`/`git pull` of the entire data repo, plus manual spin-rotation part-reassembly) left over in `__init__.py` from before the lazy pooch-fetch migration above, so a fresh import could hang indefinitely instead of being near-instant as documented. Removed the dead code; `datasets.py`'s lazy fetching already covers everything it did.
- Several error paths in `compute_overlap_with_atlases.py` called `sys.exit(1)` without `sys` ever being imported, raising a `NameError` instead of a clear error message. These now raise `ValueError`/`FileNotFoundError` instead.
- Removed dead, unreachable code in `load_example.py` (a top-level `def __init__(self):` that was never attached to any class).

## [0.3.3] - 2025-02-14
### Fixed
- we fixed a bug relating to re-sorting projected files when computing pvalues

## [0.3.2] - 2025-02-12
### Fixed
- We fixed a bug where the input data will be binarized again when computing pvalue for soft or metric data in MNI 2mm

## [0.3.1] - 2024-10-14
### Fixed
- Fixed a minor bug for spin rotation number

## [0.3.0] - 2024-10-10
### Fixed
- We fixed bugs in interpolation method for hard parcellation

### Added
- release atlases from Gordon2024, Woodward2024, Du2024

## [0.2.1] - 2024-06-11
### Added
- Upgrade pakcgaes to latest versions
- release python env file

## [0.2.0] - 2024-04-04
### Added
- Usage example and readme files are available

## [0.1.4] - 2024-03-30
### Added
- Now we can use a single command to generate everything
- Usage example and readme files are available

## [0.1.3] - 2024-03-26
### Fixed
- We fixed bugs in visualization scripts

## [0.1.2] - 2024-03-15
### Added
- we add a new feature to perform permutation test for overlap
- we allow for perform network correspondence with a given input data
- a better circular chart visualization and spider plot and summary table are provided

## [0.1.1] - 2023-05-22
### Added
- we added yan2023 homotopic atlases
- add webpage for documentation

## [0.1.0] - 2023-04-18
### Added
- Examples and description for how to use the toolbox

## [0.0.11] - 2023-04-18
### Fixed
- import error for debugging purpose and add templates

## [0.0.10] - 2023-04-18
### Added
- visualization scripts to visualize overlap data in surface or volume
- documentation in readme

### Fixed
- import error

## [0.0.9] - 2023-04-17
### Added
- we add the scripts to compute network correspondence between atlases
- we provide scripts to output the overlap data

### Fixed
- we fix and release the changelog

[Unreleased]: https://github.com/rubykong/cbig_network_correspondence/compare/0.4.0...master
[0.4.0]: https://github.com/rubykong/cbig_network_correspondence/compare/0.3.3...0.4.0
[0.3.3]: https://github.com/rubykong/cbig_network_correspondence/compare/0.3.2...0.3.3
[0.3.2]: https://github.com/rubykong/cbig_network_correspondence/compare/0.3.1...0.3.2
[0.3.1]: https://github.com/rubykong/cbig_network_correspondence/compare/0.3.0...0.3.1
[0.3.0]: https://github.com/rubykong/cbig_network_correspondence/compare/0.2.1...0.3.0
[0.2.1]: https://github.com/rubykong/cbig_network_correspondence/compare/0.2.0...0.2.1
[0.2.0]: https://github.com/rubykong/cbig_network_correspondence/compare/0.1.4...0.2.0
[0.1.4]: https://github.com/rubykong/cbig_network_correspondence/compare/0.1.3...0.1.4
[0.1.3]: https://github.com/rubykong/cbig_network_correspondence/compare/0.1.2...0.1.3
[0.1.2]: https://github.com/rubykong/cbig_network_correspondence/compare/0.1.1...0.1.2
[0.1.1]: https://github.com/rubykong/cbig_network_correspondence/compare/0.1.0...0.1.1
[0.1.0]: https://github.com/rubykong/cbig_network_correspondence/compare/0.0.11...0.1.0
[0.0.11]: https://github.com/rubykong/cbig_network_correspondence/compare/0.0.10...0.0.11
[0.0.10]: https://github.com/rubykong/cbig_network_correspondence/compare/0.0.9...0.0.10
[0.0.9]: https://github.com/rubykong/cbig_network_correspondence/tree/0.0.9
