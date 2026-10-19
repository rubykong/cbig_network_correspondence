"""Regression guard for `import cbig_network_correspondence` itself.

A prior bug had __init__.py doing a full `git clone` of the whole data repo as
an import-time side effect (left over from before the lazy pooch-fetch
migration) -- a fresh import could hang indefinitely instead of completing
near-instantly. Run in a subprocess so this measures a genuinely cold import,
not a cached re-import from an earlier test in the same session.
"""

import subprocess
import sys


def test_import_is_fast_and_network_free():
    result = subprocess.run(
        [
            sys.executable,
            "-c",
            "import time; t = time.time(); import cbig_network_correspondence; "
            "print(time.time() - t)",
        ],
        capture_output=True,
        text=True,
        timeout=30,
        check=True,
    )
    elapsed = float(result.stdout.strip())
    assert (
        elapsed < 10
    ), f"import took {elapsed:.2f}s -- should do no network I/O at import time"


def test_exposes_expected_submodules():
    import cbig_network_correspondence as cnc

    for name in (
        "datasets",
        "compute_overlap_with_atlases",
        "visualize_overlap_lib",
        "visualize_report_lib",
        "load_example",
    ):
        assert hasattr(cnc, name), f"cbig_network_correspondence.{name} missing"
