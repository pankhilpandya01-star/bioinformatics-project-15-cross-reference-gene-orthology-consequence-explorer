"""Complete-or-absent output publication."""

from __future__ import annotations

import shutil
import tempfile
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator


@contextmanager
def atomic_output_directory(destination: Path) -> Iterator[Path]:
    """Publish a sibling staging directory only after every check succeeds."""

    if destination.exists():
        raise FileExistsError(f"output directory already exists: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(
        tempfile.mkdtemp(prefix=f".{destination.name}.staging-", dir=destination.parent)
    )
    try:
        yield staging
        # A same-parent rename prevents readers from observing partial outputs.
        staging.rename(destination)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise

