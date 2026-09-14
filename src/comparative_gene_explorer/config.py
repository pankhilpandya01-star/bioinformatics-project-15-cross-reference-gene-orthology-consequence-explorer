"""Configuration for the comparative gene workflow."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True, slots=True)
class WorkflowConfig:
    reference_manifest: Path
    output_dir: Path
    baseline_accession: str = "GCF_000005845.2"
    sample_id: str = "SRR13921545"
    minimum_gene_callable_fraction: float = 0.90
    maximum_database_error_rate: float = 0.02
    upstream_downstream_length: int = 0
    threads: int = 1

