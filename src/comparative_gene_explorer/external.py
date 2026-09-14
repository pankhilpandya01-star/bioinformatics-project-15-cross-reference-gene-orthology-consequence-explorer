"""Native-tool execution, version gates, and sanitized logging."""

from __future__ import annotations

import json
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence


FAILURE_MARKERS = ("RuntimeException", "FATAL ERROR", "Exception in thread")
ANSI_ESCAPE = re.compile(r"\x1b\[[0-9;]*m")
DATABASE_CHECK = re.compile(
    r"(?P<name>CDS|Protein) check:\s+\S+\s+OK:\s+(?P<ok>\d+)"
    r"(?:\s+Warnings:\s+(?P<warnings>\d+))?"
    r"\s+Not found:\s+(?P<not_found>\d+)\s+Errors:\s+(?P<errors>\d+)"
    r"\s+Error percentage:\s+(?P<percentage>[\d.]+)%"
)


class ExternalToolError(RuntimeError):
    """Raised when a required native command cannot produce valid output."""


@dataclass(frozen=True, slots=True)
class CommandResult:
    command: tuple[str, ...]
    stdout: str
    stderr: str


def sanitize(text: str, replacements: Iterable[tuple[str, str]]) -> str:
    sanitized = text
    for source, replacement in sorted(replacements, key=lambda item: len(item[0]), reverse=True):
        if source:
            sanitized = sanitized.replace(source, replacement)
    return sanitized


def run_command(
    command: Sequence[str],
    *,
    cwd: Path,
    stdout_path: Path | None = None,
    stderr_path: Path | None = None,
    raw_stdout: bool = False,
    replacements: Iterable[tuple[str, str]] = (),
    reject_markers: Sequence[str] = FAILURE_MARKERS,
) -> CommandResult:
    try:
        completed = subprocess.run(
            list(command), cwd=cwd, capture_output=True, text=True, check=False
        )
    except OSError as error:
        raise ExternalToolError(f"could not start {command[0]}: {error}") from error
    combined = completed.stdout + "\n" + completed.stderr
    marker = next((value for value in reject_markers if value in combined), None)
    if completed.returncode or marker is not None:
        detail = marker or f"exit status {completed.returncode}"
        tail = sanitize(combined[-3000:], replacements).strip()
        raise ExternalToolError(f"{command[0]} failed ({detail}): {tail}")
    if stdout_path is not None:
        text = completed.stdout if raw_stdout else sanitize(completed.stdout, replacements)
        stdout_path.write_text(text, encoding="utf-8", newline="\n")
    if stderr_path is not None:
        stderr_path.write_text(
            sanitize(completed.stderr, replacements), encoding="utf-8", newline="\n"
        )
    return CommandResult(tuple(command), completed.stdout, completed.stderr)


def _conda_package_version(package: str) -> str:
    versions: list[str] = []
    for path in (Path(sys.prefix) / "conda-meta").glob(f"{package}-*.json"):
        metadata = json.loads(path.read_text(encoding="utf-8"))
        if metadata.get("name") == package:
            versions.append(str(metadata["version"]))
    if len(versions) != 1:
        raise ExternalToolError(f"expected one installed {package} package; found {versions}")
    return versions[0]


def verify_tool_versions(cwd: Path) -> list[dict[str, object]]:
    definitions = (
        ("OrthoFinder", "orthofinder", "3.1.5", ["orthofinder", "--version"], r"OrthoFinder:v([\w.]+)"),
        ("SnpEff", "snpeff", "5.4.0c", ["snpEff", "-version"], r"SnpEff\s+([\w.]+)"),
        ("BCFtools", "bcftools", "1.24", ["bcftools", "--version"], r"bcftools\s+([\w.]+)"),
        ("NCBI Datasets CLI", "ncbi-datasets-cli", "18.36.0", ["datasets", "version"], r"datasets version:\s+([\w.]+)"),
        ("OpenJDK", "openjdk", "21", ["java", "-version"], r'openjdk version "([^"]+)"'),
    )
    rows: list[dict[str, object]] = []
    for label, package, required, command, pattern in definitions:
        result = run_command(command, cwd=cwd, reject_markers=())
        text = ANSI_ESCAPE.sub("", result.stdout + result.stderr)
        match = re.search(pattern, text)
        if match is None:
            raise ExternalToolError(f"could not parse {label} version")
        installed = _conda_package_version(package)
        passed = installed == required if label != "OpenJDK" else installed.startswith("21.")
        if not passed:
            raise ExternalToolError(f"{label} package {installed} does not satisfy {required}")
        rows.append(
            {
                "tool": label,
                "required": required,
                "installed_package": installed,
                "cli_reported": match.group(1),
                "passed": True,
            }
        )

    import matplotlib

    python_version = f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"
    for label, package, required, reported, passed in (
        ("Python", "python", "3.12", python_version, python_version.startswith("3.12.")),
        ("Matplotlib", "matplotlib", "3.11.1", matplotlib.__version__, matplotlib.__version__ == "3.11.1"),
    ):
        installed = _conda_package_version(package)
        if not passed:
            raise ExternalToolError(f"{label} version {reported} does not satisfy {required}")
        rows.append(
            {
                "tool": label,
                "required": required,
                "installed_package": installed,
                "cli_reported": reported,
                "passed": True,
            }
        )
    return rows


def parse_database_check(text: str, expected_name: str) -> dict[str, object]:
    matches = [
        match
        for match in DATABASE_CHECK.finditer(text)
        if match.group("name").lower() == expected_name
    ]
    if len(matches) != 1:
        raise ExternalToolError(f"expected exactly one {expected_name} validation summary")
    match = matches[0]
    return {
        "check": expected_name,
        "ok": int(match.group("ok")),
        "warnings": int(match.group("warnings") or 0),
        "not_found": int(match.group("not_found")),
        "errors": int(match.group("errors")),
        "error_percentage": float(match.group("percentage")),
    }

