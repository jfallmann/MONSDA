"""Pytest wrapper around tests/tool_audit/audit_tool_flags.py.

Skips (with a visible reason) unless run on a machine that has the
`apptainer` binary on PATH and the expected container directory - this
audit needs the real external-tool containers and cannot run in a bare CI
checkout. See tests/tool_audit/README.md for the manual findings this
automates and how to run it directly:

    conda activate apptainer
    python tests/tool_audit/audit_tool_flags.py
"""

import shutil
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "tool_audit"))

CONTAINER_DIR = Path.home() / "Work2" / "Container" / "MONSDA"

pytestmark = pytest.mark.skipif(
    shutil.which("apptainer") is None or not CONTAINER_DIR.exists(),
    reason=(
        "requires 'apptainer' on PATH (conda activate apptainer) and "
        f"containers under {CONTAINER_DIR}"
    ),
)


def test_tool_flags_audit():
    import audit_tool_flags

    rc = audit_tool_flags.main([])
    assert rc == 0, (
        "tool CLI flag audit found drift/gaps - see printed table above "
        "and tests/tool_audit/README.md"
    )
