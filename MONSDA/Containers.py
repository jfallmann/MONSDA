import re

from MONSDA import _version


def normalize_container_version(version=None):
    version = version or _version.get_versions()["version"]
    version = re.sub(r"^v(?=\d)", "", version)
    version = re.sub(r"(?:[.+-]dirty)$", "", version)
    version = re.sub(r"[^A-Za-z0-9_.-]+", "-", version).strip(".-")
    if not version:
        raise ValueError("MONSDA version cannot be converted to an OCI tag")
    return version
