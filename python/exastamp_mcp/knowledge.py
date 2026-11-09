"""Curated exaStamp knowledge shipped with the MCP server."""

import os
import re

import yaml

__all__ = ["SCHEMAS", "DOC_ROOT", "doc_url"]

_HERE = os.path.dirname(os.path.abspath(__file__))

with open(os.path.join(_HERE, "schemas.yaml")) as _f:
    SCHEMAS = yaml.safe_load(_f)

# documentation pages are referenced relative to the published site
DOC_ROOT = os.environ.get("EXASTAMP_DOC_ROOT", "https://collab4exanbody.github.io/doc_exaStamp")


def doc_url(page):
    """Published URL of a documentation page given as its markdown path."""
    page = re.sub(r"(index)?\.md$", "", page)
    return f"{DOC_ROOT}/{page}"
