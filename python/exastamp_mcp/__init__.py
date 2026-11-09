"""
MCP server for exaStamp: helps an LLM agent write input files that run the
first time.

Pure Python: the server process never loads exaStamp. Anything that needs
exaStamp (operator catalog, build and smoke checks) runs in a child process
through ``onika_mcp.worker``, with ``pyexastamp`` as the application package.
"""

PACKAGE = "pyexastamp"
