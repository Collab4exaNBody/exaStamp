"""
Evaluation of the exaStamp MCP server on the requests of cases.yaml.

    python run_eval.py --reference         # no LLM: every reference spec must pass
    python run_eval.py --claude-code       # Claude Code (claude -p) drives the server: subscription
    python run_eval.py --claude-code --only ta_bcc_eam_nvt,ar_lj_nve --out results.json
    python run_eval.py --api-model claude-opus-5-5   # Anthropic API: billed per token

Each case runs in its own workspace (default: a temporary directory). The agent
must save its input as input.msp with write_input; ``grade(case, workspace)``
then validates it (smoke run) and checks the expectations. Another agent
framework (any model) can reuse ``grade`` on the workspaces it produced.

Cases outside build_input give a hand-written ``reference_input`` (in
references/) instead of a SimulationSpec, and may list ``workspace_files``
(fixtures/ or the installed data) copied into the workspace before the agent
starts. A ``refusal`` case passes when the agent writes no input and its reply
names the missing feature.

Reported per case: pass, first try (one build/validate call before the write),
number of tool calls, failed expectations.

--claude-code runs the Claude Code CLI non-interactively, signed in with your
Claude account: it uses the subscription's usage limits, not API billing. API
key variables are removed from its environment so that it cannot fall back to
API billing. Only the exastamp MCP tools are available to it.
--api-model needs ``pip install "anthropic[mcp]"`` and an API key: every
request is billed per token.
"""

import argparse
import asyncio
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time

import yaml

from onika_mcp.catalog import Catalog
from onika_mcp.units import parse_quantity
from exastamp_mcp.render import render_input
from exastamp_mcp.validate import validate_input

HERE = os.path.dirname(os.path.abspath(__file__))
INPUT_NAME = "input.msp"

SYSTEM = (
    "You write exaStamp molecular dynamics input files with the exastamp MCP tools. "
    "Follow the workflow the server describes (get_workflow). Never guess file names or parameters: "
    f"look them up with the tools. When the input validates, save it as {INPUT_NAME} with write_input, "
    "then reply with one short paragraph summarizing the simulation."
)


def _names(items):
    out = []
    for it in items if isinstance(items, list) else [items]:
        if isinstance(it, str):
            out.append(it)
        elif isinstance(it, dict) and len(it) == 1:
            out.append(next(iter(it)))
    return out


def _walk(node):
    """(key, value) of every mapping entry and every string, at any depth."""
    if isinstance(node, dict):
        for k, v in node.items():
            yield k, v
            yield from _walk(v)
    elif isinstance(node, list):
        for v in node:
            if isinstance(v, str):
                yield None, v
            yield from _walk(v)


def _words(doc):
    """Keys and string values of the input: operator names, files, region names..."""
    out = set()
    for k, v in _walk(doc):
        if isinstance(k, str):
            out.add(k.lstrip("+"))
        if isinstance(v, str):
            out.add(v)
    return out


_TEMPERATURE_KEYS = ("T", "Tstart", "Tstop", "Tend")


def _kelvins(doc):
    temps = [v for k, v in _walk(doc) if k in _TEMPERATURE_KEYS]
    out = []
    for t in temps:
        try:
            out.append(parse_quantity(t)[0])
        except Exception:
            pass
    return out


def grade(case, workspace, catalog, reply=None):
    """{"pass", "failed": [...], "validation": report} for the input the agent wrote.
    ``reply`` (the agent's final text) is needed by refusal cases only."""
    path = os.path.join(workspace, INPUT_NAME)
    e = case["expect"]
    if e.get("refusal"):
        failed = [f"{INPUT_NAME} written instead of a refusal"] if os.path.isfile(path) else []
        words = [w for w in e.get("mentions", []) if w.lower() not in (reply or "").lower()]
        if words:
            failed.append(f"reply does not mention {words}")
        return {"pass": not failed, "failed": failed}
    if not os.path.isfile(path):
        return {"pass": False, "failed": [f"no {INPUT_NAME} written"]}
    with open(path) as f:
        text = f.read()
    report = validate_input(text, catalog, workspace)
    failed = [] if report["ok"] else ["validation: " + "; ".join(i["message"][:120] for i in report["issues"] if i["level"] == "error")]
    doc = yaml.safe_load(text) or {}
    # potential: one force operator name, or a list of acceptable ones
    potentials = e["potential"] if isinstance(e.get("potential"), list) else [e.get("potential")]
    if "potential" in e and not set(potentials) & set(_names(doc.get("compute_force"))):
        failed.append(f"potential {e['potential']} not in compute_force")
    # operators: names that must appear as a key or a value anywhere ("a|b": either)
    words = _words(doc)
    for alts in e.get("operators", []):
        if not set(alts.split("|")) & words:
            failed.append(f"{alts} not used")
    # files: data files the input must reference (by name)
    for name in e.get("files", []):
        if not any(w.endswith(name) for w in words):
            failed.append(f"file {name} not referenced")
    if "scheme" in e and doc.get("numerical_scheme", "verlet_nve") != e["scheme"]:
        failed.append(f"numerical_scheme {doc.get('numerical_scheme')} != {e['scheme']}")
    if "species" in e and set(e["species"]) - set(_names(doc.get("species"))):
        failed.append(f"species {e['species']} not all declared")
    if "structure" in e:
        structures = [i["lattice"].get("structure") for i in doc.get("setup_system") or []
                      if isinstance(i, dict) and isinstance(i.get("lattice"), dict)]
        if e["structure"] not in structures:
            failed.append(f"lattice structure {structures} != {e['structure']}")
    if "atoms" in e and (report.get("thermo_after_smoke") or {}).get("particle_count") != e["atoms"]:
        failed.append(f"atoms {(report.get('thermo_after_smoke') or {}).get('particle_count')} != {e['atoms']}")
    if "temperature" in e and not any(abs(t - e["temperature"]) < 1e-6 for t in _kelvins(doc)):
        failed.append(f"temperature {e['temperature']} K not set")
    pe = (report.get("thermo_after_smoke") or {}).get("potential_energy")
    if "max_potential_energy" in e and (pe is None or pe > e["max_potential_energy"]):
        failed.append(f"potential energy {pe} > {e['max_potential_energy']} eV/atom")
    if "steps" in e and (doc.get("global") or {}).get("max_iteration") != e["steps"]:
        failed.append(f"max_iteration {(doc.get('global') or {}).get('max_iteration')} != {e['steps']}")
    if "snapshot" in e and doc.get("write_snapshot") != e["snapshot"]:
        failed.append(f"write_snapshot {doc.get('write_snapshot')} != {e['snapshot']}")
    return {"pass": not failed, "failed": failed, "validation_ok": report["ok"]}


# --------------------------------------------------------------------------- agents

def run_reference(case, workspace, catalog):
    if case["expect"].get("refusal"):
        missing = [op for op in case.get("missing_operators", []) if op in catalog.names()]
        if missing:
            raise RuntimeError(f"refusal case, but {missing} exist in this install")
        return {"tool_calls": 0, "first_try": True, "reply": " ".join(case["expect"].get("mentions", []))}
    if "reference_input" in case:
        with open(os.path.join(HERE, case["reference_input"])) as f:
            text = f.read()
    else:
        text, _ = render_input(case["reference"], catalog.data_dirs, catalog.find_data_file)
    with open(os.path.join(workspace, INPUT_NAME), "w") as f:
        f.write(text)
    return {"tool_calls": 0, "first_try": True}


def copy_workspace_files(case, workspace, catalog):
    """Files the user 'already has': from fixtures/, else from the installed data."""
    for name in case.get("workspace_files", []):
        src = os.path.join(HERE, "fixtures", name)
        if not os.path.isfile(src):
            src = next((os.path.join(d, name) for d in catalog.data_dirs if os.path.isfile(os.path.join(d, name))), None)
        if src is None:
            raise RuntimeError(f"workspace file {name} not found in fixtures/ or the data directories")
        shutil.copy(src, os.path.join(workspace, os.path.basename(name)))


# environment variables that would make the Claude Code CLI bill an API account
_API_AUTH_VARS = ("ANTHROPIC_API_KEY", "ANTHROPIC_AUTH_TOKEN", "ANTHROPIC_BASE_URL",
                  "CLAUDE_CODE_USE_BEDROCK", "CLAUDE_CODE_USE_VERTEX", "CLAUDE_CODE_USE_FOUNDRY")


def run_claude_code(case, workspace, model=None, max_turns=30, timeout=1800):
    claude = shutil.which("claude")
    if not claude:
        raise RuntimeError("the claude CLI is not in PATH (install Claude Code and sign in with 'claude')")
    mcp_config = {"mcpServers": {"exastamp": {
        "command": sys.executable, "args": ["-m", "exastamp_mcp", "--workspace", workspace]}}}
    argv = [claude, "-p", case["prompt"],
            "--mcp-config", json.dumps(mcp_config), "--strict-mcp-config",
            "--tools", "",                          # no built-in tools: only the MCP server
            "--allowedTools", "mcp__exastamp",      # every exastamp tool, without prompts
            "--append-system-prompt", SYSTEM,
            "--max-turns", str(max_turns),
            "--output-format", "stream-json", "--verbose", "--no-session-persistence"]
    if model:
        argv += ["--model", model]
    env = {k: v for k, v in os.environ.items() if k not in _API_AUTH_VARS}
    proc = subprocess.run(argv, cwd=workspace, env=env, capture_output=True, text=True, timeout=timeout)
    calls, result = [], {}
    for line in proc.stdout.splitlines():
        try:
            event = json.loads(line)
        except json.JSONDecodeError:
            continue
        if event.get("type") == "assistant":
            calls += [b["name"].removeprefix("mcp__exastamp__") for b in event["message"].get("content", [])
                      if b.get("type") == "tool_use"]
        elif event.get("type") == "result":
            result = event
    if proc.returncode != 0 and not result:
        raise RuntimeError(f"claude exited with code {proc.returncode}: {(proc.stderr or proc.stdout)[-500:]}")
    if result.get("is_error"):
        raise RuntimeError(f"claude: {result.get('result') or result.get('subtype')}")
    checks = [c for c in calls if c in ("build_input", "validate_input", "write_input")]
    return {"tool_calls": len(calls), "calls": calls, "first_try": len(checks) <= 2,
            "turns": result.get("num_turns"), "reply": result.get("result") or ""}


async def run_claude(case, workspace, model, max_iterations):
    from anthropic import AsyncAnthropic
    from anthropic.lib.tools.mcp import async_mcp_tool
    from mcp import ClientSession
    from mcp.client.stdio import StdioServerParameters, stdio_client

    server = StdioServerParameters(command=sys.executable, args=["-m", "exastamp_mcp", "--workspace", workspace],
                                   env=dict(os.environ))
    calls, reply = [], ""
    async with stdio_client(server) as (read, write):
        async with ClientSession(read, write) as session:
            await session.initialize()
            tools = (await session.list_tools()).tools
            runner = AsyncAnthropic().beta.messages.tool_runner(
                model=model, max_tokens=16000, max_iterations=max_iterations,
                system=SYSTEM, messages=[{"role": "user", "content": case["prompt"]}],
                tools=[async_mcp_tool(t, session) for t in tools],
                betas=["server-side-fallback-2026-07-01"], fallbacks="default",
            )
            async for message in runner:
                calls += [b.name for b in message.content if b.type == "tool_use"]
                reply = "".join(b.text for b in message.content if b.type == "text") or reply
    # first try: one build/validate call, then the write (no fix-and-retry round)
    checks = [c for c in calls if c in ("build_input", "validate_input", "write_input")]
    first = len(checks) <= 2
    return {"tool_calls": len(calls), "calls": calls, "first_try": first, "reply": reply}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--reference", action="store_true", help="no model: grade the reference specs")
    mode.add_argument("--claude-code", action="store_true",
                      help="Claude Code CLI with your Claude account (subscription, no API billing)")
    mode.add_argument("--api-model", metavar="MODEL", help="Anthropic API, billed per token")
    parser.add_argument("--model", help="with --claude-code: model alias or name (default: Claude Code's)")
    parser.add_argument("--cases", default=os.path.join(HERE, "cases.yaml"))
    parser.add_argument("--only")
    parser.add_argument("--workspaces", help="keep each case's workspace under this directory")
    parser.add_argument("--max-iterations", type=int, default=30)
    parser.add_argument("--out")
    args = parser.parse_args(argv)

    with open(args.cases) as f:
        cases = yaml.safe_load(f)
    if args.only:
        keep = set(args.only.split(","))
        cases = [c for c in cases if c["id"] in keep]
    catalog = Catalog("pyexastamp")
    # absolute: the agent and the server run with the workspace as current directory
    root = os.path.abspath(args.workspaces) if args.workspaces else tempfile.mkdtemp(prefix="exastamp-mcp-eval-")
    results = []
    for case in cases:
        ws = os.path.join(root, case["id"])
        os.makedirs(ws, exist_ok=True)
        start = time.monotonic()
        try:
            copy_workspace_files(case, ws, catalog)
            if args.reference:
                run = run_reference(case, ws, catalog)
            elif args.claude_code:
                run = run_claude_code(case, ws, args.model, args.max_iterations)
            else:
                run = asyncio.run(run_claude(case, ws, args.api_model, args.max_iterations))
            g = grade(case, ws, catalog, run.get("reply"))
        except Exception as e:  # report and go on with the other cases
            run, g = {}, {"pass": False, "failed": [f"{type(e).__name__}: {e}"]}
        r = {"id": case["id"], **g, **run, "seconds": round(time.monotonic() - start, 1)}
        results.append(r)
        print(f"{'PASS' if r['pass'] else 'FAIL'} {case['id']:22s} calls={r.get('tool_calls', '-'):>3} "
              f"first_try={r.get('first_try', '-')} {'; '.join(r['failed'])}", flush=True)
    n = len(results)
    passed = sum(r["pass"] for r in results)
    first = sum(bool(r["pass"] and r.get("first_try")) for r in results)
    print(f"\n{passed}/{n} pass, {first}/{n} pass first try, workspaces in {root}")
    if args.out:
        with open(args.out, "w") as f:
            json.dump(results, f, indent=1)
    return 0 if passed == n else 1


if __name__ == "__main__":
    sys.exit(main())
