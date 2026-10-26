"""
Searchable index of the installed example inputs and of the documentation.

Examples (``<XSP_DATA_DIR>/{samples,tutorial,regression_new}/**.msp``) are
tagged from their YAML: potentials, numerical scheme, lattice structures,
species, readers, outputs, plus the test flags of their first line
(``# no-seq``, ``# enable-if-...``).  Documentation pages (markdown of the
doc_exaStamp site, found through ``EXASTAMP_DOC_DIR``) are split by heading.

Search is by keywords and tags (no embeddings): entries matching more query
words come first, path and tags weigh more than body text.

Some shipped examples do not build with this install (plugin not compiled,
stale parameters): ``index_build_status`` builds every example once in strict
mode and caches the result, which search results then report as ``builds``.
"""

import json
import os
import re
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from fnmatch import fnmatch

import yaml

from .knowledge import doc_url as _doc_url

__all__ = ["Corpus"]

EXAMPLE_CATEGORIES = ("samples", "tutorial", "regression_new")
_FORCE_HINT = re.compile(r"(force|potential)$")

# shipped but obsolete inputs (operators or data removed upstream, no
# replacement): kept in the repository, never offered as examples
EXCLUDED_EXAMPLES = (
    "samples/molecule_archive/*",
    "tutorial/sample[123].msp",
    "samples/task/task_*.msp",
    "samples/task/microjet_meam_insitu.msp",
    "samples/tools/stampv3_to_grid_vtklegacy.msp",
    "regression_new/flexible_molecule/molstampv4/PB_rf_large_intramolecular_compute*.msp",
    "regression_new/flexible_molecule/molxyz/in_exaStamp_polymere.msp",
    "samples/droplet_tin_meam.msp",
)


def _first_line_flags(text):
    first = text.split("\n", 1)[0].strip()
    if not first.startswith("#"):
        return []
    return [t for t in first.lstrip("#").split() if t.startswith(("no-", "enable-if-", "disable-if-"))]


def _names(items):
    out = []
    for it in items if isinstance(items, list) else []:
        if isinstance(it, str):
            out.append(it)
        elif isinstance(it, dict) and len(it) == 1:
            out.append(next(iter(it)))
    return out


def _tags(doc):
    tags = {}
    if not isinstance(doc, dict):
        return tags
    forces = set()
    cf = doc.get("compute_force")
    for name in _names(cf if isinstance(cf, list) else [cf]):
        forces.add(name)
    for k, v in doc.items():
        if _FORCE_HINT.search(str(k)) and isinstance(v, dict):
            forces.add(k)
    if forces:
        tags["potential"] = sorted(forces)
    for key in ("numerical_scheme", "write_snapshot", "write_restart", "grid_flavor"):
        if isinstance(doc.get(key), str):
            tags[key] = [doc[key]]
    species = _names(doc.get("species"))
    if species:
        tags["species"] = species
    setup = doc.get("setup_system")
    ops = _names(setup)
    if ops:
        tags["setup"] = ops
        structures = [p["lattice"].get("structure") for p in setup
                      if isinstance(p, dict) and isinstance(p.get("lattice"), dict) and p["lattice"].get("structure")]
        if structures:
            tags["structure"] = structures
    return tags


def _tokens(text):
    """Words and identifiers of a text: 'eam_alloy_force' gives itself, 'eam', 'alloy', 'force'."""
    text = text.lower()
    idents = re.findall(r"[a-z0-9_.\-]+", text)
    return Counter(idents + [w for i in idents if not i.isalnum() for w in re.findall(r"[a-z0-9]+", i)])


# query word -> words it also matches
_SYNONYMS = {
    "lennard": ["lj"], "jones": ["lj"], "lennard-jones": ["lj"], "water": ["h2o"], "h2o": ["water"],
    "argon": ["ar"], "tantalum": ["ta"], "copper": ["cu"], "tin": ["sn"], "aluminium": ["al"],
    "aluminum": ["al"], "nickel": ["ni"], "iron": ["fe"], "carbon": ["c"],
    "langevin": ["verlet_lnvt", "langevin_thermostat"], "berendsen": ["verlet_bnvt"],
    "nose-hoover": ["verlet_nhnvt", "verlet_nhnpt"], "nvt": ["verlet_lnvt", "verlet_bnvt", "verlet_nhnvt"],
    "npt": ["verlet_nhnpt"], "nve": ["verlet_nve"], "ace": ["pace", "pace_force"], "eam": ["eam_alloy_force"],
    "minimization": ["fire", "conjugate_gradient"], "restart": ["read_dump_atoms", "write_restart_atoms"],
}


def _score(words, head, body, head_weight, body_cap):
    """Score of an entry, or None if no query word is found: entries matching
    more of the query words come first."""
    score, matched = 0, 0
    for w in words:
        best = 0
        for alt in [w] + _SYNONYMS.get(w, []):
            best = max(best, head_weight * head[alt] + min(body[alt], body_cap))
        if best:
            matched += 1
            score += best
    return matched * 1000 + score if matched else None




class Corpus:
    def __init__(self, data_dir, doc_dir=None):
        self.data_dir = os.path.abspath(data_dir)
        self.doc_dir = os.path.abspath(doc_dir) if doc_dir and os.path.isdir(doc_dir) else None
        self._examples = None
        self._docs = None
        self._status = {}
        self.status_file = None

    # ------------------------------------------------------------------ build status
    def load_build_status(self, status_file):
        """Use (and later update) the build status cached in ``status_file``."""
        self.status_file = status_file
        if os.path.isfile(status_file):
            with open(status_file) as f:
                self._status = json.load(f)

    def _mtime(self, rel):
        return int(os.stat(os.path.join(self.data_dir, rel)).st_mtime)

    def build_status(self, rel):
        st = self._status.get(rel)
        if st and st.get("mtime") == self._mtime(rel):
            return st
        return None

    def index_build_status(self, build, workers=4):
        """Build every example whose status is unknown: ``build(text, directory)``
        returns a worker result. Returns (built, failed) counts."""
        todo = [e["path"] for e in self.examples if self.build_status(e["path"]) is None]

        def one(rel):
            path = os.path.join(self.data_dir, rel)
            with open(path, errors="replace") as f:
                r = build(f.read(), os.path.dirname(path))
            return rel, {"mtime": self._mtime(rel), "builds": bool(r.get("ok")),
                         "error": (r.get("error") or "")[:300]}

        with ThreadPoolExecutor(workers) as ex:
            for rel, st in ex.map(one, todo):
                self._status[rel] = st
        if self.status_file:
            os.makedirs(os.path.dirname(self.status_file), exist_ok=True)
            with open(self.status_file, "w") as f:
                json.dump(self._status, f)
        failed = sum(1 for e in self.examples if not (self.build_status(e["path"]) or {}).get("builds", True))
        return len(self.examples) - failed, failed

    def _public(self, e):
        out = {k: v for k, v in e.items() if not k.startswith("_")}
        st = self.build_status(e["path"])
        if st is not None:
            out["builds"] = st["builds"]
            if not st["builds"]:
                out["build_error"] = st["error"]
        return out

    # ------------------------------------------------------------------ examples
    @property
    def examples(self):
        if self._examples is None:
            self._examples = []
            for cat in EXAMPLE_CATEGORIES:
                root = os.path.join(self.data_dir, cat)
                for dirpath, _, files in os.walk(root):
                    for f in sorted(files):
                        rel = os.path.relpath(os.path.join(dirpath, f), self.data_dir)
                        if f.endswith(".msp") and not any(fnmatch(rel, p) for p in EXCLUDED_EXAMPLES):
                            self._examples.append(self._index_example(os.path.join(dirpath, f)))
            self._examples.sort(key=lambda e: e["path"])
        return self._examples

    def _index_example(self, path):
        with open(path, errors="replace") as f:
            text = f.read()
        try:
            doc = yaml.safe_load(text)
        except yaml.YAMLError:
            doc = None
        rel = os.path.relpath(path, self.data_dir)
        comments = [l.lstrip("# ").strip() for l in text.splitlines()[:15] if l.startswith("#")]
        entry = {"path": rel, "category": rel.split(os.sep, 1)[0], "flags": _first_line_flags(text),
                 "tags": _tags(doc), "comment": " ".join(c for c in comments if c and not c.startswith(("no-", "enable-if")))[:300],
                 "_text": text, "_body": _tokens(text)}
        entry["_head"] = _tokens(rel + " " + " ".join(v for vs in entry["tags"].values() for v in vs) + " " + entry["comment"])
        return entry

    def search_examples(self, query="", tags=None, limit=10):
        words = [w.lower() for w in query.split()]
        hits = []
        for e in self.examples:
            if not words and not tags:
                break
            if tags and not all(any(str(v).lower() == str(want).lower() for vs in e["tags"].values() for v in vs)
                                for want in tags):
                continue
            score = _score(words, e["_head"], e["_body"], 3, 3) if words else 0
            if score is None:
                continue
            # prefer short, self-contained inputs that run everywhere
            score -= len(e["flags"]) + e["_text"].count("\n") / 200
            if not (self.build_status(e["path"]) or {}).get("builds", True):
                score -= 10000
            hits.append((score, e))
        hits.sort(key=lambda h: -h[0])
        return [self._public(e) for _, e in hits[:limit]]

    def get_example(self, path):
        full = os.path.normpath(os.path.join(self.data_dir, path))
        if not full.startswith(self.data_dir + os.sep):
            raise ValueError(f"'{path}' is outside the example directory")
        for e in self.examples:
            if e["path"] == path or e["path"].endswith(os.sep + path):
                out = self._public(e)
                out["text"] = e["_text"]
                files = sorted(f for f in os.listdir(os.path.dirname(os.path.join(self.data_dir, e["path"])))
                               if not f.endswith(".msp") and f != "CMakeLists.txt")
                if files:
                    out["files_next_to_it"] = files[:30]
                return out
        raise FileNotFoundError(f"no example '{path}' (use search_examples)")

    def tag_values(self, key):
        values = {}
        for e in self.examples:
            for v in e["tags"].get(key, []):
                values[v] = values.get(v, 0) + 1
        return dict(sorted(values.items(), key=lambda kv: -kv[1]))

    # ------------------------------------------------------------------ docs
    @property
    def docs(self):
        if self._docs is None:
            self._docs = []
            if self.doc_dir:
                for dirpath, _, files in os.walk(self.doc_dir):
                    for f in sorted(files):
                        if f.endswith(".md"):
                            self._docs += self._index_page(os.path.join(dirpath, f))
        return self._docs

    def _index_page(self, path):
        rel = os.path.relpath(path, self.doc_dir)
        with open(path, errors="replace") as f:
            text = f.read()
        text = re.sub(r"\A---\n.*?\n---\n", "", text, flags=re.S)   # front matter
        sections, titles, buf = [], [], []
        for line in text.splitlines():
            m = re.match(r"^(#{1,3})\s+(.*)", line)
            # a heading directly followed by another one (a page title) stays with the next section
            if m and any(l.strip() and not re.match(r"^#{1,6}\s", l) for l in buf):
                sections.append((" > ".join(titles) or rel, "\n".join(buf)))
                titles, buf = [], []
            if m:
                titles.append(m.group(2).strip().strip("*`").strip())
            buf.append(line)
        if buf:
            sections.append((" > ".join(titles) or rel, "\n".join(buf)))
        outline = [t for t, _ in sections]
        return [{"page": rel, "section": t, "url": _doc_url(rel), "outline": outline, "_text": b,
                 "_head": _tokens(rel + " " + t), "_body": _tokens(b)} for t, b in sections if b.strip()]

    def search_docs(self, query, limit=5):
        if not self.doc_dir:
            return []
        words = [w.lower() for w in query.split()]
        hits = []
        for s in self.docs:
            score = _score(words, s["_head"], s["_body"], 5, 5)
            if score is not None:
                hits.append((score, s))
        hits.sort(key=lambda h: -h[0])
        out = []
        for _, s in hits[:limit]:
            hit = {"page": s["page"], "section": s["section"], "url": s["url"], "text": s["_text"][:2500],
                   "page_sections": s["outline"]}
            if len(s["_text"]) > 2500 or len(s["outline"]) > 1:
                hit["more"] = ("section truncated; " if len(s["_text"]) > 2500 else "") + \
                    f"get_doc('{s['page']}') returns the whole page (syntax, parameters, examples)"
            out.append(hit)
        return out

    def get_doc(self, page):
        if not self.doc_dir:
            raise FileNotFoundError("documentation not available (set EXASTAMP_DOC_DIR to doc_exaStamp/docs)")
        full = os.path.normpath(os.path.join(self.doc_dir, page))
        if not full.startswith(self.doc_dir + os.sep) or not os.path.isfile(full):
            raise FileNotFoundError(f"no documentation page '{page}'")
        with open(full, errors="replace") as f:
            return {"page": page, "url": _doc_url(page), "text": f.read()}
