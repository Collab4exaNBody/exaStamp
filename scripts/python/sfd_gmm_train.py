#
# Licensed to the Apache Software Foundation (ASF) under one
# or more contributor license agreements. See the NOTICE file
# distributed with this work for additional information
# regarding copyright ownership. The ASF licenses this file
# to you under the Apache License, Version 2.0 (the
# "License"); you may not use this file except in compliance
# with the License. You may obtain a copy of the License at
#   http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on an
# "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied. See the License for the
# specific language governing permissions and limitations
# under the License.
#
"""Train a Gaussian mixture model on strain functional descriptors for sfd_gmm_classify.

Reads one or more extended xyz files written by exaStamp's write_xyz with sfd_<column> fields
(see strain_functional_fields), fits scikit-learn GaussianMixture models (full covariances)
for each requested number of components, prints log-likelihood, BIC, Calinski-Harabasz and
Davies-Bouldin scores, and writes the selected model as YAML (model_file of sfd_gmm_classify).
Components are sorted by decreasing weight, as in sfd_gmm_fit.

Example:
  python sfd_gmm_train.py snap_*.xyz -k 2 3 4 5 6 --select chi -o gmm_model.yaml

Requires numpy, scikit-learn and PyYAML.
"""
import argparse
import re
import sys

import numpy as np
import yaml
from sklearn.metrics import calinski_harabasz_score, davies_bouldin_score
from sklearn.mixture import GaussianMixture

DEFAULT_FEATURES = ["P1_I0", "P2_I0", "P2_I1", "P2_I2",
                    "P3_I0", "P3_I1", "P3_I2", "P3_I3", "P3_I4", "P3_I5", "P3_I6",
                    "P4_I0", "P4_I1", "P4_I2", "P4_I3", "P4_I4", "P4_I5", "P4_I6", "P4_I7",
                    "P4_I8", "P4_I9", "P4_I10", "P4_I11"]


def read_xyz_features(file_name, features, prefix):
    """Feature matrix of all frames of an extended xyz file."""
    rows = []
    with open(file_name) as f:
        lines = f.read().splitlines()
    pos = 0
    while pos < len(lines) and lines[pos].strip():
        n = int(lines[pos])
        m = re.search(r'Properties=(\S+)', lines[pos + 1])
        if m is None:
            sys.exit(f"{file_name}: no Properties= in the comment line")
        props = m.group(1).split(':')
        cols = []
        for i in range(0, len(props), 3):
            name, count = props[i], int(props[i + 2])
            cols += [name] if count == 1 else [f"{name}{j}" for j in range(count)]
        index = {c: k for k, c in enumerate(cols)}
        missing = [f for f in features if prefix + f not in index]
        if missing:
            sys.exit(f"{file_name}: missing columns {', '.join(prefix + f for f in missing)}")
        sel = [index[prefix + f] for f in features]
        for line in lines[pos + 2:pos + 2 + n]:
            v = line.split()
            rows.append([float(v[k]) for k in sel])
        pos += 2 + n
    return np.array(rows)


def sorted_model(gm):
    order = np.argsort(-gm.weights_, kind='stable')
    return gm.weights_[order], gm.means_[order], gm.covariances_[order], order


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("xyz", nargs="+", help="exaStamp extended xyz files with sfd_<column> fields")
    p.add_argument("-k", "--n-components", type=int, nargs="+", default=[2], help="numbers of components to fit")
    p.add_argument("--features", nargs="+", default=DEFAULT_FEATURES, help="descriptor columns (default: the 23 invariants P1_I0 ... P4_I11)")
    p.add_argument("--prefix", default="sfd_", help="field name prefix used by strain_functional_fields")
    p.add_argument("--standardize", action="store_true", help="center and scale each feature before fitting")
    p.add_argument("--select", choices=["dbi", "chi", "bic"], default="dbi", help="model selection criterion (default: lowest Davies-Bouldin)")
    p.add_argument("--reg-covar", type=float, default=1e-6)
    p.add_argument("--tol", type=float, default=1e-3)
    p.add_argument("--max-iter", type=int, default=100)
    p.add_argument("--n-init", type=int, default=1)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--max-samples", type=int, default=0, help="random subsample of the particles used for fitting (0: all)")
    p.add_argument("-o", "--output", default="sfd_gmm_model.yaml", help="output YAML model")
    a = p.parse_args()

    X = np.concatenate([read_xyz_features(f, a.features, a.prefix) for f in a.xyz])
    rng = np.random.default_rng(a.seed)
    if 0 < a.max_samples < len(X):
        X = X[rng.choice(len(X), a.max_samples, replace=False)]
    D = X.shape[1]
    center, scale = np.zeros(D), np.ones(D)
    if a.standardize:
        center = X.mean(axis=0)
        scale = X.std(axis=0)
        scale[scale == 0.0] = 1.0
    Xs = (X - center) / scale

    print(f"{len(X)} particles, {D} features" + (" (standardized)" if a.standardize else ""))
    print("  n_components  iterations  converged  log-likelihood            BIC        CHI        DBI")
    best, best_score = None, None
    for k in a.n_components:
        gm = GaussianMixture(k, covariance_type="full", tol=a.tol, reg_covar=a.reg_covar,
                             max_iter=a.max_iter, n_init=a.n_init, random_state=a.seed).fit(Xs)
        labels = gm.predict(Xs)
        ll, bic = gm.score(Xs), gm.bic(Xs)
        ok = len(np.unique(labels)) > 1
        chi = calinski_harabasz_score(Xs, labels) if ok else 0.0
        dbi = davies_bouldin_score(Xs, labels) if ok else 0.0
        print(f"  {k:12d}{gm.n_iter_:12d}{'yes' if gm.converged_ else 'no':>11s}{ll:16.8g}{bic:15.8g}{chi:11.5g}{dbi:11.5g}")
        score = {"bic": bic, "chi": -chi, "dbi": dbi}[a.select]
        if best is None or score < best_score:
            best, best_score = gm, score

    w, mu, cov, _ = sorted_model(best)
    model = {
        "features": list(a.features),
        "center": center.tolist(),
        "scale": scale.tolist(),
        "weights": w.tolist(),
        "means": mu.tolist(),
        "covariances": cov.tolist(),
    }
    with open(a.output, "w") as f:
        f.write("# Gaussian mixture model on strain functional descriptors (sfd_gmm_fit / sfd_gmm_train.py)\n")
        yaml.safe_dump(model, f, default_flow_style=None, sort_keys=False)
    print(f"selected n_components = {best.n_components}, weights: " + " ".join(f"{x:.5g}" for x in w))
    print(f"model written to {a.output}")


if __name__ == "__main__":
    main()
