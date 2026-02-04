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
#!/bin/bash -eu

# param covers 2 examples:
# 1. naive example (all parameters are ones)
# 2. /mnt/cesimat/exastamp/SNAP/CodeSaclaySNAP example
#    "3.5 2" (CodeSaclaySNAP example) <=> param "7. 1" (same thing but with ntype = 1)

for param in "1.0 1" "3.5 2"
do
  echo ""
  echo "param: $param"
  IFS=" " read -r jmax ntype <<< "$param"
  echo "jmax: $jmax"
  echo "ntype: $ntype"

  # Run fortran reference code => produce cg.ref

  echo ""
  ./snap-compute-cg-f90 "$jmax" "$ntype"

  # Run cpp code (+ diff with cg.ref)

  echo ""
  ./snap-compute-cg-cpp "$jmax" "$ntype"

  echo ""
  echo "============================================="



done
