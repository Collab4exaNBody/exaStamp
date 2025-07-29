# licensed to the apache software foundation (asf) under one
# or more contributor license agreements.  see the notice file
# distributed with this work for additional information
# regarding copyright ownership.  the asf licenses this file
# to you under the apache license, version 2.0 (the
# "license"); you may not use this file except in compliance
# with the license.  you may obtain a copy of the license at
# 
#   http://www.apache.org/licenses/license-2.0
# 
# unless required by applicable law or agreed to in writing,
# software distributed under the license is distributed on an
# "as is" basis, without warranties or conditions of any
# kind, either express or implied.  see the license for the
# specific language governing permissions and limitations
# under the license.

#! /usr/bin/env python
import random

data_number = int(float(input("How many lattices do you want to test? Please type an integer... ")))

print("OK, generating a .dat file containing lattice parameters for %d different lattices"%(data_number))

f_output = open("lengths_angles_input.dat", "w")
f_output.write("# a0 \tb0 \tc0 \talpha \tbeta \tgamma\n")

for i in range(0,data_number):
    a0 = random.uniform(3.0,6.0)
    b0 = random.uniform(3.0,6.0)
    c0 = random.uniform(3.0,6.0)
    alpha = random.uniform(60.,120.)
    beta  = random.uniform(60.,120.)
    gamma = random.uniform(60.,120.)
    f_output.write("%5.4f \t%5.4f \t%5.4f \t%5.4f \t%5.4f \t%5.4f\n"%(a0,b0,c0,alpha,beta,gamma))
    
f_output.close()    
