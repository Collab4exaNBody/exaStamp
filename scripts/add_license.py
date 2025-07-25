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
import os, sys

# Define the license headers
license_header_c_style = """Apache Software Foundation (ASF)"""
license_header_c_style_full = """/*
Licensed to the Apache Software Foundation (ASF) under one
or more contributor license agreements. See the NOTICE file
distributed with this work for additional information
regarding copyright ownership. The ASF licenses this file
to you under the Apache License, Version 2.0 (the
"License"); you may not use this file except in compliance
with the License. You may obtain a copy of the License at
  http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing,
software distributed under the License is distributed on an
"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
KIND, either express or implied. See the License for the
specific language governing permissions and limitations
under the License.
*/

"""

license_header_python_style = """Apache Software Foundation (ASF)"""
license_header_python_style_full = """#
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
"""

# Specify the root directory containing the files
root_directory = sys.argv[1]

# Walk through the directory tree
for foldername, subfolders, filenames in os.walk(root_directory):
    for filename in filenames:
        file_path = os.path.join(foldername, filename)

        # Determine the file extension
        _, extension = os.path.splitext(filename)

        # Open the file and read its contents
        with open(file_path, 'r+') as file:
            content = file.read()

            # Check if the header is already present
            if extension in ['.cpp', '.h', '.cxx', '.hxx', '.cu']:
                if license_header_c_style.strip() not in content:
                    # Prepend the license header if not present
                    file.seek(0, 0)
                    file.write(license_header_c_style_full + content)
            elif extension in [ '.txt', '.sh', '.py' ]:                    
                if license_header_python_style.strip() not in content:
                    # Prepend the license header if not present
                    file.seek(0, 0)
                    file.write(license_header_python_style_full + content)

print("License headers added to all relevant files in the directory tree, if not already present.")
