# exaStamp

exaStamp is a high performance molecular dynamics simulation code, originated from CEA/DIF.
 
## Important note for developpers working with legacy repository history
Please check instructions in doc/transfer_from_legacy_repo.txt

## Installation

### Minimal Requirements

To proceed with the installation, your system must meet the minimum prerequisites. The first step involves the installation of exaNBody:

```
git clone https://github.com/Collab4exaNBody/exaNBody.git
mkdir build-exaNBody/ && cd build-exaNBody/
cmake ../exaNBody/ -DCMAKE_INSTALL_PREFIX=path_to_install
make install
export exaNBody_DIR=path_to_install
```

The next step involves the installation of yaml-cpp and MPI, which can be achieved using either the spack package manager or cmake:

```
sudo apt install libyaml-cpp-dev
sudo apt install mpi-default-dev
```

### Optional Dependencies

Before proceeding further, you have the option to consider the following dependencies:

- Cuda
- Zoltan and parmetis

### exaStamp installation

To install exaStamp, follow these steps:

Set the exaNBody_DIR environment variable to the installation path. Clone the exaStamp repository using the command:

```
git clone https://github.com/Collab4exaNBody/exaStamp.git
```

Create a directory named build and navigate into it:

```
mkdir build && cd build
```

Run CMake to configure the exaStamp build, specifying that CUDA support should be turned on (or off):

```
cmake ../exaStamp -DXSTAMP_BUILD_CUDA=OFF 
```

Build ExaDEM using the make command with a specified number of parallel jobs (e.g., -j 4 for 4 parallel jobs):

```
make -j 4
```

Build Plugins

```
make UpdatePluginDataBase
```
