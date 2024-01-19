# exaStamp

![](doc/img/exaStamp-logo.png)
 
exaStamp is a high performance molecular dynamics simulation code, originated from CEA/DIF.

## Important note for developpers working with legacy repository history
Please check instructions in doc/transfer_from_legacy_repo.txt

## Installation

### Minimal Requirements

To proceed with the installation, your system must meet the minimum prerequisites. The first step involves the installation of exaNBody:

```
git clone https://github.com/Collab4exaNBody/exaNBody.git
export exaNBody_DIR=`pwd`/exaNBody
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
- Cuda, see https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html
- HIP, see https://rocm.docs.amd.com/projects/install-on-linux/en/latest/

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

Run CMake to configure the exaStamp build, with CUDA support

```
ccmake -DXNB_PRODUCT_VARIANT=all \
       -DXNB_BUILD_CUDA=ON \
       -DXSTAMP_CUDA_ARCH=86 \
       -DONIKA_HAVE_OPENMP_DETACH=OFF \
       -DONIKA_HAVE_OPENMP_TOOLS=OFF \
       -DCMAKE_BUILD_TYPE=Release \
       ../exaStamp
``

OR with HIP support

```
ccmake -DXNB_PRODUCT_VARIANT=rigidmol \
       -DCMAKE_C_COMPILER=hipcc \
       -DCMAKE_CXX_COMPILER=hipcc \
       -DCMAKE_CXX_FLAGS="-Wpass-failed" \
       -DXNB_BUILD_CUDA=ON \
       -DXNB_ENABLE_HIP=ON \
       -DCMAKE_HIP_PLATFORM=amd \
       -DONIKA_HAVE_OPENMP_DETACH=OFF \
       -DONIKA_HAVE_OPENMP_TOOLS=OFF \
       -DCMAKE_BUILD_TYPE=Release \
       ../exaStamp
```

Build exaStamp using the make command with a specified number of parallel jobs (e.g., -j 4 for 4 parallel jobs):

```
make -j 4
```

Build Plugins

```
make UpdatePluginDataBase
```
