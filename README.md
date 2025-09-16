[![CMake](https://github.com/Collab4exaNBody/exaStamp/actions/workflows/cmake.yml/badge.svg)](https://github.com/Collab4exaNBody/exaStamp/actions/workflows/cmake.yml)
[![Spack](https://github.com/Collab4exaNBody/exaStamp/actions/workflows/spack.yml/badge.svg)](https://github.com/Collab4exaNBody/exaStamp/actions/workflows/spack.yml)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

<p align="center">
  <img src="doc/img/exaStamp-logo.png" alt="Example Stamp Logo" width="200"/>
</p>

This is the **exaStamp** software package, a high performance molecular dynamics simulation code, originated at CEA/DAM. **exaStamp** stands for **S**imulations **T**emporelles **A**tomistiques et **M**oléculaires **P**arallèles à l'**exa**scale (in French) or **exa**scale **S**imulations of **T**ime-dependent **A**tomistic and **M**olecular systems in **P**arallel (in English). This software is distributed under the Apache 2.0 public license.

-----------------------------------------------------------------------------------------------------------

**exaStamp** is a high performance classical molecular dynamics simulation code designed to run efficiently on supercomputers as well as laptops or workstations. It takes advantage of hybrid parallelism with the ability to run using MPI + X where X is either OpenMP or Cuda/HIP. **exaStamp** is the result of a long-time effort at CEA/DAM/DIF, France. It is an open-source code, distributed freely under the terms of the Apache Public License version 2.0.

-----------------------------------------------------------------------------------------------------------

## Documentation

**exaStamp** documentation can be accessed using the following QR code and URls (Website and Documentation git repository.

<div align="center">
  <table>
    <tr>
      <td><img src="doc/qr_doc_exastamp.png" width="200"/></td>
      <td><a href="https://collab4exanbody.github.io/doc_exaStamp"> exaStamp Website</a><br><a href="https://github.com/Collab4exaNBody/doc_exaStamp.git"> exaStamp Documentation </a></td>
    </tr>
  </table>
</div>

Main sections of the documentation can be accessed directly using the following URLs.

| URL | Short description |
|-----|-------------|
| https://collab4exanbody.github.io/doc_exaStamp/Background/index.html | **exaStamp** background |
| https://collab4exanbody.github.io/doc_exaStamp/BuildInstall/index.html | Build and installation instructions |
| https://collab4exanbody.github.io/doc_exaStamp/User/index.html | User guide |
| https://collab4exanbody.github.io/doc_exaStamp/Beginner/index.html | Beginner's guide |
| https://collab4exanbody.github.io/doc_exaStamp/microStamp/index.html | microStamp miniApp |
| https://collab4exanbody.github.io/doc_exaStamp/References/index.html | References |
| https://collab4exanbody.github.io/doc_exaStamp/About/index.html | About |
-----------------------------------------------------------------------------------------------------------

## Community Guidelines

For more details, see `CONTRIBUTING.md`. Main guidelines are:

- For any bug, please create an issue and add the label “bug”. We welcome all feedback to make **exaStamp** as robust as possible.
- If you would like to participate and add functionality to **exaStamp**, you can find instructions for coding style, tests and pull request process in CONTRIBUTING.md.
- If you have any support-related / collaboration questions, please contact the team at pauk.lafourcade@cea.fr. If you are a `CEA` member, please request access to the group : "exaNBody & Co. (exaStamp, exaDEM, exaSPH)", an external access can also be provided.

-----------------------------------------------------------------------------------------------------------

## Authors

### Main developers

- Paul Lafourcade (CEA/DAM) (paul.lafourcade@cea.fr)
- Thierry Carrard (CEA/DAM)

### Other Developers

- Raphaël Prat (CEA/DES)
- Claire Lemarchand (CEA/DAM)

-----------------------------------------------------------------------------------------------------------

## License

See `LICENSE.txt`
