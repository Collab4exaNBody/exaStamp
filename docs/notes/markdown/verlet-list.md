# Neighbors list

Tutorial on Neighbors list

## Methods

To speed up the force computation phase, most md simulations maintain a list of neighbours for each atom. To build such lists efficiently, two typical methods are used and are often combined. The first is based on Verlet lists. A radius of Verlet (rVerlet) is added to rcut. As a consequence, the list of neighbours contains ``useless'' atoms which are too far from the considered atom. However, as long as no atom has moved from its original position (The original position is the one captured when updating the list of neighbours.) by more than 1/2 rverlet, there is not need to recompute the lists. This method is therefore the most efficient when atoms are moving very slowly, because the Verlet lists are not frequently updated. This method have a complexity of O(N^2), and can easily be performed in parallel.

The other one is the linked cells method. The domain is divided into rcut sized cells on a grid. For all atoms, the neighbourhood of an atom is included in the (27) neighbouring cells. This method avoids the cost of maintaining lists of neighbours per atom. However, it comes at the cost of considering a significant amount of atoms which are beyond the rcut distance. This method have a complexity of O(N). The strategy of building neighbour lists can influence the simulation performance because it requires a lot of memory access.

## Linked cell method

Thèse Emmanuel

## Verlet List method


### Files modified

To minimize the duplication of code while waiting for a recast, a boolean, to know if we use the lists of verlet, has been added in the referenceMap (itself included in Global):

Récapitulatif des problèmes:

Problem in integrations schemas: Neighbor lists are clear () after potential calculation. (Possibility to duplicate scheme).

Modification in nodeTimeLoops: We need to calculate whether or not Verlet's lists need to be updated.

Problem in potentials: we need to exclude particles in the distance is greater than the cut-off radius (addition of a correctionVerlet function).
--> Currently the code is duplicated.

Modification for the functions corresponding to the construction of neighbour lists (at cellList level), we use the getRverlet function which corresponds in reality to the maximum between the cut-off radius and the Verlet radius. Note that the Verlet radius is zero if the linked cell method is used. 

SingleSpecGrid problem: between iterations, the order of particles in a cellList is not changed. However, when we update the ghost cells, 
this is filled in parallel (thread). Currently to avoid inconvenience, they are filled sequentially.

Principal modified files

| file                                                          | Description                                                                                    |
| ------------------------------------------------------------- | ---------------------------------------------------------------------------------------------- |
| src/time/verletSchemes.cpp                                    | Add condition                                                                                  |
| src/time/langevinSchemes.cpp                                  | Add condition                                                                                  |
| include/grid/singleSpecGrid/singleSpecGrid_force.hxx          | Duplication of the code in include/grid/singleSpecGrid/singleSpecGrid_force_Verlet.hxx         |  
| include/cellList/particle/cellListParticle_force_*            | Code duplication (PAIR/EAM/MEAM) with buffer correction                                        |
| include/cellList/particle/cellListParticle_neighbors.hxx      | Add condition to store the position of atoms when updating + getRcut replaced by getRverlet    |
| src/parallel/nodeTimeLoop.cpp                                 | Add condition + checkVerlet function (with MPI reduction)                                      |
| include/grid/singleSpecGrid/singleSpecGrid_communications.hxx | Add condition when filling ghosts                                                              |
| include/grid/singleSpecGrid/singleSpecGrid.hpp                | Add function for filling ghosts (no thread parallelization) to keep order between iterations.  |
| src/main.cpp                                                  | Need to fill the reference quickly before initializing nodes/grids/cells                       |
| src/referenceMap.cpp/hpp                                      | Add arrays to store the verlet radius                                                          |

+ Inputs + .hpp
