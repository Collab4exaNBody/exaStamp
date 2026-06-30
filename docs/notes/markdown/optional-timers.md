# Optional Timers

Tutorial on optional timer. 

## How to add an optional timer in ExaStamp ?

All timers are included in the metrics class into node class. An additional metricsDetails class has been added to save time for each iteration.

###Into code

In order to facilitate the addition of additional timer, a generic timer has been implemented so that it can be duplicated by the user.

###Into input file

In the optional timer section, add a Chrono_XXX field. Knowing that this one is invalidated if the Chrono field is in false.

## How to instrument your code to take measure ?

We use an extern pointer ptM on the metrics of the node. This solution has chosen because the metrics being in the node, it is impossible to access the node from the domain.


| code                           | Description                                                                    
| ------------------------------ | ------------------------------ 
| ptM->data.tic(Metrics::XXX)    | Begin timer
| ptM->data.toc(Metrics::XXX)    | End timer     
| ptM->data.store(Metrics::XXX)  | Store timer into metrics.                                


## Timers currently available.

| Timer           | Description                                                                    
| --------------- | ------------------------------ 
Chrono_potential  | Corresponding to the doComputeForces() in the domain.
Chrono_neighbours | Corresponding to makeNeighborLists function in domain.
Chrono_ghost      | Corresponding to updateGhost and collectGhost() function in domain.
Chrono_refine     | Corresponding to nothing, it will be used for the adaptive mesh refinement. 
Chrono_generic    | Do nothing.
