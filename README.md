# main
 
The main library contains the MatLab code for the sensitivity analysis and the parametric model order reduction for the modified Goodwin oscillator model of the circadian clock.
 
 ## Folders
 
 The sensitivity analysis (sa) and the model order reduction (pod) are two sequential but separate steps in this code. The information obtained through 'sa' should be employed in the 'pod' to update the procedure appropriately.
 
### sys
the system : in this case it concerns the modified Goodwin oscillator and contains all details and functions concerning the simulation of the model, as well as the functions common to both other folders.

### sa
the sensitivity analysis : contains all functions necessary for a thorough Sobol and VARS analysis of the influential parameters for the system studied

### pod
the proper orthogonal decomposition : contains all functions necessary for the complete model order reduction through POD. It computes both the standard POD and the parametric POD including the influential parameters described in the sys folder, as well as the standard DEIM.

! Beware, the update of the parameters in the simulation of the model for the parametric snapshot matrix is hard-coded ; if the influential parameters change, an update is necessary in the concerned functions (simulate_sample and simulate_sample_grid).
