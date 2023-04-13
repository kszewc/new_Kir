# new_Kir
This repository contains all the data and scripts necessary to replicate the molecular dynamics (MD) simulations and analyze the results for different types of Kir6 channels. The data includes input files for running the simulations and analysis scripts to process the data.
The simulations were performed using GROMACS, systems were prepared using CHARMM-GUI. The resulting trajectories were analyzed using custom Python scripts.

Coel_Kir6x_open - all files needed to run MD simulations of coelacanth open Kir6x (x=1,2,3) homology model (prod1.gro - starting configuration, md.gro-final configuration, md.tpr - portable binary run input file for GROMACS 2019 simulation )


Coel_Kir6x_closed - the same, for the closed Kir6.x conformation

Unbiased_SURy_Kir6x - structures and input files for Kir6x-SURy unbiased simulations (x=1,2,3; y=1,2). eq6.gro - starting configuration, unb.gro final configuration. Additional parameters for SMD simulations are in SMD.mdp (the SMD simulation starts from unb.gro)




