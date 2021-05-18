# ActivePuller
Files to run active puller simulations of a molecule applying a constant tangent force to a polymer as it translocates along the chain.

## Installation Directions

1. Download LAMMPS:
 * git clone -b stable https://github.com/lammps/lammps.git lammps
2. Download ActivePuller Files
 * git clone https://github.com/peter-mm-williams/ActivePuller.git active_puller
3. Move custom source files to src directory in lammps
 * mv active_puller/*pp lammps/src/
4. Build lammps
 * cd lammps/src
 * make yes-molecule
 * make mpi mode=shlib
 * make install-python
 * cd ~

## Example Command to Run Simulation
