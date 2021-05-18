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
mpirun -np 1 python active_puller/PolymerPol_Bare.py -m active_puller/molfiles/N100b0.8r1.txt -x "insert_output_path/"N100b0.8r1F0.125s1N1pf0.001BC0.xyz -S 9778 -t 200000000 -d 0.001 -k 1000 -b 0 -p 164 -f 0.001 -s 1  -F 0.125 -N 1 -B 0

## Module dependencies:
Below are the modules that must be loaded to run the code on the two clusters to which the O'Hern group has access
#### Grace:
* module load Langs/Python
* module load Libs/MPI4PY
* module load Libs/NUMPY

#### Farnam:
* module load pip
