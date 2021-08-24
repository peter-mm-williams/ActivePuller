#!/bin/bash

#to be called from ~/chromatin directory

# directories with code
chromadir=~/chromatin
srcdir=$chromadir/active_puller
maindir=$chromadir/lammps

# directory for all output for chromatin simulations
outputdir=/gpfs/loomis/project/fas/ohern/at965/chromatin

# directory for simulations specific to active pullers
simtypedir=$outputdir/active_puller

# make directories, unless they already exist
mkdir -p $outputdir
mkdir -p $simtypedir
mkdir -p bin
mkdir -p tasks
mkdir -p slurm
mkdir -p out

# inputs
nProc=$1
molFile=$2 #active_puller/molfiles/N100b0.8r1.txt
outputFile=$3 #N100b0.8r1F0.125s1N1pf0.001BC0.xyz
startSeed=$4
numTimesteps=$5
timestep=$6
polAtt=$7
bondAngleBool=$8
packingFraction=$9
extrusionForce="${10}"
numPols="${11}"
boundaryCondition="${12}" #edit readme to not ask me to use stuff like BC and Bool, then remove those variables
partition="${13}"
time="${14}"

# name strings
basestr=ba"$bondAngleBool"_phi"$packinFraction"_force"$extrusionForce"_N"$numPols$"
runstr="$basestr"_startseed"$startSeed"

# make directory specific for this simulation
simdatadir=$simtypedir/$basestr
mkdir -p $simdatadir

# get number of jobs to submit to each array
let arraynum=$fcount
echo -- total number of array runs = $arraynum

# setup slurm files
slurmf=slurm/"$runstr".slurm
job_name="$runstr"
runout=out/"$runstr"-%a.out
rm -f $slurmf

# echo about time
echo -- running time = $time for $partition

echo -- PRINTING SLURM FILE...
echo \#\!/bin/bash >> $slurmf
echo \#SBATCH -n $nProc >> $slurmf
echo \#SBATCH -p $partition >> $slurmf
echo \#SBATCH -J $job_name >> $slurmf
echo \#SBATCH --mail-type=END,FAIL >> $slurmf
echo \#SBATCH --mail-user=andrewtondata@gmail.com >> $slurmf
echo \#SBATCH -o $runout >> $slurmf
echo sed -n \"\$\{SLURM_ARRAY_TASK_ID\}p\" "$taskf" \| /bin/bash >> $slurmf
cat $slurmf

# run sbatch file
echo -- running on slurm in partition $partition
sbatch -t $time $slurmf


# ====================
#       INPUTS
# ====================
# 1. NCELLS
# 2. NV
# 3. calA0
# 4. phiMin
# 5. phiMax
# 6. kl
# 7. kb
# 8. att
# 9. partition
# 10. time
# 11. number of runs (number of array entries, i.e. arraynum)
# 12. start seed (end seed determined by number of runs)
