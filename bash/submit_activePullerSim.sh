#!/bin/bash
# directories with code

chromatindir=/home/at965/chromatin/active_puller
inputdir=$chromatindir/molfiles
outputdir=/gpfs/ysm/scratch60/ohern/at965/chromatin/output
# make directories, unless they already exist
mkdir -p $outputdir
mkdir -p bin
mkdir -p tasks
mkdir -p slurm
mkdir -p out

molInFile=$1
xyzOutFile=$2
seed=$3
nSteps=$4
dt=$5
kPol=$6
bondAngles=$7
phiMin=$8
extForce=$9
nPols="${10}"
BC="${11}"
boolShortDat="${12}"
partition="${13}"
time="${14}"

basestr=nPols"$nPols"extF"$extForce"nt"$nsteps"dt"$dt"k"$kPol"
runstr="$basestr"_seed"$seed"

simdatadir=$outputdir/$basestr
mkdir -p $simdatadir

echo Running active_puller simulations with parameters:
echo molInFile = "$molInFile"
echo xyzOutFile = "$xyzOutFile"
echo seed = "$seed"
echo nSteps = "$nSteps"
echo dt = "$dt"
echo kPol = "$kPol"
echo phiMin = "$phiMin"
echo extForce = "$extForce"
echo nPols = "$nPols"
echo BC = "$BC"
echo boolShortDat = "$boolShortDat"

# create task file
taskf=tasks/"$runstr".task
rm -f $taskf

runString="cd `pwd`"

filestr="$basestr"_seed"$seed"

inf=$inputdir/$molInFile
xyzf=$simdatadir/$filestr.$xyzOutFile
pyf=$chromatindir/PolymerPol_Bare.py

runString="$runString ; mpirun -np 1 python $pyf -m $inf -x $xyzf -S $seed -t $nSteps -d $dt -k $kPol -b $bondAngles -f $phiMin -F $extForce -N $nPols -B $BC -O $boolShortDat"

runString="$runString ;"

echo "$runString" >> $taskf

# test if task file was created
if [[ ! -f "$taskf" ]]
then
    echo task file not created, ending before job submission
    exit 1
fi

# setup slurm files
slurmf=slurm/"$runstr".slurm
job_name="$runstr"
runout=out/"$runstr"-%a.out
rm -f $slurmf

# echo about time
echo -- running time = $time for $partition

echo -- PRINTING SLURM FILE...
echo \#\!/bin/bash >> $slurmf
echo \#SBATCH --cpus-per-task=1 >> $slurmf
echo \#SBATCH -n 1 >> $slurmf
echo \#SBATCH -p $partition >> $slurmf
echo \#SBATCH -J $job_name >> $slurmf
echo \#SBATCH --mail-type=END,FAIL >> $slurmf
echo \#SBATCH --mail-user=andrewtondata@gmail.com >> $slurmf
echo \#SBATCH -o $runout >> $slurmf
echo module load OpenMPI >> $slurmf
echo module load SciPy-bundle/2021.02-fosscuda-2020b >> $slurmf
echo sed -n \"\$\{SLURM_ARRAY_TASK_ID\}p\" "$taskf" \| /bin/bash >> $slurmf
cat $slurmf

# run sbatch file
echo -- running on slurm in partition $partition
sbatch -t $time $slurmf
