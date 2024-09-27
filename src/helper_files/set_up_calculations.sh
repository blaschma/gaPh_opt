#!/bin/bash
#this script sets up the calculations.
# $1: path to dir where turbo calc should be initialized. Must be prepared with coord and limits -> genome_to_molecule.py
# $2: path to config file
# $3: generation
# $4: individual

generation=$3
individual=$4

origin=$(pwd)
#load config file
config_file_GA_phonons=$2
source $config_file_GA_phonons

#set path
path=$1

#set cpus per task
echo $cpus_per_task

echo $queuing
if [[ "$queuing" == "SLURM" ]]; then
  echo "Using Slurm"
	sbatch --job-name=gen$3id$4 --mem-per-cpu=$mem_per_cpu --partition=$partition --time=$max_time --ntasks=1 --cpus-per-task=$cpus_per_task --signal=B:SIGUSR1@$kill_time --output $path/slurm_output.txt $helper_files/run_calculations.sh $path $config_file_GA_phonons $generation $individual
elif [[ "$queuing" == "SLURM_LOCAL" ]]; then
  echo "Using Slurm local"
	sbatch --job-name=gen$3id$4 --mem-per-cpu=$mem_per_cpu --partition=$partition --time=$max_time --ntasks=1 --cpus-per-task=$cpus_per_task --signal=B:SIGUSR1@$kill_time --output $path/slurm_output.txt $helper_files/slurm_local_job.sh $path $config_file_GA_phonons $generation $individual
elif [[ "$queuing" == "GE" ]]; then
    echo "Using Grid engine (GE)"
    qsub -N gen$3id$4p -cwd -q scc -pe openmp $cpus_per_task -l h_vmem=$mem_per_cpu $helper_files/run_calculations.sh $path $config_file_GA_phonons $generation $individual
elif [[ "$queuing" == "None" ]]; then
    echo "No Queuing system"
    . $helper_files/run_calculations.sh $path $config_file_GA_phonons $generation $individual
else
    echo "Unknown queuing system"
    return -1
fi



