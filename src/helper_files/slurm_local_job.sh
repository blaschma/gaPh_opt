#!/bin/bash

#This script copies the necessary files to the local storage of the compute node and executes the calculations.
# The results are copied back after the calculations are finished. Implemented for Slurm.

path=$1 # directory of individual
config_file_GA_phonons=$2 # path to config file
generation=$3
individual=$4
source $config_file_GA_phonons


# Print some Job info
env | egrep "(SLURM_TASKS_PER_NODE|SLURM_JOB_NODELIST|SLURM_MEM|SLURM_JOBID|SLURM_JOB_PARTITION)" | sort

# Use local Temp dir to avoid slowdowns due to disk-i/o
export CALCTMPDIR=/ltmp/${SLURM_JOBID}.calc

# Increase stack size
ulimit -s unlimited


# Create TMP directory on each node
${DIROPS_PREFIX} mkdir -p ${CALCTMPDIR}

# only the main node needs the input files, hence no DIROPS_PREFIX is needed
echo -e "\nCopying files ..."
shopt -s extglob
cp -pv "$path"/!(*.out) ${CALCTMPDIR}
cd ${CALCTMPDIR}
echo

# Actual job:

#catch signal
trap  ". $helper_files/action_before_death.sh ${CALCTMPDIR} $config_file $generation $individual"  SIGUSR1

. $helper_files/run_calculations.sh ${CALCTMPDIR} $config_file_GA_phonons $generation $individual $path

echo -e "\nCopying files back ..."
cp -rpuv ${CALCTMPDIR}/* "$path"
cd "$path"

echo -e "\nCleaning up leftover files and folders ..."
${DIROPS_PREFIX} rm -rv /ltmp/${SLURM_JOBID}.calc*

touch ../${generation}_${individual}_DONE

  #check if calculations of all individuals are ready but only if kill signal has not been set
  num_finished=$(ls ../ -1f | grep _DONE | wc -l)
  if [ "$num_finished" -eq "$population_size" ]; then
      file=../KILL_SIGNAL_SET
      if test -f "$file"; then
          echo "kill signal has been set"
      else
          echo "Everybody seems to be ready"
          #eval fitness
          python3 $helper_files/eval_fitness.py $calculation_path/generation_data/$generation/ $config_file

          #invoke next generation
          python3 $genetic_algorithm_path/src/genetic/invoke_next_generation.py $config_file $calculation_path

          #plot fitness values
          python3 $helper_files/plot_fitness_values.py $calculation_path $config_file $calculation_path

          #plot block frequency
          python3 $helper_files/analyze_block_freq.py $config_file

      fi
  else
      echo "Not everybody is ready"
  fi

echo -e "\nAll done!"