#!/bin/bash
origin_path=$(pwd)
echo "action_before_death"
dir=$1
config_file=$2
generation=$3
individual=$4
echo $config_file
source $config_file




cd $dir

#set kill notification
touch KILL_SIGNAL_SET
cp $helper_files/error_files/kappa.dat $dir/kappa.dat
touch ../${generation}_${individual}_DONE


#check if calculations of all individuals are ready
num_finished=$(ls ../ -1f | grep _DONE | wc -l)
echo $num_finished
if [ "$num_finished" -eq "$population_size" ]; then
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
