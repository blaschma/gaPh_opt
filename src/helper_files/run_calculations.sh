#!/bin/bash

calc_path=$1 # directory of individual
config_file=$2 # path to config file
generation=$3
individual=$4
path_in_file_system=$5
source $config_file

#load modules
source $helper_files/load_modules.sh

#echo git commit and branch
commit=$(git --git-dir=$genetic_algorithm_path/.git rev-parse HEAD)
branch=$(git --git-dir=$genetic_algorithm_path/.git symbolic-ref --short HEAD)
echo "Git: On branch $branch. Commit $commit"

#catch signal
trap  ". $helper_files/action_before_death.sh $calc_path $config_file $generation $individual"  SIGUSR1
trap  "echo trap"  SIGSEGV

ulimit -s unlimited
export OMP_STACKSIZE=$omp_stacksize
export OMP_NUM_THREADS=$cpus_per_task,1
export OMP_MAX_ACTIVE_LEVELS=1
export MKL_NUM_THREADS=$cpus_per_task


cd $calc_path

function relax_and_align () {
  #handle anchors for au s anchors
  if [ "$anchor_align" == "T" ];
    then
     echo "Anchor align"    
     if [ "$anchor" -eq "0" ];
     then
      python3 $helper_files/handle_anchors.py ./coord.xyz $config_file ./coord.xyz 
    else
      echo "Warning. Check Anchor align option. Might be mixed up. Anchors are NOT aligned"
     fi
   
  fi

  input_flag=$calc_path/xtb.inp
  if test -f "$input_flag"; then
    echo "Handle input"
    #relaxation
    xtb ./coord.xyz --opt $xtb_level --gfn $gfn --input $calc_path/xtb.inp --verbose > xtb_relaxation.out
    curr_dir=$(pwd)
    coord_path=$curr_dir/xtbopt.xyz
    output_coord=$curr_dir/${generation}_${individual}.xyz
    cp $coord_path $output_coord
    #for no alignment for fixed case
    if grep -q "fix" "$input_flag"; then
      echo "Skipping alignment due to fixed atoms"
    else
      #align gold along x again (and update limits)
      python3 $helper_files/align_gold.py $coord_path $config_file $output_coord > $curr_dir/align.out
    fi

    #calculate hessian with xtb
    xtb --gfn $gfn $output_coord --hess --input $calc_path/xtb.inp --verbose > hessian.out


 
  else
    echo "No input"
    #relaxation
    xtb ./coord.xyz --opt $xtb_level --gfn $gfn > xtb_relaxation.out

    #align gold along x again (and update limits)
    curr_dir=$(pwd)
    coord_path=$curr_dir/xtbopt.xyz
    output_coord=$curr_dir/${generation}_${individual}.xyz
    python3 $helper_files/align_gold.py $coord_path $config_file $output_coord > $curr_dir/align.out

    #calculate hessian with xtb
    xtb --gfn $gfn ${generation}_${individual}.xyz --hess > hessian.out

  fi

  cp ./g98.out ./modes.g98

  python3 $helper_files/x2t.py $output_coord $curr_dir/coord

}


relax_and_align

#check if xtbhess has been created -> more relaxation necessary
xtb_flag=$curr_dir/xtbhess.xyz
if test -f "$xtb_flag"; then
  echo "Relaxation was not sufficient. Relaxing one more time"
  touch .rerelax
  rm -rf coord.xyz
  cp xtbhess.xyz coord.xyz
  rm -rf xtbhess.xyz
  relax_and_align

  #check if xtbhess has been created -> relaxation not succesfull
  xtb_flag=$curr_dir/xtbhess.xyz
  if test -f "$xtb_flag"; then
    echo "Relaxation not succesfull"
    cp $helper_files/error_files/kappa.dat $curr_dir/kappa.dat
    touch RELAXATION_NOT_CONVERGED
  fi

fi

#check if calculation of hessian was successful
hess_flag=$calc_path/hessian
if test -f "$hess_flag";
  then
    #if input file with mass modification exists
    input_flag=$calc_path/xtb.inp
    if test -f "$input_flag"; then
      #check if mass modification is in xtb.inp -> what happens for mass mod and fixing
      if grep -q "modify mass" "$input_flag"; then
        #mass modification in coord file
        python3 $helper_files/mass_mod_coord.py $calc_path/coord $calc_path/xtb.inp $calc_path/coord_MM
        #setup config file for phonon transport calculation
        python3 $helper_files/set_up_phonon_transport_calculation.py phonon_config $curr_dir hessian coord_MM $M_L $M_C $gamma $E_D $N $in_plane $T_min $T_max $kappa_grid_points $anchor ${generation}_${individual}.xyz
      else
        #setup config file for phonon transport calculation
        python3 $helper_files/set_up_phonon_transport_calculation.py phonon_config $curr_dir hessian coord $M_L $M_C $gamma $E_D $N $in_plane $T_min $T_max $kappa_grid_points $anchor ${generation}_${individual}.xyz
      fi

      #start phonon transport calculation
      python3 $phonon_transport_programm/phonon_transport.py $curr_dir/phonon_config
    else 
      #setup config file for phonon transport calculation
      python3 $helper_files/set_up_phonon_transport_calculation.py phonon_config $curr_dir hessian coord $M_L $M_C $gamma $E_D $N $in_plane $T_min $T_max $kappa_grid_points $anchor ${generation}_${individual}.xyz
      #start phonon transport calculation
      python3 $phonon_transport_programm/phonon_transport.py $curr_dir/phonon_config
    fi
  #extended ana
  if [ "$extended_ana" == "T" ]; then
    python3 $helper_files/eval_participation_ratio.py $config_file $curr_dir
    python3 $helper_files/eval_propagator.py $curr_dir/phonon_config $curr_dir
    python3 $phonon_transport_programm/calculate_kappa.py $curr_dir/phonon_config
  fi


  else
  echo "Calculation of hessian was not succesfull"
  cp $helper_files/error_files/kappa.dat $curr_dir/kappa.dat
  touch HESSIAN_NOT_CONVERGED

fi




if [[ "$queuing" != "SLURM_LOCAL" ]]; then

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
fi



