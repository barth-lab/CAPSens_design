#!/bin/bash

#SBATCH --job-name ISM
#SBATCH -n 1                    # number of cores
#SBATCH --partition=serial       # fidis serial partition
#SBATCH -t 0-04:00         # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out      # STDOUT
#SBATCH -e slurm.%N.%j.err      # STDERR

touch STARTED

Rosetta_Intel/bin/rosetta.intel -s EnsembleState1_WT_rpk.pdb -nstruct 200 -resfile input/EnsembleState1_WT_rpk_N11-Y19-S259_5.0A.resfile -pose1 -cst_mode -cst_design -fa_input -pose_memb -memb_solv -memb_env -Wmbenv 0.55 -memb_hb -thickness 12.5 -steepness 10 -spanfile input/EnsembleState1_WT_rpk.span -cst_min -chi_move -bb_move -paths input/paths.txt -ex1 -ex1aro -ex2 -ex3 -ex4 -extrachi_cutoff 0 -seed_offset $$

touch FINISHED
