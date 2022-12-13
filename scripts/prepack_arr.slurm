#!/bin/bash
#SBATCH --time 01:00:00 -n 1
#SBATCH --array=1-4
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=36
#SBATCH --output slurm.out
#SBATCH --error slurm.err

R="/work/upbarth/Programs/rosetta"

declare -a commands

S=""
if [ -f ../span ]; then
	S="-in:file:spanfile ../span -membrane:Membed_init -membrane::Mhbond_depth -score:weights membrane_highres"
fi

for i in {1..4}
do
	commands[$i]=""
	for j in {1..36};do
		n=$(((i-1)*36+(j-1)))
		printf -v sn "%03d" $n
		commands[$i]="${commands[$i]} 
			$R/source/bin/FlexPepDocking.linuxgccrelease $S \
			-database $R/database/  \
			-s ${sn}.pdb \
			-flexpep_prepack\
			-ex1 -ex2aro \
			-nstruct 1 \
			-suffix "_pp" \
			-no_nstruct_label \
			> log_$sn \
			2> err_$sn &"
	done
	commands[$i]="${commands[$i]}
		wait"
done

bash -c "${commands[${SLURM_ARRAY_TASK_ID}]}"

