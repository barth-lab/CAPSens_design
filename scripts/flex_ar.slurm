#!/bin/bash
#SBATCH --time 12:00:00 -n 1
#SBATCH --array=1-4
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=36
#SBATCH --output slurm.out
#SBATCH --error slurm.err

R="/work/upbarth/Programs/rosetta"

I=input
O=output_flex

if [ ${SLURM_ARRAY_TASK_ID} -eq 1 ]
then
	D=`date +"%s"`
	mv $O $O.$D
	mkdir $O
else
	sleep 10
fi

S=""
if [ -f $I/CXCR4_SDF1.span ]; then
	S="-in:file:spanfile $I/span -membrane:Membed_init -membrane::Mhbond_depth -score:weights membrane_highres"
fi

declare -a commands

for i in {1..4}
do
	commands[$i]=""
	for j in {1..36};do
		n=$(((i-1)*36+(j-1)))
		printf -v sn "%03d" $n
		commands[$i]="${commands[$i]}
			$R/source/bin/FlexPepDocking.linuxgccrelease  $S \
			-database $R/database/  \
			-s $I/grid/${sn}_pp.pdb \
			-pep_refine \
			-ex1 -ex2aro -use_input_sc \
			-constant_seed \
			-seed_offset $n \
			-nstruct 100 \
			-out:path:all $O \
			-out:suffix _$sn \
			-out:file:silent out_${sn}.silent \
			> $O/log_$sn \
			2> $O/err_$sn &"
	done
	commands[$i]="${commands[$i]}
		wait"
done

bash -c "${commands[${SLURM_ARRAY_TASK_ID}]}"
