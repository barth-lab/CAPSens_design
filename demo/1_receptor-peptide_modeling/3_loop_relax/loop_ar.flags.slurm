#!/bin/bash
#SBATCH --time 48:00:00 -n 1
#SBATCH --array=1-6
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=36
#SBATCH --output slurm.out
#SBATCH --error slurm.err

R="/work/upbarth/Programs/rosetta"

#O=output_loop_diverse_relax_movemap
O=output_loop
I=input

if [ ${SLURM_ARRAY_TASK_ID} -eq 1 ]
then
	if [ -d "$O" ]; then
		D=`date -r $O +"%F_%T"`
		mv $O $O.$D
	fi
	mkdir $O
else
	sleep 10
fi

S=""
if [ -f $I/*.span ]; then
	S="-in:file:spanfile $I/CXCR4_SDF1.span -membrane:Membed_init -membrane::Mhbond_depth -score:weights membrane_highres"
fi

declare -a commands

for i in {1..6}
do
	commands[$i]=""
	for j in {1..36};do
		n=$(((i-1)*36+(j-1)))
		printf -v sn "%03d" $n
		commands[$i]="${commands[$i]}
		$R/source/bin/loopmodel.linuxgccrelease $S \
		-database $R/database/ \
		-in:file:s $I/diverse/diverse_$sn.pdb \
		-seed_offset $$ \
		-nstruct 100 \
		\@flags \
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
