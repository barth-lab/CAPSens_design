#!/bin/bash

printf "#%-7s %-8s %-8s %-8s %-8s %-8s %-8s %-8s %-8s %-8s %-8s %-8s %-8s %-8s %-8s %-8s %-8s %-8s %-8s %-s\n" "I_sc" "pep_sc" "I+pep" \
	"com1" "com2" "com3" "w1" "w2" "w3" \
	"v11" "v12" "v13" "v21" "v22" "v23" "v31" "v32" "v33" \
	"rmsd_na" "tag"
if [ ! -z $1 ]; then
	cd $1
fi
for f in *.pdb; do
	t=${f%.pdb}
	is=$(grep "I_sc " $f | awk '{print $NF}')
	ps=$(grep "pep_sc " $f | awk '{print $NF}')
	ip=$(echo "$is + $ps" | bc -l)
	rms=$(grep "rmsBB_allIF " $f | awk '{print $NF}')
	g=$(extract -cB $f | geometry)
	printf "%8.2f %8.2f %8.2f " $is $ps $ip
	printf "%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f " $g
	printf "%8.3f %s\n" $rms $t
done

