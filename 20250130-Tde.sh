#!/bin/bash
declare -a PotentialArray=("M2R" "M3R" "BMD192R" "M2" "M3" "BMD192") # "M2R" "M3R" "BMD192R" "M2" "M3" "BMD192" "TabGap2" "TabGap1")

for P in ${PotentialArray[@]}
	do
		echo $P
		cd $P
		for N in {1..41} #50
			do
				#sleep 2
				mkdir $N
				cd $N
				sed -e 's/PotentialTemp/'$P'/g; s/TryTemp/'$N'/g' ../../20250130-Tde.jobtemp > 20250130-Tde.job
				sed -e 's/PotentialTemp/'$P'/g; s/TryTemp/'$N'/g' ../../20250130-Tde.lammpstemp > 20250130-Tde.lammpstemp
				sbatch 20250130-Tde.job
				#chmod +rwx 20250130-Tde.job
				#./20250130-Tde.job
				cd ..
			done
		cd ..
	done