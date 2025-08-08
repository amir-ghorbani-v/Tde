#!/bin/bash
declare -a PotentialArray=("M2R") # "M3R" "BMD192R" "M2" "M3" "BMD192") # "TabGap2" "TabGap1")
declare -a SizeArray=(1)

for P in ${PotentialArray[@]}
	do
		echo $P
		#mkdir $P
		cd $P

		for S in ${SizeArray[@]}
			do
				echo $S
				#mkdir $S
				cd $S
				sed -e 's/PotentialTemp/'$P'/g; s/SizeTemp/'$S'/g' ../../20250130-Tde.jobtemp > 20250130-Tde.job
				sed -e 's/PotentialTemp/'$P'/g' ../../20250130-Tde.lammpstemp > 20250130-Tde.lammpstemp
				sbatch 20250130-Tde.job
				#chmod +rwx 20250130-Tde.job
				#./20250130-Tde.job
				cd ..
			done
		cd ..
	done
