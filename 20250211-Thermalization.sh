#!/bin/bash
declare -a PotentialArray=("BMD192") # "TabGap2" "TabGap1")

for P in ${PotentialArray[@]}
	do
		cd $P
		#sleep 2
		sed -e 's/PotentialTemp/'$P'/g' ../20250211-Thermalization.jobtemp > 20250211-Thermalization.job
		sed -e 's/PotentialTemp/'$P'/g' ../20250211-Thermalization.lammpstemp > 20250211-Thermalization.lammpsin
		sbatch 20250211-Thermalization.job
		#chmod +rwx 20250211-Thermalization.job
		#./20250211-Thermalization.job
		cd ..
	done

