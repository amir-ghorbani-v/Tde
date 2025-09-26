#!/bin/bash
declare -a TemperatureArray=("10")
declare -a PotentialArray=("BMD192") # "M2R" "M3R" "BMD192R" "M2" "M3" "BMD192" "TabGap1" "TabGap2" "MtpPd")
Tries=41
read -r -d '' directions_block <<'EOF'
directions=(
"0 0 1"
"0 1 1"
"0 3 1"
"0 1 0"
"-1 1 0"
"-1 1 1"
"-1 1 2"
"1 2 0"
"1 2 1"
"1 2 2"

"-2 3 3"
"-1 4 1"
"-1 3 1"
"-1 2 2"
"0 2 1"
"-1 2 0"
"-1 2 3"
"-4 4 1"
"-3 4 1"
"-2 4 1"
"-1 4 0"
"-1 2 1"
"0 2 3"
"-1 4 2"
"-2 2 3"
"-3 4 3"
"-2 3 2"
"-2 2 1"
"-1 3 0"
"-1 3 2"
"-3 4 2"
"-1 4 3"
"0 4 1"
"-3 4 0"
"-2 3 0"
"0 1 2"
"-2 3 1"
"-1 3 3"

"1 0 0"
"2 0 1"
"4 1 3"
)
EOF

########################################################################

host=${SLURMD_NODENAME:-$(hostname)}
ClustersAliancePrivileged=("rorqual")
ClustersAlianceNormal=("narval" "nibi" "rorqual")
ClustersToronto=("trillium" "niagara")
ClusterFrontenac=("frontenac")
Host=$(scontrol show config | awk -F= '/ClusterName/ {print $2}' | xargs)

sinfo_out=$(sinfo -o "%P %c" -h | sort -u)
Ntasks=64
if grep -q "^cpularge_interac" <<< "$sinfo_out"; then
    val=$(awk '$1 ~ /^cpularge_interac/ {print $2}' <<< "$sinfo_out")
    if [[ "$val" =~ ^[0-9]+$ ]]; then
        Ntasks=$val
    fi
elif grep -q "^cpubase_interac" <<< "$sinfo_out"; then
    val=$(awk '$1 ~ /^cpubase_interac/ {print $2}' <<< "$sinfo_out")
    if [[ "$val" =~ ^[0-9]+$ ]]; then
        Ntasks=$val
    fi
fi
echo 

PotentialsEam=("M2" "M3" "BMD192" "M2R" "M3R" "BMD192R")
PotentialsTabGap=("TabGap1" "TabGap2")
PotentialsMtp=("MtpPd")

if [[ " ${ClustersAliancePrivileged[*]} " =~ " ${Host} " ]]; then
    PotDir="/home/veshand/Zr/PotentialBank"
    LmpEam="srun lmp"
    LmpTapGap="srun /home/veshand/Software/Lammps-tabGAP/lmp4_CC"
    LmpMpt="mpirun -np ${Ntasks} /home/veshand/Software/lammps-mtp/build/lmp"
    Account="rrg-belandl1"
    Tasks="#SBATCH --ntasks=${Ntasks}"
    MemPerCpu="#SBATCH --mem-per-cpu=400M"
    Time="#SBATCH --time=6-23:59"

elif [[ " ${ClustersAlianceNormal[*]} " =~ " ${Host} " ]]; then
    PotDir="/home/veshand/Zr/PotentialBank"
    LmpEam="srun lmp"
    LmpTapGap="srun /home/veshand/Software/Lammps-tabGAP/lmp4_CC"
    LmpMpt="mpirun -np ${Ntasks} /home/veshand/Software/lammps-mtp/build/lmp"
    Account="#SBATCH --account=def-belandl1"
    Tasks="#SBATCH --ntasks=${Ntasks}"
    MemPerCpu="#SBATCH --mem-per-cpu=400M"
    Time="#SBATCH --time=6-23:59"

elif [[ " ${ClustersToronto[*]} " =~ " ${Host} " ]]; then
    PotDir="/home/veshand/Zr/PotentialBank"
    LmpEam="srun lmp"
    LmpTapGap="srun /home/veshand/Software/Lammps-tabGAP/lmp4_CC"
    LmpMpt="mpirun -np ${Ntasks} /home/veshand/Software/lammps-mtp/build/lmp"
    Account="#SBATCH --account=def-belandl1"
    Tasks="#SBATCH --ntasks=${Ntasks}"
    MemPerCpu="#SBATCH --mem-per-cpu=400M"
    Time="#SBATCH --time=6-23:59"

elif [[ " ${ClusterFrontenac[*]} " =~ " ${Host} " ]]; then
    PotDir="/global/home/hpc4974/Zr/PotentialBank"
    LmpEam="srun lmp"
    LmpTapGap="srun /home/veshand/Software/Lammps-tabGAP/lmp4_CC"
    LmpMpt="mpirun -np ${Ntasks} /global/home/hpc4974/Setup/lammps-mtp/build/lmp"
    Account="#SBATCH --account=def-hpcg1725"
    Tasks="#SBATCH --ntasks=${Ntasks}"
    MemPerCpu="#SBATCH --mem-per-cpu=400M"
    Time="#SBATCH --time=6-23:59"

else
    echo "Unknown host ($Host). Update detection logic."
    exit 1
fi

for T in "${TemperatureArray[@]}"; do
    echo "$T"
    cd "$T" || exit
    for P in "${PotentialArray[@]}"; do
        echo "$P"
        cd "$P" || exit
        for N in $(seq 1 $Tries); do
            mkdir "$N"
            cd "$N" || exit
            sed -e "s|AccountTemp|${Account}|g; \
                    s|TasksTemp|${Tasks}|g; \
                    s|MemPerCpuTemp|${MemPerCpu}|g; \
                    s|TimeTemp|${Time}|g; \
                    s|TemperatureTemp|${T}|g; \
                    s|PotentialTemp|${P}|g; \
                    s|TryTemp|${N}|g; \
                    s|PotDirTemp|${PotDir}|g; \
                    s|LmpEamTemp|${LmpEam}|g; \
                    s|LmpTapGapTemp|${LmpTapGap}|g; \
                    s|LmpMptTemp|${LmpMpt}|g" \
                -e "/^directions=\(\)/r /dev/fd/3" \
                -e "/^directions=\(\)/,+0d" \
                ../../../20250926-Tde.jobtemp 3<<<"$directions_block" > 20250926-Tde.job

            if [[ " ${PotentialsEam[*]} " =~ " ${P} " ]]; then
                sed -i 's/^echo __________TabGap__________$/#&/' 20250926-Tde.job
                sed -i 's/^module load nixpkgs\/16.09  StdEnv\/2020  intel\/2020.1.217  openmpi\/4.0.3 lammps-omp\/20201029$/#&/' 20250926-Tde.job
                sed -i 's/^module load hdf5-mpi\/1.10.6 voro++\/0.4.6 eigen\/3.3.7$/#&/' 20250926-Tde.job
                sed -i 's/^module load mlip\/2.0 fftw-mpi\/3.3.8 tbb\/2020.2$/#&/' 20250926-Tde.job
                sed -i 's/^echo __________Mtp__________$/#&/' 20250926-Tde.job
                sed -i 's/^module load StdEnv\/2023$/#&/' 20250926-Tde.job

                if [[ "$Host" == "niagara" ]]; then
                    sed -i 's/^module load StdEnv\/2023  intel\/2023.2.1  openmpi\/4.1.5 lammps-omp\/20230802$/#&/' 20250926-Tde.job
                else
                    sed -i 's/^module load intel\/2019u3  intelmpi\/2019u3 lammps\/29Mar2019 # for Niagara$/#&/' 20250926-Tde.job
                fi

            elif [[ " ${PotentialsTabGap[*]} " =~ " ${P} " ]]; then
                sed -i 's/^echo __________Lammps__________$/#&/' 20250926-Tde.job
                sed -i 's/^module load StdEnv\/2023  intel\/2023.2.1  openmpi\/4.1.5 lammps-omp\/20230802$/#&/' 20250926-Tde.job
                sed -i 's/^module load intel\/2019u3  intelmpi\/2019u3 lammps\/29Mar2019 # for Niagara$/#&/' 20250926-Tde.job
                sed -i 's/^echo __________Mtp__________$/#&/' 20250926-Tde.job
                sed -i 's/^module load StdEnv\/2023$/#&/' 20250926-Tde.job

            elif [[ " ${PotentialsMtp[*]} " =~ " ${P} " ]]; then
                sed -i 's/^echo __________Lammps__________$/#&/' 20250926-Tde.job
                sed -i 's/^module load StdEnv\/2023  intel\/2023.2.1  openmpi\/4.1.5 lammps-omp\/20230802$/#&/' 20250926-Tde.job
                sed -i 's/^module load intel\/2019u3  intelmpi\/2019u3 lammps\/29Mar2019 # for Niagara$/#&/' 20250926-Tde.job
                sed -i 's/^echo __________TabGap__________$/#&/' 20250926-Tde.job
                sed -i 's/^module load nixpkgs\/16.09  StdEnv\/2020  intel\/2020.1.217  openmpi\/4.0.3 lammps-omp\/20201029$/#&/' 20250926-Tde.job
                sed -i 's/^module load hdf5-mpi\/1.10.6 voro++\/0.4.6 eigen\/3.3.7$/#&/' 20250926-Tde.job
                sed -i 's/^module load mlip\/2.0 fftw-mpi\/3.3.8 tbb\/2020.2$/#&/' 20250926-Tde.job
            fi

            #sbatch 20250926-Tde.job
            cd ..
        done
        cd ..
    done
    cd ..
done