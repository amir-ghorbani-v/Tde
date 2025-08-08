# LAMMPS Zr Potential Workflow

This directory contains the scripts, job files, and templates for a multi-step LAMMPS workflow used to simulate and analyze Zr potentials.  
Each step in the workflow is stored in its own set of files, and the process is designed to run on both Windows (local) and Linux (HPC / bash) environments.

---

## Workflow Steps

### Step 1: **BoxSize**
- **Goal:** Optimize the simulation box size before running equilibrium calculations.
- **Files:**
  - `20250130-Equilibrium-BoxSize.job` — SLURM job submission file for HPC.
  - `20250130-Equilibrium-BoxSize.lammpstemp` — LAMMPS input template for box size optimization.
  - `20250130-Equilibrium-BoxSize.sh` — Shell script to run the job.

---

### Step 2: **Equilibrium**
- **Goal:** Run equilibrium simulations at the optimized box size.
- **Files:**
  - `20240228-Equilibrium.csv` — Output data from equilibrium simulations.
  - `20250130-Equilibrium.job` — SLURM job submission file.
  - `20250130-Equilibrium.lammpstemp` — LAMMPS input template for equilibrium.
  
---

### Step 3: **Thermalization**
- **Goal:** Prepare the system at the target temperature.  
*(Files for this step may be part of the Equilibrium stage or stored separately.)*

---

### Step 4: **Tde** (Threshold Displacement Energy)
- **Goal:** Compute threshold displacement energy and analyze defect creation using detailed molecular dynamics simulations.
- **Description:**  
  The Tde step uses LAMMPS to run atomic displacement simulations where a Primary Knock-on Atom (PKA) is given a velocity in a specified direction and magnitude. The system is evolved under an NVE ensemble with a thermalized environment, tracking kinetic, potential, and total energies as well as neighbor statistics for defect analysis.

- **Key features in the `20250130-Tde.lammpstemp` file:**
  - Initialization with metal units, periodic boundaries, and atomic style.
  - Reading pre-thermalized atomic configuration (`Thermalization.data`).
  - EAM potential specification for Zr.
  - Calculation of kinetic and potential energy per atom.
  - Identification of neighbors and minimum atomic distances.
  - Velocity assignment to a specific PKA atom using input variables `V1Temp`, `V2Temp`, `V3Temp`.
  - Running dynamics with the NVE ensemble and a custom timestep reset fix.
  - Thermo output customized to include energies, pressure, volume, neighbor count, and minimum interatomic distances.
  - Dumping atomic data and exporting simulation data files before and after the run.
  - Logging energy per atom both before and after the run for detailed energy balance analysis.

- **Files:**
  - `20250130-Tde-Defect-Dist.py` — Python script to analyze defect distribution from Tde runs.
  - `20250130-Tde-Defect-Ws.pytemp` — Wigner-Seitz defect analysis template.
  - `20250130-Tde.jobtemp` — Job submission template for Tde.
  - `20250130-Tde.lammpstemp` — LAMMPS input template detailed above.
  - `20250130-Tde.sh` — Main shell script for running Tde simulations.

---

### Step 5: **Minimization**
- **Goal:** Perform energy minimization for relaxed structures.
- **Files:**
  - `20250212-Minimization.job` — SLURM job submission file.
  - `20250212-Minimization.lammpstemp` — LAMMPS input template for minimization.

---

### Step 6: **Analysis**
- **Goal:** Post-process and extract relevant quantities from the simulations.
- **Files:**  
  *(Analysis scripts and notebooks may be stored separately.)*

---

## Example Potentials
An **`example`** folder is provided containing the potential **`M3`** as a reference for how to structure potential files in the workflow.

---

## Running on HPC (Linux)
1. Load required modules:
   ```bash
   module load gcc
   module load openmpi
   module load lammps
