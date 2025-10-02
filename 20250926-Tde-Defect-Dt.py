import numpy as np
from ase import Atoms
from scipy.spatial import cKDTree
import random
import re
import os

def read_lammps_data(filename, PKA_ID=None, return_atoms=True):
    """Read LAMMPS data.
       If PKA_ID is given -> return only that atom's position (with fallback).
       If return_atoms=True -> return ASE Atoms object, else just positions list.
    """
    with open(filename, "r") as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        if "atoms" in line:
            num_atoms = int(line.split()[0])
        if "xlo xhi" in line:
            xlo, xhi = map(float, line.split()[:2])
        if "ylo yhi" in line:
            ylo, yhi = map(float, line.split()[:2])
        if "zlo zhi" in line:
            zlo, zhi = map(float, line.split()[:2])
        if "Atoms" in line:
            atoms_start = i + 2
            break

    atom_positions, atom_types, atom_index = [], [], []
    for line in lines[atoms_start : atoms_start + num_atoms]:
        atom_data = line.split()
        atom_index.append(int(atom_data[0]))
        atom_types.append(40)  # hardcoded type
        atom_positions.append([float(x) for x in atom_data[2:5]])

    if PKA_ID is not None:  # return only one atom position
        try:
            return np.array(atom_positions[atom_index.index(PKA_ID)])
        except ValueError:
            # fallback to random or first atom
            fallback_id = random.randint(1, 5)
            if fallback_id in atom_index:
                return np.array(atom_positions[atom_index.index(fallback_id)])
            return np.array(atom_positions[0])

    cell = np.array([[xhi - xlo, 0, 0], [0, yhi - ylo, 0], [0, 0, zhi - zlo]])
    positions = np.array(atom_positions)
    return Atoms(numbers=atom_types, positions=positions, cell=cell, pbc=True) if return_atoms else positions


def extract_pka_value(file_path, search_string="PKA ID is:"):
    """Extract PKA ID from LAMMPS log using regex."""
    with open(file_path, 'r') as file:
        for line in file:
            match = re.search(rf"{re.escape(search_string)}\s*(\d+)", line)
            if match:
                return int(match.group(1))
    return None


def get_minimal_distance(atoms):
    """Get minimal nearest-neighbor distance in structure."""
    positions = atoms.positions
    tree = cKDTree(positions)
    nearest_distances, _ = tree.query(positions, k=2)
    return np.min(nearest_distances[:, 1])


def detect_vacancy(perfect_structure, defective_structure):
    """Detect vacancy by comparing minimal neighbor distances."""
    cwd = os.getcwd()
    Temperature = os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(cwd)))))
    Potential = os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(cwd))))

    cutoff_dict = {
        "10": {"M2R": 0.45, "M3R": 0.30, "BMD192R": 0.45, "M2": 0.45, "M3": 0.45, "BMD192": 0.45, "TabGap1": 0.45, "TabGap2": 0.45, "MtpPd": 0.15,  "default": 0.45},
        "300": {"M2R": 0.30, "M3R": 0.15, "BMD192R": 0.30, "M2": 0.30, "M3": 0.30, "BMD192": 0.30, "TabGap1": 0.30, "TabGap2": 0.30, "MtpPd": 0.15,  "default": 0.30},
    }

    cutoff = cutoff_dict.get(Temperature, {}).get(Potential, cutoff_dict.get(Temperature, {}).get("default", 0.45))

    min_distance_perfect = get_minimal_distance(perfect_structure)
    min_distance_defect = get_minimal_distance(defective_structure)

    return 1 if min_distance_defect < min_distance_perfect - cutoff else 0


# ---------------- MAIN SCRIPT ----------------

filename1 = "Tde-Thermalization.data"
filename2 = "Tde.data"

perfect_structure = read_lammps_data(filename1)
defective_structure = read_lammps_data(filename2)

file_path = "Tde.lammpslog"
current_dir = os.getcwd()
path_parts = current_dir.split(os.sep)
if "BoxSize" in path_parts:
    SearchFor = "group		PKA id"
else:
    SearchFor = "EperimentNo:"
PKA_ID = extract_pka_value(file_path, search_string=SearchFor)

#print(PKA_ID)
# PKA positions before/after irradiation
atom_begin_position = read_lammps_data(filename1, PKA_ID=PKA_ID, return_atoms=False)
atom_end_position   = read_lammps_data(filename2, PKA_ID=PKA_ID, return_atoms=False)
#print(atom_begin_position)
#print(atom_end_position)

# displacement (with minimum image convention)
new_atoms = Atoms(
    symbols='Zr2',
    positions=[atom_begin_position, atom_end_position],
    cell=perfect_structure.cell.array,
    pbc=True)
dist = new_atoms.get_distance(0, 1, mic=True)

# vacancy detection
vacancy_result = detect_vacancy(perfect_structure, defective_structure)
print(vacancy_result)
