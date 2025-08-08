import numpy as np
from ase import Atoms
from ase.build import bulk
from ase.neighborlist import NeighborList, natural_cutoffs
from scipy.spatial import cKDTree
import re

def read_lammps_data(filename):
    with open(filename, "r") as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        if "atoms" in line:
            num_atoms = int(line.split()[0])
        if "atom types" in line:
            num_atom_types = int(line.split()[0])
        if "xlo xhi" in line:
            xlo, xhi = map(float, line.split()[:2])
        if "ylo yhi" in line:
            ylo, yhi = map(float, line.split()[:2])
        if "zlo zhi" in line:
            zlo, zhi = map(float, line.split()[:2])
        if "Atoms" in line:
            atoms_start = i + 2
            break
    atom_positions = []
    atom_types = []
    for line in lines[atoms_start : atoms_start + num_atoms]:
        atom_data = line.split()
        atom_type = 40 #int(atom_data[1])
        position = [float(x) for x in atom_data[2:5]]
        atom_positions.append(position)
        atom_types.append(atom_type)
    cell = np.array([[xhi - xlo, 0, 0], [0, yhi - ylo, 0], [0, 0, zhi - zlo]])
    positions = np.array(atom_positions)
    atoms = Atoms(numbers=atom_types, positions=positions, cell=cell, pbc=True)
    return atoms

#filename1 = input("Enter the address for the benchmark structure: ")
filename1 = "20250130-Tde-Thermalization.data"
filename2 = "20250130-Tde.data"
perfect_structure = read_lammps_data(filename1)
defective_structure = read_lammps_data(filename2)

def extract_pka_value(file_path, search_string="PKA ID is:"):
    extracted_value = None
    with open(file_path, 'r') as file:
        for line in file:
            match = re.search(rf"{re.escape(search_string)}\s*(\d+)", line)
            if match:
                extracted_value = int(match.group(1))
                break
    return extracted_value

def extract_variable_value(file_path, variable_name="No"):
    pattern = rf"variable\s+{re.escape(variable_name)}\s+equal\s+(\d+)"
    extracted_value = None
    with open(file_path, 'r') as file:
        for line in file:
            match = re.search(pattern, line)
            if match:
                extracted_value = int(match.group(1))
                break
    return extracted_value

file_path = "20250130-Tde.lammpslog"
# PKA_ID = extract_pka_value(file_path, search_string="ExperimentNo:")
PKA_ID = extract_variable_value(file_path, variable_name="No")

def read_initial_config_data(filename,PKA_ID):
    with open(filename, "r") as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        if "atoms" in line:
            num_atoms = int(line.split()[0])
        if "atom types" in line:
            num_atom_types = int(line.split()[0])
        if "xlo xhi" in line:
            xlo, xhi = map(float, line.split()[:2])
        if "ylo yhi" in line:
            ylo, yhi = map(float, line.split()[:2])
        if "zlo zhi" in line:
            zlo, zhi = map(float, line.split()[:2])
        if "Atoms" in line:
            atoms_start = i + 2
            break
    atom_positions = []
    atom_index = []
    for line in lines[atoms_start : atoms_start + num_atoms]:
        atom_data = line.split()
        atom_index.append(int(atom_data[0]))
        position = [float(x) for x in atom_data[2:5]]
        atom_positions.append(position)
    try:
        positions = np.array(atom_positions[atom_index.index(PKA_ID)])
    except ValueError:
        # print(f"⚠️ PKA_ID {PKA_ID} not found. Assigning a random ID from 1 to 5.")
        fallback_id = random.randint(1, 5)
        if fallback_id in atom_index:
            positions = np.array(atom_positions[atom_index.index(fallback_id)])
        else:
            # fallback_id also not found: pick first available atom
            # print(f"⚠️ Random fallback ID {fallback_id} not in atom list. Using first available atom.")
            positions = np.array(atom_positions[0])
    return positions

def read_end_config_data(filename,PKA_ID):
    with open(filename, "r") as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        if "atoms" in line:
            num_atoms = int(line.split()[0])
        if "atom types" in line:
            num_atom_types = int(line.split()[0])
        if "xlo xhi" in line:
            xlo, xhi = map(float, line.split()[:2])
        if "ylo yhi" in line:
            ylo, yhi = map(float, line.split()[:2])
        if "zlo zhi" in line:
            zlo, zhi = map(float, line.split()[:2])
        if "Atoms" in line:
            atoms_start = i + 2
            break
    atom_positions = []
    atom_index = []
    for line in lines[atoms_start : atoms_start + num_atoms]:
        atom_data = line.split()
        atom_index.append(int(atom_data[0]))
        position = [float(x) for x in atom_data[2:5]]
        atom_positions.append(position)
    try:
        positions = np.array(atom_positions[atom_index.index(PKA_ID)])
    except ValueError:
        # print(f"⚠️ PKA_ID {PKA_ID} not found. Assigning a random ID from 1 to 5.")
        fallback_id = random.randint(1, 5)
        if fallback_id in atom_index:
            positions = np.array(atom_positions[atom_index.index(fallback_id)])
        else:
            # fallback_id also not found: pick first available atom
            # print(f"⚠️ Random fallback ID {fallback_id} not in atom list. Using first available atom.")
            positions = np.array(atom_positions[0])
    return positions

filename3 = "20250130-Tde-Thermalization.data"
filename4 = "20250130-Tde.data"
atom_begin_position = read_initial_config_data(filename3,PKA_ID)
atom_end_position = read_end_config_data(filename4,PKA_ID)
dist=np.linalg.norm(atom_begin_position - atom_end_position)

def read_initial_config_data(filename,PKA_ID):
    with open(filename, "r") as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        if "atoms" in line:
            num_atoms = int(line.split()[0])
        if "atom types" in line:
            num_atom_types = int(line.split()[0])
        if "xlo xhi" in line:
            xlo, xhi = map(float, line.split()[:2])
        if "ylo yhi" in line:
            ylo, yhi = map(float, line.split()[:2])
        if "zlo zhi" in line:
            zlo, zhi = map(float, line.split()[:2])
        if "Atoms" in line:
            atoms_start = i + 2
            break
    atom_positions = []
    atom_index = []
    for line in lines[atoms_start : atoms_start + num_atoms]:
        atom_data = line.split()
        atom_index.append(int(atom_data[0]))
        position = [float(x) for x in atom_data[2:5]]
        atom_positions.append(position)
    positions = np.array(atom_positions[atom_index.index(PKA_ID)])
    return positions

def read_end_config_data(filename,PKA_ID):
    with open(filename, "r") as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        if "atoms" in line:
            num_atoms = int(line.split()[0])
        if "atom types" in line:
            num_atom_types = int(line.split()[0])
        if "xlo xhi" in line:
            xlo, xhi = map(float, line.split()[:2])
        if "ylo yhi" in line:
            ylo, yhi = map(float, line.split()[:2])
        if "zlo zhi" in line:
            zlo, zhi = map(float, line.split()[:2])
        if "Atoms" in line:
            atoms_start = i + 2
            break
    atom_positions = []
    atom_index = []
    for line in lines[atoms_start : atoms_start + num_atoms]:
        atom_data = line.split()
        atom_index.append(int(atom_data[0]))
        position = [float(x) for x in atom_data[2:5]]
        atom_positions.append(position)
    positions = np.array(atom_positions[atom_index.index(PKA_ID)])
    return positions

filename3 = "20250130-Tde-Thermalization.data"
filename4 = "20250130-Tde.data"
atom_begin_position = read_initial_config_data(filename3,PKA_ID)
atom_end_position = read_end_config_data(filename4,PKA_ID)

new_atoms = Atoms(
    symbols='Zr2',
    positions=[atom_begin_position,atom_end_position],
    cell=perfect_structure.cell.array,
    pbc=True)
dist = new_atoms.get_distance(0, 1, mic=True)

def get_minimal_distance(atoms):
    positions = atoms.positions
    tree = cKDTree(positions)
    nearest_distances, _ = tree.query(positions, k=2)  # k=2 to include self-distance
    min_distance = np.min(nearest_distances[:, 1])  # [:, 1] to exclude self-distance (0)
    return min_distance

def detect_vacancy(perfect_structure, defective_structure, dist):
    min_distance_perfect = get_minimal_distance(perfect_structure)
    min_distance_defect = get_minimal_distance(defective_structure)
    if min_distance_defect <min_distance_perfect-0.45: # * lattice_constant:
        return 1  
    return 0    

vacancy_result = detect_vacancy(perfect_structure, defective_structure, dist)
print(vacancy_result)