# <editor-fold desc="######################################## Log">
print("######################################## LOG") # %%
import sys
PythonName = "20250610-Tde-Analysis"
print("######################################## ")
print("20250610-Tde-Analysis.py: START")
# sys.stdout = open(PythonName + '-Reparameterization.txt', 'w')
# </editor-fold>

# <editor-fold desc="######################################## OPTIONS"> # %%
print("######################################## OPTIONS") # %%

# <editor-fold desc="**********  Library">
print("**********  Library")
# import concurrent.futures.process
# import cv2
# import datetime
# import fileinput
from itertools import product
import mplstereonet
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker
import getpass
from orix.vector import Vector3d
from orix.quaternion import Orientation
from orix.quaternion.symmetry import (
    C1, Ci, C2, C2h, D2, D2h,
    C4, C4h, D4, D4h,
    C3, D3, D3d,
    C6, D6, D6h,
    T, Th, O, Td, Oh
)
from orix import plot
import orix
import re
from matplotlib.lines import Line2D
from matplotlib.ticker import FormatStrFormatter
# import glob
# import math
import matplotlib
import matplotlib.patches as mpatches
import matplotlib.patches as patches
import matplotlib as mpl
import matplotlib.pyplot as plt  # visualization
import matplotlib
import math
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import linregress
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
# import multiprocessing as mp
# import ntpath
import os  # I/O operations
# import pathlib
import platform
# import pickle
# import shutil
import sys
import seaborn as sns
import subprocess
# from scipy import odr
# from scipy.interpolate import CubicSpline
# from scipy.interpolate import PchipInterpolator as PchipFunc
# from scipy.interpolate import pchip_interpolate as PchipValue
# from scipy.interpolate import splrep, splev
# from scipy.interpolate import UnivariateSpline
# from scipy.interpolate import interp1d
# from scipy.misc import derivative
# from optimparallel import minimize_parallel
# from scipy.optimize import minimize
# import scipy.constants as const
# import time
# from pathlib import Path
import pandas as pd  # data process; python-version SQL
pd.options.mode.chained_assignment = None
# from ase.build import bulk
# from ase.calculators.eam import EAM


from sklearn.neighbors import KernelDensity
from orix.vector import Vector3d
from mpl_toolkits.axes_grid1 import make_axes_locatable


from orix.vector import Vector3d
from orix.quaternion import Orientation
from orix.quaternion.symmetry import Oh, Td, D4h, D3d, C4h, C1  # Correct import for your ORIX version
from orix import plot
# </editor-fold>

# <editor-fold desc="**********  Environment">
print("**********  Environment")
print(os.path.dirname(sys.executable))
print(getpass.getuser())
print(platform.system())

# SimulationEnvironment = "MyPc"
# SimulationEnvironment = "ComputeCanada"

if platform.system() == "Windows":
    SimulationEnvironment = "MyPc"

    CurrentDirectory = os.getcwd()
    # OriginalPotentialAddress = "D:/Queens_University/Project/Zr/PotentialBank/Eam/Mendelev/" + PotentialName + ".eampot"
    LammpsPythonAddress = "C:/Users/19ag40/AppData/Local/LAMMPS 64-bit 2Aug2023 with Python/Python/"
    LammpsTempAddress = "D:/Queens_University/Project/Zr/PotentialDevelopement/LammpsTemp"
    CollectionAddress = "D:/Queens_University/MME/Project/Zr/Tde/Run"
    GroupRun = False
    GroupPlot = False
    PlottingSingle = True
    RunLammpsEph = False
    MeasureRhoEos = True
    MeasureRhoQsd = False
    MeasureFDerEos = True
    MeasureFDerQsd = False
    # LammpsPythonAddress = "D:/Setups/LAMMPS/python"

elif platform.system() == "Linux":
    SimulationEnvironment = "ComputeCanada"

    CurrentDirectory = os.getcwd()
    # OriginalPotentialAddress = "/home/veshand/Zr/PotentialBank/Eam/Mendelev/" + PotentialName + ".eampot"
    LammpsPythonAddress = "/home/veshand/.local/easybuild/software/2020/avx512/MPI/intel2020/openmpi4/lammps-eph/20220623"
    LammpsTempAddress = "/home/veshand/Zr/PotentialDevelopement/LammpsTemp"

    GroupRun = True
    GroupPlot = False
    PlottingSingle = False
    RunLammpsEph = True
    MeasureRhoEos = False
    MeasureRhoQsd = False
    MeasureFDerEos = False
    MeasureFDerQsd = False
# </editor-fold>

# <editor-fold desc="**********  Function">
print("**********  Function")


# <editor-fold desc="LatexFriendlyFunc">
def LatexFriendlyFunc(direction):
    # Remove 'dir_' prefix
    direction = direction.replace('dir_', '')

    # Split the direction into individual components
    components = direction.split('_')

    # For each component, wrap it in \text{}, using \bar{\text{...}} for negatives
    new_components = []
    for comp in components:
        if comp.startswith('-'):
            # Remove the '-' sign and use \bar
            new_components.append(r"\bar{\text{" + comp[1:] + r"}}")
        else:
            new_components.append(r"\text{" + comp + r"}")

    # Join the components with spaces
    transformed = " ".join(new_components)

    # Wrap the entire label in math mode with square brackets
    return f"$[{transformed}]$"


# </editor-fold>

def miller_to_bravais(h, k, l):
    i = -h - k
    return [h, k, i, l]


# Function to format LaTeX string
def format_bravais_label(vec):
    hkil = miller_to_bravais(*vec)
    return r'$[' + ''.join(
        f'\\overline{{{abs(v)}}}' if v < 0 else str(v)
        for v in hkil
    ) + ']$'

def extract_info(filename):
    data = {}
    with open(filename, 'r') as file:
        lines = file.readlines()

        current_Try = None
        for line in lines:
            line = line.strip()
            if line.startswith("Try:"):
                current_Try = line.split(" ")[-1]
                data[current_Try] = {}
            elif line.startswith("dir_"):
                dir_info, value = line.split(":")
                data[current_Try][dir_info] = int(value)

    return data

def new_column_name(old_name):
    if not old_name.startswith("dir"):
        return old_name
    parts = old_name.split('_')[1:]
    new_parts = [str(part) if int(part) >= 0 else r'$\overline{' + str(abs(int(part))) + '}$' for part in parts]
    return ''.join(new_parts)

def DirectionFunc(Direction_V):
    Parts = Direction_V.replace('dir_', '').split('_')
    return np.array([int(x) for x in Parts], dtype=float)
# </editor-fold>

# <editor-fold desc="**********  Variables">
print("**********  Variables")
font = {"family": "serif",  # serif" "sans-serif" "cursive" "fantasy" "monospace"
        "weight": "bold",
        "size": 10}
plt.rcParams['font.family'] = 'DeJavu Serif'
plt.rcParams['font.serif'] = ['Times New Roman']
matplotlib.rc("font", **{'family': 'serif'})
matplotlib.rcParams["lines.linewidth"] = 3
matplotlib.rcParams["figure.figsize"] = (5, 5)
matplotlib.rcParams["axes.spines.right"] = True
matplotlib.rcParams["axes.spines.top"] = True
plt.rc('font', size=15)          # controls default text size
plt.rc('axes', titlesize=18)     # fontsize of the title
plt.rc('axes', labelsize=18)     # fontsize of the x and y labels
plt.rc('xtick', labelsize=18)    # fontsize of the x-axis ticks
plt.rc('ytick', labelsize=18)    # fontsize of the y-axis ticks
plt.rc('legend', fontsize=15)    # fontsize of the legend
colors = [
    "#1f77b4",  # Blue
    "#8c564b",  # Brown
    "#d62728",  # Red
    "#9467bd",  # Purple
    "#2ca02c",  # Green
    "#e377c2",  # Pink
    "#7f7f7f",  # Gray
    "#bcbd22",  # Olive
    "#17becf",  # Cyan
    "#ff7f0e",  # Orange
    "#1a55e5",  # Dark Blue
    "#ffbb78",  # Light Orange
    "#98df8a",  # Light Green
    "#ff9896",  # Light Red
    "#c5b0d5",  # Light Purple
    "#c49c94",  # Light Brown
    "#f7b6d2",  # Light Pink
    "#c7c7c7",  # Light Gray
    "#dbdb8d",  # Light Olive
    "#9edae5",  # Light Cyan
]
Colors = ["b", "cyan", "darkred", "red", "green", "lime", "w", "#FF5733", "#33FF57", "#5733FF", "#FF5733", "#33FF57",
          "#5733FF", "#FF5733", "#33FF57", "#5733FF", "#FF5733", "#33FF57", "#5733FF"]
Colors = ["b", "cyan", "white", "darkred", "red", "white", "green", "lime", "w", "#FF5733", "#33FF57", "#5733FF",
          "#FF5733", "#33FF57", "#5733FF", "#FF5733", "#33FF57", "#5733FF", "#FF5733", "#33FF57", "#5733FF"]
Colors = ["b", "red", "green", "purple", "darkgoldenrod"]

Mapping = {
    'dir_1_2_1': 'dir_1_2_-3_1',
    'dir_0_1_0': 'dir_0_1_-1_0',
    'dir_0_1_1': 'dir_0_1_-1_1',
    'dir_0_3_1': 'dir_0_3_-3_1',
    'dir_-1_1_2': 'dir_-1_1_0_2',
    'dir_-1_1_0': 'dir_-1_1_0_0',
    'dir_-1_1_1': 'dir_-1_1_0_1',
    'dir_1_2_2': 'dir_1_2_-3_2',
    'dir_0_0_1': 'dir_0_0_0_1',
    'dir_1_2_0': 'dir_1_2_-3_0',

    'dir_-2_3_3': 'dir_-2_3_-1_3',
    'dir_-1_4_1': 'dir_-1_4_-3_1',
    'dir_-1_3_1': 'dir_-1_3_-2_1',
    'dir_-1_2_2': 'dir_-1_2_-1_2',
    'dir_0_2_1': 'dir_0_2_-2_1',
    'dir_-1_2_0': 'dir_-1_2_-1_0',
    'dir_-1_2_3': 'dir_-1_2_-1_3',
    'dir_-4_4_1': 'dir_-4_4_0_1',
    'dir_-3_4_1': 'dir_-3_4_-1_1',
    'dir_-2_4_1': 'dir_-2_4_-2_1',
    'dir_-1_4_0': 'dir_-1_4_-3_0',
    'dir_-1_2_1': 'dir_-1_2_-1_1',
    'dir_0_2_3': 'dir_0_2_-2_3',
    'dir_-1_4_2': 'dir_-1_4_-3_2',
    'dir_-2_2_3': 'dir_-2_2_0_3',
    'dir_-3_4_3': 'dir_-3_4_-1_3',
    'dir_-2_3_2': 'dir_-2_3_-1_2',
    'dir_-2_2_1': 'dir_-2_2_0_1',
    'dir_-1_3_0': 'dir_-1_3_-2_0',
    'dir_-1_3_2': 'dir_-1_3_-2_2',
    'dir_-3_4_2': 'dir_-3_4_-1_2',
    'dir_-1_4_3': 'dir_-1_4_-3_3',
    'dir_0_4_1': 'dir_0_4_-4_1',
    'dir_-3_4_0': 'dir_-3_4_-1_0',
    'dir_-2_3_0': 'dir_-2_3_-1_0',
    'dir_0_1_2': 'dir_0_1_-1_2',
    'dir_-2_3_1': 'dir_-2_3_-1_1',
    'dir_-1_3_3': 'dir_-1_3_-2_3',

    'dir_4_1_3': 'dir_4_1_-5_3',
    'dir_1_0_0': 'dir_1_0_-1_0',
    'dir_2_0_1': 'dir_2_0_-2_1',

}

DirectionList = list(Mapping.keys())
# </editor-fold>

# </editor-fold>

# <editor-fold desc="######################################## Read">
print("######################################## Read")# %%

# <editor-fold desc="***** Defect">
print("***** Defect")

# <editor-fold desc="^^^ Dataframe">
print("^^^ Dataframe")
Active=False
if Active:

    # <editor-fold desc="Deleting">
    DeleteList = [
        '20431216-Nonstop.csv',
        '20431216-ErrorDf.csv',
        # '20431216-Df.csv',
        '20431216-BlankDf.csv',
        'exec_custom-Blank-One.sh',
        'exec_custom-Nonstop-One.sh'
        'exec_custom-Error-All.sh'
    ]
    Active = False
    if Active:
        for File in DeleteList:
            # Check if the file exists
            if os.path.exists(File):
                try:
                    # Attempt to delete the file
                    os.remove(File)
                    print(f"Deleted: {File}")
                except Exception as e:
                    # Print an error message if there's an issue
                    print(f"Error deleting {File}: {e}")
            else:
                print(f"File does not exist: {File}")
    # </editor-fold>

    # <editor-fold desc="Renaming">
    print("Renaming")
    Renaming = False
    if Renaming:
        Filename = "Tde-Defect.csv"
        TemperatureList = [10,300]
        PotentialList = ["M3R", "M2R", "BMD192R", "M2", "M3", "BMD192", "TabGap1", "TabGap2", "MtpPd"]
        TryList = range(1, 41)
        DirectionList = [
            "dir_0_0_1",
            "dir_0_1_1",
            "dir_0_3_1",
            "dir_0_1_0",
            "dir_-1_1_0",
            "dir_-1_1_1",
            "dir_-1_1_2",
            "dir_1_2_0",
            "dir_1_2_1",
            "dir_1_2_2",

            "dir_-2_3_3",
            "dir_-1_4_1",
            "dir_-1_3_1",
            "dir_-1_2_2",
            "dir_0_2_1",
            "dir_-1_2_0",
            "dir_-1_2_3",
            "dir_-4_4_1",
            "dir_-3_4_1",
            "dir_-2_4_1",
            "dir_-1_4_0",
            "dir_-1_2_1",
            "dir_0_2_3",
            "dir_-1_4_2",
            "dir_-2_2_3",
            "dir_-3_4_3",
            "dir_-2_3_2",
            "dir_-2_2_1",
            "dir_-1_3_0",
            "dir_-1_3_2",
            "dir_-3_4_2",
            "dir_-1_4_3",
            "dir_0_4_1",
            "dir_-3_4_0",
            "dir_-2_3_0",
            "dir_0_1_2",
            "dir_-2_3_1",
            "dir_-1_3_3",

            "dir_1_0_0",
            "dir_2_0_1",
            "dir_4_1_3",
        ]
        if SimulationEnvironment == "MyPc":
            for Temperature, Potential, Try, Direction in product(TemperatureList, PotentialList, TryList, DirectionList):
                BaseDirectory = os.path.join(CollectionAddress, str(Temperature), Potential, str(Try), Direction)
                if not os.path.isdir(BaseDirectory):
                    continue
                for Item in os.scandir(BaseDirectory):
                    if not Item.is_dir() or not Item.name.isdigit():
                        continue
                    Energy = Item.name

                    Path = os.path.join(Item.path, Filename)
                    if os.path.exists(Path):
                        continue

                    Path20250926 = os.path.join(Item.path, "20250926-Tde-Defect.csv")
                    Path20250130 = os.path.join(Item.path, "20250130-Tde-Defect.csv")

                    if any(os.path.exists(p) for p in (Path20250130, Path20250926)):
                        print(Temperature, Potential, Try, Direction, Energy)

                        if os.path.exists(Path20250926):
                            os.rename(Path20250926, Path)
                            if os.path.exists(Path20250130):
                                os.remove(Path20250130)
                        elif os.path.exists(Path20250130):
                            os.rename(Path20250130, Path)
    # </editor-fold>

    # <editor-fold desc="Extraction">
    print("Extraction")
    Method ="Collection"
    if Method == "Log":
        # <editor-fold desc="Reading">
        filename = "max_value_summary.txt"
        # DirectionList = ["dir_-1_1_0", "dir_-1_1_1", "dir_-1_1_2", "dir_0_0_1", "dir_0_1_0", "dir_0_1_1", "dir_0_3_1", "dir_1_2_0", "dir_1_2_1", "dir_1_2_2"]
        # DirectionList = ['dir_0_1_-1_1', 'dir_-1_2_-1_0', 'dir_-1_2_-1_3', 'dir_-1_2_-1_1', 'dir_-1_1_0_2', 'dir_-1_1_0_0', 'dir_-1_1_0_1', 'dir_0_1_-1_2', 'dir_0_0_0_1', 'dir_0_1_-1_0']
        # print(DirectionList)
        Df = pd.DataFrame(columns=DirectionList)

        if SimulationEnvironment == "MyPc":
            for Root, Dirs, Files in os.walk(CurrentDirectory, topdown=False):  # open the files
                # print("root is: " + str(Root))
                # print("files are: " + str(Files))
                # print("dirs is: " + str(Dirs))
                for PythonName in Files:
                    # print("File PythonName is: " + str(PythonName))
                    if filename in PythonName:
                        Address = os.path.join(Root, PythonName)
                        # print(Address)
                        FolderName = str(Root.split("\\")[-1])
                        # FolderName = str(FolderName.split("/")[-1])
                        print(FolderName)
                        # os.system("pause")
                        Read = extract_info(Address)
                        DfNew = pd.DataFrame(Read).T.reset_index()
                        DfNew = DfNew.rename(columns={'index': 'Try'})
                        # print(DfNew)
                        # os.system("pause")

                        DfNew["Potential"] = FolderName
                        Df = pd.concat([Df, DfNew], ignore_index=True, axis=0)

        elif SimulationEnvironment == "ComputeCanada":
            Address = CurrentDirectory + "/" + filename
            # FolderName = str(Address.split("\\")[-1])
            FolderName = str(Address.split("/")[-2])
            print(FolderName)
            # os.system("pause")
            Read = extract_info(Address)
            DfNew = pd.DataFrame(Read).T.reset_index()
            DfNew = DfNew.rename(columns={'index': 'Try'})
            # print(DfNew)
            # os.system("pause")

            DfNew["Potential"] = FolderName
            Df = pd.concat([Df, DfNew], ignore_index=True, axis=0)

        for Direction in DirectionList:
            Df[Direction] = Df[Direction].astype('int')
        # Df["Try"] = Df["Try"].astype('int')

        Df.to_csv(PythonName + "-Df.csv", sep=',', index=False)
        # </editor-fold>
        # <editor-fold desc="Separation">
        Tde = Df[(Df["File"] == "Tde") & (Df["LookupWord"] == "EnergyPerAtomAfter")]
        DfEquilibrium = Df[Df["File"] == "Equilibrium"]
        DfThermalization = Df[Df["File"] == "Thermalization"]
        DfMinimization = Df[Df["File"] == "Minimization"]
        # </editor-fold>
    elif Method == "Out":
        Tde = pd.read_csv("20250506-Collection-Output.csv")
        Thermal = pd.read_csv("20250506-Collection-Thermal.csv")
    elif Method == "Csv-All":
        # <editor-fold desc="Reading">
        Filename = "Tde-Defect.csv"
        TemperatureList = [10,300]
        PotentialList = ["M3R", "M2R", "BMD192R", "M2", "M3", "BMD192", "TabGap1", "TabGap2", "MtpPd"]
        DirectionList = [
            "dir_0_0_1",
            "dir_0_1_1",
            "dir_0_3_1",
            "dir_0_1_0",
            "dir_-1_1_0",
            "dir_-1_1_1",
            "dir_-1_1_2",
            "dir_1_2_0",
            "dir_1_2_1",
            "dir_1_2_2",

            "dir_-2_3_3",
            "dir_-1_4_1",
            "dir_-1_3_1",
            "dir_-1_2_2",
            "dir_0_2_1",
            "dir_-1_2_0",
            "dir_-1_2_3",
            "dir_-4_4_1",
            "dir_-3_4_1",
            "dir_-2_4_1",
            "dir_-1_4_0",
            "dir_-1_2_1",
            "dir_0_2_3",
            "dir_-1_4_2",
            "dir_-2_2_3",
            "dir_-3_4_3",
            "dir_-2_3_2",
            "dir_-2_2_1",
            "dir_-1_3_0",
            "dir_-1_3_2",
            "dir_-3_4_2",
            "dir_-1_4_3",
            "dir_0_4_1",
            "dir_-3_4_0",
            "dir_-2_3_0",
            "dir_0_1_2",
            "dir_-2_3_1",
            "dir_-1_3_3",

            "dir_1_0_0",
            "dir_2_0_1",
            "dir_4_1_3",
        ]

        Header = ["Temperature","Potential","Try","Direction","DefectDt","DefectWs"]
        Tde = pd.DataFrame(columns=Header)

        if SimulationEnvironment == "MyPc":
            for Temperature in TemperatureList:
                # print(Temperature)
                for Potential in PotentialList:
                    # print(Potential)
                    for Try in np.arange(1, 41):
                        # print(Try)
                        for Direction in DirectionList:
                            # print(Direction)
                            for Energy in np.arange(0, 200):
                                # print(Energy)
                                Path = os.path.join(CollectionAddress, str(Temperature), Potential, str(Try), Direction, str(Energy), Filename)
                                # print(Path)

                                # print(Path)
                                if os.path.exists(Path):
                                    try:
                                        print(Temperature, Potential, Try, Direction, Energy)
                                        New = pd.read_csv(Path)
                                        New["Temperature"] = Temperature
                                        New["Potential"] = Potential
                                        New["Try"] = Try
                                        New["Direction"] = Direction
                                        New["Energy"] = Energy
                                        Tde = pd.concat([Tde, New], ignore_index=True, axis=0)
                                    except Exception as e:
                                        # print(Temperature, Potential, Try, Direction, Energy, "Not Found")
                                        continue



        Tde.to_csv(PythonName + "-Collection-Tde.csv", sep=',', index=False)
        # </editor-fold>
    elif Method == "Csv-Single":
        # <editor-fold desc="Reading">
        Filename = "20250130-Tde.csv"
        PotentialList = ["M3R"]  # "M2R", "M3R", "BMD192R", "M2", "M3", "BMD192", "TabGap1", "TabGap2"]
        DirectionList = [
            # "dir_0_1_2",
            # "dir_-1_2_1",
            # "dir_-1_3_2",
            # "dir_-1_4_2",
            # "dir_-2_3_3",
            # "dir_-4_4_1",
            # "dir_-2_3_0",
            # "dir_-1_2_3",
            # "dir_-3_4_0",
            "dir_-1_2_0",
            # "dir_-3_4_3",
            # "dir_-1_4_3",
            # "dir_-1_3_3",
            # "dir_-1_2_2",
            # "dir_0_2_1",
            # "dir_-1_4_0",
            # "dir_-2_3_1",
            # "dir_-2_3_2",
            # "dir_-2_2_3",
            # "dir_0_2_3",
            # "dir_-2_4_1",
            # "dir_-2_2_1",
            # "dir_-3_4_1",
            # "dir_-1_3_1",
            # "dir_-1_4_1",
            # "dir_-3_4_2",
            # "dir_-1_3_0",
            # "dir_0_4_1"
        ]

        Header = ["Temperature", "Potential", "Try", "Direction", "DefectDt", "DefectWs"]
        Tde = pd.DataFrame(columns=Header)

        if SimulationEnvironment == "MyPc":
            for Temperature in TemperatureList:
                # print(Temperature)
                for Potential in PotentialList:
                    # print(Potential)
                    for Try in np.arange(1, 1):
                        # print(Try)
                        for Direction in DirectionList:
                            # print(Direction)
                            for Energy in np.arange(0, 200):
                                # print(Energy)
                                Path = os.path.join(CurrentDirectory, Potential, str(Try), Direction, str(Energy),
                                                    Filename)
                                # print(Path)
                                if os.path.exists(Path):
                                    try:
                                        print(Potential, Try, Direction)
                                        New = pd.read_csv(Path)
                                        New["Temperature"] = Temperature
                                        New["Potential"] = Potential
                                        New["Try"] = Try
                                        New["Direction"] = Direction
                                        New["Energy"] = Energy
                                    except Exception as e:
                                        continue

                                Tde = pd.concat([Tde, New], ignore_index=True, axis=0)

        elif SimulationEnvironment == "ComputeCanada":
            Address = CurrentDirectory + "/" + filename
            # FolderName = str(Address.split("\\")[-1])
            FolderName = str(Address.split("/")[-2])
            print(FolderName)
            # os.system("pause")
            Read = extract_info(Address)
            DfNew = pd.DataFrame(Read).T.reset_index()
            DfNew = DfNew.rename(columns={'index': 'Try'})
            # print(DfNew)
            # os.system("pause")

            DfNew["Potential"] = FolderName
            Df = pd.concat([Df, DfNew], ignore_index=True, axis=0)

        Tde.to_csv(PythonName + "-Tde" + ".csv", sep=',', index=False)
        # </editor-fold>
    elif Method == "Collection":
        Tde = pd.read_csv("D:/Queens_University/MME/Project/Zr/Tde/Run/20250810-Collection-Defect.csv")
    else:
        Tde = pd.read_csv(PythonName + "-Tde" + ".csv")

    # </editor-fold>

    # <editor-fold desc="Miller to Miller-Bravais">
    print("Miller to Miller-Bravais")
    # Tde = Tde[Tde["File"] == "Tde"].copy()
    Tde['DirectionMb'] = Tde['Direction'].replace(Mapping)
    # print(Tde)
    # </editor-fold>

    # <editor-fold desc="Reparameterized">
    print("Reparameterized")
    Tde['Reparameterization'] = Tde['Potential'].apply(lambda x: 'R' if x.endswith('R') else 'O')
    # print(DfBravaisMelt)
    # </editor-fold>

    # <editor-fold desc="Polar">
    print("Polar")

    Tde[['X', 'Y', 'Z']] = Tde['Direction'].apply(DirectionFunc).apply(pd.Series)

    Tde['R'] = np.sqrt(Tde['X'] ** 2 + Tde['Y'] ** 2 + Tde['Z'] ** 2)
    Tde['UnitX'] = Tde['X'] / Tde['R']
    Tde['UnitY'] = Tde['Y'] / Tde['R']
    Tde['UnitZ'] = Tde['Z'] / Tde['R']

    Tde['Theta'] = np.degrees(np.arctan2(Tde['UnitY'], Tde['UnitX'])) % 360
    Tde['Phi'] = np.degrees(np.arccos(Tde['UnitZ']))

    Tde[['a1', 'a2', 'a3', 'h']] = Tde['DirectionMb'].apply(DirectionFunc).apply(pd.Series)
    # </editor-fold>

    # <editor-fold desc="Type Division">
    print("Type Division")
    # print(Tde)
    TdeDt = Tde.drop(columns=['DefectWs'])
    TdeWs = Tde.drop(columns=['DefectDt'])
    # </editor-fold>

    # <editor-fold desc="True Defect">
    print("True Defect")
    # print(Tde)
    TdeDt = Tde[Tde["DefectDt"] == 1]
    TdeWs = Tde[Tde["DefectWs"] == 1]
    # </editor-fold>

    # <editor-fold desc="Group">
    print("Group")
    # print(Tde)
    TdeGroupedPotential = Tde.groupby(['Potential'])
    TdeGroupedPotentialTemperature = Tde.groupby(['Temperature', 'Potential'])
    TdeGroupedPotentialDirection = Tde.groupby(['Potential', 'Direction', 'DirectionMb', 'Reparameterization', 'X', 'Y', 'Z', 'R', 'UnitX', 'UnitY', 'UnitZ', 'Theta', 'Phi', 'a1', 'a2', 'a3', 'h'])
    TdeGroupedPotentialDirectionTemperature = Tde.groupby(['Temperature', 'Potential', 'Direction', 'DirectionMb', 'Reparameterization', 'X', 'Y', 'Z', 'R', 'UnitX', 'UnitY', 'UnitZ', 'Theta', 'Phi', 'a1', 'a2', 'a3', 'h'])

    TdeDtGroupedPotential = TdeDt.groupby(['Potential'])
    TdeDtGroupedPotentialTemperature = TdeDt.groupby(['Temperature', 'Potential'])
    TdeDtGroupedPotentialDirection = TdeDt.groupby(['Potential', 'Direction', 'DirectionMb', 'Reparameterization', 'X', 'Y', 'Z', 'R', 'UnitX', 'UnitY', 'UnitZ', 'Theta', 'Phi', 'a1', 'a2', 'a3', 'h'])
    TdeDtGroupedPotentialDirectionTemperature = TdeDt.groupby(['Temperature', 'Potential', 'Direction', 'DirectionMb', 'Reparameterization', 'X', 'Y', 'Z', 'R', 'UnitX', 'UnitY', 'UnitZ', 'Theta', 'Phi', 'a1', 'a2', 'a3', 'h'])

    TdeWsGroupedPotential = TdeWs.groupby(['Potential'])
    TdeWsGroupedPotentialTemperature = TdeWs.groupby(['Temperature', 'Potential'])
    TdeWsGroupedPotentialDirection = TdeWs.groupby(['Potential', 'Direction', 'DirectionMb', 'Reparameterization', 'X', 'Y', 'Z', 'R', 'UnitX', 'UnitY', 'UnitZ', 'Theta', 'Phi', 'a1', 'a2', 'a3', 'h'])
    TdeWsGroupedPotentialDirectionTemperature = TdeWs.groupby(['Temperature', 'Potential', 'Direction', 'DirectionMb', 'Reparameterization', 'X', 'Y', 'Z', 'R', 'UnitX', 'UnitY', 'UnitZ', 'Theta', 'Phi', 'a1', 'a2', 'a3', 'h'])

    # for group in TdeGroupedPotentialDirection:
    #     print(group)
    # print(Tde)
    # os.system("pause")
    # </editor-fold>

    # <editor-fold desc="Describe">
    print("Describe")
    TdeGroupedPotentialStd = TdeGroupedPotential.std(numeric_only=True).reset_index()
    TdeGroupedPotentialDescribe = TdeGroupedPotential["Energy"].describe().reset_index()
    TdeDtGroupedPotentialStd = TdeDtGroupedPotential.std(numeric_only=True).reset_index()
    TdeDtGroupedPotentialDescribe = TdeDtGroupedPotential["Energy"].describe().reset_index()
    TdeWsGroupedPotentialStd = TdeWsGroupedPotential.std(numeric_only=True).reset_index()
    TdeWsGroupedPotentialDescribe = TdeWsGroupedPotential["Energy"].describe().reset_index()

    TdeGroupedPotentialTemperatureStd = TdeGroupedPotentialTemperature.std(numeric_only=True).reset_index()
    TdeGroupedPotentialTemperatureDescribe = TdeGroupedPotentialTemperature["Energy"].describe().reset_index()
    TdeDtGroupedPotentialTemperatureStd = TdeDtGroupedPotentialTemperature.std(numeric_only=True).reset_index()
    TdeDtGroupedPotentialTemperatureDescribe = TdeDtGroupedPotentialTemperature["Energy"].describe().reset_index()
    TdeWsGroupedPotentialTemperatureStd = TdeWsGroupedPotentialTemperature.std(numeric_only=True).reset_index()
    TdeWsGroupedPotentialTemperatureDescribe = TdeWsGroupedPotentialTemperature["Energy"].describe().reset_index()

    TdeGroupedPotentialDirectionStd = TdeGroupedPotentialDirection.std(numeric_only=True).reset_index()
    TdeGroupedPotentialDirectionDescribe = TdeGroupedPotentialDirection["Energy"].describe().reset_index()
    TdeDtGroupedPotentialDirectionStd = TdeDtGroupedPotentialDirection.std(numeric_only=True).reset_index()
    TdeDtGroupedPotentialDirectionDescribe = TdeDtGroupedPotentialDirection["Energy"].describe().reset_index()
    TdeWsGroupedPotentialDirectionStd = TdeWsGroupedPotentialDirection.std(numeric_only=True).reset_index()
    TdeWsGroupedPotentialDirectionDescribe = TdeWsGroupedPotentialDirection["Energy"].describe().reset_index()

    TdeGroupedPotentialDirectionTemperatureStd = TdeGroupedPotentialDirectionTemperature.std(numeric_only=True).reset_index()
    TdeGroupedPotentialDirectionTemperatureDescribe = TdeGroupedPotentialDirectionTemperature["Energy"].describe().reset_index()
    TdeDtGroupedPotentialDirectionTemperatureStd = TdeDtGroupedPotentialDirectionTemperature.std(numeric_only=True).reset_index()
    TdeDtGroupedPotentialDirectionTemperatureDescribe = TdeDtGroupedPotentialDirectionTemperature["Energy"].describe().reset_index()
    TdeWsGroupedPotentialDirectionTemperatureStd = TdeWsGroupedPotentialDirectionTemperature.std(numeric_only=True).reset_index()
    TdeWsGroupedPotentialDirectionTemperatureDescribe = TdeWsGroupedPotentialDirectionTemperature["Energy"].describe().reset_index()
    # </editor-fold>

    # <editor-fold desc="Export">
    Tde.to_csv(PythonName + "-Tde" + ".csv", index=False)
    TdeGroupedPotentialDescribe.to_csv(PythonName + "-" + "TdeGroupedPotentialDescribe.csv", sep=',', index=False)
    TdeGroupedPotentialTemperatureDescribe.to_csv(PythonName + "-" + "TdeGroupedPotentialTemperatureDescribe.csv", sep=',', index=False)
    TdeGroupedPotentialDirectionDescribe.to_csv(PythonName + "-" + "TdeGroupedPotentialDirectionDescribe.csv", sep=',', index=False)
    TdeGroupedPotentialDirectionTemperatureDescribe.to_csv(PythonName + "-" + "TdeGroupedPotentialDirectionTemperatureDescribe.csv", sep=',', index=False)

    TdeDt.to_csv(PythonName + "-TdeDt" + ".csv", index=False)
    TdeDtGroupedPotentialDescribe.to_csv(PythonName + "-" + "TdeDtGroupedPotentialDescribe.csv", sep=',', index=False)
    TdeDtGroupedPotentialTemperatureDescribe.to_csv(PythonName + "-" + "TdeDtGroupedPotentialTemperatureDescribe.csv", sep=',', index=False)
    TdeDtGroupedPotentialDirectionDescribe.to_csv(PythonName + "-" + "TdeDtGroupedPotentialDirectionDescribe.csv", sep=',', index=False)
    TdeDtGroupedPotentialDirectionTemperatureDescribe.to_csv(PythonName + "-" + "TdeDtGroupedPotentialDirectionTemperatureDescribe.csv", sep=',', index=False)

    TdeWs.to_csv(PythonName + "-TdeWs" + ".csv", index=False)
    TdeWsGroupedPotentialDescribe.to_csv(PythonName + "-" + "TdeWsGroupedPotentialDescribe.csv", sep=',', index=False)
    TdeWsGroupedPotentialTemperatureDescribe.to_csv(PythonName + "-" + "TdeWsGroupedPotentialTemperatureDescribe.csv", sep=',', index=False)
    TdeWsGroupedPotentialDirectionDescribe.to_csv(PythonName + "-" + "TdeWsGroupedPotentialDirectionDescribe.csv", sep=',', index=False)
    TdeWsGroupedPotentialDirectionTemperatureDescribe.to_csv(PythonName + "-" + "TdeWsGroupedPotentialDirectionTemperatureDescribe.csv", sep=',', index=False)
    # </editor-fold>

else:
    Tde = pd.read_csv(PythonName + "-Tde" + ".csv")
    TdeGroupedPotentialDescribe = pd.read_csv(PythonName + "-" + "TdeGroupedPotentialDescribe.csv")
    TdeGroupedPotentialTemperatureDescribe = pd.read_csv(PythonName + "-" + "TdeGroupedPotentialTemperatureDescribe.csv")
    TdeGroupedPotentialDirectionDescribe = pd.read_csv(PythonName + "-" + "TdeGroupedPotentialDirectionDescribe.csv")
    TdeGroupedPotentialDirectionTemperatureDescribe = pd.read_csv(PythonName + "-" + "TdeGroupedPotentialDirectionTemperatureDescribe.csv")

    TdeDt = pd.read_csv(PythonName + "-TdeDt" + ".csv")
    TdeDtGroupedPotentialDescribe = pd.read_csv(PythonName + "-" + "TdeDtGroupedPotentialDescribe.csv")
    TdeDtGroupedPotentialTemperatureDescribe = pd.read_csv(PythonName + "-" + "TdeDtGroupedPotentialTemperatureDescribe.csv")
    TdeDtGroupedPotentialDirectionDescribe = pd.read_csv(PythonName + "-" + "TdeDtGroupedPotentialDirectionDescribe.csv")
    TdeDtGroupedPotentialDirectionTemperatureDescribe = pd.read_csv(PythonName + "-" + "TdeDtGroupedPotentialDirectionTemperatureDescribe.csv")

    TdeWs = pd.read_csv(PythonName + "-TdeWs" + ".csv")
    TdeWsGroupedPotentialDescribe = pd.read_csv(PythonName + "-" + "TdeWsGroupedPotentialDescribe.csv")
    TdeWsGroupedPotentialTemperatureDescribe = pd.read_csv(PythonName + "-" + "TdeWsGroupedPotentialTemperatureDescribe.csv")
    TdeWsGroupedPotentialDirectionDescribe = pd.read_csv(PythonName + "-" + "TdeWsGroupedPotentialDirectionDescribe.csv")
    TdeWsGroupedPotentialDirectionTemperatureDescribe = pd.read_csv(PythonName + "-" + "TdeWsGroupedPotentialDirectionTemperatureDescribe.csv")


# </editor-fold>

# <editor-fold desc="^^^ Group">
print("^^^ Group")
# print(Tde)
TdeGroupedPotential = Tde.groupby(['Potential'])
TdeGroupedPotentialTemperature = Tde.groupby(['Temperature', 'Potential'])
TdeGroupedPotentialDirection = Tde.groupby(['Potential', 'Direction', 'DirectionMb', 'Reparameterization', 'X', 'Y', 'Z', 'R', 'UnitX', 'UnitY', 'UnitZ','Theta', 'Phi', 'a1', 'a2', 'a3', 'h'])
TdeGroupedPotentialDirectionTemperature = Tde.groupby(['Temperature', 'Potential', 'Direction', 'DirectionMb', 'Reparameterization', 'X', 'Y', 'Z', 'R', 'UnitX', 'UnitY',
     'UnitZ', 'Theta', 'Phi', 'a1', 'a2', 'a3', 'h'])

TdeDtGroupedPotential = TdeDt.groupby(['Potential'])
TdeDtGroupedPotentialTemperature = TdeDt.groupby(['Temperature', 'Potential'])
TdeDtGroupedPotentialDirection = TdeDt.groupby(['Potential', 'Direction', 'DirectionMb', 'Reparameterization', 'X', 'Y', 'Z', 'R', 'UnitX', 'UnitY', 'UnitZ',
     'Theta', 'Phi', 'a1', 'a2', 'a3', 'h'])
TdeDtGroupedPotentialDirectionTemperature = TdeDt.groupby(['Temperature', 'Potential', 'Direction', 'DirectionMb', 'Reparameterization', 'X', 'Y', 'Z', 'R', 'UnitX', 'UnitY',
     'UnitZ', 'Theta', 'Phi', 'a1', 'a2', 'a3', 'h'])

TdeWsGroupedPotential = TdeWs.groupby(['Potential'])
TdeWsGroupedPotentialTemperature = TdeWs.groupby(['Temperature', 'Potential'])
TdeWsGroupedPotentialDirection = TdeWs.groupby(['Potential', 'Direction', 'DirectionMb', 'Reparameterization', 'X', 'Y', 'Z', 'R', 'UnitX', 'UnitY', 'UnitZ',
     'Theta', 'Phi', 'a1', 'a2', 'a3', 'h'])
TdeWsGroupedPotentialDirectionTemperature = TdeWs.groupby(['Temperature', 'Potential', 'Direction', 'DirectionMb', 'Reparameterization', 'X', 'Y', 'Z', 'R', 'UnitX', 'UnitY',
     'UnitZ', 'Theta', 'Phi', 'a1', 'a2', 'a3', 'h'])

# </editor-fold>

# <editor-fold desc="^^^ Dictionary">
print("^^^ Dictionary")
TdeDefectGroupedPotentialDic = {
    'Ws': TdeWsGroupedPotential,
    'Dist': TdeDtGroupedPotential
}

TdeDefectGroupedPotentialDirectionDescribeDic = {
    'Ws': TdeWsGroupedPotentialDirectionDescribe,
    'Dist': TdeDtGroupedPotentialDirectionDescribe
}

# </editor-fold>

# </editor-fold>

# <editor-fold desc="***** Report">
print("***** Report")

# <editor-fold desc="^^^ Dataframe">
print("^^^ Dataframe")
ReadBatch = False
Method = "ExtractIndividually"
if Method == "ExtractBatch":

    # <editor-fold desc="Deleting">
    DeleteList = [
        '20431216-Nonstop.csv',
        '20431216-ErrorDf.csv',
        # '20431216-Df.csv',
        '20431216-BlankDf.csv',
        'exec_custom-Blank-One.sh',
        'exec_custom-Nonstop-One.sh'
        'exec_custom-Error-All.sh'
    ]

    Active = False
    if Active:
        for File in DeleteList:
            # Check if the file exists
            if os.path.exists(File):
                try:
                    # Attempt to delete the file
                    os.remove(File)
                    print(f"Deleted: {File}")
                except Exception as e:
                    # Print an error message if there's an issue
                    print(f"Error deleting {File}: {e}")
            else:
                print(f"File does not exist: {File}")
    # </editor-fold>

    # <editor-fold desc="Renaming">
    print("Renaming")
    Renaming = True
    if Renaming:
        Filename = "Tde.csv"
        TemperatureList = [10,300]
        PotentialList = ["M3R", "M2R", "BMD192R", "M2", "M3", "BMD192", "TabGap1", "TabGap2", "MtpPd"]
        TryList = range(1, 41)
        DirectionList = [
            "dir_0_0_1",
            "dir_0_1_1",
            "dir_0_3_1",
            "dir_0_1_0",
            "dir_-1_1_0",
            "dir_-1_1_1",
            "dir_-1_1_2",
            "dir_1_2_0",
            "dir_1_2_1",
            "dir_1_2_2",

            "dir_-2_3_3",
            "dir_-1_4_1",
            "dir_-1_3_1",
            "dir_-1_2_2",
            "dir_0_2_1",
            "dir_-1_2_0",
            "dir_-1_2_3",
            "dir_-4_4_1",
            "dir_-3_4_1",
            "dir_-2_4_1",
            "dir_-1_4_0",
            "dir_-1_2_1",
            "dir_0_2_3",
            "dir_-1_4_2",
            "dir_-2_2_3",
            "dir_-3_4_3",
            "dir_-2_3_2",
            "dir_-2_2_1",
            "dir_-1_3_0",
            "dir_-1_3_2",
            "dir_-3_4_2",
            "dir_-1_4_3",
            "dir_0_4_1",
            "dir_-3_4_0",
            "dir_-2_3_0",
            "dir_0_1_2",
            "dir_-2_3_1",
            "dir_-1_3_3",

            "dir_1_0_0",
            "dir_2_0_1",
            "dir_4_1_3",
        ]
        if SimulationEnvironment == "MyPc":
            for Temperature, Potential, Try, Direction in product(TemperatureList, PotentialList, TryList, DirectionList):
                BaseDirectory = os.path.join(CollectionAddress, str(Temperature), Potential, str(Try), Direction)
                if not os.path.isdir(BaseDirectory):
                    continue
                for Item in os.scandir(BaseDirectory):
                    if not Item.is_dir() or not Item.name.isdigit():
                        continue
                    Energy = Item.name

                    Path = os.path.join(Item.path, Filename)
                    if os.path.exists(Path):
                        continue

                    Path20250926 = os.path.join(Item.path, "20250926-Tde.csv")
                    Path20250130 = os.path.join(Item.path, "20250130-Tde.csv")

                    if any(os.path.exists(p) for p in (Path20250130, Path20250926)):
                        print(Temperature, Potential, Try, Direction, Energy)

                        if os.path.exists(Path20250926):
                            os.rename(Path20250926, Path)
                            if os.path.exists(Path20250130):
                                os.remove(Path20250130)
                        elif os.path.exists(Path20250130):
                            os.rename(Path20250130, Path)
    # </editor-fold>

    # <editor-fold desc="Extraction">
    print("Extraction")
    Method ="Tde-Defected"
    if Method == "Tde-All":
        # <editor-fold desc="Reading">
        Filename = "20250130-Tde.csv"
        # DirectionList = ["dir_0_1_2",
        #                  "dir_-1_2_1",
        #                  "dir_-1_3_2",
        #                  "dir_-1_4_2",
        #                  "dir_-2_3_3",
        #                  "dir_-4_4_1",
        #                  "dir_-2_3_0",
        #                  "dir_-1_2_3",
        #                  "dir_-3_4_0",
        #                  "dir_-1_2_0",
        #                  "dir_-3_4_3",
        #                  "dir_-1_4_3",
        #                  "dir_-1_3_3",
        #                  "dir_-1_2_2",
        #                  "dir_0_2_1",
        #                  "dir_-1_4_0",
        #                  "dir_-2_3_1",
        #                  "dir_-2_3_2",
        #                  "dir_-2_2_3",
        #                  "dir_0_2_3",
        #                  "dir_-2_4_1",
        #                  "dir_-2_2_1",
        #                  "dir_-3_4_1",
        #                  "dir_-1_3_1",
        #                  "dir_-1_4_1",
        #                  "dir_-3_4_2",
        #                  "dir_-1_3_0",
        #                  "dir_0_4_1"
        #                  ]
        PotentialList = ["M2R", "M3R", "BMD192R", "M2", "M3", "BMD192", "TabGap1", "TabGap2"]
        Header = ["Step","Time","Temp","Pressure","Volume","EnergyKinetic","EnergyPotential","EnergyTotal","EnergyTotalPerAtom","NoOfNeighbors","MinDist"]
        ReportDf = pd.DataFrame(columns=Header)

        if SimulationEnvironment == "MyPc":
            for Potential in PotentialList:
                print(Potential)
                for Try in np.arange(1,41):
                    # print(Try)
                    for Direction in DirectionList:
                        # print(Direction)
                        for Energy in np.arange(0,200):
                            # print(Energy)
                            Path = os.path.join(CurrentDirectory,Potential,str(Try), Direction,str(Energy),Filename)
                            if os.path.exists(Path):
                                try:
                                    # print(Potential, Try, Direction)
                                    New = pd.read_csv(Path)
                                except Exception as e:
                                    continue

                            ReportDf = pd.concat([ReportDf, New], ignore_index=True, axis=0)

        elif SimulationEnvironment == "ComputeCanada":
            Address = CurrentDirectory + "/" + filename
            # FolderName = str(Address.split("\\")[-1])
            FolderName = str(Address.split("/")[-2])
            print(FolderName)
            # os.system("pause")
            Read = extract_info(Address)
            DfNew = pd.DataFrame(Read).T.reset_index()
            DfNew = DfNew.rename(columns={'index': 'Try'})
            # print(DfNew)
            # os.system("pause")

            DfNew["Potential"] = FolderName
            Df = pd.concat([Df, DfNew], ignore_index=True, axis=0)

        ReportDf.to_csv(PythonName + "-Collection-ReportDf.csv", sep=',', index=False)
        # </editor-fold>
    elif Method == "Tde-Defected":
        # <editor-fold desc="Reading">
        PotentialList = ["M2R", "M3R", "BMD192R", "M2", "M3", "BMD192", "TabGap1", "TabGap2"]
        # DirectionList = ["dir_0_1_2",
        #                  "dir_-1_2_1",
        #                  "dir_-1_3_2",
        #                  "dir_-1_4_2",
        #                  "dir_-2_3_3",
        #                  "dir_-4_4_1",
        #                  "dir_-2_3_0",
        #                  "dir_-1_2_3",
        #                  "dir_-3_4_0",
        #                  "dir_-1_2_0",
        #                  "dir_-3_4_3",
        #                  "dir_-1_4_3",
        #                  "dir_-1_3_3",
        #                  "dir_-1_2_2",
        #                  "dir_0_2_1",
        #                  "dir_-1_4_0",
        #                  "dir_-2_3_1",
        #                  "dir_-2_3_2",
        #                  "dir_-2_2_3",
        #                  "dir_0_2_3",
        #                  "dir_-2_4_1",
        #                  "dir_-2_2_1",
        #                  "dir_-3_4_1",
        #                  "dir_-1_3_1",
        #                  "dir_-1_4_1",
        #                  "dir_-3_4_2",
        #                  "dir_-1_3_0",
        #                  "dir_0_4_1"
        #                  ]

        Header = ["Potential","Try","Direction","Energy","Step","Time","Temp","Pressure","Volume","EnergyKinetic","EnergyPotential","EnergyTotal","EnergyTotalPerAtom","NoOfNeighbors","MinDist"]
        ReportDf = pd.DataFrame(columns=Header)

        if SimulationEnvironment == "MyPc":
            for Potential in PotentialList:
                # print(Potential)
                for Try in np.arange(1,41):
                    # print(Try)
                    for Direction in DirectionList:
                        # print(Direction)
                        for Energy in np.arange(0,200):
                            # print(Energy)
                            PathMinDist = os.path.join(CurrentDirectory,Potential,str(Try), Direction,str(Energy),"20250130-Tde.csv")
                            PathDefect = os.path.join(CurrentDirectory,Potential,str(Try), Direction,str(Energy),"20250130-Tde-Defect.csv")
                            if os.path.exists(PathDefect):
                                with open(PathDefect, 'r') as File:
                                    Lines = File.readlines()
                                    if len(Lines) >= 2 and Lines[1].strip() == "1,1":
                                        try:
                                            # print(Potential, Try, Direction)
                                            New = pd.read_csv(PathMinDist)
                                            New["Potential"] = Potential
                                            New["Try"] = Try
                                            New["Direction"] = Direction
                                            New["Energy"] = Energy
                                            ReportDf = pd.concat([ReportDf, New], ignore_index=True, axis=0)
                                        except Exception as e:
                                            continue

        elif SimulationEnvironment == "ComputeCanada":
            Address = CurrentDirectory + "/" + filename
            # FolderName = str(Address.split("\\")[-1])
            FolderName = str(Address.split("/")[-2])
            print(FolderName)
            # os.system("pause")
            Read = extract_info(Address)
            DfNew = pd.DataFrame(Read).T.reset_index()
            DfNew = DfNew.rename(columns={'index': 'Try'})
            # print(DfNew)
            # os.system("pause")

            DfNew["Potential"] = FolderName
            Df = pd.concat([Df, DfNew], ignore_index=True, axis=0)

        ReportDf.to_csv(PythonName + "-Collection-ReportDf.csv", sep=',', index=False)
        # </editor-fold>
    elif Method == "Tde-Single":
        # <editor-fold desc="Reading">
        Filename = "20250130-Tde.csv"
        PotentialList = ["M3R"] #"M2R", "M3R", "BMD192R", "M2", "M3", "BMD192", "TabGap1", "TabGap2"]
        DirectionList = [
            #"dir_0_1_2",
             # "dir_-1_2_1",
             # "dir_-1_3_2",
             # "dir_-1_4_2",
             # "dir_-2_3_3",
             # "dir_-4_4_1",
             # "dir_-2_3_0",
             # "dir_-1_2_3",
             # "dir_-3_4_0",
             "dir_-1_2_0",
             # "dir_-3_4_3",
             # "dir_-1_4_3",
             # "dir_-1_3_3",
             # "dir_-1_2_2",
             # "dir_0_2_1",
             # "dir_-1_4_0",
             # "dir_-2_3_1",
             # "dir_-2_3_2",
             # "dir_-2_2_3",
             # "dir_0_2_3",
             # "dir_-2_4_1",
             # "dir_-2_2_1",
             # "dir_-3_4_1",
             # "dir_-1_3_1",
             # "dir_-1_4_1",
             # "dir_-3_4_2",
             # "dir_-1_3_0",
             # "dir_0_4_1"
             ]

        Header = ["Step","Time","Temp","Pressure","Volume","EnergyKinetic","EnergyPotential","EnergyTotal","EnergyTotalPerAtom","NoOfNeighbors","MinDist"]
        ReportDf = pd.DataFrame(columns=Header)

        if SimulationEnvironment == "MyPc":
            for Potential in PotentialList:
                print(Potential)
                for Try in np.arange(1,1):
                    print(Try)
                    for Direction in DirectionList:
                        print(Direction)
                        for Energy in np.arange(0,200):
                            print(Energy)
                            Path = os.path.join(CurrentDirectory,Potential,str(Try), Direction,str(Energy),Filename)
                            print(Path)
                            if os.path.exists(Path):
                                try:
                                    print(Potential, Try, Direction)
                                    New = pd.read_csv(Path)
                                except Exception as e:
                                    continue

                            ReportDf = pd.concat([ReportDf, New], ignore_index=True, axis=0)

        elif SimulationEnvironment == "ComputeCanada":
            Address = CurrentDirectory + "/" + filename
            # FolderName = str(Address.split("\\")[-1])
            FolderName = str(Address.split("/")[-2])
            print(FolderName)
            # os.system("pause")
            Read = extract_info(Address)
            DfNew = pd.DataFrame(Read).T.reset_index()
            DfNew = DfNew.rename(columns={'index': 'Try'})
            # print(DfNew)
            # os.system("pause")

            DfNew["Potential"] = FolderName
            Df = pd.concat([Df, DfNew], ignore_index=True, axis=0)

        ReportDf.to_csv(PythonName + "-Collection-ReportDf.csv", sep=',', index=False)
        # </editor-fold>
    # </editor-fold>

    # <editor-fold desc="Miller to Miller-Bravais">
    # Tde = Tde[Tde["File"] == "Tde"].copy()
    print(ReportDf)
    ReportDf['DirectionMb'] = ReportDf['Direction'].replace(Mapping)
    # print(Tde)
    # </editor-fold>

    # <editor-fold desc="Reparameterized">
    print("Reparameterized")
    ReportDf['Reparameterization'] = ReportDf['Potential'].apply(lambda x: 'R' if x.endswith('R') else 'O')
    # print(DfBravaisMelt)
    # </editor-fold>

    # <editor-fold desc="Polar">
    print("Polar")

    ReportDf[['X', 'Y', 'Z']] = ReportDf['Direction'].apply(DirectionFunc).apply(pd.Series)

    ReportDf['R'] = np.sqrt(ReportDf['X'] ** 2 + ReportDf['Y'] ** 2 + ReportDf['Z'] ** 2)
    ReportDf['UnitX'] = ReportDf['X'] / ReportDf['R']
    ReportDf['UnitY'] = ReportDf['Y'] / ReportDf['R']
    ReportDf['UnitZ'] = ReportDf['Z'] / ReportDf['R']

    ReportDf['Theta'] = np.degrees(np.arctan2(ReportDf['UnitY'], ReportDf['UnitX'])) % 360
    ReportDf['Phi'] = np.degrees(np.arccos(ReportDf['UnitZ']))

    ReportDf[['a1', 'a2', 'a3', 'h']] = ReportDf['DirectionMb'].apply(DirectionFunc).apply(pd.Series)
    # </editor-fold>

    # <editor-fold desc="Group">
    print("Group")
    # print(ReportDf)
    ReportDfGroupedPotential = ReportDf.groupby(['Potential'])
    ReportDfGroupedPotentialDirection = ReportDf.groupby(['Potential', 'Direction', 'DirectionMb', 'Reparameterization', 'X', 'Y', 'Z', 'R', 'UnitX', 'UnitY', 'UnitZ', 'Theta', 'Phi', 'a1', 'a2', 'a3', 'h'])
    # </editor-fold>

    # <editor-fold desc="Describe">
    print("Describe")
    # print(list(ReportDf.columns.values))
    # print(ReportDf)
    ReportDfGroupedPotentialDirectionStd = ReportDfGroupedPotentialDirection.std(numeric_only=True).reset_index()
    ReportDfGroupedPotentialDirectionDescribe = ReportDfGroupedPotentialDirection["Energy"].describe().reset_index()
    # </editor-fold>

    # <editor-fold desc="Export">
    ReportDf.to_csv(PythonName + "-ReportDf" + ".csv", index=False)
    ReportDfGroupedPotentialDirectionDescribe.to_csv(PythonName + "-" + "ReportDfGroupedPotentialDirectionDescribe.csv", sep=',', index=False)
    # </editor-fold>


elif Method == "ReadBatch":
    ReportDf = pd.read_csv(PythonName + "-ReportDf" + ".csv") # commented our because of the size, activate it later
    ReportDfGroupedPotentialDirectionDescribe = pd.read_csv(PythonName + "-" + "ReportDfGroupedPotentialDirectionDescribe.csv")

    # <editor-fold desc="Group">
    print("Group")
    ReportDfGroupedPotential = ReportDf.groupby(['Potential']) # commented our because of the size, activate it later
    ReportDfGroupedPotentialDirection = ReportDf.groupby(
        ['Potential', 'Direction', 'DirectionMb', 'Reparameterization', 'X', 'Y', 'Z', 'R', 'UnitX', 'UnitY', 'UnitZ',
         'Theta', 'Phi', 'a1', 'a2', 'a3', 'h'])
    # </editor-fold>

elif Method == "ExtractIndividually":
    print("ExtractIndividually")


# </editor-fold>

# </editor-fold>

# <editor-fold desc="***** Report-Single">
print("***** Report-Single")

# <editor-fold desc="Dataframe">
print("Dataframe")

Active=False
if Active:

    # <editor-fold desc="Deleting">
    DeleteList = [
        '20431216-Nonstop.csv',
        '20431216-ErrorDf.csv',
        # '20431216-Df.csv',
        '20431216-BlankDf.csv',
        'exec_custom-Blank-One.sh',
        'exec_custom-Nonstop-One.sh'
        'exec_custom-Error-All.sh'
    ]

    Active = False
    if Active:
        for File in DeleteList:
            # Check if the file exists
            if os.path.exists(File):
                try:
                    # Attempt to delete the file
                    os.remove(File)
                    print(f"Deleted: {File}")
                except Exception as e:
                    # Print an error message if there's an issue
                    print(f"Error deleting {File}: {e}")
            else:
                print(f"File does not exist: {File}")
    # </editor-fold>

    # <editor-fold desc="Extraction">
    print("Extraction")
    Method ="Tde-Single"
    if Method == "Tde-Single":
        # <editor-fold desc="Reading">
        Filename = "20250130-Tde.csv"
        DirectionList = [
            #"dir_0_1_2",
             # "dir_-1_2_1",
             # "dir_-1_3_2",
             # "dir_-1_4_2",
             # "dir_-2_3_3",
             # "dir_-4_4_1",
             # "dir_-2_3_0",
             # "dir_-1_2_3",
             # "dir_-3_4_0",
             "dir_-1_2_0",
             # "dir_-3_4_3",
             # "dir_-1_4_3",
             # "dir_-1_3_3",
             # "dir_-1_2_2",
             # "dir_0_2_1",
             # "dir_-1_4_0",
             # "dir_-2_3_1",
             # "dir_-2_3_2",
             # "dir_-2_2_3",
             # "dir_0_2_3",
             # "dir_-2_4_1",
             # "dir_-2_2_1",
             # "dir_-3_4_1",
             # "dir_-1_3_1",
             # "dir_-1_4_1",
             # "dir_-3_4_2",
             # "dir_-1_3_0",
             # "dir_0_4_1"
             ]
        PotentialList = ["M3R"] #"M2R", "M3R", "BMD192R", "M2", "M3", "BMD192", "TabGap1", "TabGap2"]
        Header = ["Potential", "Try", "Direction", "Energy",
                  "Step", "Time", "Temp", "Pressure", "Volume", "EnergyKinetic", "EnergyPotential", "EnergyTotal", "EnergyTotalPerAtom", "NoOfNeighbors", "MinDist"]
        ReportSingleDf = pd.DataFrame(columns=Header)

        if SimulationEnvironment == "MyPc":
            for Potential in PotentialList:
                # print(Potential)
                for Try in np.arange(1,2):
                    # print(Try)
                    for Direction in DirectionList:
                        # print(Direction)
                        # for Energy in np.arange(0,200):
                        for Energy in [26,28,30]:
                            # print(Energy)
                            Path = os.path.join(CurrentDirectory,Potential,str(Try), Direction,str(Energy),Filename)
                            # print(Path)
                            if os.path.exists(Path):
                                try:
                                    # print(Potential, Try, Direction)
                                    New = pd.read_csv(Path)
                                    # print(New)

                                    FilterSize = False
                                    if FilterSize:
                                        print("FilterSize")
                                        New = New.iloc[::10].reset_index(drop=True)

                                    New["Potential"] = Potential
                                    New["Try"] = Try
                                    New["Direction"] = Direction
                                    New["Energy"] = Energy
                                except Exception as e:
                                    continue

                            ReportSingleDf = pd.concat([ReportSingleDf, New], ignore_index=True, axis=0)

        elif SimulationEnvironment == "ComputeCanada":
            Address = CurrentDirectory + "/" + filename
            # FolderName = str(Address.split("\\")[-1])
            FolderName = str(Address.split("/")[-2])
            print(FolderName)
            # os.system("pause")
            Read = extract_info(Address)
            DfNew = pd.DataFrame(Read).T.reset_index()
            DfNew = DfNew.rename(columns={'index': 'Try'})
            # print(DfNew)
            # os.system("pause")

            DfNew["Potential"] = FolderName
            Df = pd.concat([Df, DfNew], ignore_index=True, axis=0)

        ReportSingleDf.to_csv(PythonName + "-Collection-ReportSingleDf.csv", sep=',', index=False)
        # </editor-fold>

    # </editor-fold>

    # <editor-fold desc="Miller to Miller-Bravais">
    print("Miller to Miller-Bravais")
    # Tde = Tde[Tde["File"] == "Tde"].copy()
    # print(ReportSingleDf)
    ReportSingleDf['DirectionMb'] = ReportSingleDf['Direction'].replace(Mapping)
    # print(Tde)
    # </editor-fold>

    # <editor-fold desc="Reparameterized">
    print("Reparameterized")
    ReportSingleDf['Reparameterization'] = ReportSingleDf['Potential'].apply(lambda x: 'R' if x.endswith('R') else 'O')
    # print(DfBravaisMelt)
    # </editor-fold>

    # <editor-fold desc="Polar">
    print("Polar")
    # print(ReportSingleDf)
    ReportSingleDf[['X', 'Y', 'Z']] = ReportSingleDf['Direction'].apply(DirectionFunc).apply(pd.Series)

    ReportSingleDf['R'] = np.sqrt(ReportSingleDf['X'] ** 2 + ReportSingleDf['Y'] ** 2 + ReportSingleDf['Z'] ** 2)
    ReportSingleDf['UnitX'] = ReportSingleDf['X'] / ReportSingleDf['R']
    ReportSingleDf['UnitY'] = ReportSingleDf['Y'] / ReportSingleDf['R']
    ReportSingleDf['UnitZ'] = ReportSingleDf['Z'] / ReportSingleDf['R']

    ReportSingleDf['Theta'] = np.degrees(np.arctan2(ReportSingleDf['UnitY'], ReportSingleDf['UnitX'])) % 360
    ReportSingleDf['Phi'] = np.degrees(np.arccos(ReportSingleDf['UnitZ']))

    ReportSingleDf[['a1', 'a2', 'a3', 'h']] = ReportSingleDf['DirectionMb'].apply(DirectionFunc).apply(pd.Series)
    # </editor-fold>

    # <editor-fold desc="Group">
    print("Group")
    # print(ReportSingleDf)
    ReportSingleDfGroupedPotential = ReportSingleDf.groupby(['Potential'])
    ReportSingleDfGroupedPotentialDirection = ReportSingleDf.groupby(['Potential', 'Direction', 'DirectionMb', 'Reparameterization', 'X', 'Y', 'Z', 'R', 'UnitX', 'UnitY', 'UnitZ', 'Theta', 'Phi', 'a1', 'a2', 'a3', 'h'])
    # </editor-fold>

    # <editor-fold desc="Describe">
    print("Describe")
    # print(list(ReportSingleDf.columns.values))
    # print(ReportSingleDf)
    ReportSingleDfGroupedPotentialDirectionStd = ReportSingleDfGroupedPotentialDirection.std(numeric_only=True).reset_index()
    ReportSingleDfGroupedPotentialDirectionDescribe = ReportSingleDfGroupedPotentialDirection["Energy"].describe().reset_index()
    # </editor-fold>

    # <editor-fold desc="Export">
    ReportSingleDf.to_csv(PythonName + "-ReportSingleDf" + ".csv", index=False)
    ReportSingleDfGroupedPotentialDirectionDescribe.to_csv(PythonName + "-" + "ReportSingleDfGroupedPotentialDirectionDescribe.csv", sep=',', index=False)
    # </editor-fold>

else:
    # ReportSingleDf = pd.read_csv(PythonName + "-ReportSingleDf" + ".csv") # commented our because of the size, activate it later
    # ReportSingleDfGroupedPotentialDirectionDescribe = pd.read_csv(PythonName + "-" + "ReportSingleDfGroupedPotentialDirectionDescribe.csv")

    # <editor-fold desc="Group">
    print("Group")
    # print(Tde)
    # ReportSingleDfGroupedPotential = ReportSingleDf.groupby(['Potential']) # commented our because of the size, activate it later
    # ReportSingleDfGroupedPotentialDirection = ReportSingleDf.groupby(
    #     ['Potential', 'Direction', 'DirectionMb', 'Reparameterization', 'X', 'Y', 'Z', 'R', 'UnitX', 'UnitY', 'UnitZ',
    #      'Theta', 'Phi', 'a1', 'a2', 'a3', 'h'])
    # </editor-fold>

# </editor-fold>

# </editor-fold>

# <editor-fold desc="***** Thermal">
print("***** Thermal")

# <editor-fold desc="Read">
Thermal10 = pd.read_csv("D:/Queens_University/MME/Project/Zr/Tde/Run/10/20250506-Collection-Thermal.csv")
Thermal300 = pd.read_csv("D:/Queens_University/MME/Project/Zr/Tde/Run/300/20250506-Collection-Thermal.csv")
Thermal10["Temperature"] = 10
Thermal300["Temperature"] = 300
Thermal = pd.concat([Thermal10, Thermal300], ignore_index=True, axis=0)
# </editor-fold>

# <editor-fold desc="Plot">
Plotting = True
if Plotting:
    files = Thermal['File'].unique()

    for f in files:
        subset = Thermal[Thermal['File'] == f]
        plt.figure(figsize=(6, 6))
        sns.barplot(data=subset, x='Potential', y='EnergyPerAtom', hue='LookupWord')
        plt.title(f'Energy per Atom vs Potential ({f})')
        plt.xlabel('Potential')
        plt.ylabel('Energy per Atom')
        plt.tight_layout()
        plt.show()

# </editor-fold>


# </editor-fold>

# <editor-fold desc="######################################## Ws vs. Distance">
print("######################################## Ws vs. Distance") # %%

Active = False
if Active:
    print("hi")

# </editor-fold>

# <editor-fold desc="######################################## Tde Stats">
print("######################################## Tde Stats") # %%

# <editor-fold desc="***** Stat">
print("***** Stat")
# <editor-fold desc="Dt">
# print(TdeDtGroupedPotentialDirection)
# os.system("pause")
TdeDtGroupedPotentialDirectionTemperatureMin = TdeDtGroupedPotentialDirectionTemperature["Energy"].min().reset_index()
TdeDtGroupedPotentialDirectionTemperatureMin.to_csv(PythonName + "-" + "TdeDtGroupedPotentialDirectionTemperatureMin.csv", index=False)

TdeDtGroupedPotentialDirectionTemperatureMean = TdeDtGroupedPotentialDirectionTemperature["Energy"].mean().reset_index()
TdeDtGroupedPotentialDirectionTemperatureMean.to_csv(PythonName + "-" + "TdeDtGroupedPotentialDirectionTemperatureMean.csv", index=False)

TdeDtGroupedPotentialDirectionTemperatureMedian = (TdeDtGroupedPotentialDirectionTemperature["Energy"].median().reset_index())
TdeDtGroupedPotentialDirectionTemperatureMedian.to_csv(PythonName + "-" + "TdeDtGroupedPotentialDirectionTemperatureMedian.csv", index=False)

TdeDtGroupedPotentialDirectionTemperatureMode = (TdeDtGroupedPotentialDirectionTemperature["Energy"].apply(lambda x: x.mode().iloc[0]).reset_index())
TdeDtGroupedPotentialDirectionTemperatureMode.to_csv(PythonName + "-" + "TdeDtGroupedPotentialDirectionTemperatureMode.csv", index=False)

TdeDtGroupedPotentialDirectionTemperatureVariance = (TdeDtGroupedPotentialDirectionTemperature["Energy"].var().reset_index())
TdeDtGroupedPotentialDirectionTemperatureVariance.to_csv(PythonName + "-" + "TdeDtGroupedPotentialDirectionTemperatureVar.csv", index=False)
# </editor-fold>

# <editor-fold desc="Ws">
TdeWsGroupedPotentialDirectionTemperatureMin = TdeWsGroupedPotentialDirectionTemperature["Energy"].min().reset_index()
TdeWsGroupedPotentialDirectionTemperatureMin.to_csv(PythonName + "-" + "TdeWsGroupedPotentialDirectionTemperatureMin.csv", index=False)

TdeWsGroupedPotentialDirectionTemperatureMean = TdeWsGroupedPotentialDirectionTemperature["Energy"].mean().reset_index()
TdeWsGroupedPotentialDirectionTemperatureMean.to_csv(PythonName + "-" + "TdeWsGroupedPotentialDirectionTemperatureMean.csv", index=False)

TdeWsGroupedPotentialDirectionTemperatureMedian = (TdeWsGroupedPotentialDirectionTemperature["Energy"].median().reset_index())
TdeWsGroupedPotentialDirectionTemperatureMedian.to_csv(PythonName + "-" + "TdeWsGroupedPotentialDirectionTemperatureMedian.csv", index=False)

TdeWsGroupedPotentialDirectionTemperatureMode = (TdeWsGroupedPotentialDirectionTemperature["Energy"].apply(lambda x: x.mode().iloc[0]).reset_index())
TdeWsGroupedPotentialDirectionTemperatureMode.to_csv(PythonName + "-" + "TdeWsGroupedPotentialDirectionTemperatureMode.csv", index=False)

TdeWsGroupedPotentialDirectionTemperatureVariance = (TdeWsGroupedPotentialDirectionTemperature["Energy"].var().reset_index())
TdeWsGroupedPotentialDirectionTemperatureVariance.to_csv(PythonName + "-" + "TdeWsGroupedPotentialDirectionTemperatureVar.csv", index=False)
# </editor-fold>

# </editor-fold>

# <editor-fold desc="***** Literature">
print("***** Literature")
TdeDtUniqueDirections = TdeDt["Direction"].unique()
# print(TdeDt["Temperature"].dtype)
TdeDtUniqueTemperatures = TdeDt["Temperature"].unique()
# print(TdeDtUniqueTemperatures)
# os.system("pause")

Combination = pd.MultiIndex.from_product(
    [TdeDtUniqueDirections, TdeDtUniqueTemperatures],
    names=['Direction', 'Temperature']
)

CombinationDf = Combination.to_frame(index=False)
# print(CombinationDf)
# os.system("pause")
LiteratureBiget = CombinationDf.copy()
LiteratureBiget["Potential"] = "Biget"
LiteratureBiget["Energy"] = 21

LiteratureNeely = CombinationDf.copy()
LiteratureNeely["Potential"] = "Neely"
LiteratureNeely["Energy"] = 24

LiteratureBiget['DirectionMb'] = LiteratureBiget['Direction'].replace(Mapping)
LiteratureNeely['DirectionMb'] = LiteratureNeely['Direction'].replace(Mapping)

TdeLiterature = pd.concat([LiteratureBiget, LiteratureNeely], ignore_index=True)
# print(TdeLiterature)
# TdeLiterature.to_csv("test.csv", sep=',', index=False)
# os.system("pause")
TdeDtGroupedPotentialDirectionTemperatureMin = pd.concat([TdeDtGroupedPotentialDirectionTemperatureMin, TdeLiterature], ignore_index=True)
TdeDtGroupedPotentialDirectionTemperatureMean = pd.concat([TdeDtGroupedPotentialDirectionTemperatureMean, TdeLiterature], ignore_index=True)

TdeWsGroupedPotentialDirectionTemperatureMin = pd.concat([TdeWsGroupedPotentialDirectionTemperatureMin, TdeLiterature], ignore_index=True)
TdeWsGroupedPotentialDirectionTemperatureMean = pd.concat([TdeWsGroupedPotentialDirectionTemperatureMean, TdeLiterature], ignore_index=True)
# </editor-fold>

# <editor-fold desc="***** Dic">
print("***** Dic")
TdeDtGroupedPotentialDirectionTemperatureStatDic = {
    "Min": TdeDtGroupedPotentialDirectionTemperatureMin,
    "Mean": TdeDtGroupedPotentialDirectionTemperatureMean,
    "Median": TdeDtGroupedPotentialDirectionTemperatureMedian,
    "Mode": TdeDtGroupedPotentialDirectionTemperatureMode,
    "Variance": TdeDtGroupedPotentialDirectionTemperatureVariance,
}

TdeWsGroupedPotentialDirectionTemperatureStatDic = {
    "Min": TdeWsGroupedPotentialDirectionTemperatureMin,
    "Mean": TdeWsGroupedPotentialDirectionTemperatureMean,
    "Median": TdeWsGroupedPotentialDirectionTemperatureMedian,
    "Mode": TdeWsGroupedPotentialDirectionTemperatureMode,
    "Variance": TdeWsGroupedPotentialDirectionTemperatureVariance,
}

TdeGroupedPotentialDirectionTemperatureStatDic = {
    "Ws": TdeWsGroupedPotentialDirectionTemperatureStatDic,
    "Dt": TdeDtGroupedPotentialDirectionTemperatureStatDic,
}
# </editor-fold>

# <editor-fold desc="***** Plot">
print("***** Plot")
Plotting = False
if Plotting:
    for Analysis, Dic in TdeGroupedPotentialDirectionTemperatureStatDic.items():
        # print(Analysis)
        for Stat, Df in TdeDtGroupedPotentialDirectionTemperatureStatDic.items():
            # print(Stat)
            # print(Df)
            # Df.to_csv("test.csv", sep=',', index=False)
            for Temperature in Df['Temperature'].unique():
                # print(Temperature)
                print(Analysis,Stat,Temperature)
                HeatmapDf = Df[Df['Temperature'] == Temperature]

                DirectionListRemove = ["dir_4_1_3","dir_1_0_0","dir_2_0_1"]  # e.g., ["dir_0_1_2"] if you want to exclude some
                DirectionListRemove = []
                HeatmapDf = HeatmapDf[~HeatmapDf["Direction"].isin(DirectionListRemove)].copy()

                HeatmapDf["DirectionMb"] = HeatmapDf["DirectionMb"].apply(LatexFriendlyFunc)

                potential_order = [
                    "M2", "M3", "BMD192", "M2R", "M3R", "BMD192R",
                    "TabGap1", "TabGap2", "Biget", "Neely"
                ]
                HeatmapDf["Potential"] = pd.Categorical(
                    HeatmapDf["Potential"], categories=potential_order, ordered=True
                )

                Title = PythonName + "-TdeGroupedPotentialDirectionTemperatureStatDic" + "-" + str(Analysis) + "-" + str(Stat) + "-" + str(Temperature)

                HeatmapPivot = HeatmapDf.pivot_table(
                    index="Potential",
                    columns="DirectionMb",
                    values="Energy",
                    aggfunc="mean"  # safe even if duplicates exist
                )

                # Compute average across each row (Potential)
                row_means = HeatmapPivot.mean(axis=1, skipna=True).round(2)

                # Append averages to the row labels
                HeatmapPivot.index = [
                    f"{idx} ({row_means.loc[idx]})" for idx in HeatmapPivot.index
                ]

                plt.figure(figsize=(18, 4))
                sns.heatmap(
                    HeatmapPivot,
                    annot=False,
                    cmap="coolwarm",
                    cbar_kws={'label': str(Stat)}
                )
                plt.xlabel("")
                plt.ylabel("")
                plt.tight_layout()
                plt.savefig(f"{Title}.png", dpi=600)
                plt.show()
                plt.close()

# </editor-fold>
# </editor-fold>

# <editor-fold desc="######################################## Missing">
print("######################################## Missing") # %%

Active = False
if Active:
    # <editor-fold desc="Deleting">
    DeleteList = [
        '20431216-Nonstop.csv',
        '20431216-ErrorDf.csv',
        # '20431216-Tde.csv',
        '20431216-BlankDf.csv',
        'exec_custom-Blank-One.sh',
        'exec_custom-Nonstop-One.sh'
        'exec_custom-Error-All.sh'
    ]
    Active = False
    if Active:
        for File in DeleteList:
            # Check if the file exists
            if os.path.exists(File):
                try:
                    # Attempt to delete the file
                    os.remove(File)
                    print(f"Deleted: {File}")
                except Exception as e:
                    # Print an error message if there's an issue
                    print(f"Error deleting {File}: {e}")
            else:
                print(f"File does not exist: {File}")
    # </editor-fold>

    # <editor-fold desc="Checking">
    BlankList = []
    NonstopList = []

    # Iterate through each row in the dataframe
    for index, row in Tde.iterrows():
        potential = row['Potential']
        Try = row['Try']
        # print(row)
        for column in Tde.columns[:-2]:  # Exclude 'Try' and 'Potential' columns
            print(row[column])
            if pd.isnull(row[column]):
                # Append data for BlankDf
                BlankList.append([potential, Try, column])
            if row[column] == 200:
                # Append data for ErrorDf
                NonstopList.append([potential, Try, column])

    # Create dataframes
    BlankDf = pd.DataFrame(BlankList, columns=['Potential', 'Try', 'Directory'])
    BlankDf["ErrorType"] = "Blank"
    NonstopDf = pd.DataFrame(NonstopList, columns=['Potential', 'Try', 'Directory'])
    NonstopDf["ErrorType"] = "Nonstop"
    ErrorDf = pd.concat([BlankDf, NonstopDf], ignore_index=True)

    # BlankDf.to_csv(PythonName + "-BlankDf.csv", sep=',',index=False)
    # NonstopDf.to_csv(PythonName + "-Nonstop.csv", sep=',',index=False)
    ErrorDf.to_csv(PythonName + "-ErrorDf.csv", sep=',', index=False)
    # </editor-fold>

    # <editor-fold desc="***** All">
    print("***** All")
    # <editor-fold desc="************ Blank">
    print("************ Blank")
    BlankDfGrouped = BlankDf.groupby(['Potential', 'Try']).filter(lambda x: x['Directory'].nunique() > 1)
    BlankDfGroupedUnique = BlankDfGrouped[['Potential', 'Try']].drop_duplicates()
    BlankDfGroupedUniqueGrouped = BlankDfGroupedUnique.groupby('Potential')['Try'].apply(list)
    # print(BlankDfGroupedUniqueGrouped)
    # os.system("pause")
    for potential, Trys in BlankDfGroupedUniqueGrouped.items():
        Trys_str = ' '.join(map(str, Trys))
        print(f"{potential}: {' '.join(map(str, Trys))}")
        # if SimulationEnvironment == "MyPc":
        #     FileName = f"{potential}-Blank-All.sh"
        # elif SimulationEnvironment == "ComputeCanada":
        #     FileName = "exec_custom-Blank-All.sh"
        #
        # with open(FileName, "w") as file:
        #     file.write(f"echo {potential}:\n + "
        #                f"declare -a ReRunList=({Trys_str})\n" +
        #                'for k in "${ReRunList[@]}";  do\n'
        #                ' echo $k\n'
        #                ' cd $k\n'
        #                ' find -mindepth 1 -maxdepth 1 -type d -exec rm -r {} \;\n'
        #                ' sbatch run_active2\n'
        #                ' cd ..\n'
        #                'done')
    # </editor-fold>

    # <editor-fold desc="************ Nonstop">
    print("************ Nonstop")
    NonstopDfGrouped = NonstopDf.groupby(['Potential', 'Try']).filter(lambda x: x['Directory'].nunique() > 1)
    NonstopDfGroupedUnique = NonstopDfGrouped[['Potential', 'Try']].drop_duplicates()
    NonstopDfGroupedUniqueGrouped = NonstopDfGroupedUnique.groupby('Potential')['Try'].apply(list)
    # print(NonstopDfGroupedUniqueGrouped)

    for potential, Trys in NonstopDfGroupedUniqueGrouped.items():
        Trys_str = ' '.join(map(str, Trys))
        print(f"{potential}: {' '.join(map(str, Trys))}")
        # if SimulationEnvironment == "MyPc":
        #     FileName = f"{potential}-Nonstop-All.sh"
        # elif SimulationEnvironment == "ComputeCanada":
        #     FileName = "exec_custom-Nonstop-All.sh"
        #
        # with open(FileName, "w") as file:
        #     file.write(f"echo {potential}:\n" +
        #                f"declare -a ReRunList=({Trys_str})\n" +
        #                'for k in "${ReRunList[@]}";  do\n'
        #                ' echo $k\n'
        #                ' cd $k\n'
        #                ' find -mindepth 1 -maxdepth 1 -type d -exec rm -r {} \;\n'
        #                ' sbatch run_active2\n'
        #                ' cd ..\n'
        #                'done')
    # </editor-fold>

    # <editor-fold desc="************ Merging">
    print("************ Merging")
    ErrorAll = pd.concat([BlankDfGroupedUnique, NonstopDfGroupedUnique], axis=0)
    ErrorAllUnique = ErrorAll[['Potential', 'Try']].drop_duplicates()
    ErrorAllUniqueGrouped = ErrorAllUnique.groupby('Potential')['Try'].apply(list)
    # print(ErrorAllUniqueGrouped)

    for potential, Trys in ErrorAllUniqueGrouped.items():
        Trys_str = ' '.join(map(str, Trys))
        # print(f"{potential}: {' '.join(map(str, Trys))}")
        if SimulationEnvironment == "MyPc":
            FileName = f"{potential}-Error-All.sh"  # if needed, shift the writing part to write bash files in windows too
        elif SimulationEnvironment == "ComputeCanada":
            FileName = "exec_custom-Error-All.sh"

            with open(FileName, "w") as file:
                file.write("cluster_name=$(scontrol show config | grep ClusterName | awk '{print $3}')\n"
                           'if [ "$cluster_name" == "beluga" ]; then\n'
                           '    NewAccount="#SBATCH --account=rrg-belandl1"\n'
                           '    NewNtask="#SBATCH --ntasks=32"\n'
                           '    NewMemPerCpu="#SBATCH --mem-per-cpu=1000M"\n'
                           '    NewTime="#SBATCH --time=0-23:59"\n'
                           'elif [ "$cluster_name" == "graham" ]; then\n'
                           '    NewAccount="#SBATCH --account=def-belandl1"\n'
                           '    NewNtask="#SBATCH --ntasks=32"\n'
                           '    NewMemPerCpu="#SBATCH --mem-per-cpu=1000M"\n'
                           '    NewTime="#SBATCH --time=0-23:59"\n'
                           '    NewScratch="/scratch/"\n'
                           '    NewSoftware="/Software/"\n'
                           'elif [ "$cluster_name" == "cedar" ]; then\n'
                           '    NewAccount="#SBATCH --account=def-belandl1"\n'
                           '    NewNtask="#SBATCH --ntasks=32"\n'
                           '    NewMemPerCpu="#SBATCH --mem-per-cpu=1000M"\n'
                           '    NewTime="#SBATCH --time=0-23:59"\n'
                           '    NewScratch="/scratch/"\n'
                           '    NewSoftware="/Software/"\n'
                           'fi\n'
    
                           f"echo {potential}:\n" +
                           f"declare -a ReRunList=({Trys_str})\n" +
                           'for k in "${ReRunList[@]}";  do\n'
                           '    echo $k\n'
                           '    cd $k\n'
                           '    rm -rf ./*/\n'
    
                           '    OldKeyword="--account"\n'
                           '    sed -i "/$OldKeyword/c' + fr"\\{''}" + '$NewAccount" "run_active2"\n'
                                                                       '    OldKeyword="--ntasks"\n'
                                                                       '    sed -i "/$OldKeyword/c' + fr"\\{''}" + '$NewNtask" "run_active2"\n'
                                                                                                                   '    OldKeyword="--mem-per-cpu"\n'
                                                                                                                   '    sed -i "/$OldKeyword/c' + fr"\\{''}" + '$NewMemPerCpu" "run_active2"\n'
                                                                                                                                                               '    OldKeyword="--time"\n'
                                                                                                                                                               '    sed -i "/$OldKeyword/c' + fr"\\{''}" + '$NewTime" "run_active2"\n'
                                                                                                                                                                                                           '    OldKeyword="/lustre04/scratch/"\n'
                                                                                                                                                                                                           '    sed -i "s|$OldKeyword|$NewScratch|g" "run_active2"\n'
                                                                                                                                                                                                           '    OldKeyword="/software/"\n'
                                                                                                                                                                                                           '    sed -i "s|$OldKeyword|$NewSoftware|Ig" "run_active2"\n'
                                                                                                                                                                                                           '    sbatch run_active2\n'
                                                                                                                                                                                                           '    cd ..\n'
                                                                                                                                                                                                           'done\n')
    # </editor-fold>
    # print(merged_unique)
    # </editor-fold>

    # <editor-fold desc="***** One by One">
    print("***** One by One")
    if BlankDfGroupedUniqueGrouped.empty and NonstopDfGroupedUniqueGrouped.empty:

        # <editor-fold desc="************ Blank">
        print("************ Blank")
        BlankDfGrouped = BlankDf.groupby(['Potential', 'Try']).filter(lambda x: x['Directory'].nunique() == 1)
        BlankDfGroupedUnique = BlankDfGrouped[['Potential', 'Try']].drop_duplicates()
        for index, row in BlankDfGroupedUnique.iterrows():
            Potential = row['Potential']
            Try = row['Try']
            # print(Try)
            Direction = BlankDfGrouped[
                (BlankDfGrouped['Potential'] == row['Potential']) & (BlankDfGrouped['Try'] == row['Try'])][
                'Directory'].iloc[0]
            DirectionX = str(Direction.split("_")[-3])
            DirectionY = str(Direction.split("_")[-2])
            DirectionZ = str(Direction.split("_")[-1])
            DirectionLine = "  " + '"' + DirectionX + " " + DirectionY + " " + DirectionZ + '"'
            # print(Direction + ":" + DirectionX + "," + DirectionY + "," + DirectionZ)
            # print(DirectionLine)
            # print(f"{row['Potential']}: {row['Try']}: {BlankDfGrouped[(BlankDfGrouped['Potential'] == row['Potential']) & (BlankDfGrouped['Try'] == row['Try'])]['Directory'].iloc[0]}")

            FilePath = CurrentDirectory + "/" + Try + "/run_active2"  # Update with your file path
            # print(FilePath)

            DirectionLineList = ['  "0 0 1"', '  "0 1 1"', '  "0 3 1"', '  "0 1 0"', '  "-1 1 0"', '  "-1 1 1"',
                                 '  "-1 1 2"', '  "1 2 0"', '  "1 2 1"', '  "1 2 2"']
            for i in DirectionLineList:
                if SimulationEnvironment == "ComputeCanada":
                    try:
                        with open(FilePath) as r:
                            text = r.read().replace(i, '#' + i)
                        with open(FilePath, "w") as w:
                            w.write(text)
                    except:
                        print(FilePath + " Not Exists")

            try:
                with open(FilePath) as r:
                    text = r.read().replace('#' + DirectionLine, DirectionLine)
                with open(FilePath, "w") as w:
                    w.write(text)
            except:
                print(FilePath + " Not Exists")

            if DirectionLine == '  "1 2 2"':  # to make sure it wasn't because of run_out_of_time error
                try:
                    with open(FilePath) as r:
                        text = r.read().replace('#  "1 2 1"', '  "1 2 1"')
                    with open(FilePath, "w") as w:
                        w.write(text)
                except:
                    print(FilePath + " Not Exists")

        BlankDfGroupedUniqueGrouped = BlankDfGroupedUnique.groupby('Potential')['Try'].apply(list)
        # print(BlankDfGroupedUniqueGrouped)

        for potential, Trys in BlankDfGroupedUniqueGrouped.items():
            # print(potential)
            # print(Trys)
            # os.system("pause")
            Trys_str = ' '.join(map(str, Trys))
            # print(f"{potential}:\nReRunList = ({Trys_str})\n" +
            #       'for k in "${ReRunList[@]}";  do\n'
            #       ' cd $k\n'
            #       ' find -mindepth 1 -maxdepth 1 -type d -exec rm -r {} \;\n'
            #       ' sbatch run_active2\n'
            #       ' cd ..\n'
            #       'done')
            if SimulationEnvironment == "MyPc":
                FileName = f"{potential}-Blank-One.sh"
            elif SimulationEnvironment == "ComputeCanada":
                FileName = "exec_custom-Blank-One.sh"

                # SubText = ""
                # for Sub in Trys:
                #     SubText = SubText + "cd " + Sub + "\n" + "sbatch run_active2" + "\n" + "cd .."  + "\n"

                with open(FileName, "w") as file:
                    file.write("cluster_name=$(scontrol show config | grep ClusterName | awk '{print $3}')\n"
                               'if [ "$cluster_name" == "beluga" ]; then\n'
                               '    NewAccount="#SBATCH --account=rrg-belandl1"\n'
                               '    NewNtask="#SBATCH --ntasks=32"\n'
                               '    NewMemPerCpu="#SBATCH --mem-per-cpu=1000M"\n'
                               '    NewTime="#SBATCH --time=0-23:59"\n'
                               'elif [ "$cluster_name" == "graham" ]; then\n'
                               '    NewAccount="#SBATCH --account=def-belandl1"\n'
                               '    NewNtask="#SBATCH --ntasks=32"\n'
                               '    NewMemPerCpu="#SBATCH --mem-per-cpu=1000M"\n'
                               '    NewTime="#SBATCH --time=0-23:59"\n'
                               '    NewScratch="/scratch/"\n'
                               '    NewSoftware="/Software/"\n'
                               'elif [ "$cluster_name" == "cedar" ]; then\n'
                               '    NewAccount="#SBATCH --account=def-belandl1"\n'
                               '    NewNtask="#SBATCH --ntasks=32"\n'
                               '    NewMemPerCpu="#SBATCH --mem-per-cpu=1000M"\n'
                               '    NewTime="#SBATCH --time=0-23:59"\n'
                               '    NewScratch="/scratch/"\n'
                               '    NewSoftware="/Software/"\n'
                               'fi\n'
    
                               f"echo {potential}:\nReRunList=({Trys_str})\n" +
                               'for dir in "${ReRunList[@]}"; do\n'
                               '    cd "$dir" || { echo "Directory $dir not found"; exit 1; }\n'
    
                               '    OldKeyword="--account"\n'
                               '    sed -i "/$OldKeyword/c' + fr"\\{''}" + '$NewAccount" "run_active2"\n'
                                                                           '    OldKeyword="--ntasks"\n'
                                                                           '    sed -i "/$OldKeyword/c' + fr"\\{''}" + '$NewNtask" "run_active2"\n'
                                                                                                                       '    OldKeyword="--mem-per-cpu"\n'
                                                                                                                       '    sed -i "/$OldKeyword/c' + fr"\\{''}" + '$NewMemPerCpu" "run_active2"\n'
                                                                                                                                                                   '    OldKeyword="--time"\n'
                                                                                                                                                                   '    sed -i "/$OldKeyword/c' + fr"\\{''}" + '$NewTime" "run_active2"\n'
                                                                                                                                                                                                               '    OldKeyword="/lustre04/scratch/"\n'
                                                                                                                                                                                                               '    sed -i "s|$OldKeyword|$NewScratch|g" "run_active2"\n'
                                                                                                                                                                                                               '    OldKeyword="/software/"\n'
                                                                                                                                                                                                               '    sed -i "s|$OldKeyword|$NewSoftware|Ig" "run_active2"\n'
                                                                                                                                                                                                               '    sbatch run_active2\n'
                                                                                                                                                                                                               '    cd ..\n'
                                                                                                                                                                                                               'done\n'
                               )
        # </editor-fold>

        # <editor-fold desc="************ Nonstop">
        print("************ Nonstop")
        NonstopDfGrouped = NonstopDf.groupby(['Potential', 'Try']).filter(lambda x: x['Directory'].nunique() == 1)
        NonstopDfGroupedUnique = NonstopDfGrouped[['Potential', 'Try']].drop_duplicates()
        for index, row in NonstopDfGroupedUnique.iterrows():
            Potential = row['Potential']
            Try = row['Try']
            # print(Try)
            Direction = NonstopDfGrouped[
                (NonstopDfGrouped['Potential'] == row['Potential']) & (NonstopDfGrouped['Try'] == row['Try'])][
                'Directory'].iloc[0]
            DirectionX = str(Direction.split("_")[-3])
            DirectionY = str(Direction.split("_")[-2])
            DirectionZ = str(Direction.split("_")[-1])
            DirectionLine = "  " + '"' + DirectionX + " " + DirectionY + " " + DirectionZ + '"'
            # print(Direction + ":" + DirectionX + "," + DirectionY + "," + DirectionZ)
            # print(DirectionLine)
            # print(f"{row['Potential']}: {row['Try']}: {NonstopDfGrouped[(NonstopDfGrouped['Potential'] == row['Potential']) & (NonstopDfGrouped['Try'] == row['Try'])]['Directory'].iloc[0]}")

            FilePath = CurrentDirectory + "/" + str(Try) + "/run_active2"  # Update with your file path
            # print(FilePath)

            DirectionLineList = ['  "0 0 1"', '  "0 1 1"', '  "0 3 1"', '  "0 1 0"', '  "-1 1 0"', '  "-1 1 1"',
                                 '  "-1 1 2"', '  "1 2 0"', '  "1 2 1"', '  "1 2 2"']
            for i in DirectionLineList:
                try:
                    with open(FilePath) as r:
                        text = r.read().replace(i, '#' + i)
                    with open(FilePath, "w") as w:
                        w.write(text)
                except:
                    print(FilePath + " Not Exists")

            try:
                with open(FilePath) as r:
                    text = r.read().replace('#' + DirectionLine, DirectionLine)
                with open(FilePath, "w") as w:
                    w.write(text)
            except:
                print(FilePath + " Not Exists")

            if DirectionLine == '  "1 2 2"':  # to make sure it wasn't because of run_out_of_time error
                try:
                    with open(FilePath) as r:
                        text = r.read().replace('#  "1 2 1"', '  "1 2 1"')
                    with open(FilePath, "w") as w:
                        w.write(text)
                except:
                    print(FilePath + " Not Exists")

        NonstopDfGroupedUniqueGrouped = NonstopDfGroupedUnique.groupby('Potential')['Try'].apply(list)
        # print(NonstopDfGroupedUniqueGrouped)

        for potential, Trys in NonstopDfGroupedUniqueGrouped.items():
            # print(potential)
            # print(Trys)
            # os.system("pause")
            Trys_str = ' '.join(map(str, Trys))
            # print(f"{potential}:\nReRunList = ({Trys_str})\n" +
            #       'for k in "${ReRunList[@]}";  do\n'
            #       ' cd $k\n'
            #       ' find -mindepth 1 -maxdepth 1 -type d -exec rm -r {} \;\n'
            #       ' sbatch run_active2\n'
            #       ' cd ..\n'
            #       'done')
            if SimulationEnvironment == "MyPc":
                FileName = f"{potential}-Nonstop-One.sh"
            elif SimulationEnvironment == "ComputeCanada":
                FileName = "exec_custom-Nonstop-One.sh"

                # SubText = ""
                # for Sub in Trys:
                #     SubText = SubText + "cd " + Sub + "\n" + "sbatch run_active2" + "\n" + "cd .."  + "\n"

                with open(FileName, "w") as file:
                    file.write("cluster_name=$(scontrol show config | grep ClusterName | awk '{print $3}')\n"
                               'if [ "$cluster_name" == "beluga" ]; then\n'
                               '    NewAccount="#SBATCH --account=rrg-belandl1"\n'
                               '    NewNtask="#SBATCH --ntasks=32"\n'
                               '    NewMemPerCpu="#SBATCH --mem-per-cpu=1000M"\n'
                               '    NewTime="#SBATCH --time=0-23:59"\n'
                               'elif [ "$cluster_name" == "graham" ]; then\n'
                               '    NewAccount="#SBATCH --account=def-belandl1"\n'
                               '    NewNtask="#SBATCH --ntasks=32"\n'
                               '    NewMemPerCpu="#SBATCH --mem-per-cpu=1000M"\n'
                               '    NewTime="#SBATCH --time=0-23:59"\n'
                               '    NewScratch="/scratch/"\n'
                               '    NewSoftware="/Software/"\n'
                               'elif [ "$cluster_name" == "cedar" ]; then\n'
                               '    NewAccount="#SBATCH --account=def-belandl1"\n'
                               '    NewNtask="#SBATCH --ntasks=32"\n'
                               '    NewMemPerCpu="#SBATCH --mem-per-cpu=1000M"\n'
                               '    NewTime="#SBATCH --time=0-23:59"\n'
                               '    NewScratch="/scratch/"\n'
                               '    NewSoftware="/Software/"\n'
                               'fi\n'
    
                               f"echo {potential}:\nReRunList=({Trys_str})\n" +
                               'for dir in "${ReRunList[@]}"; do\n'
                               '    cd "$dir" || { echo "Directory $dir not found"; exit 1; }\n'
    
                               '    OldKeyword="--account"\n'
                               '    sed -i "/$OldKeyword/c' + fr"\\{''}" + '$NewAccount" "run_active2"\n'
                                                                           '    OldKeyword="--ntasks"\n'
                                                                           '    sed -i "/$OldKeyword/c' + fr"\\{''}" + '$NewNtask" "run_active2"\n'
                                                                                                                       '    OldKeyword="--mem-per-cpu"\n'
                                                                                                                       '    sed -i "/$OldKeyword/c' + fr"\\{''}" + '$NewMemPerCpu" "run_active2"\n'
                                                                                                                                                                   '    OldKeyword="--time"\n'
                                                                                                                                                                   '    sed -i "/$OldKeyword/c' + fr"\\{''}" + '$NewTime" "run_active2"\n'
                                                                                                                                                                                                               '    OldKeyword="/lustre04/scratch/"\n'
                                                                                                                                                                                                               '    sed -i "s|$OldKeyword|$NewScratch|g" "run_active2"\n'
                                                                                                                                                                                                               '    OldKeyword="/software/"\n'
                                                                                                                                                                                                               '    sed -i "s|$OldKeyword|$NewSoftware|Ig" "run_active2"\n'
                                                                                                                                                                                                               '    sbatch run_active2\n'
                                                                                                                                                                                                               '    cd ..\n'
                                                                                                                                                                                                               'done\n'
                               )
        # </editor-fold>
    else:
        print("!-Blank or Nonstop is not Empty")
    # </editor-fold>

# </editor-fold>

# <editor-fold desc="######################################## Energy per atom">
print("######################################## Energy per atom") # %%

# <editor-fold desc="***** Plot">
print("***** Plot")
Active = True
if Active:
    for Potential in Tde["Potential"].unique():
        DfFilteredPotential = Tde[Tde["Potential"] == Potential]
        Equilibrium = Thermal[(Thermal["Potential"] == Potential) &
                         (Thermal["File"] == "Equilibrium") &
                         (Thermal["LookupWord"] == "EnergyPerAtomAfter")]["EnergyPerAtom"].iloc[0]
        Thermalization = Thermal[(Thermal["Potential"] == Potential) &
                            (Thermal["File"] == "Thermalization") &
                            (Thermal["LookupWord"] == "EnergyPerAtomNve")]["EnergyPerAtom"].iloc[0]

        for File in Tde["File"].unique():
            if File in ["Tde", "Minimization"]:
                fig, axs = plt.subplots(figsize=(4.5, 4))
                Title = Potential + "-" + File
                print(Title)

                DfFilteredPotentialFile = DfFilteredPotential[DfFilteredPotential["File"] == File]

                hue_mapping = {
                    "EnergyPerAtomBefore": "Before",
                    "EnergyPerAtomAfter": "After"
                }
                DfFilteredPotentialFile = DfFilteredPotentialFile.copy()
                DfFilteredPotentialFile["LookupWord"] = DfFilteredPotentialFile["LookupWord"].map(hue_mapping)

                # Set seaborn style
                sns.set(style="whitegrid")

                # Create the plot
                g = sns.catplot(
                    data=DfFilteredPotentialFile,
                    estimator='mean',
                    x="File",
                    y="EnergyPerAtom",
                    hue="LookupWord",
                    col='DirectionMb',
                    col_wrap=5,
                    kind="bar",
                    height=2,
                    aspect=1.2,
                    legend=False
                )

                g.set_ylabels("")
                g.set_xlabels("")
                g.set_xticklabels("")

                # g._legend.set_title("")

                g.fig.text(-0.00, 0.5, "Energy Per Atom (eV)", va='center', rotation='vertical', fontsize=12)

                for ax in g.axes.flat:
                    ymin = Equilibrium-0.005
                    ymax = DfFilteredPotentialFile["EnergyPerAtom"].max() + 0.005
                    ax.set_ylim(ymin, ymax)
                    ax.axhline(y=Equilibrium, color='b', linestyle='-', label="Equilibrium")
                    ax.axhline(y=Thermalization, color='r', linestyle='-', label="Thermalization")
                    ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))

                g.set_titles("Direction: {col_name}")

                legend_elements = [
                    mpatches.Patch(color="b", label="Before"),
                    mpatches.Patch(color="peru", label="After"),
                    Line2D([0], [0], color='b', lw=2, label='Equilibrium'),
                    Line2D([0], [0], color='r', lw=2, label='Thermalization')
                ]

                plt.legend(handles=legend_elements, bbox_to_anchor=(-1.5, -0.1), loc='center', ncol=4, frameon=False,
                           prop={'weight': 'normal', 'size': 12})
                # plt.tight_layout()
                plt.subplots_adjust(left=0.070, bottom=None, right=None, top=None, wspace=None, hspace=None)
                plt.savefig(PythonName + "-" + Title, dpi=600, bbox_inches='tight')
                plt.show()
                plt.close()
# </editor-fold>

# </editor-fold>

# <editor-fold desc="######################################## Probability">
print("######################################## Probability") # %%
Active=False
if Active:

    # <editor-fold desc="Literature">
    Ackland = pd.read_csv("D:\Queens_University\Literature\Ackland.csv", usecols=["T", "Dir", "Complete"])
    Ackland = Ackland[Ackland['Complete'] == True]
    Griffiths = pd.read_csv("D:\Queens_University\Literature\Griffiths.csv", usecols=["T", "Dir"])
    # print(Griffiths)
    Ackland['DirectionMb'] = Ackland['Dir'].apply(new_column_name)
    Griffiths['DirectionMb'] = Griffiths['Dir'].apply(new_column_name)
    # </editor-fold>

    # <editor-fold desc="Plotting-TDE">
    Plotting = False
    if Plotting:
        Title = "DfTdeBox"
        Tde = Tde.copy()
        Tde['DirectionMb'] = Tde['DirectionMb'].apply(LatexFriendlyFunc)
        plt.figure(figsize=(10, 7))  # Your preferred figure size
        sns.boxplot(data=Tde, x='DirectionMb', y="Energy", hue="Potential")

        plt.xlabel("PKA Direction")
        plt.ylabel("PKA Energy (eV)")
        # plt.title("Boxplot of Pka Energy by Direction and Potential")

        plt.savefig(PythonName + "-" + Title, dpi=600, bbox_inches='tight')
        plt.show()
        plt.close()
    # </editor-fold>

    # <editor-fold desc="Plotting-TDE Min">
    print("Plotting-TDE Min")
    Plotting = False
    if Plotting:
        TdeGroupedPotentialDirectionDescribe = TdeGroupedPotentialDirectionDescribe.reset_index()
        TdeGroupedPotentialDirectionDescribe['DirectionMb'] = TdeGroupedPotentialDirectionDescribe['DirectionMb'].apply(LatexFriendlyFunc)

        Title = "TdeGroupedPotentialDirectionDescribe"
        plt.figure(figsize=(10, 6))
        sns.barplot(data=TdeGroupedPotentialDirectionDescribe, x='DirectionMb', y="min", hue="Potential")

        plt.xlabel('DirectionMb')
        plt.ylabel("Min PKA Energy (eV)")
        # plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(PythonName + "-" + Title, dpi=600, bbox_inches='tight')
        plt.show()
        plt.close()
    # </editor-fold>

    # <editor-fold desc="Plotting-Probability">
    Plotting = False
    if Plotting:
        for Name, Group in TdeGroupedPotentialDirection:
            print(Name)
            Potential = Name[0]
            Direction = Name[1]
            # print(Potential)
            # print(Group)
            # print(list(Group))
            GroupValues = Group.drop(columns=['Potential', 'Try', 'DirectionMb', 'File', 'LookupWord', 'EnergyPerAtom', 'Reparameterization']).values
            # print(GroupValues)
            GroupValuesFlatten = GroupValues.flatten()
            # print(GroupValuesFlatten)
            GroupValuesFlattenNp = np.asarray(GroupValuesFlatten)
            # print(GroupValuesFlattenNp)

            bar_gap = 1
            # Temperature ranges and corresponding values

            # plt.xticks(x_positions, temperature_ranges)
            temperature_ranges = ["$T_{d}$<21 eV", "21<=$T_{d}$<30 eV", "30<=$T_{d}$<50 eV", "50<=$T_{d}$<110 eV"]
            values = [0, 0.16, 0.55, 0.61]
            x_positions = np.arange(len(temperature_ranges)) * (1 + bar_gap)

            sns.kdeplot(data=GroupValuesFlattenNp,
                        cumulative=True, lw=3.5,
                        color='r',
                        label=Potential)

            plt.step([0, 21, 30, 50, 110, 120], [0, 0.16, 0.55, 0.61, 1, 1], 'black', linewidth=3.5, where='post',
                     label="Biget")
            plt.xticks([21, 30, 50, 110], ["21", "30", "50", "110"])

            # ["$T_{d}$<21 eV", "21<=$T_{d}$<30 eV", "30<=$T_{d}$<50 eV", "50<=$T_{d}$<110 eV"])
            # Create a bar plot
            # plt.bar(x_positions, values, color='blue', width=0.5, edgecolor='black')

            plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:.0%}'))
            # Add labels and title
            # plt.title("Probability Distribution (tabGAP1)")
            plt.xlabel(r"$T_{d} (eV)$")
            plt.ylabel("Cumulative distribution probability")
            plt.ylim(0, 1)  # Set y-axis limits
            # plt.xlim(0, 110)  # Set x-axis limits
            plt.legend()
            plt.grid(False)
            plt.tight_layout()
            plt.savefig(PythonName + "-" + str(Potential) + "-" + str(Direction) + "-" + "Prob.jpg", bbox_inches='tight')
            plt.show()
            plt.close()
    # </editor-fold>

    # <editor-fold desc="Plotting-Probability-PotentialDirection">
    Plotting = False
    if Plotting:
        Iteration = 0
        for Potential, Group in TdeGroupedPotentialDirection:
            Iteration += 1
            # print(Potential)
            # print(Group)
            GroupValues = Group.drop(columns=['Potential', 'Try', 'DirectionMb', 'File', 'LookupWord', 'EnergyPerAtom', 'Reparameterization']).values
            # print(GroupValues)
            GroupValuesFlatten = GroupValues.flatten()
            # print(GroupValuesFlatten)
            GroupValuesFlattenNp = np.asarray(GroupValuesFlatten)
            # print(GroupValuesFlattenNp)

            sns.kdeplot(data=GroupValuesFlattenNp,
                        cumulative=True, lw=3.5,
                        # color=colors[Iteration],
                        label=Potential)

        bar_gap = 1
        # Temperature ranges and corresponding values

        # plt.xticks(x_positions, temperature_ranges)
        temperature_ranges = ["$T_{d}$<21 eV", "21<=$T_{d}$<30 eV", "30<=$T_{d}$<50 eV", "50<=$T_{d}$<110 eV"]
        values = [0, 0.16, 0.55, 0.61]
        x_positions = np.arange(len(temperature_ranges)) * (1 + bar_gap)
        plt.xticks([21, 30, 50, 110], ["21", "30", "50", "110"])

        # ["$T_{d}$<21 eV", "21<=$T_{d}$<30 eV", "30<=$T_{d}$<50 eV", "50<=$T_{d}$<110 eV"])
        # Create a bar plot
        # plt.bar(x_positions, values, color='blue', width=0.5, edgecolor='black')

        plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:.0%}'))
        # Add labels and title
        # plt.title("Probability Distribution (tabGAP1)")
        plt.xlabel(r"$T_{d} (eV)$")
        plt.ylabel("Cumulative distribution probability")
        plt.ylim(0, 1)  # Set y-axis limits
        # plt.xlim(0, 110)  # Set x-axis limits
        plt.legend()
        plt.grid(False)
        plt.step([0, 21, 30, 50, 110, 120], [0, 0.16, 0.55, 0.61, 1, 1], 'black', linewidth=3.5, where='post',
                 label="Biget")
        plt.savefig(PythonName + "-" + "All" + "-Prob.jpg", bbox_inches='tight')
        plt.show()
        plt.close()
    # </editor-fold>

    # <editor-fold desc="Plotting-Probability-PotentialListAll">
    print("Plotting-Probability-PotentialListAll")
    Plotting = False
    if Plotting:
        Iteration = 0
        PotentialColorMap = {}  # Mapping of potentials without 'R' to colors
        PotentialList = ["M2R", "M3R", "BMD192R", "M2", "M3", "BMD192", "TabGap1", "TabGap2"]
        # PotentialList = ["M2R", "M3R", "BMD192R", "TabGap1", "TabGap2"]
        for Name, Group in TdeGroupedPotential:
            print(Name)
            Potential = Name
            # Direction = Name[1]
            if Potential in PotentialList:
                # print(Name)
                # print(Group)

                potential_key = Potential.rstrip('R')  # Remove 'R' from potential
                if potential_key not in PotentialColorMap:
                    PotentialColorMap[potential_key] = colors[Iteration]

                GroupValues = Group.drop(columns=['Potential', 'Try', 'DirectionMb', 'File', 'LookupWord', 'EnergyPerAtom', 'Reparameterization']).values
                # print(GroupValues)
                GroupValuesFlatten = GroupValues.flatten()
                # print(GroupValuesFlatten)
                GroupValuesFlattenNp = np.asarray(GroupValuesFlatten)
                # print(GroupValuesFlattenNp)

                LineStyle = '-' if Name[-1] != 'R' else 'dotted'

                sns.kdeplot(data=GroupValuesFlattenNp,
                            cumulative=True, lw=2.5,
                            color=PotentialColorMap[potential_key],
                            label=Name,
                            linestyle=LineStyle,
                            )
                Iteration += 1
        plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:.0%}'))

        plt.step([0, 21, 30, 50, 110, 120], [0, 0.16, 0.55, 0.61, 1, 1], 'black', linewidth=3.5, where='post', label="Biget")
        bar_gap = 1
        # Temperature ranges and corresponding values

        # plt.xticks(x_positions, temperature_ranges)
        temperature_ranges = ["$T_{d}$<21 eV", "21<=$T_{d}$<30 eV", "30<=$T_{d}$<50 eV", "50<=$T_{d}$<110 eV"]

        FontSize = 8
        SpaceH = 2

        plt.text(21 + SpaceH, 0, 'P$_\mathregular{d}$(T) = 0', verticalalignment='bottom', horizontalalignment='left',
                 fontsize=FontSize, fontstyle='italic')
        plt.text(30 + SpaceH, 0.16, 'P$_\mathregular{d}$(T) = 0.16', verticalalignment='bottom', horizontalalignment='left',
                 fontsize=FontSize, fontstyle='italic')
        plt.text(50 + SpaceH, 0.55, 'P$_\mathregular{d}$(T) = 0.55', verticalalignment='bottom', horizontalalignment='left',
                 fontsize=FontSize, fontstyle='italic')
        plt.text(110 + SpaceH, 0.61, 'P$_\mathregular{d}$(T) = 0.61', verticalalignment='bottom',
                 horizontalalignment='left', fontsize=FontSize, fontstyle='italic')

        # Add labels and title
        # plt.title("Probability Distribution (tabGAP1)")
        plt.xlabel(r"$T_{d} (eV)$")
        plt.ylabel("Cumulative distribution probability")
        plt.ylim(0, 1)  # Set y-axis limits
        plt.xlim(0, 200)  # Set x-axis limits
        plt.legend()
        plt.grid(False)
        plt.tight_layout()
        plt.savefig(PythonName + "-" + "All" + "-Prob.jpg", bbox_inches='tight')
        plt.show()
        plt.close()
    # </editor-fold>

    # <editor-fold desc="Plotting-Probability-PotentialListReparamterized">
    print("Plotting-Probability-PotentialListReparamterized")
    Plotting = False
    if Plotting:
        Iteration = 0
        PotentialList = ["M2R", "M3R", "BMD192R", "M2", "M3", "BMD192", "TabGap1", "TabGap2"]
        # PotentialList = ["M2R", "M3R", "BMD192R", "TabGap1", "TabGap2"]
        ColorsDic = {
            "M2": Colors[0],
            "M2R": Colors[0],
            "M3": Colors[1],
            "M3R": Colors[1],
            "BMD192": Colors[2],
            "BMD192R": Colors[2],
            "TabGap1": Colors[3],
            "TabGap2": Colors[4],
        }
        for Name, Group in TdeGroupedPotential:
            # print(Name)
            # print(Group.columns.tolist())
            Potential = Name
            # Direction = Name[1]
            if Potential in PotentialList:
                print(Name)
                print(Group.columns.tolist())
                # print(Group)
                GroupDropped = Group.drop(columns=['Potential', 'Try', 'DirectionMb', 'Direction',
                                                  'File', 'LookupWord', 'EnergyPerAtom', 'Reparameterization',
                                                  'X', 'Y', 'Z', 'R',
                                                  'UnitX', 'UnitY', 'UnitZ',
                                                  'Theta', 'Phi',
                                                  'a1', 'a2', 'a3', 'h'])
                print(GroupDropped.columns.tolist())
                GroupValues = GroupDropped.values
                # print(GroupValues)

                GroupValuesFlatten = GroupValues.flatten()
                # print(GroupValuesFlatten)
                GroupValuesFlattenNp = np.asarray(GroupValuesFlatten)
                # print(GroupValuesFlattenNp)

                LineStyle = '--' if Name[-1] != 'R' else '-'

                sns.kdeplot(data=GroupValuesFlattenNp,
                            cumulative=True, lw=2,
                            color= ColorsDic[Potential],
                            label=Name,
                            linestyle=LineStyle,
                            )
                Iteration += 1
        plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:.0%}'))

        plt.step([0, 21, 30, 50, 110, 150], [0, 0.16, 0.55, 0.61, 1, 1], 'black', linewidth=3.5, where='post', label="Biget")
        bar_gap = 1
        # Temperature ranges and corresponding values

        # plt.xticks(x_positions, temperature_ranges)
        temperature_ranges = ["$T_{d}$<21 eV", "21<=$T_{d}$<30 eV", "30<=$T_{d}$<50 eV", "50<=$T_{d}$<110 eV"]

        FontSize = 18
        SpaceH = 2

        # plt.text(21 + SpaceH, 0, 'P$_\mathregular{d}$(T) = 0', verticalalignment='bottom', horizontalalignment='left',
        #          fontsize=FontSize, fontstyle='italic')
        # plt.text(30 + SpaceH, 0.16, 'P$_\mathregular{d}$(T) = 0.16', verticalalignment='bottom', horizontalalignment='left',
        #          fontsize=FontSize, fontstyle='italic')
        # plt.text(50 + SpaceH, 0.55, 'P$_\mathregular{d}$(T) = 0.55', verticalalignment='bottom', horizontalalignment='left',
        #          fontsize=FontSize, fontstyle='italic')
        # plt.text(110 + SpaceH, 0.61, 'P$_\mathregular{d}$(T) = 0.61', verticalalignment='bottom',
        #          horizontalalignment='left', fontsize=FontSize, fontstyle='italic')

        plt.text(20 + SpaceH, 0, '%0', verticalalignment='bottom', horizontalalignment='left',
                 fontsize=FontSize, fontstyle='italic',color='black')
        plt.text(30 + SpaceH, 0.16, '%16', verticalalignment='bottom', horizontalalignment='left',
                 fontsize=FontSize, fontstyle='italic')
        plt.text(50 + SpaceH, 0.53, '%55', verticalalignment='bottom', horizontalalignment='left',
                 fontsize=FontSize, fontstyle='italic')
        plt.text(110 + SpaceH, 0.61, '%61', verticalalignment='bottom',
                 horizontalalignment='left', fontsize=FontSize, fontstyle='italic')

        # Add labels and title
        # plt.title("Probability Distribution (tabGAP1)")
        # plt.xlabel(r"$T_{d} (eV)$")
        plt.xlabel("TDE (eV)")
        plt.ylabel("Cumulative distribution probability")
        plt.ylim(0, 1)  # Set y-axis limits
        plt.xlim(0, 200)  # Set x-axis limits
        # plt.xticks([0,21,30,50,110])
        handles, labels = plt.gca().get_legend_handles_labels()
        # order = [1, 2, 0, 3, 4, 5] # when no original version is added
        order = [0, 1, 2, 3, 4, 5, 6, 7] # when original versions are added
        plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order],frameon=False, fontsize=15, loc='center',bbox_to_anchor=(0.75, 0.27),labelspacing=0.2 )
        plt.grid(False)
        plt.tight_layout()
        plt.savefig(PythonName + "-" + "All" + "-Prob.png", bbox_inches='tight')
        plt.show()
        plt.close()
    # </editor-fold>

    # <editor-fold desc="Plotting-Probability-ecdfplot">
    Plotting = False
    if Plotting:
        sns.ecdfplot(data=DfBravaisMelt, x="Tde",
                     # cumulative=True,
                     hue="Potential",
                     lw=3.5,
                     legend=True,
                     # color=colors[Iteration],
                     # label="Potential"
                     )

        plt.step([0, 21, 30, 50, 110, 120], [0, 0.16, 0.55, 0.61, 1, 1], 'black', linewidth=3.5, where='post',
                 label="Biget")

        # bar_gap = 1
        # Temperature ranges and corresponding values

        # plt.xticks(x_positions, temperature_ranges)
        # temperature_ranges = ["$T_{d}$<21 eV", "21<=$T_{d}$<30 eV", "30<=$T_{d}$<50 eV", "50<=$T_{d}$<110 eV"]
        # values = [0, 0.16, 0.55, 0.61]
        # x_positions = np.arange(len(temperature_ranges)) * (1 + bar_gap)
        # plt.xticks([21, 30, 50, 110], ["21", "30", "50", "110"])

        # ["$T_{d}$<21 eV", "21<=$T_{d}$<30 eV", "30<=$T_{d}$<50 eV", "50<=$T_{d}$<110 eV"])
        # Create a bar plot
        # plt.bar(x_positions, values, color='blue', width=0.5, edgecolor='black')

        plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:.0%}'))
        # Add labels and title
        # plt.title("Probability Distribution (tabGAP1)")
        plt.xlabel(r"$T_{d} (eV)$")
        plt.ylabel("Cumulative distribution probability")
        plt.ylim(0, 1)  # Set y-axis limits
        # plt.xlim(0, 110)  # Set x-axis limits
        plt.legend()
        plt.grid(False)

        plt.savefig(PythonName + "-" + "All" + "-Prob.jpg", bbox_inches='tight')
        plt.show()
        plt.close()
    # </editor-fold>

    # <editor-fold desc="Plotting-Probability-ecdf">
    Plotting = False
    if Plotting:
        sns.ecdfplot(data=DfBravaisMelt, x="Tde",
                     # cumulative=True,
                     hue="Potential",
                     lw=3.5,
                     legend=True,
                     # color=colors[Iteration],
                     # label=Potential
                     )

        plt.step([0, 21, 30, 50, 110, 120], [0, 0.16, 0.55, 0.61, 1, 1], 'black', linewidth=3.5, where='post',
                 label="Biget")

        bar_gap = 1
        # Temperature ranges and corresponding values

        # plt.xticks(x_positions, temperature_ranges)
        # temperature_ranges = ["$T_{d}$<21 eV", "21<=$T_{d}$<30 eV", "30<=$T_{d}$<50 eV", "50<=$T_{d}$<110 eV"]
        # values = [0, 0.16, 0.55, 0.61]
        # x_positions = np.arange(len(temperature_ranges)) * (1 + bar_gap)
        # plt.xticks([21, 30, 50, 110], ["21", "30", "50", "110"])

        # ["$T_{d}$<21 eV", "21<=$T_{d}$<30 eV", "30<=$T_{d}$<50 eV", "50<=$T_{d}$<110 eV"])
        # Create a bar plot
        # plt.bar(x_positions, values, color='blue', width=0.5, edgecolor='black')

        plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:.0%}'))
        # Add labels and title
        # plt.title("Probability Distribution (tabGAP1)")
        plt.xlabel(r"$T_{d} (eV)$")
        plt.ylabel("Cumulative distribution probability")
        plt.ylim(0, 1)  # Set y-axis limits
        # plt.xlim(0, 110)  # Set x-axis limits
        plt.legend()
        plt.grid(False)

        plt.savefig(PythonName + "-" + "All" + "-Prob.jpg", bbox_inches='tight')
        plt.show()
        plt.close()
    # </editor-fold>

    # <editor-fold desc="Table">
    Active = False
    if Active:
        data = pd.DataFrame(columns=["Potential", 'DirectionMb', "Value"])
        for column in Tde.columns[:-1]:
            # essentially just the new_column_name function from your code
            indices = column.split("_")[1:]
            indSubStrs = [r"$\bar{" + index[1:] + "}$" if "-" in index else index for index in indices]
            indStr = "".join(indSubStrs)

            pots = Tde["Potential"]
            vals = Tde[column]
            lenVals = len(vals)

            dataNew = pd.DataFrame({
                "Potential": pots,
                'DirectionMb': pd.Series([indStr] * lenVals),
                "Value": vals
            })
            data = pd.concat([data, dataNew], ignore_index=True, axis=0)
        print(data)
        # mostly for debugging you can get rid of this
        # data.to_csv("20240516-transform.csv")

        # %% plot
        sns.set_theme(style="whitegrid", rc={"figure.figsize": (25, 5), "font.sans-serif": "DejaVu Sans"})
        axes = sns.boxplot(data=data, x='DirectionMb', y="Value", hue="Potential", gap=.1)
        sns.move_legend(axes, "lower center", bbox_to_anchor=(.5, 1), ncol=len(Tde.columns[:-1]))

        fig = axes.get_figure()
        fig.savefig("complete_maxval_plot.png")

        # %% stats
        statDf = pd.DataFrame(columns=["Potential",
                                       "T_d,avg", "T_d,avg Error",
                                       "Tavg_d,avg", "Tavg_d,avg Error",
                                       "Tmed_d,avg", "Tmed_d,avg Error"])

        for pot in data["Potential"].unique():
            potData = data[data["Potential"] == pot]
            cumMean = np.nanmean(potData["Value"])
            cumStErr = stats.sem(potData["Value"], nan_policy="omit")

            pdMeans = np.array([])
            pdMeds = np.array([])

            for dir in potData['DirectionMb'].unique():
                pdData = potData[potData['DirectionMb'] == dir]
                pdMean = np.nanmean(pdData["Value"])
                pdMed = np.nanmedian(pdData["Value"])

                pdMeans = np.append(pdMeans, pdMean)
                pdMeds = np.append(pdMeds, pdMed)

            meanOfMeans = np.nanmean(pdMeans)
            mMeanStErr = stats.sem(pdMeans, nan_policy="omit")

            meanOfMedians = np.nanmean(pdMeds)
            mMedStErr = stats.sem(pdMeds, nan_policy="omit")

            dfTemp = pd.DataFrame({
                "Potential": [pot],
                "T_d,avg": [cumMean],
                "T_d,avg Error": [cumStErr],
                "Tavg_d,avg": [meanOfMeans],
                "Tavg_d,avg Error": [mMeanStErr],
                "Tmed_d,avg": [meanOfMedians],
                "Tmed_d,avg Error": [mMedStErr]
            })

            statDf = pd.concat([statDf, dfTemp], ignore_index=True, axis=0)

        statDf.to_csv("20240516-stats.csv")

    # </editor-fold>
# </editor-fold>

# <editor-fold desc="######################################## Stereograph">
print("######################################## Stereograph") # %%

# <editor-fold desc="***** Plot-stereograph-contour-TdeGroupedPotential">
print("***** Plot-stereograph-contour-TdeGroupedPotential")
Active = False
if Active:
    PotentialList = ["M2R", "M3R", "BMD192R", "M2", "M3", "BMD192", "TabGap1", "TabGap2"]

    for Key, Content in TdeDefectGroupedPotentialDic.items():
        print(Key)
        for Name, Group in Content:
            Potential = Name[0]
            if Potential in PotentialList:
                print(Potential)
                # print(Group)
                # print(Group["Direction"])
                Title = "TdeDefectGroupedPotentialDic" + "-" + Key + "-" + Potential
                fig, ax = mplstereonet.subplots(figsize=(10, 8))

                Phi = Group["Phi"] + 90
                Theta = Group["Theta"]
                X = Group["X"]
                Y = Group["Y"]
                Z = Group["Z"]

                GroupGroupedCount = Group.groupby('Direction')["Energy"].describe().reset_index()["min"]
                # print(GroupGroupedCount)
                IntensityMin = GroupGroupedCount.min()
                IntensityMax = GroupGroupedCount.max()
                # print(IntensityMax, IntensityMin)

                cax = ax.density_contourf(
                    Phi, Theta,
                    measurement='poles',
                    cmap='jet',
                    # vmin=IntensityMin, vmax=IntensityMax  # set color range here
                    levels=Levels  # <--- ADD THIS LINE
                )

                # ---- Scatter plot the points ----
                Xs, Ys = [], []
                for phi, theta in zip(Phi, Theta):
                    x_point, y_point = mplstereonet.stereonet_math.pole(phi, theta)
                    Xs.append(x_point)
                    Ys.append(y_point)
                ax.scatter(Xs, Ys, color='black', s=20, alpha=0.7, label='_nolegend_')

                # ---- Optional: annotate points ----
                for x_point, y_point, x, y, z in zip(Xs, Ys, X, Y, Z):
                    ax.text(x_point, y_point+0.1, f"[{x:.0f},{y:.0f},{z:.0f}]",
                            fontsize=20, ha='center', va='center')

                ax.set_azimuth_ticks([])
                ax.grid(False)

                cbar = fig.colorbar(cax)
                # cbar.set_label("Pole Density")

                plt.savefig(PythonName + "-" + Title, dpi=600, bbox_inches='tight')
                plt.show()
                plt.close()
# </editor-fold>
# os.system("pause")
# <editor-fold desc="***** Plot-stereograph-points-TdeGroupedPotentialDirectionDescribe">
print("***** Plot-stereograph-points-TdeGroupedPotentialDirectionDescribe")
Active = False
if Active:
    for Key, Content in TdeDefectGroupedPotentialDirectionDescribeDic.items():
        print(Key)
        for Potential, group in Content.groupby('Potential'):
            print(Potential, group[["Direction","min"]])
            Title = "Stereograph" + "-" + Key + "-" + Potential + "-" + str(Based)
            r = np.sqrt(1 / (1 + group['UnitZ']))
            x_proj = r * group['UnitX']
            y_proj = r * group['UnitY']

            plt.figure(figsize=(8,8))
            sc = plt.scatter(x_proj, y_proj, c=group['mean'], cmap='viridis', s=100, edgecolor='k')

            circle = plt.Circle((0, 0), 1, color='black', fill=False)
            plt.gca().add_artist(circle)

            plt.colorbar(sc, label='Mean Intensity')
            plt.xlim([-1, 1])
            plt.ylim([-1, 1])
            plt.gca().set_aspect('equal')
            plt.xlabel('X')
            plt.ylabel('Y')
            # plt.title(f'Polar Stereographic Projection - {Potential}')
            plt.grid(True)
            plt.savefig(PythonName + "-" + Title, dpi=600, bbox_inches='tight')
            plt.show()
# </editor-fold>

# <editor-fold desc="***** Repeat">
print("***** Repeat")

Active = False
if Active:
    BasedList = ["min","mean"]

    for Based in BasedList:
        # <editor-fold desc="Repeat">
        print(Based)

        for Key, Content in TdeDefectGroupedPotentialDirectionDescribeDic.items():
            print(Key)
            # print(Content)
            Content['Repeating'] = Content[Based].round().astype(int)

            # print(Content['Repeating'])
            # os.system("pause")

            DfTdeRepeated = Content.loc[Content.index.repeat(Content['Repeating'])].drop(columns=['Repeating','count','std','25%','50%','75%']).reset_index(drop=True)
            # print(DfTdeRepeated)
            # </editor-fold>

            # <editor-fold desc="Export">
            print("Export")
            # DfTdeRepeated.to_csv(PythonName + "-DfTdeRepeated-" + str(Based) + ".csv", index=False)
            # DfTdeRepeated = pd.read_csv(PythonName + "-DfTdeRepeated" + ".csv")
            # </editor-fold>

            # <editor-fold desc="Definition">
            print("Definition")
            def stereographic_projection(v):
                """Project a 3D vector onto a 2D stereographic plane."""
                v = np.array(v)
                v = v / np.linalg.norm(v)  # Normalize
                x, y, z = v
                X = x / (1 + z)
                Y = y / (1 + z)
                return X, Y

            DirectionList = [
                [1,2,1],
                [0,1,0],
                [0,1,1],
                [0,3,1],
               [-1,1,2],
                [-1,1,0],
                [-1,1,1],
                [1,2,2],
                [0,0,1],
                [1,2,0],
            ]
            # </editor-fold>

            # <editor-fold desc="Plot-Repeat">
            print("Plot-Repeat")
            Active = False
            if Active:
                LetterLabels = list("abcdefghijklmnopqrstuvwxyz")
                LabelCounter = 0
                # for potential, group in DfTdeRepeated.groupby('Potential'):
                for Potential, group in DfTdeRepeated.groupby('Potential'):
                    Title = "ReversePoleFigure" + "-" + Potential + "-" + str(Based)
                    # print(Potential, group[["Direction","DirectionMb","UnitX","UnitY","UnitZ"]])

                    # Create Vector3d objects for UnitX, UnitY, UnitZ
                    Vectors = Vector3d(group[['UnitX', 'UnitY', 'UnitZ']].to_numpy())
                    # print(VectorsNumpy)

                    # Set the crystal symmetry
                    symmetry = C1

                    # Set the reference direction (Z-axis)
                    direction = Vector3d.zvector()


                    fig, ax = plt.subplots(figsize=(9, 4), subplot_kw={
                        "projection": "ipf",
                        "symmetry": symmetry,
                        "direction": direction
                    })

                    ax.scatter(Vectors, s=50, alpha=0.7, color="black")

                    # for idx, vector in enumerate(DirectionList):
                    #     X, Y = stereographic_projection(vector)
                    #     ax.scatter(X, Y, s=50, alpha=0.7, color="red")
                    #     ax.annotate(str(vector), (X, Y), textcoords="offset points", xytext=(5, 5), ha='left')

                    GroupGroupedCount = group.groupby('Direction')['Direction'].describe().reset_index()["count"]
                    # print(GroupGroupedCount)
                    IntensityMin = GroupGroupedCount.min()
                    IntensityMax = GroupGroupedCount.max()
                    # print(IntensityMin,IntensityMax)
                    IntensityMin = 0
                    IntensityMax = 1

                    group['GroupCount'] = group.groupby('Direction')['Direction'].transform('count')

                    # Then:
                    weights = group['GroupCount'].values

                    ax.pole_density_function(
                        Vectors,
                        weights=weights,
                        cmap="viridis",
                        vmin=IntensityMin,
                        vmax=IntensityMax,
                    )
                    fig.suptitle("(" + str(LetterLabels[LabelCounter]) + "): " + str(Potential), fontsize=24, y=1.00)
                    ax.set_title(Potential)
                    plt.tight_layout()
                    plt.show()
            # </editor-fold>

            # <editor-fold desc="Plot-Repeat-Stereograph">
            print("Plot-Repeat-Stereograph")
            Active = True
            if Active:
                PotentialList = ["M2R", "M3R", "BMD192R", "M2", "M3", "BMD192", "TabGap1", "TabGap2"]
                for Name, Group in DfTdeRepeated.groupby('Potential'):
                    print(Name)
                    # print(Group["Direction"])
                    Potential = Name

                    LetterLabels = list("abcdefghijklmnopqrstuvwxyz")
                    LabelCounter = 0
                    if Potential in PotentialList:
                        Title = "Stereograph" + "-" + Potential + "-" + str(Based)
                        fig, ax = mplstereonet.subplots(figsize=(10, 8))

                        Phi = Group["Phi"] + 90
                        Theta = Group["Theta"]
                        X = Group["X"]
                        Y = Group["Y"]
                        Z = Group["Z"]

                        GroupAgg = Group[Based].agg(['min', 'max', 'mean'])
                        # print(GroupAgg)
                        GroupAggMin = GroupAgg.loc["min"]
                        # GroupAggMin = 0
                        GroupAggMax = GroupAgg.loc["max"]+20
                        # GroupAggMax = np.percentile(Group[Based], 99)  # Use 99th percentile to exclude extreme peaks
                        # GroupAggMean = Group[Based].mean()
                        # GroupAggStd = Group[Based].std()
                        # GroupAggMax = max(GroupAgg.loc["max"], GroupAggMean + 2 * GroupAggStd)
                        GroupAggMean = GroupAgg.loc["mean"]
                        # print(GroupAggMin, GroupAggMax, GroupAggMean)

                        import matplotlib.colors as mcolors
                        NormalColors = mcolors.Normalize(vmin=GroupAggMin, vmax=GroupAggMax)
                        Levels = np.linspace(GroupAggMin, GroupAggMax, 200)

                        cax = ax.density_contourf(
                            Phi, Theta,
                            measurement='poles',
                            cmap='rainbow',
                            # vmin=GroupAggMin, vmax=GroupAggMax,
                            norm=NormalColors,
                            levels=Levels  # <--- ADD THIS LINE
                        )

                        Xs, Ys = [], []
                        for phi, theta in zip(Phi, Theta):
                            x_point, y_point = mplstereonet.stereonet_math.pole(phi, theta)
                            Xs.append(x_point)
                            Ys.append(y_point)

                        ax.scatter(
                            Xs, Ys,
                            edgecolors='black',  # Note: plural 'edgecolors'
                            facecolors='none',  # Empty center
                            s=20,
                            alpha=0.7,
                            label='_nolegend_'
                        )

                        for x_point, y_point, x, y, z in zip(Xs, Ys, X, Y, Z):
                            ax.text(x_point, y_point+0.1, f"[{x:.0f},{y:.0f},{z:.0f}]",
                                    fontsize=10, ha='center', va='center')

                        ax.set_azimuth_ticks([])
                        ax.grid(False)

                        cbar = fig.colorbar(cax)
                        cbar.locator = mticker.MaxNLocator(nbins=5)
                        cbar.update_ticks()
                        cbar.ax.yaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f'))
                        cbar.set_label(str(Based) + " of PKA energy (eV)")


                        # plt.tight_layout()
                        fig.suptitle("(" + str(LetterLabels[LabelCounter]) + "): " +  str(Potential), fontsize=24, y=1.00)
                        plt.savefig(PythonName + "-" + Title, dpi=600, bbox_inches='tight')
                        # plt.show()
                        plt.close()
                        LabelCounter +=1
            # </editor-fold>

# </editor-fold>

# </editor-fold>

# <editor-fold desc="######################################## Inverse Pole Figure">
print("######################################## Inverse Pole Figure") # %%

# <editor-fold desc="***** Plot-inverse pole figures-example-TdeGroupedPotential">
print("***** Plot-inverse pole figures-example-TdeGroupedPotential")
Active = False
if Active:
    from orix import data, plot
    from orix.vector import Vector3d
    for potential, group in TdeGroupedPotential:
        print(potential, group)


        xmap = data.sdss_ferrite_austenite(allow_download=True)
        print(xmap)

        # Extract orientations, O
        pg_m3m = xmap.phases[1].point_group.laue
        print(pg_m3m)
        O_fe = xmap["ferrite"].orientations
        # O_fe = group[['UnitX','UnitY','UnitZ','min']].to_numpy()
        print(O_fe)
        with open('O_fe', 'w') as outfile:
            outfile.write('\n'.join(str(i) for i in O_fe))


        # Some sample direction, v
        v = Vector3d([0, 0, 1])
        v_title = "Z"

        # Rotate sample direction v into every crystal orientation O
        t_fe = O_fe * v
        with open('t_fe', 'w') as outfile:
            outfile.write('\n'.join(str(i) for i in t_fe))

        # Set IPDF range
        vmin, vmax = (0, 3)

        subplot_kw = {"projection": "ipf", "symmetry": pg_m3m, "direction": v}
        fig = plt.figure(figsize=(9, 8))

        ax0 = fig.add_subplot(221, **subplot_kw)
        ax0.scatter(O_fe, alpha=0.05)
        ax2 = fig.add_subplot(223, **subplot_kw)
        ax2.pole_density_function(t_fe, vmin=vmin, vmax=vmax)


        plt.show()

# </editor-fold>

# <editor-fold desc="***** Plot-inverse pole figures-contour and points-TdeGroupedPotentialDirectionDescribe">
Active = False
if Active:
    print("***** Plot-inverse pole figures-contour and points-TdeGroupedPotentialDirectionDescribe")
    BasedList = ["min", "mean"]
    for Based in BasedList:
        print("---" + str(Based))

        for Key, Content in TdeDefectGroupedPotentialDirectionDescribeDic.items():
            print("--" + str(Key))
            # print(Content)

            # <editor-fold desc="Repeat">
            # print(TdeGroupedPotentialDirectionDescribe)
            Content['Repeating'] = Content[Based].round().astype(int)

            DfTdeRepeated = Content.loc[Content.index.repeat(Content['Repeating'])
            ].drop(columns=['Repeating', 'count', 'std', '25%', '50%', '75%']).reset_index(drop=True)
            # print(DfTdeRepeated)
            # </editor-fold>

            # <editor-fold desc="Export">
            # print("Export")
            # DfTdeRepeated.to_csv(PythonName + "-DfTdeRepeated-" + str(Based) + ".csv", index=False)
            # DfTdeRepeated = pd.read_csv(PythonName + "-DfTdeRepeated" + ".csv")
            # </editor-fold>

            # <editor-fold desc="Plot-inverse pole figures-contour and points-TdeGroupedPotentialDirectionDescribe-Pole Density">
            Active = False
            if Active:
                print("Plot-inverse pole figures-contour and points-TdeGroupedPotentialDirectionDescribe-Pole Density")
                LetterLabels = list("abcdefghijklmnopqrstuvwxyz")
                LabelCounter = 0

                # for Potential, Group in DfTdeRepeated.groupby('Potential'):
                for Potential, Group in Content.groupby('Potential'):
                    Title = "PoleFigureDensity" + "-" + Potential + "-" + str(Based)

                    # <editor-fold desc="Reflection">
                    Reflection = True
                    if Reflection:
                        MirrorMatrix = np.array([
                            [-1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]
                        ])

                        OriginalStack = np.column_stack([
                            Group['UnitX'].values,
                            Group['UnitY'].values,
                            Group['UnitZ'].values
                        ])

                        MirroredStack = np.array([
                            MirrorMatrix @ vec if vec[0] > 0 else vec
                            for vec in OriginalStack
                        ])

                        Group["TransformedX"] = MirroredStack[:, 0]
                        Group["TransformedY"] = MirroredStack[:, 1]
                        Group["TransformedZ"] = MirroredStack[:, 2]
                    else:
                        Group["TransformedX"] = Group['UnitX']
                        Group["TransformedY"] = Group['UnitY']
                        Group["TransformedZ"] = Group['UnitZ']
                    # </editor-fold>

                    # <editor-fold desc="Rotation">
                    Rotating = True
                    if Rotating:
                        RotationDeg = -90
                        RotationRad = np.deg2rad(RotationDeg)
                        # RotationRad = np.deg2rad(-1.1)
                        print(RotationRad)

                        RotationMatrix = np.array([
                            [np.cos(RotationRad), -np.sin(RotationRad), 0],
                            [np.sin(RotationRad), np.cos(RotationRad), 0],
                            [0, 0, 1]
                        ])

                        OriginalStack = np.column_stack([
                            Group['TransformedX'].values,
                            Group['TransformedY'].values,
                            Group['TransformedZ'].values
                        ])
                        RotatedStack = OriginalStack @ RotationMatrix.T
                        Group["TransformedX"] = RotatedStack[:, 0]
                        Group["TransformedY"] = RotatedStack[:, 1]
                        Group["TransformedZ"] = RotatedStack[:, 2]
                    else:
                        RotationMatrix = np.array([
                            [1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]
                        ])

                    UnitStack = np.column_stack([
                        Group['TransformedX'].values,
                        Group['TransformedY'].values,
                        Group['TransformedZ'].values
                    ])
                    Vectors = Vector3d(UnitStack)
                        # print(Potential, Group[["Direction", "DirectionMb", "X", "Y", "Z", "UnitX", "UnitY", "UnitZ"]])

                    # print(Potential, Group[
                    #     ["Direction", "DirectionMb", "X", "Y", "Z", "UnitX", "UnitY", "UnitZ", 'RotatedX', 'RotatedY',
                    #      'RotatedZ']])
                    print(Potential, Group[[
                        "Direction", "DirectionMb" ,"min" ,"UnitX" ,"UnitY" ,"UnitZ"
                    ]])
                    # </editor-fold>

                    SymmetryList = [
                        C1,
                        # Ci, C2, C2h, D2, D2h,
                        # C4, C4h, D4,
                        D4h,
                        # C3, D3, D3d,
                        # C6, D6,
                        D6h,
                        # T, Th, O, Td, Oh
                    ]
                    for Symmerty in SymmetryList:
                        print(Symmerty.name)
                        # Set default font color and size
                        plt.rcParams['text.color'] = 'black'
                        # plt.rcParams['axes.labelcolor'] = 'black'
                        # plt.rcParams['axes.titlesize'] = 20  # Optional: global title size
                        # plt.rcParams['axes.titlecolor'] = 'black'
                        # plt.rcParams['xtick.color'] = 'white'
                        # plt.rcParams['ytick.color'] = 'white'
                        # plt.rcParams['xtick.labelsize'] = 1
                        # plt.rcParams['ytick.labelsize'] = 1
                        # plt.rcParams['legend.edgecolor'] = 'white'
                        # plt.rcParams['legend.labelcolor'] = 'white'
                        # # plt.rcParams['figure.facecolor'] = 'black'  # Optional: dark background
                        # # plt.rcParams['axes.facecolor'] = 'black'  # Optional: dark plot area
                        # # plt.rcParams['savefig.facecolor'] = 'black'  # Optional: dark saved image
                        # plt.rcParams['font.size'] = 1

                        # print(UnitStack)
                        # print(Vectors)

                        # <editor-fold desc="Exporting Vectors">
                        Exporting = False
                        if Exporting:
                            Vectors_array = Vectors.data  # shape: (N, 3)
                            filename = f"Vectors_{Potential}.txt"
                            np.savetxt(filename, Vectors_array, fmt="%.6f", header="X Y Z", comments='')
                        # </editor-fold>

                        RefDirection = Vector3d.zvector()

                        fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={
                            "projection": "ipf",
                            "symmetry": Symmerty,  # D4h
                            "direction": RefDirection
                        })

                        # print(Vectors)
                        ax.scatter(Vectors, s=50, alpha=0.7)

                        # <editor-fold desc="Range">
                        RangeType = "Manual"
                        if RangeType == "Agg":
                            GroupAgg = Group[Based].agg(['min', 'max', 'mean'])
                            # print(GroupAgg)
                            RangeMin = GroupAgg.loc["min"]
                            RangeMax = GroupAgg.loc["max"]
                            GroupAggMean = GroupAgg.loc["mean"]
                            # print(Group[Based].values)
                        elif RangeType == "Direct":
                            RangeMin = Group["min"].min()
                            RangeMax = Group["min"].max()
                        elif RangeType == "Manual":
                            RangeMin = 0.5
                            RangeMax = 10
                        print(RangeMin, RangeMax)
                        NormalColors = mcolors.Normalize(vmin=RangeMin, vmax=RangeMax)
                        Levels = np.linspace(RangeMin, RangeMax, 100)
                        # </editor-fold>

                        # print(Group["min"].values)
                        # print(Group["Direction"].values)
                        ax.pole_density_function(
                            Vectors,
                            weights=Group["min"].values,
                            cmap='rainbow',
                            # vmin=RangeMin,
                            # vmax=RangeMax,
                            # log=True
                            norm=NormalColors,
                        )

                        # for vec, label in zip(Vectors, Group["DirectionMb"]):  # or any other column like group.index
                        #     xy = ax._projection.vector2xy(vec)
                        #     ax.annotate(str(label), xy=xy, fontsize=10, ha='center', va='center')

                        FontSize = 10
                        Color = "black"
                        FontWeight = "normal"
                        ZOrder = 10
                        Offset = 0.0

                        cbar = ax.figure.axes[-1]  # Typically the last axis is the colorbar
                        cbar.set_ylabel("TDE (eV)", color=Color, fontsize=FontSize)  # Change text, color, size
                        cbar.tick_params(labelsize=FontSize)

                        # <editor-fold desc="Annotation">
                        Annotation = True
                        if Annotation:
                            Annotation = [
                                [-1 + Offset, 1 + Offset, 0 + Offset],
                                [-1 + Offset, 1 + Offset, 1 + Offset],
                                [-1 + Offset, 1 + Offset, 2 + Offset],
                                [0 + Offset, 0 + Offset, 1 + Offset],
                                [0 + Offset, 1 + Offset, 0 + Offset],
                                [0 + Offset, 1 + Offset, 1 + Offset],
                                [0 + Offset, 3 + Offset, 1 + Offset],
                                [1 + Offset, 2 + Offset, 0 + Offset],
                                [1 + Offset, 2 + Offset, 1 + Offset],
                                [1 + Offset, 2 + Offset, 2 + Offset],
                            ]

                            Annotation = [
                                MirrorMatrix @ np.array(vec) if vec[0] > 0 else np.array(vec)
                                for vec in Annotation
                            ]

                            Annotation = [RotationMatrix @ vec for vec in Annotation]

                            print(Annotation)

                            for vec in Annotation:
                                vec = np.array(vec)
                                unit_vec = vec / np.linalg.norm(vec)
                                rounded = np.round(vec).astype(int)

                                if str(Symmerty.name) in ["1", "4/mmm"]:  # Miller
                                    label = r'$[' + ''.join(
                                        f'\\overline{{{abs(i)}}}' if i < 0 else str(i)
                                        for i in rounded
                                    ) + ']$'
                                elif str(Symmerty.name) == "6/mmm":  # MillerBravis
                                    label = format_bravais_label(rounded)
                                    print(vec,label)


                                xy = ax._projection.vector2xy(Vector3d(unit_vec))
                                ax.annotate(label, xy=xy, fontsize=FontSize, color=Color,
                                            ha='center', va='center', fontweight='bold', zorder=ZOrder)
                        # </editor-fold>

                        # ax.annotate('[001]', xy=ax._projection.vector2xy(Vector3d([0+Offset,0+Offset,1+Offset])), fontsize=FontSize, color= Color, ha='center', va='center',fontweight='bold', zorder=ZOrder)
                        # ax.annotate('[100]', xy=ax._projection.vector2xy(Vector3d([1+Offset,0+Offset,0+Offset])), fontsize=FontSize, color= Color, ha='center', va='center',fontweight='bold', zorder=ZOrder)
                        # ax.annotate('[110]', xy=ax._projection.vector2xy(Vector3d([1+Offset,1+Offset,0+Offset])), fontsize=FontSize, color= Color, ha='center', va='center',fontweight='bold', zorder=ZOrder)
                        # ax.annotate('[210]', xy=ax._projection.vector2xy(Vector3d([2+Offset,1+Offset,0+Offset])), fontsize=FontSize, color= Color, ha='center', va='center',fontweight='bold', zorder=ZOrder)
                        # ax.annotate('[111]', xy=ax._projection.vector2xy(Vector3d([1+Offset,1+Offset,1+Offset])), fontsize=FontSize, color= Color, ha='center', va='center',fontweight='bold', zorder=ZOrder)
                        # ax.annotate('[112]', xy=ax._projection.vector2xy(Vector3d([1+Offset,1+Offset,2+Offset])), fontsize=FontSize, color= Color, ha='center', va='center',fontweight='bold', zorder=ZOrder)
                        # ax.annotate('[101]', xy=ax._projection.vector2xy(Vector3d([1+Offset,0+Offset,1+Offset])), fontsize=FontSize, color= Color, ha='center', va='center',fontweight='bold', zorder=ZOrder)
                        # ax.annotate('[301]', xy=ax._projection.vector2xy(Vector3d([3+Offset,0+Offset,1+Offset])), fontsize=FontSize, color= Color, ha='center', va='center',fontweight='bold', zorder=ZOrder)
                        # ax.annotate('[212]', xy=ax._projection.vector2xy(Vector3d([2+Offset,1+Offset,2+Offset])), fontsize=FontSize, color= Color, ha='center', va='center',fontweight='bold', zorder=ZOrder)
                        # ax.annotate('[211]', xy=ax._projection.vector2xy(Vector3d([2+Offset,1+Offset,1+Offset])), fontsize=FontSize, color= Color, ha='center', va='center',fontweight='bold', zorder=ZOrder)

                        Offset = 0
                        # [0,0,1] [0,0,1]
                        # [0,1,0] [1,0,0]
                        # [-1,1,0] [1,1,0]
                        # [-1,2,0] [2,1,0]
                        # [-1,1,1] [1,1,1]
                        # [-1,1,2] [1,1,2]
                        # [0,1,1] [1,0,1]
                        # [0,3,1] [3,0,1]
                        # [-1,2,2] [2,1,2]
                        # [-1,2,1] [2,1,1]

                        # ax.annotate('[001]', xy=ax._projection.vector2xy(Vector3d([0,0,1])), fontsize=FontSize, color= Color, ha='center', va='center',fontweight='bold', zorder=ZOrder)
                        # ax.annotate('[010]', xy=ax._projection.vector2xy(Vector3d([1,0,0])), fontsize=FontSize, color= Color, ha='center', va='center',fontweight='bold', zorder=ZOrder)
                        # ax.annotate('[-110]', xy=ax._projection.vector2xy(Vector3d([1,1,0])), fontsize=FontSize, color= Color, ha='center', va='center',fontweight='bold', zorder=ZOrder)
                        # ax.annotate('[-120]', xy=ax._projection.vector2xy(Vector3d([2,1,0])), fontsize=FontSize, color= Color, ha='center', va='center',fontweight='bold', zorder=ZOrder)
                        # ax.annotate('[-111]', xy=ax._projection.vector2xy(Vector3d([1,1,1])), fontsize=FontSize, color= Color, ha='center', va='center',fontweight='bold', zorder=ZOrder)
                        # ax.annotate('[-112]', xy=ax._projection.vector2xy(Vector3d([1,1,2])), fontsize=FontSize, color= Color, ha='center', va='center',fontweight='bold', zorder=ZOrder)
                        # ax.annotate('[011]', xy=ax._projection.vector2xy(Vector3d([1,0,1])), fontsize=FontSize, color= Color, ha='center', va='center',fontweight='bold', zorder=ZOrder)
                        # ax.annotate('[031]', xy=ax._projection.vector2xy(Vector3d([3,0,1])), fontsize=FontSize, color= Color, ha='center', va='center',fontweight='bold', zorder=ZOrder)
                        # ax.annotate('[-122]', xy=ax._projection.vector2xy(Vector3d([2,1,2])), fontsize=FontSize, color= Color, ha='center', va='center',fontweight='bold', zorder=ZOrder)
                        # ax.annotate('[-121]', xy=ax._projection.vector2xy(Vector3d([2,1,1])), fontsize=FontSize, color= Color, ha='center', va='center',fontweight='bold', zorder=ZOrder)

                        # ax.annotate(r'$[001]$', xy=ax._projection.vector2xy(Vector3d([0, 0, 1])), fontsize=FontSize, color= Color, ha='center', va='center',fontweight='bold', zorder=ZOrder)
                        # ax.annotate(r'$[010]$', xy=ax._projection.vector2xy(Vector3d([1, 0, 0])), fontsize=FontSize, color= Color, ha='center', va='center',fontweight='bold', zorder=ZOrder)
                        # ax.annotate(r'$[\overline{1}10]$', xy=ax._projection.vector2xy(Vector3d([1, 1, 0])), fontsize=FontSize, color= Color, ha='center', va='center',fontweight='bold', zorder=ZOrder)
                        # ax.annotate(r'$[\overline{1}20]$', xy=ax._projection.vector2xy(Vector3d([2, 1, 0])), fontsize=FontSize, color= Color, ha='center', va='center',fontweight='bold', zorder=ZOrder)
                        # ax.annotate(r'$[\overline{1}11]$', xy=ax._projection.vector2xy(Vector3d([1, 1, 1])), fontsize=FontSize, color= Color, ha='center', va='center',fontweight='bold', zorder=ZOrder)
                        # ax.annotate(r'$[\overline{1}12]$', xy=ax._projection.vector2xy(Vector3d([1, 1, 2])), fontsize=FontSize, color= Color, ha='center', va='center',fontweight='bold', zorder=ZOrder)
                        # ax.annotate(r'$[011]$', xy=ax._projection.vector2xy(Vector3d([1, 0, 1])), fontsize=FontSize, color= Color, ha='center', va='center',fontweight='bold', zorder=ZOrder)
                        # ax.annotate(r'$[031]$', xy=ax._projection.vector2xy(Vector3d([3, 0, 1])), fontsize=FontSize, color= Color, ha='center', va='center',fontweight='bold', zorder=ZOrder)
                        # ax.annotate(r'$[\overline{1}22]$', xy=ax._projection.vector2xy(Vector3d([2, 1, 2])), fontsize=FontSize, color= Color, ha='center', va='center',fontweight='bold', zorder=ZOrder)
                        # ax.annotate(r'$[\overline{1}21]$', xy=ax._projection.vector2xy(Vector3d([2, 1, 1])), fontsize=FontSize, color= Color, ha='center', va='center',fontweight='bold', zorder=ZOrder)

                        # ax.annotate(r'$[0001]$', xy=ax._projection.vector2xy(Vector3d([0, 0, 1])), fontsize=FontSize, color= Color, ha='center', va='center',fontweight='bold', zorder=ZOrder)
                        # ax.annotate(r'$[\overline{1}2\overline{1}0]$', xy=ax._projection.vector2xy(Vector3d([1, 0, 0])), fontsize=FontSize, color= Color, ha='center', va='center',fontweight='bold', zorder=ZOrder)
                        # ax.annotate(r'$[\overline{1}10 0]$', xy=ax._projection.vector2xy(Vector3d([1, 1, 0])), fontsize=FontSize, color= Color, ha='center', va='center',fontweight='bold', zorder=ZOrder)
                        # ax.annotate(r'$[01\overline{1}0]$', xy=ax._projection.vector2xy(Vector3d([2, 1, 0])), fontsize=FontSize, color= Color, ha='center', va='center',fontweight='bold', zorder=ZOrder)
                        # ax.annotate(r'$[\overline{1}10 1]$', xy=ax._projection.vector2xy(Vector3d([1, 1, 1])), fontsize=FontSize, color= Color, ha='center', va='center',fontweight='bold', zorder=ZOrder)
                        # ax.annotate(r'$[\overline{1}10 2]$', xy=ax._projection.vector2xy(Vector3d([1, 1, 2])), fontsize=FontSize, color= Color, ha='center', va='center',fontweight='bold', zorder=ZOrder)
                        # ax.annotate(r'$[\overline{1}2\overline{1}3]$', xy=ax._projection.vector2xy(Vector3d([1, 0, 1])), fontsize=FontSize, color= Color, ha='center', va='center',fontweight='bold', zorder=ZOrder)
                        # ax.annotate(r'$[\overline{1}2\overline{1}1]$', xy=ax._projection.vector2xy(Vector3d([3, 0, 1])), fontsize=FontSize, color= Color, ha='center', va='center',fontweight='bold', zorder=ZOrder)
                        # ax.annotate(r'$[01\overline{1}2]$', xy=ax._projection.vector2xy(Vector3d([2, 1, 2])), fontsize=FontSize, color= Color, ha='center', va='center',fontweight='bold', zorder=ZOrder)
                        # ax.annotate(r'$[01\overline{1}1]$', xy=ax._projection.vector2xy(Vector3d([2, 1, 1])), fontsize=FontSize, color= Color, ha='center', va='center',fontweight='bold', zorder=ZOrder)

                        fig.suptitle("(" + str(LetterLabels[LabelCounter]) + "): " + str(Potential), fontsize=24, y=1.00,
                                     color=Color)
                        # ax.set_title(Potential)
                        plt.tight_layout()
                        plt.show()
                        plt.savefig(PythonName + "-" + Title, dpi=600, bbox_inches='tight')
                        plt.close()
                        LabelCounter += 1
            # </editor-fold>

            # <editor-fold desc="Plot-inverse pole figures-contour and points-TdeGroupedPotentialDirectionDescribe-Fitting">
            Active = True
            if Active:
                # print("Plot-inverse pole figures-contour and points-TdeGroupedPotentialDirectionDescribe-Fitting")
                LetterLabels = list("abcdefghijklmnopqrstuvwxyz")
                LabelCounter = 0

                PotentialList = ["M2R", "M3R", "BMD192R", "TabGap1", "TabGap2"]

                for Potential in PotentialList:
                    if Potential in Content['Potential'].values:
                        Group = Content[Content['Potential'] == Potential]
                        print("-" + Potential)

                        # print(Group.columns)
                        # print(Group[["Direction","Repeating"]])

                        Group["TransformedX"] = Group['UnitX']
                        Group["TransformedY"] = Group['UnitY']
                        Group["TransformedZ"] = Group['UnitZ']

                        # <editor-fold desc="Reflection">
                        Reflection = False
                        if Reflection:
                            MirrorMatrix = np.array([
                                [-1, 0, 0],
                                [0, 1, 0],
                                [0, 0, 1]
                            ])

                            OriginalStack = np.column_stack([
                                Group['UnitX'].values,
                                Group['UnitY'].values,
                                Group['UnitZ'].values
                            ])

                            MirroredStack = np.array([
                                MirrorMatrix @ vec if vec[0] > 0 else vec
                                for vec in OriginalStack
                            ])

                            Group["TransformedX"] = MirroredStack[:, 0]
                            Group["TransformedY"] = MirroredStack[:, 1]
                            Group["TransformedZ"] = MirroredStack[:, 2]

                        # print(Group.columns)
                        # print(Group[["Direction", "Repeating", 'TransformedX', 'TransformedY', 'TransformedZ']])

                        # </editor-fold>

                        # <editor-fold desc="Rotation">
                        Rotating = True
                        if Rotating:
                            RotationDeg = -90
                            RotationRad = np.deg2rad(RotationDeg)
                            # RotationRad = np.deg2rad(-1.1)
                            # print(RotationRad)

                            RotationMatrix = np.array([
                                [np.cos(RotationRad), -np.sin(RotationRad), 0],
                                [np.sin(RotationRad), np.cos(RotationRad), 0],
                                [0, 0, 1]
                            ])

                            OriginalStack = np.column_stack([
                                Group['X'].values,
                                Group['Y'].values,
                                Group['Z'].values
                            ])
                            RotatedStack = OriginalStack @ RotationMatrix.T # this is just to save the rotation. not used anywhere
                            Group["RotateX"] = RotatedStack[:, 0].round().astype(int)
                            Group["RotateY"] = RotatedStack[:, 1].round().astype(int)
                            Group["RotateZ"] = RotatedStack[:, 2].round().astype(int)

                            OriginalStack = np.column_stack([
                                Group['TransformedX'].values,
                                Group['TransformedY'].values,
                                Group['TransformedZ'].values
                            ])
                            RotatedStack = OriginalStack @ RotationMatrix.T
                            Group["TransformedX"] = RotatedStack[:, 0]
                            Group["TransformedY"] = RotatedStack[:, 1]
                            Group["TransformedZ"] = RotatedStack[:, 2]

                            # print(Group[["Direction","RotateX","RotateY","RotateZ"]])
                            # os.system("pause")
                        else:
                            RotationMatrix = np.array([
                                [1, 0, 0],
                                [0, 1, 0],
                                [0, 0, 1]
                            ])

                        UnitStack = np.column_stack([
                            Group['TransformedX'].values,
                            Group['TransformedY'].values,
                            Group['TransformedZ'].values
                        ])
                        Vectors = Vector3d(UnitStack)
                        # print(Vectors)
                        # </editor-fold>

                        # <editor-fold desc="Filter">
                        Filtering = True
                        if Filtering:
                            # Filter out directions with positive x-component
                            OriginalStack = np.column_stack([
                                Group['TransformedX'].values,
                                Group['TransformedY'].values,
                                Group['TransformedZ'].values
                            ])
                            # print(OriginalStack)

                            Mask = OriginalStack[:, 1] >= 0

                            # Apply mask to keep only valid vectors
                            FilteredStack = OriginalStack[Mask]

                            # Also filter the corresponding rows in Group
                            Group = Group[Mask].copy()

                            # Assign transformed values (same as original in this case)
                            Group["TransformedX"] = FilteredStack[:, 0]
                            Group["TransformedY"] = FilteredStack[:, 1]
                            Group["TransformedZ"] = FilteredStack[:, 2]
                        # </editor-fold>

                        # <editor-fold desc="Gridding">
                        n_grid = 500#5000
                        SectionMinDeg = 0
                        SectionMaxDeg = 45
                        SectionMinRad = np.deg2rad(SectionMinDeg)
                        SectionMaxRad = np.deg2rad(SectionMaxDeg)
                        GridSpacePhi = np.linspace(SectionMinRad, SectionMaxRad, n_grid)
                        GridSpaceTheta = np.linspace(0, np.pi / 2, n_grid)
                        MeshPhi, MeshTheta = np.meshgrid(GridSpacePhi, GridSpaceTheta)
                        GridX = np.sin(MeshTheta) * np.cos(MeshPhi)
                        GridY = np.sin(MeshTheta) * np.sin(MeshPhi)
                        GridZ = np.cos(MeshTheta)
                        X_proj = GridX / (1 + GridZ)
                        Y_proj = GridY / (1 + GridZ)
                        GridPoints = np.column_stack([GridX.ravel(), GridY.ravel(), GridZ.ravel()])
                        # print(GridPoints)
                        Weights = Group[Based].values
                        # </editor-fold>

                        # <editor-fold desc="Fitting">
                        # print(list(Group))
                        # print(Group[['Direction',"RotateX","RotateY","RotateZ", Based]])
                        FittingList = ["nearest","rbf"]
                        for Fitting in FittingList:
                            print("Fitting")
                            if Fitting == "Kde":
                                KdeModel = KernelDensity(kernel='gaussian', bandwidth=0.1)
                                KdeModel.fit(Vectors.data, sample_weight=Weights)
                                log_density = KdeModel.score_samples(GridPoints)
                                # print(log_density)
                                Density = np.exp(log_density).reshape(GridX.shape)
                                # print("Min Density:", np.min(Density))
                                # print("Max Density:", np.max(Density))
                                # print(Density)
                            elif Fitting == "nearest":
                                from scipy.interpolate import griddata
                                GridPoints = np.column_stack([GridX.ravel(), GridY.ravel(), GridZ.ravel()])
                                Points = np.column_stack([Group['TransformedX'], Group['TransformedY'], Group['TransformedZ']])
                                Density = griddata(Points, Weights, GridPoints, method='nearest')
                                Density = Density.reshape(GridX.shape)
                            elif Fitting == "rbf":
                                from scipy.interpolate import RBFInterpolator
                                Points = np.column_stack([Group['TransformedX'], Group['TransformedY'], Group['TransformedZ']])
                                rbf = RBFInterpolator(Points, Weights, kernel='multiquadric', epsilon=10000)
                                Density = rbf(GridPoints)
                                Density = Density.reshape(GridX.shape)
                            elif Fitting == "griddata":
                                from scipy.interpolate import griddata
                                Points2D = np.column_stack([Group['TransformedX'], Group['TransformedY']])
                                GridPoints2D = np.column_stack([GridX.ravel(), GridY.ravel()])

                                Density = griddata(Points2D, Weights, GridPoints2D, method='cubic')
                                Density = Density.reshape(GridX.shape)
                            # </editor-fold>

                            # print(X_proj, Y_proj, Density)

                            fig, ax = plt.subplots(figsize=(8, 6))
                            Title = "PoleFigureFitting" + "-" + Key + "-" + Potential + "-" + str(Based) + "-" + Fitting
                            sc = ax.contourf(X_proj, Y_proj, Density, levels=100, cmap="jet")
                            ax.scatter(Group["TransformedX"] / (1 + Group["TransformedZ"]),
                                       Group["TransformedY"] / (1 + Group["TransformedZ"]),
                                       color='k', s=20, alpha=1, zorder=10)

                            divider = make_axes_locatable(ax)
                            cax = divider.append_axes("right", size="5%", pad=0.5)
                            cbar = fig.colorbar(sc, cax=cax)
                            cbar.set_label("TDE (eV)")

                            cbar.locator = mticker.MaxNLocator(nbins=5, integer=True)
                            cbar.update_ticks()
                            cbar.ax.yaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f'))

                            FontSize = 10
                            Color = "black"
                            FontWeight = "normal"
                            ZOrder = 10
                            Offset = 0.0

                            # <editor-fold desc="Annotation">
                            Annotation = True
                            if Annotation:
                                AnnotationType = "Auto" # Manual
                                if AnnotationType == "Manual":
                                    Annotation = [
                                        [-4 + Offset, 4 + Offset, 1 + Offset],  # [4, 4, 1]
                                        [-3 + Offset, 4 + Offset, 2 + Offset],  # [4, 4, 1]
                                        [-2 + Offset, 2 + Offset, 1 + Offset],  # [2, 2, 1]
                                        [-2 + Offset, 3 + Offset, 0 + Offset],  # [3, 2, 0]
                                        [-2 + Offset, 3 + Offset, 1 + Offset],  # [3, 2, 1]
                                        [-2 + Offset, 4 + Offset, 1 + Offset],  # [4, 2, 1]
                                        [-1 + Offset, 1 + Offset, 0 + Offset],  # [1, 1, 0]
                                        [-1 + Offset, 1 + Offset, 1 + Offset],  # [1, 1, 0]
                                        [-1 + Offset, 1 + Offset, 2 + Offset],  # [1, 1, 2]
                                        [-1 + Offset, 2 + Offset, 0 + Offset],  # [2, 1, 0]
                                        [-1 + Offset, 2 + Offset, 3 + Offset],  # [2, 1, 3]
                                        [-1 + Offset, 3 + Offset, 1 + Offset],  # [3, 1, 1]
                                        [-1 + Offset, 3 + Offset, 2 + Offset],  # [3, 1, 2]
                                        # # [-1 + Offset, 3 + Offset, 5 + Offset],  # [3, 1, 5]
                                        [-1 + Offset, 2 + Offset, 1 + Offset],  # [2, 1, 1]
                                        [-1 + Offset, 2 + Offset, 2 + Offset],  # [2, 1, 2]
                                        [-1 + Offset, 3 + Offset, 0 + Offset],  # [3, 1, 0]
                                        [-1 + Offset, 3 + Offset, 3 + Offset],  # [-1, 2, 2]
                                        [-1 + Offset, 4 + Offset, 0 + Offset],  # [4, 1, 0]
                                        [0 + Offset, 0 + Offset, 1 + Offset],  # [0, 0, 1]
                                        [0 + Offset, 1 + Offset, 0 + Offset],  # [0, 1, 0]
                                        [0 + Offset, 1 + Offset, 1 + Offset],  # [1, 0, 1]
                                        [0 + Offset, 1 + Offset, 2 + Offset],  # [1, 0, 2]
                                        # # # [0 + Offset, 1 + Offset, 3 + Offset],  # [1, 0, 3]
                                        [0 + Offset, 2 + Offset, 1 + Offset],  # [2, 0, 1]
                                        [0 + Offset, 2 + Offset, 3 + Offset],  # [2, 0, 3]
                                        [0 + Offset, 3 + Offset, 1 + Offset],  # [3, 0, 1]
                                        [0 + Offset, 4 + Offset, 1 + Offset],  # [4, 0, 1]
                                        [1 + Offset, 2 + Offset, 0 + Offset],  # [1, 2, 0]
                                        [1 + Offset, 2 + Offset, 1 + Offset],  # [1, 2, 1]
                                        # [1 + Offset, 3 + Offset, 0 + Offset],  # [3, 1, 0]
                                        # [1 + Offset, 4 + Offset, 0 + Offset],  # [4, 1, 0]
                                        # [1 + Offset, 4 + Offset, 1 + Offset],  # [1, 4, 1]
                                        # [1 + Offset, 4 + Offset, 2 + Offset],  # [1, 4, 2]
                                        # [1 + Offset, 4 + Offset, 3 + Offset],  # [1, 4, 3]
                                        # [2 + Offset, 2 + Offset, 1 + Offset],  # [2, 2, 1]
                                        # [2 + Offset, 2 + Offset, 3 + Offset],  # [2, 2, 3]
                                        # [2 + Offset, 3 + Offset, 0 + Offset],  # [2, 3, 0]
                                        # [2 + Offset, 3 + Offset, 1 + Offset],  # [2, 3, 1]
                                        # [2 + Offset, 3 + Offset, 2 + Offset],  # [2, 3, 2]
                                        [2 + Offset, 3 + Offset, 3 + Offset],  # [2, 3, 3]
                                        # # # [3 + Offset, 0 + Offset, 1 + Offset],  # [3, 1, 0]
                                        # # # [3 + Offset, 1 + Offset, 0 + Offset],  # [3, 1, 0]
                                        # # [3 + Offset, 1 + Offset, 3 + Offset],  # [3, 4, 0]
                                        [-3 + Offset, 4 + Offset, 0 + Offset],  # [4, 3, 0]
                                        [-3 + Offset, 4 + Offset, 1 + Offset],  # [4, 3, 1]
                                        [-3 + Offset, 4 + Offset, 2 + Offset],  # [4, 3, 2]
                                        [-3 + Offset, 4 + Offset, 3 + Offset],  # [4, 3, 3]
                                        ]
                                elif AnnotationType =="Auto":
                                    def AnnotationFunc(group_df, offset=0):
                                        annotation = []

                                        for _, row in group_df.iterrows():
                                            x, y, z = int(row['X']), int(row['Y']), int(row['Z'])
                                            annotation.append([x + offset, y + offset, z + offset])

                                        return annotation


                                    Annotation = AnnotationFunc(Group, Offset)

                                if Reflection:
                                    Annotation = [
                                        MirrorMatrix @ np.array(vec) if vec[0] > 0 else np.array(vec)
                                        for vec in Annotation
                                    ]

                                if Rotating:
                                    Annotation = [RotationMatrix @ vec for vec in Annotation]
                                # print(Annotation)
                                for vec in Annotation:
                                    vec = np.array(vec)
                                    unit_vec = vec / np.linalg.norm(vec)
                                    rounded = np.round(vec).astype(int)

                                    Symmerty = "4/mmm"
                                    if str(Symmerty) in ["1", "4/mmm"]:  # Miller
                                        label = r'$[' + ''.join(
                                            f'\\overline{{{abs(i)}}}' if i < 0 else str(i)
                                            for i in rounded
                                        ) + ']$'
                                    elif str(Symmerty) == "6/mmm":  # MillerBravis
                                        label = format_bravais_label(rounded)
                                        # print(vec, label)

                                    # xy = ax._projection.vector2xy(Vector3d(unit_vec))
                                    # x, y = unit_vec[0], unit_vec[1]  # project to XY plane
                                    x = unit_vec[0] / (1 + unit_vec[2])
                                    y = unit_vec[1] / (1 + unit_vec[2])
                                    import matplotlib.patheffects as pe

                                    ax.annotate(
                                        label,
                                        xy=(x, y),
                                        fontsize=FontSize,
                                        color=Color,
                                        ha='center',
                                        va='center',
                                        fontweight='bold',
                                        zorder=ZOrder,
                                        path_effects=[
                                            pe.withStroke(linewidth=3, foreground='white')
                                        ]
                                    )
                            # </editor-fold>

                            ax.set_aspect("equal")
                            ax.axis("off")

                            if Potential == "M2R" and  Fitting == "nearest":
                                fig.suptitle("Nearest", fontsize=18, y=0.9)
                            elif Potential == "M2R" and  Fitting == "rbf":
                                fig.suptitle("RBF", fontsize=18, y=0.9)

                            fig.suptitle(f"({LetterLabels[LabelCounter]})", fontsize=18, x=0.05, y=0.775)

                            # fig.suptitle(f"({LetterLabels[LabelCounter]}): {Potential}", fontsize=18, y=0.9)
                            if Fitting == "nearest":
                                fig.text(0.02, 0.5, Potential, va='center', ha='center', rotation=90)

                            plt.tight_layout()
                            # plt.show()
                            plt.savefig(PythonName + "-" + Title, dpi=600, bbox_inches='tight', transparent=True)
                            plt.close()
                            LabelCounter += 1
            # </editor-fold>

# <editor-fold desc="***** Plot-inverse pole figures-Yu">
print("***** Plot-inverse pole figures-Yu")
Active = False
if Active:
    Z1 = [[1,1,1,1,1,1,1],[2,2,2,2,2,2,2],[3,3,3,3,3,3,3],[4,4,4,4,4,4,4],[5,5,5,5,5,5,5],[6,6,6,6,6,6,6],[7,7,7,7,7,7,7]]

    from mpl_toolkits.mplot3d import Axes3D
    from scipy.interpolate import griddata

    # plt.style.use('seaborn-white')
    angle = 30
    angle_rad = np.deg2rad(angle)

    points = 7
    r = np.linspace(0.001, 1, 7)
    theta = np.linspace(-angle_rad, angle_rad, points)
    r, theta = np.meshgrid(r, theta)
    X = r * np.cos(theta)
    Y = r * np.sin(theta)


    r = np.linspace(0.001, 1, 800)
    theta = np.linspace(-angle_rad, angle_rad, 800)
    r, theta = np.meshgrid(r, theta)
    X1 = r * np.cos(theta)
    Y1 = r * np.sin(theta)

    Z = griddata((X.flatten(), Y.flatten()), np.array(Z1).flatten(), (X1, Y1), method='cubic')

    plt.figure(figsize=(8, 6))

    # Plot the sector using a colormap
    plt.pcolormesh(X1, Y1, Z, cmap='coolwarm', vmin=0, vmax=35)

    # Add a color bar
    cbar=plt.colorbar(pad=0.08)
    plt.axis('off')
    cbar.set_label('Threshold displacement energy(eV)', fontsize=12)
    plt.tight_layout()
    plt.legend(fontsize=15)
    plt.savefig("./figure_tde_GAP_surface.jpg",dpi=600, bbox_inches='tight')
    plt.show()
# </editor-fold>
# </editor-fold>

# </editor-fold>

# <editor-fold desc="######################################## Report">
print("######################################## Report")  # %%
Active = False
if Active:
    # <editor-fold desc="***** Single">
    print("***** Single")

    # <editor-fold desc="Excess">
    print("Excess")
    ReportSingleDf["EnergyPotentialExcess"] = (ReportSingleDf["EnergyPotential"] - ReportSingleDf.groupby(["Potential", "Try", "Direction", "Energy"])["EnergyPotential"].transform("min"))
    ReportSingleDf["EnergyKineticExcess"] = (ReportSingleDf["EnergyKinetic"] - ReportSingleDf.groupby(["Potential", "Try", "Direction", "Energy"])["EnergyKinetic"].transform("min"))
    ReportSingleDf["EnergyTotalExcess"] = (ReportSingleDf["EnergyTotal"] - ReportSingleDf.groupby(["Potential", "Try", "Direction", "Energy"])["EnergyTotal"].transform("min"))

    ReportSingleDf["NoOfNeighborsAverage"] = ReportSingleDf["NoOfNeighbors"] / 10752
    # </editor-fold>

    # <editor-fold desc="Group">
    print("Group")
    # print(ReportSingleDf)
    ReportSingleDfGroupedPotentialTryDirectionEnergy = ReportSingleDf.groupby(["Potential", "Try", "Direction", "Energy"])
    # </editor-fold>

    EnergyList = [28]
    EnergyList = [26,28,30]

    # <editor-fold desc="Plotting-NoOfNeighbors-Single">
    print("Plotting-NoOfNeighbors-Single")
    Plotting = False
    if Plotting:
        Title = "ReportSingleDfNoOfNeighbors"

        # PotentialList = ["M3R"]
        # TryList = [1]
        # DirectionList = ["dir_-1_2_0"]
        # ReportDfFiltered = ReportDf[ReportDf['Potential'].isin(PotentialList)].reset_index(drop=True)
        # ReportDfFiltered = ReportDfFiltered[ReportDfFiltered['Try'].isin(TryList)]
        # ReportDfFiltered = ReportDfFiltered[ReportDfFiltered['Direction'].isin(DirectionList)].reset_index(drop=True)

        # TimeMatch = ReportDfMinDist[
        #     (ReportDfMinDist['Potential'] == "M3R") &
        #     (ReportDfMinDist['Try'] == 1) &
        #     (ReportDfMinDist['Direction'] == "dir_-1_2_0")
        # ]
        #
        # if not TimeMatch.empty:
        #     target_time = TimeMatch.iloc[0]['Time']
        # else:
        #     target_time = None
        #
        plt.figure(figsize=(12, 6))

        ReportSingleDfFiltered = ReportSingleDf[ReportSingleDf['Energy'].isin(EnergyList)]
        # ReportSingleDfFiltered = (
        #     ReportSingleDfFiltered
        #     .sort_values("Time")  # make sure ordering is stable
        #     .groupby("Energy", group_keys=False)
        #     .nth[::10]  # pick every 10th row
        #     .reset_index()
        # )
        ReportSingleDfFiltered["Energy"] = ReportSingleDf["Energy"].astype("str")

        colors = ["blue", "black", "red"]
        black_seismic = mcolors.LinearSegmentedColormap.from_list("black_seismic", colors, N=256)
        n_colors = len(ReportSingleDfFiltered["Energy"].unique())
        BlueRedPallet = [black_seismic(i) for i in np.linspace(0, 1, n_colors)]

        sns.lineplot(
            data=ReportSingleDfFiltered,
            x="Time",
            y="NoOfNeighbors",
            hue="Energy",
            alpha=0.75,
            linewidth=4,
            palette=BlueRedPallet,
            # estimator=None,
            # sort=False,
        )


        # plt.xticks(np.arange(0, 1.1, 0.5),fontsize=40)
        # plt.yticks(np.arange(1.5, 3.1, 0.5),fontsize=40)
        plt.tick_params(axis='y', labelsize=30)  # set ytick label size (optional)
        plt.tick_params(labelbottom=False)
        plt.xscale("log")
        plt.xlabel("Time (ps)", fontsize=40)
        plt.ylabel("No. of Neighbors", fontsize=40, labelpad=17)
        plt.legend(title="Energy", title_fontsize=30, fontsize=25)
        # plt.legend([], [], frameon=False)  # removes the legend completely
        plt.tight_layout()
        plt.savefig(PythonName + "-" + Title, dpi=600, bbox_inches='tight', transparent=True)
        plt.show()
        plt.close()
    # </editor-fold>

    # <editor-fold desc="Plotting-MinDist-Single">
    print("Plotting-MinDist-Single")
    Plotting = False
    if Plotting:
        import math

        Title = "ReportSingleDfMinDist"

        # PotentialList = ["M3R"]
        # TryList = [1]
        # DirectionList = ["dir_-1_2_0"]
        # ReportDfFiltered = ReportDf[ReportDf['Potential'].isin(PotentialList)].reset_index(drop=True)
        # ReportDfFiltered = ReportDfFiltered[ReportDfFiltered['Try'].isin(TryList)]
        # ReportDfFiltered = ReportDfFiltered[ReportDfFiltered['Direction'].isin(DirectionList)].reset_index(drop=True)

        # TimeMatch = ReportDfMinDist[
        #     (ReportDfMinDist['Potential'] == "M3R") &
        #     (ReportDfMinDist['Try'] == 1) &
        #     (ReportDfMinDist['Direction'] == "dir_-1_2_0")
        # ]
        #
        # if not TimeMatch.empty:
        #     target_time = TimeMatch.iloc[0]['Time']
        # else:
        #     target_time = None
        #
        plt.figure(figsize=(12, 6))
        ReportSingleDfFiltered = ReportSingleDf[ReportSingleDf['Energy'].isin(EnergyList)]
        # ReportSingleDfFiltered = (
        #     ReportSingleDfFiltered
        #     .sort_values("Time")  # make sure ordering is stable
        #     .groupby("Energy", group_keys=False)
        #     .nth[::10]  # pick every 10th row
        #     .reset_index()
        # )
        ReportSingleDfFiltered["Energy"] = ReportSingleDf["Energy"].astype("str")

        colors = ["blue", "black", "red"]
        black_seismic = mcolors.LinearSegmentedColormap.from_list("black_seismic", colors, N=256)
        n_colors = len(ReportSingleDfFiltered["Energy"].unique())
        BlueRedPallet = [black_seismic(i) for i in np.linspace(0, 1, n_colors)]

        sns.lineplot(
            data=ReportSingleDfFiltered,
            x="Time",
            y="MinDist",
            hue="Energy",
            alpha=0.75,
            linewidth=4,
            palette=BlueRedPallet,
            # estimator=None,
            # sort=False,
        )


        # plt.xticks(np.arange(0, 1.1, 0.5),fontsize=40)
        # plt.yticks(np.arange(1.5, 3.1, 0.5),fontsize=40)
        # ax.tick_params(labelbottom=False)
        plt.tick_params(axis='x', labelsize=30)  # set xtick label size
        plt.tick_params(axis='y', labelsize=30)  # set ytick label size (optional)
        plt.legend([], [], frameon=False)  # removes the legend completely
        plt.xscale("log")
        plt.xlabel("Time (ps)", fontsize=40)
        plt.ylabel("Min Dist (Å)", fontsize=40)
        plt.tight_layout()
        plt.savefig(PythonName + "-" + Title, dpi=600, bbox_inches='tight', transparent=False)
        plt.show()
        plt.close()
    # </editor-fold>

    # <editor-fold desc="Plotting-EnergyTotalExcess-Single">
    print("Plotting-EnergyTotalExcess-Single")
    Plotting = False
    if Plotting:
        Title = "ReportSingleDfEnergyTotalExcess"

        # PotentialList = ["M3R"]
        # TryList = [1]
        # DirectionList = ["dir_-1_2_0"]
        # ReportDfFiltered = ReportDf[ReportDf['Potential'].isin(PotentialList)].reset_index(drop=True)
        # ReportDfFiltered = ReportDfFiltered[ReportDfFiltered['Try'].isin(TryList)]
        # ReportDfFiltered = ReportDfFiltered[ReportDfFiltered['Direction'].isin(DirectionList)].reset_index(drop=True)

        # TimeMatch = ReportDfMinDist[
        #     (ReportDfMinDist['Potential'] == "M3R") &
        #     (ReportDfMinDist['Try'] == 1) &
        #     (ReportDfMinDist['Direction'] == "dir_-1_2_0")
        # ]
        #
        # if not TimeMatch.empty:
        #     target_time = TimeMatch.iloc[0]['Time']
        # else:
        #     target_time = None
        #
        plt.figure(figsize=(12, 6))
        ReportSingleDfFiltered = ReportSingleDf[ReportSingleDf['Energy'].isin(EnergyList)]
        # ReportSingleDfFiltered = (
        #     ReportSingleDfFiltered
        #     .sort_values("Time")  # make sure ordering is stable
        #     .groupby("Energy", group_keys=False)
        #     .nth[::10]  # pick every 10th row
        #     .reset_index()
        # )
        ReportSingleDfFiltered["Energy"] = ReportSingleDf["Energy"].astype("str")

        colors = ["blue", "black", "red"]
        black_seismic = mcolors.LinearSegmentedColormap.from_list("black_seismic", colors, N=256)
        n_colors = len(ReportSingleDfFiltered["Energy"].unique())
        BlueRedPallet = [black_seismic(i) for i in np.linspace(0, 1, n_colors)]

        sns.lineplot(
            data=ReportSingleDfFiltered,
            x="Time",
            y="EnergyTotalExcess",
            hue="Energy",
            alpha=0.75,
            linewidth=4,
            palette=BlueRedPallet,
            # estimator=None,
            # sort=False,
        )


        # plt.xticks(np.arange(0, 1.1, 0.5),fontsize=40)
        # plt.yticks(np.arange(1.5, 3.1, 0.5),fontsize=40)
        plt.tick_params(axis='y', labelsize=30)  # set ytick label size (optional)
        plt.tick_params(labelbottom=False)
        plt.xscale("log")
        plt.xlabel("Time (ps)", fontsize=40)
        plt.ylabel("Exc. Pot. E. (eV)", fontsize=40, labelpad=17)
        plt.legend([], [], frameon=False)  # removes the legend completely
        plt.tight_layout()
        plt.savefig(PythonName + "-" + Title, dpi=600, bbox_inches='tight', transparent=False)
        plt.show()
        plt.close()
    # </editor-fold>

    # <editor-fold desc="Plotting-EnergyPotentialExcess-Single">
    print("Plotting-EnergyPotentialExcess-Single")
    Plotting = False
    if Plotting:
        Title = "ReportSingleDfEnergyPotentialExcess"

        # PotentialList = ["M3R"]
        # TryList = [1]
        # DirectionList = ["dir_-1_2_0"]
        # ReportDfFiltered = ReportDf[ReportDf['Potential'].isin(PotentialList)].reset_index(drop=True)
        # ReportDfFiltered = ReportDfFiltered[ReportDfFiltered['Try'].isin(TryList)]
        # ReportDfFiltered = ReportDfFiltered[ReportDfFiltered['Direction'].isin(DirectionList)].reset_index(drop=True)

        # TimeMatch = ReportDfMinDist[
        #     (ReportDfMinDist['Potential'] == "M3R") &
        #     (ReportDfMinDist['Try'] == 1) &
        #     (ReportDfMinDist['Direction'] == "dir_-1_2_0")
        # ]
        #
        # if not TimeMatch.empty:
        #     target_time = TimeMatch.iloc[0]['Time']
        # else:
        #     target_time = None
        #
        plt.figure(figsize=(12, 6))
        ReportSingleDfFiltered = ReportSingleDf[ReportSingleDf['Energy'].isin(EnergyList)]
        # ReportSingleDfFiltered = (
        #     ReportSingleDfFiltered
        #     .sort_values("Time")  # make sure ordering is stable
        #     .groupby("Energy", group_keys=False)
        #     .nth[::10]  # pick every 10th row
        #     .reset_index()
        # )
        ReportSingleDfFiltered["Energy"] = ReportSingleDf["Energy"].astype("str")

        colors = ["blue", "black", "red"]
        black_seismic = mcolors.LinearSegmentedColormap.from_list("black_seismic", colors, N=256)
        n_colors = len(ReportSingleDfFiltered["Energy"].unique())
        BlueRedPallet = [black_seismic(i) for i in np.linspace(0, 1, n_colors)]

        sns.lineplot(
            data=ReportSingleDfFiltered,
            x="Time",
            y="EnergyPotentialExcess",
            hue="Energy",
            alpha=0.75,
            linewidth=4,
            palette=BlueRedPallet,
            # estimator=None,
            # sort=False,
        )


        # plt.xticks(np.arange(0, 1.1, 0.5),fontsize=40)
        # plt.yticks(np.arange(1.5, 3.1, 0.5),fontsize=40)
        plt.tick_params(axis='y', labelsize=30)  # set ytick label size (optional)
        plt.tick_params(labelbottom=False)
        plt.xscale("log")
        plt.xlabel("Time (ps)", fontsize=40)
        plt.ylabel("Exc. Pot. E. (eV)", fontsize=40, labelpad=17)
        plt.legend([], [], frameon=False)  # removes the legend completely
        plt.tight_layout()
        plt.savefig(PythonName + "-" + Title, dpi=600, bbox_inches='tight', transparent=False)
        plt.show()
        plt.close()
    # </editor-fold>

    # <editor-fold desc="Plotting-EnergyKineticExcess-Single">
    print("Plotting-EnergyKineticExcess-Single")
    Plotting = False
    if Plotting:
        Title = "ReportSingleDfEnergyKineticExcess"

        # PotentialList = ["M3R"]
        # TryList = [1]
        # DirectionList = ["dir_-1_2_0"]
        # ReportDfFiltered = ReportDf[ReportDf['Potential'].isin(PotentialList)].reset_index(drop=True)
        # ReportDfFiltered = ReportDfFiltered[ReportDfFiltered['Try'].isin(TryList)]
        # ReportDfFiltered = ReportDfFiltered[ReportDfFiltered['Direction'].isin(DirectionList)].reset_index(drop=True)

        # TimeMatch = ReportDfMinDist[
        #     (ReportDfMinDist['Potential'] == "M3R") &
        #     (ReportDfMinDist['Try'] == 1) &
        #     (ReportDfMinDist['Direction'] == "dir_-1_2_0")
        # ]
        #
        # if not TimeMatch.empty:
        #     target_time = TimeMatch.iloc[0]['Time']
        # else:
        #     target_time = None
        #
        plt.figure(figsize=(12, 6))
        ReportSingleDfFiltered = ReportSingleDf[ReportSingleDf['Energy'].isin(EnergyList)]
        # ReportSingleDfFiltered = (
        #     ReportSingleDfFiltered
        #     .sort_values("Time")  # make sure ordering is stable
        #     .groupby("Energy", group_keys=False)
        #     .nth[::10]  # pick every 10th row
        #     .reset_index()
        # )
        ReportSingleDfFiltered["Energy"] = ReportSingleDf["Energy"].astype("str")

        colors = ["blue", "black", "red"]
        black_seismic = mcolors.LinearSegmentedColormap.from_list("black_seismic", colors, N=256)
        n_colors = len(ReportSingleDfFiltered["Energy"].unique())
        BlueRedPallet = [black_seismic(i) for i in np.linspace(0, 1, n_colors)]

        sns.lineplot(
            data=ReportSingleDfFiltered,
            x="Time",
            y="EnergyKineticExcess",
            hue="Energy",
            alpha=0.75,
            linewidth=4,
            palette=BlueRedPallet,
            # estimator=None,
            # sort=False,
        )


        # plt.xticks(np.arange(0, 1.1, 0.5),fontsize=40)
        # plt.yticks(np.arange(1.5, 3.1, 0.5),fontsize=40)
        plt.tick_params(axis='y', labelsize=30)  # set ytick label size (optional)
        plt.tick_params(labelbottom=False)
        plt.xscale("log")
        plt.xlabel("Time (ps)", fontsize=40)
        plt.ylabel("Exc. Kin. E. (eV)", fontsize=40, labelpad=17)
        plt.legend([], [], frameon=False)  # removes the legend completely
        plt.tight_layout()
        plt.savefig(PythonName + "-" + Title, dpi=600, bbox_inches='tight', transparent=False)
        plt.show()
        plt.close()
    # </editor-fold>

    # </editor-fold>

    # <editor-fold desc="***** MinDist">
    print("***** MinDist")

    # <editor-fold desc="Min of MinDist">
    print("Min of MinDist")
    Active = True
    if Active:
        # print(ReportDf['Potential'].unique())
        ReportDfMinDistIdxmin = ReportDf.groupby(['Potential', 'Try', 'Direction', 'Energy'])['MinDist'].idxmin()
        ReportDfMinDist = ReportDf.loc[ReportDfMinDistIdxmin].reset_index(drop=True)
        ReportDfMinDist.to_csv(PythonName + "-ReportDfMinDist" + ".csv", index=False)
    else:
        ReportDfMinDist = pd.read_csv(PythonName + "-ReportDfMinDist" + ".csv")
        # print(ReportDfMinDist['Potential'].unique())
    # </editor-fold>

    # <editor-fold desc="Stat">
    print("Stat")
    from scipy import stats

    ReportDfMinDistStat = (
        ReportDfMinDist
        .groupby(['Potential', 'Direction', 'DirectionMb'], as_index=False)
        .agg(
            MinDist_mean=('MinDist', 'mean'),
            MinDist_median=('MinDist', 'median'),
            MinDist_mode=('MinDist', lambda x: stats.mode(x, keepdims=True)[0][0]),
            MinDist_var=('MinDist', 'var'),
            Energy_mean=('Energy', 'mean'),
            Energy_median=('Energy', 'median'),
            Energy_mode=('Energy', lambda x: stats.mode(x, keepdims=True)[0][0]),
            Energy_var=('Energy', 'var')
        )
    )

    ReportDfMinDistStat.to_csv(PythonName + "-ReportDfMinDistStat" + ".csv", index=False)
    # </editor-fold>

    # <editor-fold desc="Group">
    print("Group")
    # print(Tde)
    ReportDfMinDistGroupedPotential = ReportDfMinDist.groupby(['Potential'])
    ReportDfMinDistGroupedPotentialDirection = ReportDfMinDist.groupby(
        ['Potential',
         'Direction', 'DirectionMb',
         'Reparameterization',
         'X', 'Y', 'Z', 'R',
         'UnitX', 'UnitY', 'UnitZ',
         'Theta', 'Phi',
         'a1', 'a2', 'a3', 'h'])


    # os.system("pause")
    # </editor-fold>

    # <editor-fold desc="Plotting-TDE Min-bar">
    print("Plotting-TDE Min-bar")
    Plotting = False
    if Plotting:
        import math
        # print(ReportDfMinDist['Potential'].unique())
        Title = "ReportDfMinDistGroupedPotentialDirection-Bar"
        PotentialList = ["M2R", "M3R", "BMD192R", "TabGap1", "TabGap2"]
        ReportDfMinDistFiltered = ReportDfMinDist[ReportDfMinDist['Potential'].isin(PotentialList)]
        ReportDfMinDistFiltered = ReportDfMinDistFiltered.reset_index()
        ReportDfMinDistFiltered['DirectionMb'] = ReportDfMinDistFiltered['DirectionMb'].apply(LatexFriendlyFunc)
        # print(ReportDfMinDistFiltered['Potential'].unique())
        order = (
            ReportDfMinDistFiltered
            .groupby("Potential")["MinDist"]
            .mean()
            .sort_values()
            .index
        )

        # Order for Directions (hue/colors)
        hue_order = (
            ReportDfMinDistFiltered
            .groupby("DirectionMb")["MinDist"]
            .mean()
            .sort_values()
            .index
        )

        plt.figure(figsize=(7, 7))
        sns.barplot(
            data=ReportDfMinDistFiltered,
            x="Potential",
            y="MinDist",
            hue="DirectionMb",
            order=order,
            hue_order=hue_order
        )
        num_entries = len(ReportDfMinDistFiltered['DirectionMb'].unique())
        ncol_legend = math.ceil(num_entries / 8)
        plt.legend(
            fontsize=12,
            loc='lower center',
            bbox_to_anchor=(0.5, 0),  # Move outside the plot
            ncol=ncol_legend
        )
        ymin, ymax = plt.ylim()
        plt.yticks(np.arange(ymin, ymax + 0.5, 0.5), fontsize=24)
        plt.xticks(fontsize=20)
        plt.xlabel('')
        plt.ylabel("Min Dist (Å)", fontsize=24)
        # plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(PythonName + "-" + Title, dpi=600, bbox_inches='tight')
        plt.show()
        plt.close()
    # </editor-fold>

    # <editor-fold desc="Plotting-TDE Min Scatter - All">
    print("Plotting-TDE Min Scatter - All")
    Plotting = False
    if Plotting:
        import math

        Title = "ReportDfMinDistGroupedPotentialDirection-Scatter-All"

        # Filter potentials
        PotentialList = ["M2R", "M3R", "BME192R", "TabGap1", "TabGap2"]
        ReportDfMinDistFiltered = ReportDfMinDist[ReportDfMinDist['Potential'].isin(PotentialList)].reset_index()

        plt.figure(figsize=(6, 6))
        sns.scatterplot(
            data=ReportDfMinDistFiltered,
            x='MinDist',
            y='Energy',
            hue='Potential',
            alpha=0.75,
            # edgecolor='k',
            s=100,  # marker size
            # style = "Potential"
        )

        # Legend formatting
        num_entries = len(ReportDfMinDistFiltered['Potential'].unique())
        ncol_legend = math.ceil(num_entries / 2)
        # plt.legend(
        #     fontsize=16,
        #     loc='upper right',
        #     bbox_to_anchor=(1, 1),
        #     ncol=ncol_legend
        # )
        plt.legend('', frameon=False)

        plt.xticks(np.arange(1.4, 2.2, 0.2), fontsize=24)
        plt.yticks(fontsize=20)

        # Axis labels
        plt.xlabel("Min Dist (Å)", fontsize=24)
        plt.ylabel("TDE (eV)", fontsize=24)
        plt.tight_layout()
        # plt.xscale('log')
        # plt.yscale('log')
        plt.savefig(PythonName + "-" + Title, dpi=600, bbox_inches='tight')
        plt.show()
        plt.close()
    # </editor-fold>

    # <editor-fold desc="Plotting-TDE Min Scatter-One">
    print("Plotting-TDE Min Scatter-One")
    Plotting = False
    if Plotting:
        import math
        PotentialList = ["M2R", "M3R", "BME192R", "TabGap1", "TabGap2"]
        ReportDfMinDistFiltered = ReportDfMinDist[ReportDfMinDist['Potential'].isin(PotentialList)].reset_index()
        for Name, Group in ReportDfMinDistFiltered.groupby("Potential"):
            print(Name)
            Title = "ReportDfMinDistGroupedPotential-Scatter-One-" + str(Name)
            plt.figure(figsize=(6, 6))
            sns.scatterplot(
                data=Group,
                x='MinDist',
                y='Energy',
                hue='Direction',
                # alpha=0.75,
                palette="Set1",
                edgecolor='none',
                s=100,  # marker size
                # style = "Potential"
            )

            # Legend formatting
            num_entries = len(ReportDfMinDistFiltered['Potential'].unique())
            ncol_legend = math.ceil(num_entries / 2)
            # plt.legend(
            #     fontsize=16,
            #     loc='upper right',
            #     bbox_to_anchor=(1, 1),
            #     ncol=ncol_legend
            # )
            plt.legend('', frameon=False)

            plt.xticks(np.arange(1.4, 2.2, 0.2), fontsize=24)
            plt.yticks(np.arange(0, 150, 40), fontsize=20)
            plt.xlim([1.3, 2.2])
            plt.ylim([5, 150])
            plt.xlabel("Min Dist (Å)", fontsize=24)
            plt.ylabel("TDE (eV)", fontsize=24)
            plt.tight_layout()
            # plt.xscale('log')
            # plt.yscale('log')
            plt.savefig(PythonName + "-" + Title, dpi=600, bbox_inches='tight')
            plt.show()
            plt.close()
    # </editor-fold>

    # <editor-fold desc="Plotting-TDE Min Scatter-One with FacetGrid">
    print("Plotting-TDE Min Scatter-One with FacetGrid")
    Plotting = False
    if Plotting:
        import math

        PotentialList = ["M2R", "M3R", "BME192R", "TabGap1", "TabGap2"]
        ReportDfMinDistFiltered = ReportDfMinDist[ReportDfMinDist['Potential'].isin(PotentialList)].reset_index()
        ReportDfMinDistFiltered['DirectionMb'] = ReportDfMinDistFiltered['DirectionMb'].apply(LatexFriendlyFunc)
        Title = "ReportDfMinDistGroupedPotential-Scatter-One-Facet"

        # Create FacetGrid with 2 rows × 2 columns
        g = sns.FacetGrid(
            ReportDfMinDistFiltered,
            col="Potential",
            col_wrap=2,   # 2 columns
            height=4,
            sharex=True,
            sharey=True
        )

        # Scatterplot inside each facet
        g.map_dataframe(
            sns.scatterplot,
            x='MinDist',
            y='Energy',
            hue='DirectionMb',
            palette="Set1",
            edgecolor='none',
            s=100
        )

        # g.add_legend()
        # g._legend.set_title("DirectionMb")

        # Format ticks, limits, remove axis labels
        for ax, title in zip(g.axes.flatten(), g.col_names):
            ax.set_xticks(np.arange(1.4, 2.2, 0.2))
            ax.set_yticks(np.arange(0, 150, 40))
            ax.set_xlim([1.3, 2.2])
            ax.set_ylim([5, 150])
            ax.set_xlabel("")   # remove xlabel
            ax.set_ylabel("")   # remove ylabel
            ax.set_title(title, fontsize=18)  # clean title without "Potential = "

        plt.tight_layout()
        plt.savefig(PythonName + "-" + Title, dpi=600, bbox_inches='tight')
        plt.show()
        plt.close()
    # </editor-fold>

    # <editor-fold desc="Plotting-TDE Min Scatter-Direction - All">
    print("Plotting-TDE Min Scatter-Direction - All")
    Plotting = False
    if Plotting:
        import math

        Title = "ReportDfMinDistGroupedPotentialDirection-Scatter-Direction-All"

        # Filter potentials
        PotentialList = ["M2R", "M3R", "BME192R", "TabGap1", "TabGap2"]
        DirectionList = [
                         # "dir_0_1_2",
                         # "dir_-1_2_1",
                         # "dir_-1_3_2",
                         # "dir_-1_4_2",
                         # "dir_-2_3_3",
                         # "dir_-4_4_1",
                         # "dir_-2_3_0",
                         # "dir_-1_2_3",
                         # "dir_-3_4_0",
                         "dir_-1_2_0",
                         # "dir_-3_4_3",
                         # "dir_-1_4_3",
                         # "dir_-1_3_3",
                         # "dir_-1_2_2",
                         # "dir_0_2_1",
                         # "dir_-1_4_0",
                         # "dir_-2_3_1",
                         # "dir_-2_3_2",
                         # "dir_-2_2_3",
                         # "dir_0_2_3",
                         # "dir_-2_4_1",
                         # "dir_-2_2_1",
                         # "dir_-3_4_1",
                         # "dir_-1_3_1",
                         # "dir_-1_4_1",
                         # "dir_-3_4_2",
                         # "dir_-1_3_0",
                         # "dir_0_4_1"
                         ]
        ReportDfMinDistFiltered = ReportDfMinDist[ReportDfMinDist['Potential'].isin(PotentialList)].reset_index()
        ReportDfMinDistFiltered = ReportDfMinDistFiltered[ReportDfMinDistFiltered['Direction'].isin(DirectionList)].reset_index()

        plt.figure(figsize=(6, 6))
        sns.scatterplot(
            data=ReportDfMinDistFiltered,
            x='MinDist',
            y='Energy',
            hue='Potential',
            alpha=0.75,
            # edgecolor='k',
            s=100,  # marker size
            # style = "Potential"
        )

        # Legend formatting
        num_entries = len(ReportDfMinDistFiltered['Potential'].unique())
        ncol_legend = math.ceil(num_entries / 2)
        # plt.legend(
        #     fontsize=16,
        #     loc='upper right',
        #     bbox_to_anchor=(1, 1),
        #     ncol=ncol_legend
        # )
        plt.legend('', frameon=False)


        plt.xticks(np.arange(1.4, 2.2, 0.1), fontsize=24)
        plt.yticks(np.arange(15, 46, 10), fontsize=24)
        plt.xlim([1.55, 1.8])
        plt.ylim([15, 45])

        # Axis labels
        plt.xlabel("TDE (eV)", fontsize=24)
        plt.ylabel("Min Dist (Å)", fontsize=24)
        plt.tight_layout()
        # plt.xscale('log')
        # plt.yscale('log')
        plt.savefig(PythonName + "-" + Title, dpi=600, bbox_inches='tight')
        plt.show()
        plt.close()
    # </editor-fold>

    # <editor-fold desc="Plotting-TDE Min Scatter-Direction - One with Regression">
    print("Plotting-TDE Min Scatter-Direction - One with Regression")
    Plotting = False
    if Plotting:
        PotentialList = ["M2R", "M3R", "BME192R", "TabGap1", "TabGap2"]
        DirectionList = ["dir_-1_2_0"]
        # DirectionList = ["dir_0_1_2"]
        # print(ReportDfMinDist['Direction'])

        ReportDfMinDistFiltered = ReportDfMinDist[
            ReportDfMinDist['Potential'].isin(PotentialList)
        ].reset_index(drop=True)

        ReportDfMinDistFiltered = ReportDfMinDistFiltered[
            ReportDfMinDistFiltered['Direction'].isin(DirectionList)
        ].reset_index(drop=True)

        for Name, Group in ReportDfMinDistFiltered.groupby("Potential"):
            print(Name)
            Title = "ReportDfMinDistGroupedPotentialDirection-Scatter-One-" + str(Name)
            plt.figure(figsize=(6, 6))

            # Linear regression fit
            slope, intercept, r_value, p_value, std_err = linregress(Group['MinDist'], Group['Energy'])
            slope_int = int(round(slope))
            intercept_int = int(round(intercept))
            print(f"{Name}: slope={slope_int}, intercept={intercept_int}")

            # Plot regression line first (behind scatter)
            x_vals = np.linspace(Group['MinDist'].min(), Group['MinDist'].max(), 100)
            y_vals = slope * x_vals + intercept
            plt.plot(x_vals, y_vals, 'r--', zorder=1)  # regression line

            # Scatter plot (above line)
            sns.scatterplot(
                data=Group,
                x='MinDist',
                y='Energy',
                hue='Potential',
                alpha=0.75,
                s=400,
                zorder=2
            )

            # Mean point for Group
            mean_x = Group['MinDist'].mean()
            mean_y = Group['Energy'].mean()
            plt.scatter(mean_x, mean_y, color="red", edgecolor = "black", marker="*", s=500, zorder=3)

            # Annotate slope & intercept
            plt.text(
                0.3, 0.95,
                f"Slope = {slope_int}\nIntercept = {intercept_int}",
                transform=plt.gca().transAxes,
                fontsize=40,
                verticalalignment='top',
                bbox=dict(facecolor="none", edgecolor="none")
            )

            # Legend formatting
            plt.legend('', frameon=False)

            # plt.xticks(np.arange(1.4, 2.2, 0.1), fontsize=50)
            # plt.yticks(np.arange(15, 46, 10), fontsize=50)
            # plt.xlim([1.55, 1.8])
            # plt.ylim([15, 45])

            # Hide spines
            ax = plt.gca()
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)

            plt.xlabel("", fontsize=24)
            plt.ylabel("", fontsize=24)
            plt.tight_layout()

            # Save & show
            plt.savefig(PythonName + "-" + Title, dpi=600, bbox_inches='tight', transparent=True)
            plt.show()
            plt.close()
    # </editor-fold>

    # <editor-fold desc="Regression">
    print("Regression")
    ReportDfMinDistRegression = []

    for (Potential, Direction,DirectionMb), group in ReportDfMinDist.groupby(["Potential", "Direction", "DirectionMb"]):
        if len(group) >= 2:  # Need at least 2 points for regression
            slope, intercept, r_value, p_value, std_err = linregress(group["MinDist"], group["Energy"])
            ReportDfMinDistRegression.append({
                "Potential": Potential,
                "Direction": Direction,
                "DirectionMb": DirectionMb,
                "Slope": slope,
                "Intercept": intercept,
                "R_value": r_value,
                "P_value": p_value,
                "Std_err": std_err
            })

    ReportDfMinDistRegression = pd.DataFrame(ReportDfMinDistRegression)
    ReportDfMinDistRegression.to_csv(PythonName + "-ReportDfMinDistRegression.csv", index=False)
    # </editor-fold>

    # <editor-fold desc="Regression-Heatmap">
    print("Regression-Heatmap")
    Plotting = False
    if Plotting:
        for Parameter in ["Slope","Intercept"]:
            DirectionList = ["dir_0_1_2"]
            ReportDfMinDistRegressionFiltered = ReportDfMinDistRegression[
                ~ReportDfMinDistRegression['Direction'].isin(DirectionList)
            ].reset_index(drop=True)
            ReportDfMinDistRegressionFiltered['DirectionMb'] = ReportDfMinDistRegressionFiltered['DirectionMb'].apply(LatexFriendlyFunc)
            Title = "ReportDfMinDistRegressionFiltered-" + Parameter

            # direction_order = (
            #     ReportDfMinDistRegression
            #     .groupby("DirectionMb")[Parameter]
            #     .mean()
            #     .sort_values()
            #     .index
            # )

            HeatMap = ReportDfMinDistRegressionFiltered.pivot(
                index="Potential",
                columns="DirectionMb",
                values=Parameter
            )
            HeatMap = HeatMap.rename(columns=str)

            # HeatMap = HeatMap[direction_order]

            plt.figure(figsize=(15, 5))
            sns.heatmap(
                HeatMap,
                annot=False,
                # fmt=".0f",
                cmap="coolwarm",
                cbar_kws={'label': str(Parameter)}
            )

            # plt.title("Slope Heatmap (Energy vs MinDist)")
            plt.xlabel("")
            plt.ylabel("")
            plt.tight_layout()
            plt.savefig(PythonName + "-" + Title, dpi=600)
            plt.show()
            plt.close()
    # </editor-fold>

    # <editor-fold desc="Plotting-TDE Min Kde">
    print("Plotting-TDE Min Kde")
    Plotting = False
    if Plotting:
        import math

        Title = "ReportDfMinDistGroupedPotentialDirection-Kde-All"

        # Filter potentials
        PotentialList = ["M2R", "M3R", "BME192R", "TabGap1", "TabGap2"]
        ReportDfMinDistFiltered = ReportDfMinDist[ReportDfMinDist['Potential'].isin(PotentialList)].reset_index()

        plt.figure(figsize=(6, 6))
        sns.kdeplot(
            data=ReportDfMinDistFiltered,
            x="MinDist", y="Energy",
            hue="Potential",
            # fill=True,
            alpha=0.8, levels=10
        )

        # Legend formatting
        num_entries = len(ReportDfMinDistFiltered['Potential'].unique())
        ncol_legend = math.ceil(num_entries / 2)
        # plt.legend(
        #     fontsize=16,
        #     loc='upper right',
        #     bbox_to_anchor=(1, 1),
        #     ncol=ncol_legend
        # )
        plt.legend('', frameon=False)

        plt.xticks(np.arange(1.4, 2.2, 0.2), fontsize=24)
        plt.yticks(np.arange(0, 140, 20), fontsize=20)
        plt.xlim([1.3, 2.2])
        plt.ylim([5, 130])

        # Axis labels
        plt.xlabel("Min Dist (Å)", fontsize=24)
        plt.ylabel("TDE (eV)", fontsize=24)
        plt.tight_layout()
        # plt.xscale('log')
        # plt.yscale('log')
        plt.savefig(PythonName + "-" + Title, dpi=600, bbox_inches='tight')
        plt.show()
        plt.close()
    # </editor-fold>

    # <editor-fold desc="Plotting-TDE-Min-Kde-One-Mean">
    print("Plotting-TDE-Min-Kde-One-Mean")
    Plotting = False
    if Plotting:
        import math
        import numpy as np
        import seaborn as sns
        import matplotlib.pyplot as plt

        PotentialList = ["M2R", "M3R", "BME192R", "TabGap1", "TabGap2"]
        ReportDfMinDistFiltered = ReportDfMinDist[ReportDfMinDist['Potential'].isin(PotentialList)].reset_index()

        for Name, Group in ReportDfMinDistFiltered.groupby("Potential"):
            print(Name)
            Title = "ReportDfMinDistGroupedPotentialDirection-Kde-One-" + str(Name)
            plt.figure(figsize=(6, 6))

            # KDE density
            sns.kdeplot(
                data=Group,
                x="MinDist", y="Energy",
                alpha=0.8, levels=10
            )

            # Overlay stars for each Direction (from grouped means)
            GroupMean = ReportDfMinDistStat[ReportDfMinDistStat['Potential'] == Name]
            GroupMean = GroupMean.nsmallest(2, 'MinDist_var')
            DirectionList = ["dir_0_0_1", "dir_-1_2_0", "dir_1_0_0", "dir_-1_1_0",'dir_-1_2_0','dir_-1_2_3']
            GroupMean = GroupMean[GroupMean['Direction'].isin(DirectionList)].reset_index()
            GroupMean['DirectionMb'] = GroupMean['DirectionMb'].apply(LatexFriendlyFunc)
            # putting stars for each mean
            # plt.scatter(
            #     GroupMean['MinDist_mean'],
            #     GroupMean['Energy_mean'],
            #     marker="*",
            #     color="red",
            #     s=400,
            #     edgecolor="black",
            #     linewidth=1.2,
            #     zorder=3
            # )

            for i, (_, row) in enumerate(GroupMean.iterrows()):
                shift = (-1) ** i * 0.2  # alternates: +0.1, -0.1, +0.1, ...
                plt.text(
                    row['MinDist_mean'] + shift,  # apply alternating horizontal shift
                    row['Energy_mean'] + shift,
                    row['DirectionMb'],
                    fontsize=24,
                    color="black",
                    weight="bold",
                    ha="center", va="center"
                )
                plt.scatter(
                    row['MinDist_mean'],
                    row['Energy_mean'],
                    marker="*",
                    color="red",
                    s=400,
                    edgecolor="black",
                    linewidth=1.2,
                    zorder=3
                )

            # Legend formatting
            plt.legend('', frameon=False)
            plt.title(Name, fontsize=30)

            plt.xticks(np.arange(1.4, 2.2, 0.2), fontsize=30)
            plt.yticks(np.arange(0, 140, 40), fontsize=30)
            plt.xlim([1.3, 2.2])
            plt.ylim([5, 140])

            ax = plt.gca()
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)

            # plt.xlabel("Min Dist (Å)", fontsize=24)
            # plt.ylabel("TDE (eV)", fontsize=24)
            plt.xlabel("", fontsize=24)
            plt.ylabel("", fontsize=24)
            plt.tight_layout()

            # Save & show
            plt.savefig(PythonName + "-" + Title, dpi=600, bbox_inches='tight', transparent=True)
            plt.show()
            plt.close()
    # </editor-fold>

    # <editor-fold desc="Plotting-TDE Min Convex">
    print("Plotting-TDE Min Convex")
    Plotting = False
    if Plotting:
        import math
        from scipy.spatial import ConvexHull

        Title = "ReportDfMinDistGroupedPotentialDirection-Convex"

        # Filter potentials
        PotentialList = ["M2R", "M3R", "BME192R", "TabGap1", "TabGap2"]
        ReportDfMinDistFiltered = ReportDfMinDist[ReportDfMinDist['Potential'].isin(PotentialList)].reset_index()

        plt.figure(figsize=(6, 6))


        for Name, Group in ReportDfMinDistFiltered.groupby("Potential"):
            hull = ConvexHull(Group[["MinDist", "Energy"]])
            plt.fill(
                Group.iloc[hull.vertices]["MinDist"],
                Group.iloc[hull.vertices]["Energy"],
                alpha=0.5, label=Name
            )

        # Legend formatting
        num_entries = len(ReportDfMinDistFiltered['Potential'].unique())
        ncol_legend = math.ceil(num_entries / 2)
        # plt.legend(
        #     fontsize=16,
        #     loc='upper right',
        #     bbox_to_anchor=(1, 1),
        #     ncol=ncol_legend
        # )
        plt.legend('', frameon=False)

        plt.xticks(np.arange(1.4, 2.2, 0.2), fontsize=24)
        plt.yticks(fontsize=20)

        # Axis labels
        plt.xlabel("Min Dist (Å)", fontsize=24)
        plt.ylabel("TDE (eV)", fontsize=24)
        plt.tight_layout()
        # plt.xscale('log')
        # plt.yscale('log')
        plt.savefig(PythonName + "-" + Title, dpi=600, bbox_inches='tight')
        plt.show()
        plt.close()
    # </editor-fold>

    # </editor-fold>
# </editor-fold>

# <editor-fold desc="######################################## Formation Energy">
print("######################################## Formation Energy")  # %%

# <editor-fold desc="***** Read">
print("***** Read")
Header= ["Potential","Type","NoAtoms","EnergyInitial","EnergyFinal","FormationEnergy","AbsEnergyInitialPerAtom"]
FormationSiaEam = pd.read_csv(
    "D:/Queens_University/MME/Project/Zr/SiaFormation/20240819-EquiBased/20250909-FormationDf.csv",
    usecols=["Potential", "Type", "FormationEnergy"]
)
# </editor-fold>

# <editor-fold desc="***** Merge">
print("***** Merge")
# print(TdeWs)
TdeWsFormationSiaEam = pd.merge(FormationSiaEam, TdeWs[["Potential", "Energy"]], on="Potential", how="left")
TdeWsFormationSiaEam.to_csv(PythonName + "-TdeWsFormationSiaEam.csv", index=False)

TdeWsGroupedPotentialDescribeFormationSiaEam = pd.merge(FormationSiaEam, TdeWsGroupedPotentialDescribe, on="Potential", how="left")
TdeWsGroupedPotentialDescribeFormationSiaEam.to_csv(PythonName + "-TdeWsGroupedPotentialDescribeFormationSiaEam.csv", index=False)

TdeWsGroupedPotentialDirectionDescribeFormationSiaEam = pd.merge(FormationSiaEam, TdeWsGroupedPotentialDirectionDescribe, on="Potential", how="left")
TdeWsGroupedPotentialDirectionDescribeFormationSiaEam.to_csv(PythonName + "-TdeWsGroupedPotentialDirectionDescribeFormationSiaEam.csv", index=False)
# </editor-fold>

# <editor-fold desc="***** Pivot">
print("***** Pivot")
# print(TdeWsFormationSiaEam.columns)
# print(TdeWsFormationSiaEam.head())
TdeWsFormationSiaEamPivotMin = TdeWsFormationSiaEam.pivot_table(
    index='Potential',
    columns=['Type'],
    values='Energy',
    aggfunc='min'
).reset_index()
TdeWsFormationSiaEamPivotMean = TdeWsFormationSiaEam.pivot_table(
    index='Potential',
    columns=['Type'],
    values='Energy',
    aggfunc='mean'
).reset_index()
TdeWsFormationSiaEamPivotMax = TdeWsFormationSiaEam.pivot_table(
    index='Potential',
    columns=['Type'],
    values='Energy',
    aggfunc='max'
).reset_index()

TdeWsFormationSiaEamPivotMin.to_csv(PythonName + "-TdeWsFormationSiaEamPivotMin.csv", index=False)
TdeWsFormationSiaEamPivotMean.to_csv(PythonName + "-TdeWsFormationSiaEamPivotMean.csv", index=False)
TdeWsFormationSiaEamPivotMax.to_csv(PythonName + "-TdeWsFormationSiaEamPivotMax.csv", index=False)
# </editor-fold>

# <editor-fold desc="***** Plot-Tde range">
print("***** Plot-Tde range")
Plotting = False

if Plotting:
    PotentialList = ["M2R", "M3R", "BMD192R", "TabGap1", "TabGap2"]
    TypeList = [
    "BC",
    "BO",
    "C",
    "O",
    "T",
    # "BT",
        ]

    TdeWsFormationSiaEamFiltered = (TdeWsFormationSiaEam[TdeWsFormationSiaEam['Potential'].isin(PotentialList)].reset_index(drop=True))
    TdeWsFormationSiaEamFiltered = (TdeWsFormationSiaEamFiltered[TdeWsFormationSiaEamFiltered['Type'].isin(TypeList)].reset_index(drop=True))

    g = sns.catplot(
        data=TdeWsFormationSiaEamFiltered,
        kind="box",
        x="FormationEnergy",
        y="Energy",
        hue="Type",
        col="Potential",
        col_wrap=2,
        native_scale=True,
        height=4,
        aspect=1.0,
        width=500
    )

    g.set_axis_labels("Formation Energy (eV)", "TDE (eV)")
    g.set_titles("{col_name}")
    g.add_legend(title="Type")

    # Example vertical line at FormationEnergy = 25
    # for ax in g.axes.flat:
    #     ax.axvline(25, color=".3", dashes=(2, 2))

    g.savefig(PythonName + "-TdeWsFormationSiaEamFiltered", dpi=600, bbox_inches="tight", transparent=True)

    plt.show()
# </editor-fold>

# <editor-fold desc="***** Plot: Potential-Type">
print("***** Plot: Potential-Type")
Plotting = False
if Plotting:
    StatList = ["mean", "min"]
    for Stat in StatList:
        Title = "TdeWsGroupedPotentialDescribeFormationSiaEamFiltered-" + Stat
        PotentialList = ["M2R", "M3R", "BMD192R", "TabGap1", "TabGap2"]
        # PotentialList = ["M2", "M3", "BMD192", "M2R", "M3R", "BMD192R", "TabGap1", "TabGap2"]
        TypeList = [
        "BC",
        "BO",
        "C",
        "O",
        "T",
        # "BT",
            ]

        TdeWsGroupedPotentialDescribeFormationSiaEamFiltered = (TdeWsGroupedPotentialDescribeFormationSiaEam[TdeWsGroupedPotentialDescribeFormationSiaEam['Potential'].isin(PotentialList)].reset_index(drop=True))
        TdeWsGroupedPotentialDescribeFormationSiaEamFiltered = (TdeWsGroupedPotentialDescribeFormationSiaEamFiltered[TdeWsGroupedPotentialDescribeFormationSiaEamFiltered['Type'].isin(TypeList)].reset_index(drop=True))
        plt.figure(figsize=(7, 5))
        sns.scatterplot(
            data=TdeWsGroupedPotentialDescribeFormationSiaEamFiltered,
            x=Stat,
            y="FormationEnergy",
            hue="Potential",
            style="Type",
            s=200,
            legend = "full"
        )

        plt.legend(
            frameon=False,
            loc="center left",
            bbox_to_anchor=(1.02, 0.5)  # adjust as needed
        )

        # plt.xticks(np.arange(1.4, 2.2, 0.1), fontsize=50)
        # plt.yticks(np.arange(15, 46, 10), fontsize=50)
        # plt.xlim([1.55, 1.8])
        # plt.ylim([15, 45])

        # Hide spines
        # ax = plt.gca()
        # ax.spines['right'].set_visible(False)
        # ax.spines['top'].set_visible(False)

        plt.xlabel(Stat + " TDE (eV)", fontsize=24)
        plt.ylabel("Formation Energy (eV)", fontsize=24)
        plt.tight_layout()

        # Save & show
        plt.savefig(PythonName + "-" + Title, dpi=600, bbox_inches='tight', transparent=True)
        plt.show()
        plt.close()

        plt.show()
# </editor-fold>

# <editor-fold desc="***** Plot: Potential-Type-Direction">
print("***** Plot: Potential-Type-Direction")
Plotting = False
if Plotting:
    StatList = ["mean", "min"]
    for Stat in StatList:
        Title = "TdeWsGroupedPotentialDirectionDescribeFormationSiaEamFiltered-" + Stat
        PotentialList = ["M2R", "M3R", "BMD192R", "TabGap1", "TabGap2"]
        # PotentialList = ["M2", "M3", "BMD192", "M2R", "M3R", "BMD192R", "TabGap1", "TabGap2"]
        TypeList = [
        "BC",
        "BO",
        "C",
        "O",
        "T",
        # "BT",
            ]
        DirectionList = ["dir_0_0_1", "dir_-1_2_0", "dir_1_0_0", "dir_-1_1_0",'dir_-1_2_0','dir_-1_2_3']

        TdeWsGroupedPotentialDirectionDescribeFormationSiaEamFiltered = (TdeWsGroupedPotentialDirectionDescribeFormationSiaEam[TdeWsGroupedPotentialDirectionDescribeFormationSiaEam['Potential'].isin(PotentialList)].reset_index(drop=True))
        TdeWsGroupedPotentialDirectionDescribeFormationSiaEamFiltered = (TdeWsGroupedPotentialDirectionDescribeFormationSiaEamFiltered[TdeWsGroupedPotentialDirectionDescribeFormationSiaEamFiltered['Type'].isin(TypeList)].reset_index(drop=True))
        # TdeWsGroupedPotentialDirectionDescribeFormationSiaEamFiltered = (TdeWsGroupedPotentialDirectionDescribeFormationSiaEamFiltered[TdeWsGroupedPotentialDirectionDescribeFormationSiaEamFiltered['Direction'].isin(DirectionList)].reset_index(drop=True))
        TdeWsGroupedPotentialDirectionDescribeFormationSiaEamFiltered['DirectionMb'] = TdeWsGroupedPotentialDirectionDescribeFormationSiaEamFiltered['DirectionMb'].apply(LatexFriendlyFunc)
        # Create the relplot (FacetGrid) with square facets
        # plt.figure(figsize=(12, 7))  # Your preferred figure size
        g = sns.relplot(
            data=TdeWsGroupedPotentialDirectionDescribeFormationSiaEamFiltered,
            x=Stat,
            y="FormationEnergy",
            hue="Potential",
            style="Type",
            col="DirectionMb",
            col_wrap=6,
            kind="scatter",
            height=4,  # height of each facet
            aspect=0.75,  # ensures square facets
            s=200  # increase marker size
        )

        # Remove "DirectionMb = " from facet titles
        for ax in g.axes.flat:
            ax.set_title(ax.get_title().replace("DirectionMb = ", ""), fontsize=16)

        # Place legend outside further to the right
        g._legend.set_bbox_to_anchor((1.115, 0.5))  # moved further right
        g._legend.set_frame_on(False)
        for text in g._legend.get_texts():
            text.set_fontsize(18)

            # Remove individual axis labels
        g.set_axis_labels("", "")

        # Shared X label below the last row, further away
        g.fig.text(0.5, 0.05, Stat + " TDE (eV)", ha='center', fontsize=24)

        # Shared Y label left of the grid, further away
        g.fig.text(0.03, 0.5, "Formation Energy (eV)", va='center', rotation='vertical', fontsize=24)

        # Adjust layout to prevent overlap with legend and labels
        plt.tight_layout(rect=[0.05, 0.05, 0.9, 0.95])

        # Save figure
        plt.savefig(PythonName + "-" + Title, dpi=600, bbox_inches='tight', transparent=False)#, bbox_extra_artists=(g._legend,))
        plt.show()
        plt.close()

# </editor-fold>

# </editor-fold>
