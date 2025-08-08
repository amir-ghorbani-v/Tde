# <editor-fold desc="######################################## OPTIONS">
print("######################################## OPTIONS")
# <editor-fold desc="**********  Library">
# %%
print("***** Library")
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import os
import numpy as np
import pandas as pd
import shutil
import matplotlib
import seaborn as sns
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib import colormaps
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
matplotlib.use('Qt5Agg')
# </editor-fold>

# <editor-fold desc="**********  Variables">
print("**********  Variables")
PythonName = "20250709-BoxSize"
font = {"family": "serif", #serif" "sans-serif" "cursive" "fantasy" "monospace"
        "weight": "bold",
        "size": 10}
matplotlib.rc("font", **font)
matplotlib.rcParams["lines.linewidth"] = 3
matplotlib.rcParams["axes.spines.right"] = True
matplotlib.rcParams["axes.spines.top"] = True
matplotlib.rcParams['figure.figsize'] = (8, 8)
plt.rc('font', size=25)          # controls default text size
plt.rc('axes', titlesize=25)     # fontsize of the title
plt.rc('axes', labelsize=25)     # fontsize of the x and y labels
plt.rc('xtick', labelsize=20)    # fontsize of the x-axis ticks
plt.rc('ytick', labelsize=20)    # fontsize of the y-axis ticks
plt.rc('legend', fontsize=20)    # fontsize of the legend
Colors = ["b", "r", "y", "m", "g", "k", "c", "w", "#FF5733", "#33FF57", "#5733FF", "#FF5733", "#33FF57", "#5733FF", "#FF5733", "#33FF57", "#5733FF", "#FF5733", "#33FF57", "#5733FF"]
Colors = ["b", "red", "green", "olive", "m", "lime", "#FF5733", "#33FF57", "#5733FF", "#FF5733", "#33FF57", "#5733FF", "#FF5733", "#33FF57", "#5733FF", "#FF5733", "#33FF57", "#5733FF"]
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
# </editor-fold>

# </editor-fold>

# <editor-fold desc="######################################## Reading">
print("######################################## Reading")
# <editor-fold desc="***** Reading">
# %%
print("***** Reading")
SizeList = [1,2, 3, 4, 5, 10, 15, 20, 40]
PotentialList = ['M2R']
TryArray = np.arange(1, 42)
DirectionList = ['dir_0_0_1', 'dir_-1_2_0','dir_1_0_0']
BoxAtomMap = {
1:8,
2:64,
3:198,
4:512,
5:1000,
10:8000,
15:26550,
20:64000,
40:512000,
}

BoxSize = pd.DataFrame()

Active = False
if Active:
    for Potential in PotentialList:
        for Size in SizeList:
            for Try in TryArray:
                for Direction in DirectionList:
                    base_path = os.path.join(Potential, str(Size), str(Try), Direction)
                    if not os.path.isdir(base_path):
                        continue

                    for Energy in os.listdir(base_path):
                        Energy_path = os.path.join(base_path, Energy)
                        if os.path.isdir(Energy_path):

                            Active = True
                            if Active:
                                File = os.path.join(Energy_path, '20250130-Tde-Defect.csv')
                                df = pd.read_csv(File)
                                df.dropna(inplace=True)
                                df['DefectDist'] = df['DefectDist'].astype(int)
                                df['DefectWs'] = df['DefectWs'].astype(int)

                            else:
                                df = pd.read_csv(copied_csv)

                            df['Size'] = Size
                            df['Potential'] = Potential
                            df['Try'] = Try
                            df['Direction'] = Direction
                            df['Energy'] = Energy

                            BoxSize = pd.concat([BoxSize, df], ignore_index=True)

    BoxSize = BoxSize.rename(columns={'DefectDist': 'DT', 'DefectWs': 'WS'})
    BoxSize.to_csv('20250709-BoxSize.csv', index=False)
else:
    BoxSize = pd.read_csv("20250709-BoxSize.csv")
# </editor-fold>

# <editor-fold desc="***** Miller to Miller-Bravais">
print("***** Miller to Miller-Bravais")
# %%
# print(Tde)
Mapping = {
    'dir_-2_3_3': 'dir_-2_3_-1_3',
    'dir_-1_4_1': 'dir_-1_4_-3_1',
    'dir_-1_3_1': 'dir_-1_3_-2_1',
    'dir_-1_2_2': 'dir_-1_2_-1_2',
    'dir_0_2_1': 'dir_0_2_-2_1',
    'dir_-1_2_0': 'dir_-1_2_-1_0', # Dense
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
    'dir_1_0_0': 'dir_1_0_-1_0', #Dense
    'dir_0_0_1': 'dir_0_0_0_1', #Dense
        
        
}

# Tde = Tde[Tde["File"] == "Tde"].copy()

BoxSize['DirectionMb'] = BoxSize['Direction'].replace(Mapping)
# print(Tde)
# </editor-fold>

# <editor-fold desc="***** Separation">
print("***** Separation")
# %%
DefectDfMap = {
    "WS": BoxSize[["WS", "Size", "Potential", "Try", "Direction", "Energy"]],
    "DT": BoxSize[["DT", "Size", "Potential", "Try", "Direction", "Energy"]]
}
# </editor-fold>
# </editor-fold>

# <editor-fold desc="######################################## Energy Evolution">
print("######################################## Energy Evolution")

# <editor-fold desc="***** Plot-histogram">
print("***** Plot-histogram")
# %%
Active = False
if Active:
    for Name, Data in DefectDfMap.items():
        Data["Size"] = Data["Size"].astype(str)
        print(Name)

        for Potential, group_df in Data.groupby("Potential"):
            print(Potential)
            plt.figure(figsize=(6, 6))
            Title = Name + "-" + Potential

            cmap = plt.cm.get_cmap("jet", 7)
            Colors = [cmap(i) for i in range(cmap.N)]
            unique_sizes = Data["Size"].unique()
            unique_sizes_sorted = sorted(unique_sizes, key=lambda x: int(x))
            palette_dict = {k: Colors[i % len(Colors)] for i, k in enumerate(unique_sizes_sorted)}

            sns.histplot(
                data=group_df,
                x=Name,
                hue="Size",
                kde=True,
                palette=palette_dict,
                bins=10,
                stat='density',
                element='step',
                legend=False
            )

            handles = [
                Patch(facecolor=color, edgecolor='black', label=label)
                for label, color in palette_dict.items()
            ]

            plt.xticks([0, 1])  # <-- restrict x-axis ticks to 0 and 1
            plt.tight_layout()

            plt.xlabel(f'{Name} (eV)')
            plt.ylabel('Density')
            plt.legend(
                handles=handles,
                bbox_to_anchor=(0.99, 0.75),
                loc='center',
                ncol=1,
                frameon=False,
                prop={'weight': 'normal', 'size': 16}
            )
            plt.tight_layout()
            plt.savefig("20250709-" + Title + ".png", dpi=600, bbox_inches='tight')
            plt.show()
            plt.close()
# </editor-fold>
# </editor-fold>

# <editor-fold desc="######################################## Tde">
print("######################################## Tde")

# <editor-fold desc="***** Min">
print("***** Min")
# %%
Active = True
if Active:
    DfEnergyMin = []
    DfEnergyMinDic = {}


    for Key, Df in DefectDfMap.items():
        print(Key)
        Df = Df[Df[Key] == 1]
        DfEnergyMinNew = Df.groupby(["Size", "Potential", "Direction"], as_index=False)["Energy"].min()
        DfEnergyMinNew["SizeAtom"] = DfEnergyMinNew["Size"].map(BoxAtomMap)
        DfEnergyMinNew["Type"] = Key
        print(DfEnergyMinNew)
        DfEnergyMin.append(DfEnergyMinNew)
        DfEnergyMinDic[Key] = DfEnergyMinNew

        # os.system("pause")

    DfEnergyMin = pd.concat(DfEnergyMin, ignore_index=True)
    DfEnergyMin["SizeAtom"] = DfEnergyMin["SizeAtom"].astype(int)
    DfEnergyMin['DirectionMb'] = DfEnergyMin['Direction'].replace(Mapping)
    DfEnergyMin['DirectionMbLatex'] = DfEnergyMin['DirectionMb'].apply(LatexFriendlyFunc)


    # print(DfEnergyMinDic)
    # print(DfEnergyMin)
    DfEnergyMin.to_csv(PythonName + '-DfEnergyMin.csv', index=False)
else:
    DfEnergyMin = pd.read_csv(PythonName + '-DfEnergyMin.csv')

# </editor-fold>

# <editor-fold desc="***** Plot-Line-MinEnergy per direction">
print("***** Plot-Line-MinEnergy per direction")
# %%
Active = False
if Active:
    for Potential, Group in DfEnergyMin.groupby("Potential"):
        print(f"Potential: {Potential}")

        for Direction, dir_df in Group.groupby("Direction"):
            print(f"  Direction: {Direction}")

            plt.figure(figsize=(6, 6))
            Title = "MinEnergy" + "-" + Potential + "-" + Direction

            sns.lineplot(
                data=dir_df.sort_values("SizeAtom"),
                x="SizeAtom",
                y="Energy",
                hue="Type",
                marker="o",
                linewidth=4,
                markersize=20,
                palette=Colors
            )

            custom = [
                # Line2D([], [], marker='o', markersize=15, color=Colors[0], linestyle='None', label='Classic MD'),
                # Line2D([], [], marker='*', markersize=15, color=Colors[1], linestyle='None', label='UTTM'),
                Line2D([], [], color=Colors[0], linestyle='-', label='Wigner-Seitz'),
                Line2D([], [], color=Colors[1], linestyle='-', label='Displacement'),
            ]

            plt.xlabel("Number of Atoms")
            plt.ylabel("Minimum Energy (eV)")
            # plt.title(Title)
            plt.xscale("log")
            plt.legend(handles=custom, title="Type", fontsize=20)
            plt.tight_layout()
            plt.savefig(PythonName + "-" + Title + ".png", dpi=600, bbox_inches="tight")
            plt.show()
            plt.close()
# </editor-fold>

# <editor-fold desc="***** Plot-Line-MinEnergy">
print("***** Plot-Line-MinEnergy")
# %%
Active = True
if Active:
    for Potential, Group in DfEnergyMin.groupby("Potential"):
        print(f"Potential: {Potential}")

        plt.figure(figsize=(6, 6))
        Title = f"MinEnergy-{Potential}"

        ax = sns.lineplot(
            data=Group.sort_values("SizeAtom"),
            x="SizeAtom",
            y="Energy",
            hue="DirectionMbLatex",
            style="Type",
            style_order=["DT", "WS"],
            markers=True,
            linewidth=4,
            markersize=20,
            palette=Colors,
        )

        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles=handles[1:], labels=labels[1:])
        ax.yaxis.set_major_locator(MultipleLocator(5))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))
        plt.xlabel("Number of Atoms")
        plt.ylabel("min TDE (eV)")
        plt.xscale("log")
        plt.tight_layout()
        plt.savefig(f"{PythonName}-{Title}.png", dpi=600, bbox_inches="tight")
        plt.show()
        plt.close()

# </editor-fold>

# <editor-fold desc="***** Plot-Line-MinEnergy">
print("***** Plot-Line-MinEnergy")
# %%
Active = False
if Active:
    for Potential, Group in DfEnergyMin.groupby("Potential"):
        print(f"Potential: {Potential}")

        plt.figure(figsize=(6, 6))
        Title = f"MinEnergy-{Potential}"

        ax = plt.gca()

        # First plot solid lines ("WS")
        sns.lineplot(
            data=Group[Group["Type"] == "WS"].sort_values("SizeAtom"),
            x="SizeAtom",
            y="Energy",
            hue="Direction",
            linestyle="-",
            markers='x',
            linewidth=4,
            markersize=20,
            palette=Colors,
            ax=ax,
            legend=False,
        )

        # Then plot dashed lines ("DT")
        sns.lineplot(
            data=Group[Group["Type"] == "DT"].sort_values("SizeAtom"),
            x="SizeAtom",
            y="Energy",
            hue="Direction",
            linestyle="--",
            markers='o',
            linewidth=4,
            markersize=20,
            palette=Colors,
            ax=ax,
        )

        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles=handles[1:], labels=labels[1:])

        ax.yaxis.set_major_locator(MultipleLocator(5))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))
        plt.xlabel("Number of Atoms")
        plt.ylabel("min TDE (eV)")
        plt.xscale("log")
        plt.tight_layout()
        plt.savefig(f"{PythonName}-{Title}.png", dpi=600, bbox_inches="tight")
        plt.show()
        plt.close()
# </editor-fold>

# </editor-fold>
