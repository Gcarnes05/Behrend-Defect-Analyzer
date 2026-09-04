# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 13:15:19 2024
@author: evanp
Updated: April 2026 (Gcarnes05)
========================================================================================
INPUT: target_vertices.yaml, energies_correction.csv
OUTPUT: Charge Defect Plot with all defects at all specified points in .yaml file -- Optional: Single defect plots
========================================================================================
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import yaml
import argparse
import math

def create_output_folder(folder_name: str):
    if not os.path.exists(folder_name):
        os.mkdir(folder_name)


# Chemical potentials
def read_chemical_potentials(yaml_path: str):

    with open(yaml_path, "r") as file:
        data = yaml.safe_load(file)

    if not isinstance(data, dict):
        raise ValueError("target_vertices.yaml must contain a YAML dictionary.")

    phases = []
    phase_data = []

    # Find every top level phase containing a chem_pot dictionary
    for phase_name, phase_contents in data.items():

        if not isinstance(phase_contents, dict):
            continue

        if "chem_pot" not in phase_contents:
            continue

        chem_pot_block = phase_contents["chem_pot"]

        if not isinstance(chem_pot_block, dict):
            raise ValueError(
                f"Phase '{phase_name}' has an invalid chem_pot section."
            )

        phases.append(phase_name)
        phase_data.append(chem_pot_block)

    if len(phases) == 0:
        raise ValueError("No chemical potential phases were found in target_vertices.yaml.")

    # Determine all unique elements
    elements = []

    for chem_pot_block in phase_data:

        for element in chem_pot_block:

            if element not in elements:
                elements.append(element)

    if len(elements) == 0:
        raise ValueError("No elements were found in the chem_pot sections.")

    # Make sure every phase contains every element
    for phase_name, chem_pot_block in zip(
        phases,
        phase_data
    ):

        missing = [
            element
            for element in elements
            if element not in chem_pot_block
        ]

        if missing:
            raise ValueError(
                f"Phase '{phase_name}' is missing chemical potentials "
                f"for: {', '.join(missing)}"
            )

    # Arrange chemical potentials in a predictable order
    chem_pot = []

    for chem_pot_block in phase_data:

        for element in elements:

            chem_pot.append(
                float(chem_pot_block[element])
            )

    print(f"Chemical potential values: {len(chem_pot)}")

    return phases, elements, chem_pot


# POSCAR reader
def read_poscar(poscar_path: str):
    """
    Read element names and atom counts from a VASP POSCAR.
    """

    with open(poscar_path, "r") as f:
        poscar_lines = f.readlines()

    if len(poscar_lines) < 7:
        raise ValueError(
            f"POSCAR appears to be incomplete: {poscar_path}"
        )

    element_names = poscar_lines[5].split()

    try:
        atom_counts = [int(x) for x in poscar_lines[6].split()]
    except ValueError:
        raise ValueError(
            f"Could not read atom counts from POSCAR: {poscar_path}"
        )

    if len(element_names) != len(atom_counts):
        raise ValueError(
            f"Number of elements does not match number of atom counts "
            f"in {poscar_path}"
        )

    return element_names, atom_counts


# Composition comparison
def compare_compositions(bulk_elements, bulk_counts, defect_elements, defect_counts):

    bulk_composition = dict(zip(bulk_elements, bulk_counts))
    defect_composition = dict(zip(defect_elements, defect_counts))

    all_elements = sorted(
        set(bulk_composition.keys()) |
        set(defect_composition.keys())
    )

    delta_N = {}

    for element in all_elements:

        bulk_count = bulk_composition.get(element, 0)
        defect_count = defect_composition.get(element, 0)

        delta_N[element] = defect_count - bulk_count

    return delta_N
    

# get chem pots
def get_effective_chemical_potentials(elements, chem_pot, reference_elements, reference_mu, phase_index, num_phases):
 
    effective_mu = {}

    for i, element in enumerate(elements):

        # Chemical potential from YAML
        delta_mu = float(
            chem_pot[phase_index * len(elements) + i]
        )

        # Find elemental reference energy
        if element not in reference_elements:
            raise ValueError(
                f"No elemental reference energy was supplied for "
                f"{element}. Add it to -mu."
            )

        reference_index = reference_elements.index(element)

        reference_energy = float(
            reference_mu[reference_index]
        )

        effective_mu[element] = (
            reference_energy + delta_mu
        )

    return effective_mu


# Print composition information
def print_composition_information(bulk_elements, bulk_counts, defect_elements, defect_counts, delta_N):
    print("\nBulk POSCAR composition:")
    for element, count in zip(bulk_elements, bulk_counts):
        print(f"    {element:<5} {count}")

    print("\nDefect POSCAR composition:")
    for element, count in zip(defect_elements, defect_counts):
        print(f"    {element:<5} {count}")

    print("\nPOSCAR Composition change:")
    for element, change in delta_N.items():

        if change > 0:
            sign = "+"
        else:
            sign = ""

        print(f"    {element:<5} {sign}{change}")

    print()


# Determine degeneracy
def determine_degeneracy(bulk_elements, bulk_counts, delta_N):

    bulk_composition = dict(zip(bulk_elements, bulk_counts))

    # Automatically select an element that was removed
    removed_elements = [
        element
        for element, change in delta_N.items()
        if change < 0
    ]

    if len(removed_elements) == 1:

        element = removed_elements[0]

        if element in bulk_composition:
            return bulk_composition[element]

    # If no unique removed element exists, fall back to one defect
    return 1


# Format labels
def format_label(label):

    # Handle arbitrary names without underscores
    if "_" not in label:
        return label

    parts = label.split("_")

    if len(parts) == 2:
        return f"{parts[0]}$_{{{parts[1]}}}$"

    # For names such as Va_Si_N or Va_Si_O
    formatted = parts[0]

    for part in parts[1:]:
        formatted += f"$_{{{part}}}$"

    return formatted


# Validate energies_final.csv
def validate_energies_final(df: pd.DataFrame):

    required_columns = ["Defect Name", " Charge", " Bulk Energy", " Correction Energy", " Delta V", " Std Deviation"]

    for col in required_columns:
        if col not in df.columns:
            raise ValueError(f"energies_final.csv missing required column: '{col}'")

    if df[required_columns].isnull().any().any():
        raise ValueError("energies_final.csv contains NaN values")

    first = df.iloc[0]

    if first["Defect Name"].lower() != "bulk":
        raise ValueError("First row must be the bulk reference ('bulk')")

    if int(round(first[" Charge"])) != 0:
        raise ValueError("Bulk charge must be 0")

    if abs(first[" Correction Energy"]) > 1e-6:
        raise ValueError("Bulk correction energy must be 0")

    for i, q in enumerate(df[" Charge"][1:], start=2):

        if not float(q).is_integer():
            raise ValueError(f"Non-integer charge at row {i}: {q}")

    duplicates = df[["Defect Name", " Charge"]].duplicated()

    if duplicates.any():

        dup_rows = np.where(duplicates)[0] + 2

        raise ValueError(
            f"Duplicate defect/charge entries at rows {dup_rows}"
        )

    print("energies_final.csv validation passed")
    print()


# Main
def main():

    parser = argparse.ArgumentParser(description="Charge defect formation-energy plotter", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-plotsingledefect", nargs="?", type=bool, default=False, help="Generates individual formation energy vs Fermi energy plots for each defect")
    parser.add_argument("-poscar", nargs="?", default="./POSCAR", help="Path to defect POSCAR")
    parser.add_argument("-bulkposcar", nargs="?", default="../bulk/POSCAR", help="Path to perfect/bulk POSCAR")
    parser.add_argument("-correction", nargs="?", default="./energies_final.csv", help="Final correction energies file location")
    parser.add_argument("-chempot", nargs="?", default="./target_vertices.yaml", help="Desired chemical potential file location (.yaml)")
    parser.add_argument("-ymax", nargs="?", type=float, default=7, help="ymax for defect graph")
    parser.add_argument("-xmax", nargs="?", type=float, default=-3, help="xmax for defect graph")
    parser.add_argument("-ymin", nargs="?", type=float, default=-7, help="ymin for defect graph")
    parser.add_argument("-xmin", nargs="?", type=float, default=0, help="xmin for defect graph")
    parser.add_argument("-testfe", nargs="?", type=float, default=-1, help="Displays Q information and defect information at specified fermi energy")
    parser.add_argument("-kT", nargs="?", type=float, default=0.05, help="kT value")
    parser.add_argument("-printQ", nargs="?", type=bool, default=False, help="Prints Q values of all defects at intrinsic fermi level")
    parser.add_argument("-colors", nargs="+", default=["red", "green", "blue", "orange"], help="Color array for charge neutrality plot")
    parser.add_argument("-legloc", nargs="?", default=8, help="Sets the location of the legend in charge neutrality plot")
    parser.add_argument("-hse", nargs=2, type=float, help="HSE band gap and VBM")
    parser.add_argument("--save_as", nargs="?", default="combinedDefects", help="Custom filename (without extension) for the saved formation energy plot")
    parser.add_argument("-bg", type=float, required=True, help="Band gap")
    parser.add_argument("-vbm", type=float, required=True, help="VBM offset")
    parser.add_argument("-mu", nargs="+", type=float, required=True, help="Per-atom bulk energies corresponding to the elements, in target_verticies order.")
    args = parser.parse_args()
    config = vars(args)

    # Validate files
    for path_key in ["poscar","bulkposcar","correction","chempot"]:

        if not os.path.exists(config[path_key]):

            raise FileNotFoundError(f"Input file not found: {config[path_key]}")

    # Setup
    save_folder = "chargeDefectPlots"
    create_output_folder(save_folder)

    if config["plotsingledefect"]:
        single_defect_folder = os.path.join(save_folder, "singleDefects")
        create_output_folder(single_defect_folder)

    energies_final = pd.read_csv(config["correction"])

    phases, elements, chem_pot = read_chemical_potentials(config["chempot"])

    numOfElements = len(elements)
    numPhases = len(phases)

    if len(chem_pot) != numPhases * numOfElements:
        raise ValueError(
            "Chemical potential data does not contain exactly one "
            "chemical potential for every element in every phase."
        )

    
    # Read BULK and DEFECT POSCARs    
    bulk_elements, bulk_counts = read_poscar(config["bulkposcar"])
    
    defect_elements, defect_counts = read_poscar(config["poscar"])
    
    # Determine all unique elements in the system    
    all_elements = list(dict.fromkeys(bulk_elements + defect_elements))
    
    numOfElements = len(all_elements)
    
    print(f"Number of unique elements: {numOfElements}")
    print(f"Elements: {', '.join(all_elements)}")
    
    # Determine composition change    
    delta_N = compare_compositions(bulk_elements, bulk_counts, defect_elements, defect_counts)
    
    print_composition_information(bulk_elements, bulk_counts, defect_elements, defect_counts, delta_N)

    # Validate chemical potentials
    mu = config["mu"]
    
    if len(mu) == 0:
        raise ValueError("At least one elemental reference energy must be supplied.")
    
    if len(mu) != len(elements):
        raise ValueError(
            "Number of elemental reference energies supplied with -mu "
            "must match the number of elements in target_vertices.yaml."
        )

    # Determine degeneracy
    degeneracy = determine_degeneracy(bulk_elements, bulk_counts, delta_N)

    # Validate energies file
    validate_energies_final(energies_final)

    # Validate numerical arguments
    if config["bg"] <= 0:
        raise ValueError("Band gap must be positive")

    if config["kT"] <= 0:
        raise ValueError("kT must be positive")

    if config["xmin"] >= config["xmax"] and config["xmax"] != -3:
        raise ValueError("xmin must be less than xmax")

    if config["ymin"] >= config["ymax"]:
        raise ValueError("ymin must be less than ymax")

    # Band structure information
    E_f = config["vbm"]
    gap = config["bg"]

    if config["xmax"] == -3:
        config["xmax"] = gap

    if config["hse"] is not None:

        originalVBM = E_f
        originalGap = gap

        E_f = config["hse"][1]
        gap = config["hse"][0]

        config["xmax"] = gap

    # Fermi-energy grid
    stepSize = 0.0001

    iterations = gap / stepSize

    if iterations <= 0:
        raise ValueError("Invalid iteration count")

    fermiEnergies = [
        stepSize * i
        for i in range(int(iterations))]

    # Plot setup
    colors = config["colors"]

    if len(colors) == 0:
        raise ValueError("At least one color must be provided")

    lineStyles = ["solid",(0, (5, 7)),"dotted","dashdot","dashed"]

    ylimmax = config["ymax"]
    ylimmin = config["ymin"]
    xlimmax = config["xmax"]
    xlimmin = config["xmin"]

    bulkEnergy = float(energies_final.iloc[0, 2])

    # Determine which defects are present in energies_final.csv
    defect_names = list(energies_final["Defect Name"].iloc[1:].unique())

    if len(defect_names) == 0:
        raise ValueError("No defects found in energies_final.csv")

    print("Defect found in energies_final.csv:")

    for defect in defect_names:
        print(f"    {defect}")

    print()

    # Arrays used by plotting / charge neutrality
    allValues = []
    allCharges = []

    completeGraph = []
    completeMinCharge = []

    namesArray = []
    defectSpots = []

    colorName = []
    finalColorNames = []

    for p, phase_name in enumerate(phases):

        for defect_name in defect_names:

            defect_rows = energies_final[
                energies_final["Defect Name"] == defect_name
            ]

            if len(defect_rows) == 0:
                continue

            effective_mu = get_effective_chemical_potentials(elements=elements, chem_pot=chem_pot, reference_elements=elements, reference_mu=mu, phase_index=p, num_phases=numPhases)
        
            chemical_potential_term = 0.0
            
            for element, change in delta_N.items():
            
                if change == 0:
                    continue
            
                if element not in effective_mu:
                    raise ValueError(f"No chemical potential found for element {element}.")
            
                chemical_potential_term -= (change * effective_mu[element])
            
            print("Effective chemical potentials:")
            
            for element in sorted(effective_mu):
                print(f"    {element:<5} "f"{effective_mu[element]:.6f} eV")
            
            print("Chemical potential terms:")

            for element, change in delta_N.items():

                if change == 0:
                    continue

                term = -change * effective_mu[element]

                print(f"    {element}: "f"Delta_N = {change:+d}, "f"mu = {effective_mu[element]:.6f}, "f"contribution = {term:.6f} eV")

            print()

            # Charge-state formation energies
            defect_graphs = []
            defect_charges = []
            print(f"Defect Formation Energies at VBM ({E_f}) in eV:")
            for _, row in defect_rows.iterrows():

                bulkDefectEnergy = float(row[" Bulk Energy"])

                q = int(row[" Charge"])

                V = float(row[" Delta V"])

                correction = float(row[" Correction Energy"])

                # Base formation energy
                finalDefectEnergy = (bulkDefectEnergy - bulkEnergy + chemical_potential_term)

                # Charge correction
                finalDefectEnergy += (q * (E_f + V) + correction)

                print(f"{defect_name}_{q:<3} "f"{finalDefectEnergy:>12.6f} eV")

                # Formation energy vs Fermi energy
                graph = []

                charges = []

                for k in range(int(iterations)):

                    graph.append(finalDefectEnergy + q * stepSize * k)

                    charges.append(q)

                defect_graphs.append(graph)
                defect_charges.append(charges)

            old_charge = None

            for m in range(len(fermiEnergies)):

                energies_at_fe = [defect_graphs[q][m] for q in range(len(defect_graphs))]

                charges_at_fe = [defect_charges[q][m] for q in range(len(defect_charges))]

                minimum_energy = min(energies_at_fe)

                minimum_index = energies_at_fe.index(minimum_energy)

                minimum_charge = charges_at_fe[minimum_index]

                completeGraph.append(minimum_energy)

                completeMinCharge.append(minimum_charge)

                defectSpots.append(degeneracy)

                if old_charge is not None and minimum_charge != old_charge:

                    print(
                        f"Transition from {old_charge:2d} "
                        f"to {minimum_charge:2d} "
                        f"at {fermiEnergies[m]:.5f} eV"
                    )

                old_charge = minimum_charge

            namesArray.append(defect_name)

            color_group = defect_name.split("_")[0]

            colorName.append(color_group)

            if color_group not in finalColorNames:

                finalColorNames.append(color_group)

            # Individual defect plot
            if config["plotsingledefect"]:

                plt.figure(figsize=(10, 6))

                formattedTitle = format_label(str(defect_name))

                plt.title("Defect Plot of " + formattedTitle)

                plt.xlabel("Fermi Energy (eV)")

                plt.ylabel("Formation Energy (eV)")

                for graph, row in zip(defect_graphs, defect_rows.itertuples()):

                    plt.plot(fermiEnergies, graph, label=f"q = {row[2]}")

                plt.xlim(xlimmin, xlimmax)

                plt.ylim(ylimmin, ylimmax)

                plt.legend()

                saveLocation = os.path.join(single_defect_folder, str(defect_name) + ".png")

                plt.savefig(saveLocation)

                plt.show()

        # Charge neutrality calculation
        numberOfDefects = len(defect_names)

        plt.figure(figsize=(5, 7))

        plt.xlabel("Fermi Energy (eV)")

        plt.ylabel("Formation Energy (eV)")

        plt.xlim(xlimmin,xlimmax)

        plt.ylim(ylimmin,ylimmax)

        # Determine intrinsic Fermi level
        kT = config["kT"]

        e = np.exp(1)

        qValue = None

        previous_sign = None

        for i, fermi_energy in enumerate(fermiEnergies):

            Q = 0.0
            qArray = []

            for j in range(numberOfDefects):

                index = i + j * len(fermiEnergies)

                formation_energy = completeGraph[index]
                charge = int(completeMinCharge[index])
                N_i = defectSpots[index]
                Q += (N_i * charge * np.exp(-formation_energy / kT))

                # Store the cumulative Q after each defect
                qArray.append(Q)

                if (
                    config["testfe"] != -1
                    and abs(fermi_energy - config["testfe"]) < stepSize / 2
                ):
                    print("charge state of defect", j, "=", charge)
                    print("degeneracy states of defect", j, "=", N_i)
                    print("formation energy of defect", j, "=", formation_energy)
                    print()

            current_sign = Q > 0

            if previous_sign is not None and current_sign != previous_sign:

                qValue = fermi_energy

                print()
                print(
                    "Intrinsic Fermi Defect Level: "
                    f"{fermi_energy:.4f} eV"
                )
                print()

                if config["printQ"]:

                    tempQ = 0.0

                    for j in range(numberOfDefects):

                        tempQ = qArray[j] - tempQ
                        print("Q value of defect", j, "=", tempQ)

                    print()
                    print("Total Q Value =", Q)
                    print()

            if (
                config["testfe"] != -1
                and abs(fermi_energy - config["testfe"]) < stepSize / 2
            ):
                print("Q value =", Q)
                print()

            previous_sign = current_sign

        # HSE shaded regions
        if config["hse"] is not None:

            plt.fill(
                [
                    xlimmin,
                    xlimmin,
                    originalVBM - E_f,
                    originalVBM - E_f
                ],
                [
                    ylimmin,
                    ylimmax,
                    ylimmax,
                    ylimmin
                ],
                color="silver"
            )

            plt.fill(
                [
                    xlimmax,
                    xlimmax,
                    originalVBM - E_f + originalGap,
                    originalVBM - E_f + originalGap
                ],
                [
                    ylimmin,
                    ylimmax,
                    ylimmax,
                    ylimmin
                ],
                color="silver"
            )

        # Plot lowest energy charge state for every defect
        lineStyleCount = [0 for _ in range(len(finalColorNames))]

        for i in range(numberOfDefects):

            tempData = []

            start = (i * len(fermiEnergies))

            end = ((i + 1)* len(fermiEnergies))

            tempData = completeGraph[start:end]

            color_index = finalColorNames.index(colorName[i])

            plt.plot(
                fermiEnergies,
                tempData,
                label=format_label(
                    namesArray[i]
                ),
                color=colors[
                    color_index % len(colors)
                ],
                linestyle=lineStyles[
                    lineStyleCount[
                        color_index
                    ] % len(lineStyles)
                ]
            )

            lineStyleCount[color_index] += 1

        # Intrinsic Fermi level
        if qValue is not None:

            plt.axvline(
                qValue,
                color="black",
                linestyle="dashed"
            )

        # Legend and save
        plt.legend(loc=config["legloc"])

        plot_name = config["save_as"]

        saveLocation = (f"{save_folder}/{plot_name}.png")

        plt.savefig(saveLocation)

        plt.show()

# Run
if __name__ == "__main__":
    main()
