# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 13:15:19 2024
@author: evanp
Updated: September 2026 (Gcarnes05)
========================================================================================
INPUT: target_vertices.yaml, energies_correction.csv, defect_poscars.yaml
OUTPUT: Charge Defect Plot with all defects at all specified points -- Optional: Single defect plots
========================================================================================
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import yaml
import argparse
import math


# Create output folder
def create_output_folder(folder_name: str):
    if not os.path.exists(folder_name):
        os.mkdir(folder_name)


# Parse a chemical potential in the form Element=Energy
def parse_mu_entry(value):
    if "=" not in value:
        raise argparse.ArgumentTypeError(f"Invalid -mu entry '{value}'. Expected Element=Energy.")

    element, energy = (part.strip() for part in value.split("=", 1))

    if not element:
        raise argparse.ArgumentTypeError(f"Invalid -mu entry '{value}': element is empty.")

    try:
        energy = float(energy)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"Invalid energy in -mu entry '{value}'.") from exc

    return element, energy


# Store chemical potentials in a dictionary
def parse_mu_entries(entries):
    mu = {}

    for element, energy in entries:
        if element in mu:
            raise ValueError(f"Chemical potential for element '{element}' was provided more than once.")
        mu[element] = energy

    return mu


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
            raise ValueError(f"Phase '{phase_name}' has an invalid chem_pot section.")

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
    for phase_name, chem_pot_block in zip(phases, phase_data):

        missing = [
            element
            for element in elements
            if element not in chem_pot_block]

        if missing:
            raise ValueError(
                f"Phase '{phase_name}' is missing chemical potentials "
                f"for: {', '.join(missing)}"
            )

    # Arrange chemical potentials in a predictable order
    chem_pot = []

    for chem_pot_block in phase_data:

        for element in elements:

            chem_pot.append(float(chem_pot_block[element]))

    print(f"Chemical potential values: {len(chem_pot)}")

    return phases, elements, chem_pot


# POSCAR reader
def read_poscar(poscar_path: str):

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
        raise ValueError(f"Could not read atom counts from POSCAR: {poscar_path}")

    if len(element_names) != len(atom_counts):
        raise ValueError(f"Number of elements does not match number of atom counts " f"in {poscar_path}")

    return element_names, atom_counts


# Read all neutral defect compositions from defect_poscars.yaml
def read_defect_yaml(yaml_path: str):
    with open(yaml_path, "r") as file:
        data = yaml.safe_load(file)

    if not isinstance(data, dict) or "defects" not in data:
        raise ValueError("defect_poscars.yaml must contain a top-level 'defects' dictionary.")

    defects = data["defects"]

    if not isinstance(defects, dict) or not defects:
        raise ValueError("No defects were found in defect_poscars.yaml.")

    for defect_name, defect in defects.items():

        if not isinstance(defect, dict) or "elements" not in defect or "counts" not in defect:
            raise ValueError(f"Defect '{defect_name}' must contain 'elements' and 'counts'.")

        if len(defect["elements"]) != len(defect["counts"]):
            raise ValueError(f"Defect '{defect_name}' has mismatched elements and counts.")

        try:
            defect["counts"] = [int(x) for x in defect["counts"]]
        except (TypeError, ValueError) as exc:
            raise ValueError(f"Defect '{defect_name}' contains invalid atom counts.") from exc

    return defects


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

        # Positive values mean atoms were added, negative values mean atoms were removed
        delta_N[element] = defect_count - bulk_count

    return delta_N


# Get effective chemical potentials
def get_effective_chemical_potentials(
    elements, chem_pot, reference_mu, phase_index, num_phases, selected_elements=None
):

    effective_mu = {}
    selected_elements = elements if selected_elements is None else selected_elements

    for element in selected_elements:

        element_index = elements.index(element)

        # Chemical potential
        delta_mu = float(chem_pot[phase_index * len(elements) + element_index])

        # Add the elemental reference energy supplied with -mu
        if element not in reference_mu:
            raise ValueError(f"No elemental reference energy was supplied for "f"{element}. Add it to -mu.")

        effective_mu[element] = (float(reference_mu[element]) + delta_mu)

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

    required_columns = [
        "Defect Name",
        " Charge",
        " Bulk Energy",
        " Correction Energy",
        " Delta V",
        " Std Deviation"
    ]

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

    # Make sure all defect charges are integers
    for i, q in enumerate(df[" Charge"][1:], start=2):

        if not float(q).is_integer():
            raise ValueError(f"Non-integer charge at row {i}: {q}")

    # Make sure no defect/charge combinations are repeated
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

    parser = argparse.ArgumentParser(description="Charge defect formation-energy plotter",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-plotsingledefect", nargs="?", type=bool, default=False, help="Generates individual formation energy vs Fermi energy plots for each defect")
    parser.add_argument("-poscar", nargs="?", default="./POSCAR", help="Legacy single-defect POSCAR path")
    parser.add_argument("-defectyaml", nargs="?", default="./defect_poscars.yaml", help="YAML file containing all neutral defect compositions")
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
    parser.add_argument("-colors", nargs="+", default=["red", "blue", "orange", "purple", "gray", "pink", "olive", "cyan", "green"], help="Color array for charge neutrality plot")
    parser.add_argument("-legloc", nargs="?", default=8, help="Sets the location of the legend in charge neutrality plot")
    parser.add_argument("-hse", nargs=2, type=float, help="HSE band gap and VBM")
    parser.add_argument("--save_as", nargs="?", default="combinedDefects", help="Custom filename (without extension) for the saved formation energy plot")
    parser.add_argument("-bg", type=float, required=True, help="Band gap")
    parser.add_argument("-vbm", type=float, required=True, help="VBM offset")
    parser.add_argument("-mu", nargs="+", type=parse_mu_entry, required=True, metavar="ELEMENT=VALUE", help="Per-atom bulk energies, e.g. Ga=-3.20 N=-8.10",)

    args = parser.parse_args()
    config = vars(args)

    # Validate files
    for path_key in ["bulkposcar", "defectyaml", "correction", "chempot"]:

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
            "chemical potential for every element in every phase.")

    # Read BULK and all DEFECT compositions
    bulk_elements, bulk_counts = read_poscar(config["bulkposcar"])
    defects = read_defect_yaml(config["defectyaml"])

    # Validate defect elements against target_vertices.yaml
    all_defect_elements = []

    for defect_name, defect in defects.items():

        for element in defect["elements"]:

            if element not in all_defect_elements:
                all_defect_elements.append(element)

        missing_mu = sorted(set(defect["elements"]) - set(elements))

        if missing_mu:
            raise ValueError(
                f"Defect '{defect_name}' contains elements without chemical potentials "
                f"in target_vertices.yaml: {', '.join(missing_mu)}"
            )

    # Validate chemical potentials against defect_poscars.yaml
    mu = parse_mu_entries(config["mu"])

    missing = sorted(set(all_defect_elements) - set(mu))
    extra = sorted(set(mu) - set(all_defect_elements))

    if missing or extra:

        problems = []

        if missing:
            problems.append(f"missing: {', '.join(missing)}")

        if extra:
            problems.append(f"unexpected: {', '.join(extra)}")

        raise ValueError("Chemical potentials do not match defect POSCAR elements ("+ "; ".join(problems) + ").")

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

    # Replace PBE band information with HSE values when supplied
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
        for i in range(int(iterations))
    ]

    # Plot setup
    colors = config["colors"]

    if len(colors) == 0:
        raise ValueError("At least one color must be provided")

    lineStyles = ["solid", (0, (5, 7)), "dotted", "dashdot", "dashed"]

    ylimmax = config["ymax"]
    ylimmin = config["ymin"]
    xlimmax = config["xmax"]
    xlimmin = config["xmin"]

    bulkEnergy = float(energies_final.iloc[0, 2])

    # Determine which defects are processed from defect_poscars.yaml
    defect_names = list(defects.keys())

    energy_names = set(energies_final["Defect Name"].iloc[1:].unique())
    missing_energy_defects = [
        name for name in defect_names
        if name not in energy_names
    ]

    if missing_energy_defects:
        raise ValueError(
            "No energy entries were found for defect(s): "
            + ", ".join(missing_energy_defects)
        )

    print("Defects found in defect_poscars.yaml:")

    for defect in defect_names:
        print(f"    {defect}")

    print()

    # Prevent the same composition block from being printed once per phase
    printed_compositions = set()
    printed_defect_sections = 0

    # Store formation energies and transition levels for CSV output
    formation_energy_rows = []
    transition_level_rows = []

    # Process each chemical potential phase
    for p, phase_name in enumerate(phases):

        # Arrays used for plotting and charge neutrality
        completeGraph = []
        completeMinCharge = []

        namesArray = []
        defectSpots = []

        colorName = []
        finalColorNames = []

        # Process each defect
        for defect_name in defect_names:

            defect_rows = energies_final[
                energies_final["Defect Name"] == defect_name
            ]

            if len(defect_rows) == 0:
                continue

            # Print phase and defect information
            if printed_defect_sections > 0:
                print("\n" + "=" * 80 + "\n")

            print(f"Phase: {phase_name}")
            print(f"Defect: {defect_name}")
            printed_defect_sections += 1

            # Determine composition change and degeneracy
            defect = defects[defect_name]

            delta_N = compare_compositions(
                bulk_elements,
                bulk_counts,
                defect["elements"],
                defect["counts"]
            )

            degeneracy = determine_degeneracy(
                bulk_elements,
                bulk_counts,
                delta_N
            )

            if defect_name not in printed_compositions:

                print_composition_information(bulk_elements, bulk_counts, defect["elements"], defect["counts"], delta_N)

                printed_compositions.add(defect_name)

            # Only calculate chemical potentials for elements that changed
            changed_elements = [
                element
                for element, change in delta_N.items()
                if change != 0
            ]

            effective_mu = get_effective_chemical_potentials(elements=elements, chem_pot=chem_pot, reference_mu=mu, phase_index=p, num_phases=numPhases, selected_elements=changed_elements)

            # Calculate the chemical potential contribution
            chemical_potential_term = 0.0

            for element in changed_elements:

                change = delta_N[element]

                chemical_potential_term -= ( change * effective_mu[element])

            print("Effective chemical potentials:")

            for element in sorted(effective_mu):
                print(
                    f"    {element:<5} "
                    f"{effective_mu[element]:.6f} eV"
                )

            print("Chemical potential terms:")

            for element, change in delta_N.items():

                if change == 0:
                    continue

                term = -change * effective_mu[element]

                print(
                    f"    {element}: "
                    f"Delta_N = {change:+d}, "
                    f"mu = {effective_mu[element]:.6f}, "
                    f"contribution = {term:.6f} eV"
                )

            print()

            # Charge-state formation energies
            defect_graphs = []
            defect_charges = []

            print(
                f"Defect Formation Energies at VBM ({E_f}) in eV "
                f"[phase: {phase_name}]:"
            )

            for _, row in defect_rows.iterrows():

                bulkDefectEnergy = float(row[" Bulk Energy"])
                q = int(row[" Charge"])
                V = float(row[" Delta V"])
                correction = float(row[" Correction Energy"])

                # Base formation energy
                finalDefectEnergy = (bulkDefectEnergy - bulkEnergy + chemical_potential_term)

                # Add charge and correction terms
                finalDefectEnergy += (q * (E_f + V) + correction)

                print(
                    f"{defect_name}_{q:<3} "
                    f"{finalDefectEnergy:>12.6f} eV"
                )

                # Store formation energy for CSV output
                formation_energy_rows.append({
                    "Phase": phase_name,
                    "Defect": defect_name,
                    "Charge": q,
                    "Formation_Energy_eV": finalDefectEnergy
                })

                # Formation energy vs Fermi energy
                graph = []
                charges = []

                for k in range(int(iterations)):

                    graph.append(finalDefectEnergy + q * stepSize * k)

                    charges.append(q)

                defect_graphs.append(graph)
                defect_charges.append(charges)

            # Find charge state transition levels
            old_charge = None

            for m in range(len(fermiEnergies)):

                energies_at_fe = [
                    defect_graphs[q][m]
                    for q in range(len(defect_graphs))
                ]

                charges_at_fe = [
                    defect_charges[q][m]
                    for q in range(len(defect_charges))
                ]

                minimum_energy = min(energies_at_fe)
                minimum_index = energies_at_fe.index(minimum_energy)
                minimum_charge = charges_at_fe[minimum_index]

                completeGraph.append(minimum_energy)
                completeMinCharge.append(minimum_charge)
                defectSpots.append(degeneracy)

                # Detect when the lowest energy charge state changes
                if old_charge is not None and minimum_charge != old_charge:

                    transition_energy = fermiEnergies[m]

                    print(
                        f"Transition from {old_charge:2d} "
                        f"to {minimum_charge:2d} "
                        f"at {fermiEnergies[m]:.5f} eV"
                    )

                    transition_level_rows.append({
                        "Phase": phase_name,
                        "Defect": defect_name,
                        "Charge_1": old_charge,
                        "Charge_2": minimum_charge,
                        "Transition_Fermi_Energy_eV": transition_energy
                    })

                old_charge = minimum_charge

            namesArray.append(defect_name)

            # Group defects by the first part of their name for coloring
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

                # Plot the lowest energy charge state for this defect
                defect_minimum_graph = [
                    min(
                        defect_graphs[q][m]
                        for q in range(len(defect_graphs))
                    )
                    for m in range(len(fermiEnergies))
                ]

                plt.plot(fermiEnergies, defect_minimum_graph, label=format_label(str(defect_name)))

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

        plt.xlim(xlimmin, xlimmax)
        plt.ylim(ylimmin, ylimmax)

        # Determine intrinsic Fermi level
        kT = config["kT"]

        qValue = None
        previous_sign = None

        for i, fermi_energy in enumerate(fermiEnergies):

            Q = 0.0
            qArray = []

            # Calculate total charge contribution from each defect
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

            # Find where the total charge changes sign
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

            # Print Q at a user specified Fermi energy
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

            start = i * len(fermiEnergies)
            end = (i + 1) * len(fermiEnergies)

            tempData = completeGraph[start:end]

            color_index = i % len(colors)
            line_style_index = (i // len(colors)) % len(lineStyles)

            plt.plot(
                fermiEnergies,
                tempData,
                label=format_label(namesArray[i]),
                color=colors[color_index],
                linestyle=lineStyles[line_style_index]
            )

        # Intrinsic Fermi level
        if qValue is not None:

            plt.axvline(qValue, color="black", linestyle="dashed")

        # Legend and save
        plt.legend(loc=config["legloc"])

        plot_name = config["save_as"]
        saveLocation = f"{save_folder}/{plot_name}.png"

        plt.savefig(saveLocation)
        plt.show()

        # Save formation energies to CSV
        formation_energies_df = pd.DataFrame(formation_energy_rows)
        formation_energies_df.to_csv("formation_energies.csv", index=False)

        # Save transition levels to CSV
        transition_levels_df = pd.DataFrame(transition_level_rows)
        transition_levels_df.to_csv("transition_levels.csv", index=False)

        print("Formation energies saved to formation_energies.csv")
        print("Transition levels saved to transition_levels.csv")

# Run
if __name__ == "__main__":
    main()
