# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 13:15:19 2024
@author: evanp
Updated: September 2026 (Gcarnes05)
========================================================================================
Input: vAtoms_output.csv, energies_correction.csv
Output: defect_poscars.yaml, energies_final.csv -- Optional: vAtoms plots
========================================================================================
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import argparse
import yaml


# Convert command line string inputs into boolean values
def str2bool(v):
    return str(v).lower() in ("true", "1", "yes")


# Read element names and atom counts from a POSCAR
def read_poscar(poscar_path):
    with open(poscar_path) as f:
        lines = f.readlines()
    elements = lines[5].split()
    sites = list(map(int, lines[6].split()))
    return elements, sites


# Find every neutral defect POSCAR
def find_neutral_defect_poscars(root):
    defects = {}

    for current_root, dirs, _ in os.walk(root):
        for directory in dirs:
            # Only use directories ending in "_0"
            if not directory.endswith("_0"):
                continue

            poscar = os.path.join(current_root, directory, "POSCAR")
            if not os.path.isfile(poscar):
                continue

            # Remove the final "_0" to get the defect name
            defect_name = directory.rsplit("_", 1)[0]

            if defect_name in defects:
                raise ValueError(f"Multiple neutral POSCARs found for defect '{defect_name}'.")

            elements, counts = read_poscar(poscar)

            # Store information needed for formation energy calculations
            defects[defect_name] = {
                "poscar": os.path.relpath(poscar, root),
                "elements": elements,
                "counts": counts,
            }

    if not defects:
        raise FileNotFoundError(f"No <defectname>_0/POSCAR directories found under {root}")

    return defects


# Write neutral defect information to a YAML file
def write_defect_yaml(defects, yaml_path):
    with open(yaml_path, "w") as file:
        yaml.safe_dump({"defects": defects}, file, sort_keys=False)


# Compute the delta V correction
def compute_delta_v(sortedData, percent, number):
    last_sum = 0.0
    std_vals = []

    # Use a percentage based cutoff when number is not specified
    if number < 0:
        minDistance = sortedData.iloc[0, 0] * percent
        i = 0

        while i < len(sortedData) and sortedData.iloc[i, 0] > minDistance:
            last_sum += sortedData.iloc[i, 1]
            std_vals.append(sortedData.iloc[i, 1])
            i += 1

        delV = last_sum / i if i > 0 else 0.0
        cutoff_index = max(i - 1, 0)

    # Use a fixed number of atoms when specified
    else:
        i = 0
        while i < min(number, len(sortedData)):
            last_sum += sortedData.iloc[i, 1]
            std_vals.append(sortedData.iloc[i, 1])
            i += 1

        delV = last_sum / len(std_vals) if len(std_vals) > 0 else 0.0
        cutoff_index = max(len(std_vals) - 1, 0)

    return delV, np.std(std_vals), cutoff_index, i


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


# Store all chemical potentials in a dictionary
def parse_mu_entries(entries):
    mu = {}

    for element, energy in entries:
        if element in mu:
            raise ValueError(f"Chemical potential for element '{element}' was provided more than once.")
        mu[element] = energy

    return mu


# Plot vAtoms for a defect
def plot_vatoms(defect_name, c1, c2, c3, c4, sortedData, delV, cutoff_i, config, saveFolder, i):
    title = "vAtoms_for_" + defect_name
    plt.figure(figsize=(10, 6))
    plt.title(title)
    plt.xlabel("Radial Distance (bohr)")
    plt.ylabel("Energy (eV)")

    plt.scatter(c1, c2, label="V(long-range)")
    plt.scatter(c1, c3, label="V(defect)-V(ref)")
    plt.scatter(c1, c4, label="Corrected potential")

    xmin = 0 if config["vatomsxmin"] == -100 else config["vatomsxmin"]
    xmax = sortedData.iloc[0, 0] + 1 if config["vatomsxmax"] == -100 else config["vatomsxmax"]
    ymin = plt.ylim()[0] if config["vatomsymin"] == -100 else config["vatomsymin"]
    ymax = plt.ylim()[1] if config["vatomsymax"] == -100 else config["vatomsymax"]

    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)

    # Show the delta V cutoff used for the correction
    plt.plot([sortedData.iloc[cutoff_i, 0], xmax], [delV, delV], color="black", linestyle="dashed")
    plt.plot([sortedData.iloc[cutoff_i, 0], sortedData.iloc[cutoff_i, 0]], [ymin, delV], color="black", linestyle="dashed")

    plt.legend(loc="upper right")
    plt.tight_layout()
    plt.savefig(os.path.join(saveFolder, title + ".png"))
    plt.close()


def main():

    parser = argparse.ArgumentParser(description="Arguments for charged defect correction", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-plotvatoms", nargs="?", type=str2bool, default=True)
    parser.add_argument("-defectdirectory", nargs="?", default=".", help="Directory containing <defectname>_0 directories")
    parser.add_argument("-defectyaml", nargs="?", default="defect_poscars.yaml", help="Output YAML for neutral defect POSCARs")
    parser.add_argument("-vatoms", nargs="?", default="./vAtoms_output.csv", help="Path to vAtoms_output.csv (default: ./vAtoms_output.csv)")
    parser.add_argument("-correction", nargs="?", default="./energies_correction.csv", help="Path to energies_correction.csv (default: ./energies_correction.csv)")
    parser.add_argument("-vatomsymax", nargs="?", type=float, default=-100, help="Maximum y-axis for vAtoms plots")
    parser.add_argument("-vatomsxmax", nargs="?", type=float, default=-100, help="Maximum x-axis for vAtoms plots")
    parser.add_argument("-vatomsymin", nargs="?", type=float, default=-100, help="Minimum y-axis for vAtoms plots")
    parser.add_argument("-vatomsxmin", nargs="?", type=float, default=-100, help="Minimum x-axis for vAtoms plots")
    parser.add_argument("-percent", nargs="?", type=float, default=0.8, help="Fraction of the furthest atoms used to compute delta V (default: 0.8)")
    parser.add_argument("-number", nargs="?", type=int, default=-1, help="Number of furthest atoms used for delta V (default: -1)")
    parser.add_argument("-mu", nargs="+", type=parse_mu_entry, required=True, metavar="ELEMENT=VALUE", help="Per-atom bulk energies, e.g. Ga=-3.20 N=-8.10")
    config = vars(parser.parse_args())

    # Create output folder for vAtoms plots
    if not os.path.exists("vAtomsImages"):
        os.mkdir("vAtomsImages")

    # Load vAtoms and correction data
    data = pd.read_csv(config["vatoms"]).astype(str)
    finalFile = pd.read_csv(config["correction"])

    # Find neutral defect POSCARs and create defect_poscars.yaml
    defects = find_neutral_defect_poscars(config["defectdirectory"])
    write_defect_yaml(defects, config["defectyaml"])

    # Find every unique element in the defect POSCARs
    all_element_names = []

    for defect in defects.values():
        for element in defect["elements"]:
            if element not in all_element_names:
                all_element_names.append(element)

    # Read and validate chemical potentials
    mu = parse_mu_entries(config["mu"])

    missing = sorted(set(all_element_names) - set(mu))
    extra = sorted(set(mu) - set(all_element_names))

    # Make sure the supplied chemical potentials match the POSCAR elements
    if missing or extra:
        problems = []

        if missing:
            problems.append(f"missing: {', '.join(missing)}")

        if extra:
            problems.append(f"unexpected: {', '.join(extra)}")

        raise ValueError("Chemical potentials do not match POSCAR elements (" + "; ".join(problems) + ").")

    print(f"\nFound {len(defects)} neutral defect POSCAR(s).")
    print(f"{config['defectyaml']} written successfully.\n")
    print("Chemical potentials:")

    for el in all_element_names:
        print(f" μ_{el} = {mu[el]} eV")

    # Initialize storage
    print_records = []
    excelFile = [0]
    allDev = [0]
    defectNames = ["bulk"]
    charges = [0]
    column_buffers = ([], [], [], [])
    start = 0

    # Process each defect in vAtoms_output.csv
    while start <= len(data) - 2:
        c1, c2, c3, c4 = column_buffers
        j = start + 1

        # Read defect name and charge from the final "_"
        defect_name = data.iloc[start, 1].replace("/", "")
        charge_str = defect_name.rsplit("_", 1)[1]

        try:
            charge = int(charge_str)
        except ValueError:
            raise ValueError(f"Invalid charge in defect name '{defect_name}'. Expected an integer after the final '_'.")

        # Read one full vAtoms block
        while data.iloc[j, 0] != "stop":
            if data.iloc[j, 0] != "nan":
                c1.append(float(data.iloc[j, 0]))
                c2.append(float(data.iloc[j, 1]))
                c3.append(float(data.iloc[j, 2]))
                c4.append(float(data.iloc[j, 3]))
            j += 1

        # Sort atoms from furthest to closest to the defect
        sortedData = pd.DataFrame({"distance": c1, "values": c4}, dtype=float).sort_values("distance", ascending=False)

        # Compute delta V
        delV, std, cutoff_i, i = compute_delta_v(sortedData, config["percent"], config["number"])

        # Store corrections for the final CSV file
        print_records.append({"defect": defect_name, "delta_v": delV, "std": std, "charge": charge})
        excelFile.append(delV)
        allDev.append(std)

        # Create a vAtoms plot if requested
        if config["plotvatoms"]:
            plot_vatoms(defect_name, c1, c2, c3, c4, sortedData, delV, cutoff_i, config, "vAtomsImages", i)

        # Reset buffers before reading the next defect
        column_buffers = ([], [], [], [])
        start = j

    # Build defect names and charges for energies_final.csv
    for i in range(1, len(finalFile)):
        name = finalFile.iloc[i, 0].replace("/", "")

        # Everything before the final "_" is the defect name
        defect_name = name.rsplit("_", 1)[0]

        # Everything after the final "_" is the charge
        charge_str = name.rsplit("_", 1)[1]

        try:
            charge = int(charge_str)
        except ValueError:
            raise ValueError(f"Invalid defect name: {name}. Expected an integer charge after the final '_ '.")

        defectNames.append(defect_name)
        charges.append(charge)

    # Add defect names, charges, delta V, and standard deviation to the final table
    finalFile = finalFile.drop(finalFile.columns[0], axis=1)
    finalFile.insert(0, "Defect Name", defectNames)
    finalFile.insert(1, " Charge", charges)
    finalFile.insert(4, " Delta V", excelFile)
    finalFile.insert(5, " Std Deviation", allDev)

    # Save the final corrected dataset
    finalFile.to_csv("energies_final.csv", index=False)

    # Print the correction summary
    print("\nDefect summary:")

    for rec in print_records:
        defect = rec["defect"]
        delV = rec["delta_v"]
        std = rec["std"]
        charge = rec["charge"]
        qstd = round(charge * std, 5)

        print(f"{defect:<15}" f"delta V={delV:>9.5f}, " f"std={std:>9.5f}, " f"q*std={qstd:>9.5f}")

        # Warn when the charge multiplied by the standard deviation is large
        if abs(charge * std) >= 0.1:
            print(f"Warning: |q*std| greater than or equal to 0.1 for {defect}")

    print("\nenergies_final.csv written successfully.\n")


if __name__ == "__main__":
    main()
