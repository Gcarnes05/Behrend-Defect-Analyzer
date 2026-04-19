# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 13:15:19 2024
@author: evanp
Updated: April 2026 (Gcarnes05)
========================================================================================
Input: vAtoms_output.csv, energies_correction.csv
Output: energies_final.csv -- Optional: vAtoms plots
========================================================================================
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import argparse
import math

# Convert command line string inputs into boolean values
def str2bool(v):
    return str(v).lower() in ("true", "1", "yes")


# Read POSCAR
def read_poscar(poscar_path):
    with open(poscar_path) as f:
        lines = f.readlines()

    elements = lines[5].split()
    sites = list(map(int, lines[6].split()))

    factor = math.gcd(*sites)
    sites = [s // factor for s in sites]

    return elements, sites


# Compute delta V correction
def compute_delta_v(sortedData, percent, number):

    last_sum = 0.0
    std_vals = []

    # Use percent based cutoff when number is not specified
    if number < 0:
        minDistance = sortedData.iloc[0, 0] * percent
        i = 0

        while i < len(sortedData) and sortedData.iloc[i, 0] > minDistance:
            last_sum += sortedData.iloc[i, 1]
            std_vals.append(sortedData.iloc[i, 1])
            i += 1

        delV = last_sum / i if i > 0 else 0.0
        cutoff_index = max(i - 1, 0)

    # Use fixed number of atoms when specified
    else:
        for i in range(min(number, len(sortedData))):
            last_sum += sortedData.iloc[i, 1]
            std_vals.append(sortedData.iloc[i, 1])

        delV = last_sum / len(std_vals)
        cutoff_index = len(std_vals) - 1

    return delV, np.std(std_vals), cutoff_index, i


# Plot vAtoms for a defect
def plot_vatoms(defect_name, c1, c2, c3, c4, sortedData,
                delV, cutoff_i, config, saveFolder, i):

    title = "vAtoms_for_" + defect_name

    plt.figure(figsize=(10, 6))
    plt.title(title)
    plt.xlabel("Radial Distance (bohr)")
    plt.ylabel("Energy (eV)")

    # Raw electrostatic components from vAtoms output
    plt.scatter(c1, c2, label="V(long-range)")
    plt.scatter(c1, c3, label="V(defect)-V(ref)")
    plt.scatter(c1, c4, label="Corrected potential")

    # Axis limits
    xmin = 0 if config["vatomsxmin"] == -100 else config["vatomsxmin"]
    xmax = sortedData.iloc[0, 0] + 1 if config["vatomsxmax"] == -100 else config["vatomsxmax"]
    ymin = plt.ylim()[0] if config["vatomsymin"] == -100 else config["vatomsymin"]
    ymax = plt.ylim()[1] if config["vatomsymax"] == -100 else config["vatomsymax"]

    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)

    # Show delta V cutoff used for correction
    plt.plot([sortedData.iloc[i, 0], xmax], [delV, delV], color='black', linestyle="dashed")
    plt.plot([sortedData.iloc[i, 0], sortedData.iloc[i, 0]], [ymin, delV], color='black', linestyle="dashed")

    plt.legend(loc='upper right')

    saveLocation = os.path.join(saveFolder, title + ".png")
    plt.tight_layout()
    plt.savefig(saveLocation)
    plt.close()


def main():
    # Argument parsing
    parser = argparse.ArgumentParser(description="Arguments for charged defect correction", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-plotvatoms", nargs='?', type=str2bool, default=True)
    parser.add_argument("-poscar", nargs='?', default="./POSCAR", help=" Path to the POSCAR file (default: ./POSCAR)")
    parser.add_argument("-vatoms", nargs='?', default="./vAtoms_output.csv", help="Path to vAtoms_output.csv (default: ./vAtoms_output.csv)")
    parser.add_argument("-correction", nargs='?', default="./energies_correction.csv", help="Path to energies_correction.csv (default: ./energies_correction.csv)")
    parser.add_argument("-vatomsymax", nargs='?', type=float, default=-100, help="Maximum y-axis for vAtoms plots")
    parser.add_argument("-vatomsxmax", nargs='?', type=float, default=-100, help="Maximum x-axis for vAtoms plots")
    parser.add_argument("-vatomsymin", nargs='?', type=float, default=-100, help="Minimum y-axis for vAtoms plots")
    parser.add_argument("-vatomsxmin", nargs='?', type=float, default=-100, help="Minimum x-axis for vAtoms plots")
    parser.add_argument("-percent", nargs='?', type=float, default=0.8, help="Fraction of the furthest atoms used to compute delta V (default: 0.8)")
    parser.add_argument("-number", nargs='?', type=int, default=-1, help="Number of furthest atoms used for delta V (default: -1)")
    parser.add_argument("-mu", nargs="+", type=float, required=True, help="Per-atom bulk energies for each element in POSCAR order")
    config = vars(parser.parse_args())

    # Create output folder for vAtoms plots
    if not os.path.exists("vAtomsImages"):
        os.mkdir("vAtomsImages")

    # Load vAtoms and correction data
    data = pd.read_csv(config["vatoms"]).astype(str)
    finalFile = pd.read_csv(config["correction"])

    # Read POSCAR
    element_names, defectSites = read_poscar(config["poscar"])
    mu = config["mu"]

    # Ensure chemical potentials match POSCAR elements
    if len(mu) != len(element_names):
        raise ValueError("Mismatch between mu values and POSCAR elements")

    print("\nChemical potentials:")
    for el, val in zip(element_names, mu):
        print(f"  μ_{el} = {val} eV")
    print()

    # Initialize storage for delta V extraction
    vatoms_defects = []
    print_records = []
    excelFile = [0]
    allDev = [0]
    defectNames = ["bulk"]
    charges = [0]

    column_buffers = ([], [], [], [])
    start = 0

    # Loop over vAtoms blocks
    while start <= len(data) - 2:

        c1, c2, c3, c4 = column_buffers
        j = start + 1

        defect_name = data.iloc[start, 1].replace("/", "")
        base_name = "_".join(defect_name.split("_")[:2])
        vatoms_defects.append(base_name)

        # Read one full vAtoms block
        while data.iloc[j, 0] != "stop":
            if data.iloc[j, 0] != "nan":
                c1.append(float(data.iloc[j, 0]))
                c2.append(float(data.iloc[j, 1]))
                c3.append(float(data.iloc[j, 2]))
                c4.append(float(data.iloc[j, 3]))
            j += 1

        # Sort by distance for far field averaging
        sortedData = pd.DataFrame(
            {"distance": c1, "values": c4},
            dtype=float
        ).sort_values("distance", ascending=False)

        # Compute delta V from far field region
        delV, std, cutoff_i, i = compute_delta_v(
            sortedData, config["percent"], config["number"]
        )

        # Extract charge state from defect name
        charge_str = defect_name.split("_")[-1]
        try:
            charge = float(charge_str)
            charge_ok = True
        except ValueError:
            charge = None
            charge_ok = False

        print_records.append({
            "defect": defect_name,
            "delta_v": delV,
            "std": std,
            "charge": charge,
            "charge_ok": charge_ok
        })

        excelFile.append(delV)
        allDev.append(std)

        # Optional vAtoms visualization
        if config["plotvatoms"]:
            plot_vatoms(defect_name, c1, c2, c3, c4,
                        sortedData, delV, cutoff_i,
                        config, "vAtomsImages", i)

        # Reset buffers for next defect
        column_buffers = ([], [], [], [])
        start = j

    # Build final corrected dataset
    for i in range(1, len(finalFile)):
        name = finalFile.iloc[i, 0].replace("/", "")
        defectNames.append("_".join(name.split("_")[:2]))

        charge_str = name.split("_")[-1]
        try:
            charges.append(float(charge_str))
        except ValueError:
            raise ValueError(f"Invalid defect name: {name}")

    # Insert delta V corrections into final output table
    finalFile = finalFile.drop(finalFile.columns[0], axis=1)
    finalFile.insert(0, "Defect Name", defectNames)
    finalFile.insert(1, " Charge", charges)
    finalFile.insert(4, " Delta V", excelFile)
    finalFile.insert(5, " Std Deviation", allDev)

    finalFile.to_csv("energies_final.csv", index=False)
    print("energies_final.csv written successfully.\n")

    # Print summary
    print("Defect summary:\n")

    for rec in print_records:
        defect = rec["defect"]
        delV = rec["delta_v"]
        std = rec["std"]
        charge = rec["charge"]
        charge_ok = rec["charge_ok"]

        qstd = round(charge * std, 5) if charge_ok else "N/A"

        print(f"{defect:<12} delta V={round(delV,5):>6}, std={round(std,5):>6}, q*std={qstd}")

        if not charge_ok:
            print(f"Warning: Unparsed charge in '{defect}'")
        elif abs(charge * std) >= 0.1:
            print(f"Warning: |q*std| greater than or equal to 0.1 for {defect}")

    print("\nComplete")

if __name__ == "__main__":
    main()
