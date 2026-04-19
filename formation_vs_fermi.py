# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 13:15:19 2024
@author: evanp
Updated: April 2026 (Gcarnes05)
========================================================================================
Input: target_vertices.yaml, energies_correction.csv
Output: Charge Defect Plot with all defects at all specified points in .yaml file -- Optional: Single defect plots
========================================================================================
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import yaml
import argparse
import math


# Create output directory if it does not exist
def create_output_folder(folder_name: str):
    if not os.path.exists(folder_name):
        os.mkdir(folder_name)


# Read chemical potentials from YAML configuration file
def read_chemical_potentials(yaml_path: str):
    with open(yaml_path, "r") as file:
        data = yaml.safe_load(file)

    if not isinstance(data, dict):
        raise ValueError("YAML file is not properly formatted")

    elements = []
    chem_pot = []

    counter1 = 0
    counter2 = 0

    # Parse nested YAML structure
    for x in data:
        if counter1 > 0:
            for y in data[x]:
                if counter2 < 1:
                    for z in data[x][y]:
                        elements.append(z)
                        chem_pot.append(data[x][y][z])
                counter2 += 1
        counter1 += 1
        counter2 = 0

    if len(elements) == 0 or len(chem_pot) == 0:
        raise ValueError("No chemical potentials found in YAML file")

    return elements, chem_pot


# Read POSCAR
def read_poscar(poscar_path: str):
    with open(poscar_path, "r") as f:
        poscar_lines = f.readlines()

    element_names = poscar_lines[5].split()
    defect_sites = [int(x) for x in poscar_lines[6].split()]

    factor = math.gcd(*defect_sites)
    defect_sites = [x / factor for x in defect_sites]

    return element_names, defect_sites, factor
    
# Function to format labels for plotting with subscript
def format_label(label):
    base, subscript = label.split('_')
    return f"{base}$_{{{subscript}}}$"


# Validate energies final file
def validate_energies_final(df: pd.DataFrame):
    required_columns = ["Defect Name"," Charge"," Bulk Energy"," Correction Energy"," Delta V"," Std Deviation"]

    # Column check
    for col in required_columns:
        if col not in df.columns:
            raise ValueError(f"energies_final.csv missing required column: '{col}'")

    # NaN check
    if df[required_columns].isnull().any().any():
        raise ValueError("energies_final.csv contains NaN values")

    # Bulk row checks
    first = df.iloc[0]
    if first["Defect Name"].lower() != "bulk":
        raise ValueError("First row must be the bulk reference ('bulk')")

    if int(round(first[" Charge"])) != 0:
        raise ValueError("Bulk charge must be 0")

    if abs(first[" Correction Energy"]) > 1e-6:
        raise ValueError("Bulk correction energy must be 0")

    # Defect naming
    for i, name in enumerate(df["Defect Name"][1:], start=2):
        if "_" not in name:
            raise ValueError(f"Invalid defect name format at row {i}: '{name}'")

    # Charge must be integer
    for i, q in enumerate(df[" Charge"][1:], start=2):
        if not float(q).is_integer():
            raise ValueError(f"Non-integer charge at row {i}: {q}")

    # Duplicate (defect, charge) check
    duplicates = df[["Defect Name", " Charge"]].duplicated()
    if duplicates.any():
        dup_rows = np.where(duplicates)[0] + 2
        raise ValueError(f"Duplicate defect/charge entries at rows {dup_rows}")

    print("energies_final.csv validation passed")
    print()


def main():
    # Argument parsing
    parser = argparse.ArgumentParser(description="Arguments for charge defect ", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-plotsingledefect", nargs='?', type=bool, default = False, help="Generates individual formation energy vs Fermi energy plots for each defect")
    parser.add_argument("-poscar", nargs='?', default = "./POSCAR", help="poscar file location, POSCAR is used to determine number of defect sites so use a version that includes all defects you are looking at")
    parser.add_argument("-correction", nargs='?', default = "./energies_final.csv", help="final correction energies and defect names file location")
    parser.add_argument("-chempot", nargs='?', default = "./target_vertices.yaml", help="chemical potential file location (.yaml)")
    parser.add_argument("-ymax", nargs='?', type=float, default = 7, help="ymax for defect graph")
    parser.add_argument("-xmax", nargs='?', type=float, default = -3, help="xmax for defect graph")
    parser.add_argument("-ymin", nargs='?', type=float, default = -7, help="ymin for defect graph")
    parser.add_argument("-xmin", nargs='?', type=float, default = 0, help="xmin for defect graph")
    parser.add_argument("-testfe", nargs='?', type=float, default = -1, help="displayes Q information and defect information at specified fermi energy")
    parser.add_argument("-kT", nargs='?', type=float, default = 0.05, help="kT value")
    parser.add_argument("-printQ", nargs='?', type=bool, default = False, help="prints Q values of all defects at intrinsic fermi level")
    parser.add_argument("-colors", nargs='+', default=["red", "green", "blue", "orange"], help="color array for charge neutrality plot")
    parser.add_argument("-legloc", nargs='?', default = 8, help="sets the location of the legend in charge neutrality plot")
    parser.add_argument("-hse", nargs=2, type=float, help="enter in values for band gap and VBM for HSE calculation to generate PBE 'prediction'")
    parser.add_argument("--save_as", nargs="?", default="combinedDefects", help="Custom filename (without extension) for the saved formation energy (charge neutrality) plot")
    parser.add_argument("-bg", type=float, required=True, help="Band Gap")
    parser.add_argument("-vbm", type=float, required=True, help="VBM Offset")
    parser.add_argument("-mu", nargs="+", type=float, required=True, help="Per-atom bulk energies for each element in POSCAR order")
    args = parser.parse_args()
    config = vars(args)
    
    # Validate input files exist
    for path_key in ["poscar", "correction", "chempot"]:
        if not os.path.exists(config[path_key]):
            raise FileNotFoundError(f"Input file not found: {config[path_key]}")

    #Declaration of arrays for charge neutrality plot
    graphValues = [] #Values stored for plot, used as a temp to analyze lowest energy charge state
    minCharge = [] #Array that stores the minimum charge at given "fermi energy"
    fermiEnergies = [] #Array that contains "x" values, based on the number of steps set above

    # Defaults and Setup
    if config["xmax"] == -3:
        config["xmax"] = config["bg"]

    save_folder = "chargeDefectPlots"
    create_output_folder(save_folder)
    
    # Read Input Files
    energies_final = pd.read_csv(config["correction"])
    
    elements, chem_pot = read_chemical_potentials(config["chempot"])
    element_names, defect_sites, factor = read_poscar(config["poscar"])

    mu = config["mu"]

    # Validate energies final
    validate_energies_final(energies_final)

    # Validate chemical potentials
    if len(mu) != len(element_names):
        raise ValueError("Number of chemical potentials must match POSCAR elements")
    
    # Validate numeric arguments
    if config["bg"] <= 0:
        raise ValueError("Band gap (-bg) must be positive")

    if config["kT"] <= 0:
        raise ValueError("kT must be positive")

    if config["xmin"] >= config["xmax"]:
        raise ValueError("xmin must be less than xmax")

    if config["ymin"] >= config["ymax"]:
        raise ValueError("ymin must be less than ymax")

    if config["testfe"] != -1 and (config["testfe"] < config["xmin"] or config["testfe"] > config["xmax"]):
        print("Warning: testfe is outside plotting range")
    

    E_f = config['vbm'] # Fermi Energy (eV)
    gap = config['bg'] # Band Gap (eV)

    #This determines if the hse tag was set to T/F and adjusts calues accordingly
    if(config["hse"] != None):
        originalVBM = E_f
        originalGap = gap
        
        E_f = config["hse"][1]
        gap = config["hse"][0]
        
        config['xmax'] = gap

    stepSize = 0.0001

    # Calculate and validate iterations
    iterations = gap / stepSize
    if iterations <= 0:
        raise ValueError("Invalid iteration count (check band gap and step size)")

    # Generates all Fermi Energies for Plot
    fermiEnergies = [stepSize * i for i in range(int(iterations))]

    # Define and validate colors
    colors = config["colors"]
    if len(colors) == 0:
        raise ValueError("At least one color must be provided")
    
    # Plot styling setup
    lineStyles = ["solid", (0, (5, 7)), "dotted", "dashdot", "dashed"]

    ylimmax = config["ymax"]
    ylimmin = config["ymin"]
    xlimmax = config["xmax"]
    xlimmin = config["xmin"]

    bulkEnergy = float(energies_final.iloc[0, 2])

    count = 0
    
    # Gets the name of the first defect for plots
    storedName = energies_final.iloc[1,0]

    # Variable for determining number of spots
    oldElement = " " 

    #Declartion of variables/arrays for charge neutrality plot
    numOfElements = 0 #Counts the number of "unique" elements
    tempArray = [] #Used to calculate number of elements
    tempValue = 0 #Used to calculate number of elements
    allValues = [] #Stores all of the defect formation energies calculated below at each fermi level
    allCharges = [] #Stores the charge state for each fermi level and defect
    colorName = [] #Used to keep color coding for the plot
    degenArray = [] #Array that stores the number of degeneracy states in supercell/primitive cell
    finalColorNames = [] #Used to keep color coding for the plot

    #Determing the number of unique elements
    for i in range (0,len(elements)):
        for j in range (0, len(tempArray)):
            if(tempArray[j] == elements[i]):
                tempValue = tempValue + 1
        if (tempValue == 0):
            tempArray.append(elements[i])
    numOfElements = len(tempArray)
    del (i, j, tempArray, tempValue)

    #Large for loop that will calculate charge neutrality of system based on given chemical potential points given in .yaml file
    for p in range(0, int(len(elements)/numOfElements)):
        oldIndex = 0 #Used to determine when charge state switches accross fermi levels
        elementNames = []
        elementEPA = [] #Stores energy needed to add/subtract specific atom from defect
        
        lineStyleCount = [] #Used to ensure that all of the same color lines gave different line styles
        for i in range (0, len(colors)):
            lineStyleCount.append(0)
        
        for j in range(0, numOfElements):
            elementNames.append(str(elements[numOfElements*p + j]))
            elementEPA.append(float(chem_pot[numOfElements*p + j]) + mu[j])
        
        print(f"{'Elements:':<15} {elementNames}")
        print(f"{'mu:':<15} {[f'{x:.4f}' for x in mu]}")
        print(f"{'Delta mu:':<15} {[f'{chem_pot[numOfElements*p + j]:.4f}' for j in range(numOfElements)]}")
        print(f"{'Effective mu:':<15} {[f'{x:.4f}' for x in elementEPA]}")
        print(f"Defect Formation Energies at VBM ({E_f:.4f}) in eV:")

        #Declaration of arrays for charge neutrality plot
        completeGraph = [] #Stores defect energies for all of the defects across the fermi energies
        namesArray = [] #Stores the names of the defects
        completeMinCharge = [] #Stores minimum charge state for all of the defects across the fermi energies
        defectSpots = [] #Number of locations for a specific defect 
        
        for i in range (1, len(energies_final)):
            bulkDefectEnergy = float(energies_final.iloc[i,2])
            
            defectName = energies_final.iloc[i,0]
            if "_" not in defectName:
                raise ValueError(f"Invalid defect name format: {defectName}")
            
            j = 0
            firstElement = "" #stores the name of the first/added "element" in defect
            secondElement = "" #stores the name of the second/removed "second" in defect
            
            # Gets the name of the first/added element
            while(defectName[j] != "_"):
                firstElement = firstElement + defectName[j]
                j = j + 1
                
            j = j + 1
            
            # Gets the name of the second/removed element
            while(j != len(defectName)):
                secondElement = secondElement + defectName[j]
                j = j + 1
                
            if(oldElement == " "):
                oldElement = secondElement
            
            finalDefectEnergy = bulkDefectEnergy - bulkEnergy 
            
            # Calculates the defect energy at specific charge state
            for k in range(0, len(elementNames)):
                # Subtract Energy From "Added" Element
                if(firstElement == elementNames[k]):
                    finalDefectEnergy = finalDefectEnergy - elementEPA[k]
                
                # Add Energy From "Subtracted" Element
                if(secondElement == elementNames[k]):
                    finalDefectEnergy = finalDefectEnergy + elementEPA[k]
            q = int(energies_final.iloc[i, 1])
            V = float(energies_final.iloc[i, 4])
            correction = float(energies_final.iloc[i, 3])
            
            if(storedName != defectName and i != 1):
                #Temp Array Declarations all needed for the charge neutrality plot
                tempArray = [] #Temp for charge states
                tempChargeArray = [] #Temp for charge states
                forGraph = [] #Temp for defect energies, merges tempArray
                forCharge = [] #Temp for charge states, merges tempChargeArray
                
                qw = 0 #used as a temp to determine color for plot
                feColor = "" #used as a temp to determine color for plot
                while(storedName[qw] != "_"):
                    feColor = feColor + storedName[qw]
                    qw += 1
                    
                colorName.append(feColor)
                
                if(len(finalColorNames) == 0):
                    finalColorNames.append(feColor)
                
                for bb in range (0, len(finalColorNames)):
                    if(feColor == finalColorNames[bb]):
                        break
                    elif(bb == len(finalColorNames) - 1):
                        finalColorNames.append(feColor)
                
                del(feColor, qw)
                
                for n in range (0, len(graphValues)):
                    allValues.append(graphValues[n])
                    allCharges.append(minCharge[n])
                
                if count == 0:
                    raise ValueError("No charge states found for defect, check input data")
                
                # Appends minimum defect energy at each fermi energy
                for m in range (0, int(len(graphValues)/count)):
                    for n in range (0, count):
                        tempArray.append(graphValues[m + int(len(graphValues)/count)*n])
                        tempChargeArray.append(minCharge  [m + int(len(minCharge)/count)*n])
                    forGraph.append(min(tempArray))
                    
                            
                    newIndex = tempArray.index(min(tempArray))
                    
                    forCharge.append(tempChargeArray[newIndex])
                    defectSpots.append(oldElement)
                        
                    if (newIndex != oldIndex):
                        oldIndex = newIndex
                        if(m != 0):
                            print("Transition from", forCharge[m - 1], "to", forCharge[m], "at", round(fermiEnergies[m], 5), "eV")
                                                    
                    completeGraph.append(forGraph[m])
                    completeMinCharge.append(forCharge[m])
                    tempArray = []
                    tempChargeArray = []
                    
                #Plots the individual charge defect plots, contains all of the charge states for an individual defect spanned across entire band gap
                if(config["plotsingledefect"] == True): 
                    formattedTitle = format_label(str(storedName))
                    plt.figure(figsize=(10,6))
                    plt.title("Defect Plot of " + formattedTitle)
                    plt.xlabel("Fermi Energy (eV)")
                    plt.ylabel("Formation Energy (eV)")
                    plt.plot(fermiEnergies, forGraph)
                    plt.xlim(xlimmin, xlimmax)
                    plt.ylim(ylimmin, ylimmax)
                    saveLocation = save_folder + "/" + str(storedName) + ".png"
                    plt.savefig(saveLocation)
                    plt.show()
                
                namesArray.append(storedName)
                storedName = defectName
                oldElement = secondElement
                
                # Clear Graph Values
                graphValues = []
                minCharge = []
                count = 0
            
            # Account for Charge Defect and Correction Values
            finalDefectEnergy = finalDefectEnergy + q*(E_f + V) + correction
            energy = '{:<12}  {:>6}'.format(defectName + "_" + str(q), str(round(finalDefectEnergy,5)))
            print(energy)

            # Everything Below is for plotting defect energy vs fermi energy
            for k in range(0, int(iterations)):
                graphValues.append(finalDefectEnergy + q*stepSize*k) # Adds the total fermi energy multiplied by charge
                minCharge.append(q)
            
            count = count + 1

        # Erases temporary data for last graph
        tempArray = [] 
        forGraph = []
        forCharge = []
        tempChargeArray = []

        #Everything below repeats above for the last defect, until noted
        qw = 0
        feColor = ""
        while(storedName[qw] != "_"):
            feColor = feColor + storedName[qw]
            qw += 1
            
        colorName.append(feColor)
        
        if(len(finalColorNames) == 0):
            finalColorNames.append(feColor)
        
        for bb in range (0, len(finalColorNames)):
            if(feColor == finalColorNames[bb]):
                break
            elif(bb == len(finalColorNames) - 1):
                finalColorNames.append(feColor)
        
        del(feColor, qw)
        
        for n in range (0, len(graphValues)):
            allValues.append(graphValues[n])
            allCharges.append(minCharge[n])
        
        for m in range (0, int(len(graphValues)/count)):
            for n in range (0, count):
                tempArray.append(graphValues[m + int(len(graphValues)/count)*n])
                tempChargeArray.append(minCharge[m + int(len(minCharge)/count)*n])
            forGraph.append(min(tempArray))
            
            newIndex = tempArray.index(min(tempArray))
            
            forCharge.append(tempChargeArray[newIndex])
            defectSpots.append(oldElement)
            
            if (newIndex != oldIndex):
                oldIndex = newIndex
                if(m != 0):
                    print(f"Transition from {forCharge[m - 1]:>2} to {forCharge[m]:>2} at {fermiEnergies[m]:.5f} eV")

            completeGraph.append(forGraph[m])
            completeMinCharge.append(forCharge[m])
            tempArray = []
            tempChargeArray = []
        # This is the last of the repeated analysis
    
        # This plots the last individual defect
        plt.figure(figsize=(10,6))
        storedName = defectName
        namesArray.append(storedName)

        if(config["plotsingledefect"] == True): 
            formattedTitle = format_label(str(storedName))
            plt.title("Defect Plot of " + formattedTitle)
            plt.xlabel("Fermi Energy (eV)")
            plt.ylabel("Formation Energy (eV)")
            plt.plot(fermiEnergies, forGraph, label = str(q))
            plt.xlim(xlimmin, xlimmax)
            plt.ylim(ylimmin, ylimmax)
            saveLocation = save_folder + "/" + str(storedName) + ".png"
            plt.savefig(saveLocation)
            plt.show()

        numberOfDefects = int(len(completeGraph)/len(fermiEnergies))

        plt.figure(figsize=(5,7))
        # plt.title("Charge Defect Plot")
        plt.xlabel("Fermi Energy (eV)")
        plt.ylabel("Formation Energy (eV)")
        plt.xlim(xlimmin, xlimmax)
        plt.ylim(ylimmin, ylimmax)
        
        print("")

        #Below contains calculations to determine the charge neutrality of the system
        temp1 = [] #Contains the defect energy of each defect at each fermi energy
        temp2 = 0 #currect fermi energy being analyzed
        temp3 = [] #Contains the minimum charge states
        temp4 = [] #Contains the number of possible degeneracy states for the type of defect, obtained from original POSCAR
        sign1 = False #Temps used to determine if the charge of the whole system flips
        sign2 = False #Temps used to determine if the charge of the whole system flips
        Q = 0 #Temps used to determine if the charge of the whole system flips
        kT = config['kT'] #specfied boltzmann * temp value
        e = np.exp(1) #exp value
        qArray = [] #stores all Q values, and is used to see if charge of system changes
        qValue = 0 #fermi energy value stored to print when sign changes
        #Determines intrinsic fermi level of defects
        for i in range(0, int(iterations)):
            qArray = []
            temp2 = float("{:0.4f}".format(fermiEnergies[i]))
            for j in range(0, numberOfDefects):
                temp1.append("{:0.7f}".format(completeGraph[i + j*int(iterations)]))
                temp3.append(completeMinCharge[i + j*int(iterations)])
                temp4.append(defectSpots[i + j*int(iterations)])
                
                q_i = int(temp3[j])
                
                for k in range(0, len(element_names)):
                    # Subtract Energy From "Added" Element
                    if(temp4[j] == element_names[k]):
                        N_i = defect_sites[k]
                                                        
                Q = Q + N_i*q_i*(e**(-1 * float(temp1[j]) / (kT))) #Calculates total Q value at fermi energy
                
                if(float(config['testfe']) == float("{:0.4f}".format(fermiEnergies[i]))):
                    print("charge state of defect", j, "=", q_i)
                    print("degeneracy states of defect", j, "=", N_i)
                    print("formation energy of defect", j, "=", temp1[j])
                    print()
            
                qArray.append(Q)
                    
            if(Q > 0):
                sign1 = True
            else:
                sign1 = False
            
            if(i != 0 and sign2 != sign1):
                print()
                print("Intrinisc Fermi Defect Level: " + "{:0.4f}".format(temp2) + " eV")
                print()
                if(config['printQ']):
                    tempQ = 0 
                    for k in range(0, numberOfDefects):
                        tempQ = qArray[k] - tempQ 
                        print("Q value of defect", k, "=", tempQ)
                    print("")
                    print("Total Q Value =", Q)
                    print("")
                    del(tempQ)
                qValue = temp2
                        
            if(config["testfe"] == float("{:0.4f}".format(fermiEnergies[i]))):
                print("Q value =", Q)
                print()
            
            sign2 = sign1
            
            temp1 = []
            temp3 = []
            temp4 = []
            Q = 0
            
        if(config["hse"] != None):
            plt.fill([xlimmin, xlimmin, originalVBM - E_f, originalVBM - E_f], [ylimmin, ylimmax, ylimmax, ylimmin], color = "silver")
            plt.fill([xlimmax, xlimmax, originalVBM - E_f + originalGap, originalVBM - E_f + originalGap], [ylimmin, ylimmax, ylimmax, ylimmin], color = "silver")

        # Format the labels in namesArray
        formatted_labels = [format_label(label) for label in namesArray]

        for i in range(0, numberOfDefects):
            tempData = []
            for j in range (0, int(len(fermiEnergies))):
                tempData.append(completeGraph[i*len(fermiEnergies) + j])
            
            for j in range(0, len(finalColorNames)):
                if(colorName[i] == finalColorNames[j]):
            
                    plt.plot(fermiEnergies, tempData, label=formatted_labels[i], color = colors[j], linestyle = lineStyles[lineStyleCount[j] % len(lineStyles)])
                    lineStyleCount[j] = lineStyleCount[j] + 1
        plt.axvline(qValue, color="black", linestyle="dashed")
        colNum = math.ceil(numberOfDefects/7)
        plt.legend(loc = config["legloc"])
        plot_name = config["save_as"]
        saveLocation = f"{save_folder}/{plot_name}.png"
        plt.savefig(saveLocation)
        plt.show()
        
        namesArray = []
        storedName = energies_final.iloc[1,0]
        oldElement = " "
        
        colorName = []
        
if __name__ == "__main__":
    main()
