## Tutorial of Behrend Defect Analyzer (BDA)

This page explains how to use the `BDA` code.

To follow this tutorial, it is important that defect calculations have already been performed using pydefect or another method. Once these calculations are complete, charged defect formation energy analysis can begin.

The BDA is a tool used to calculate formation energies using the quantum simulation package known as [VASP](https://www.vasp.at/) [1]. It also uses [`sxdefectalign`](https://sxrepo.mpie.de/attachments/download/73/sxdefectalign-manual.pdf) to compute correction schemes [2]. As an end goal, the tool can plot formation energy as a function of fermi energy.  These plots help identify the stability of defects, defect charge, and the effects of doping (i.e., changes in fermi energy). The BDA helps to automate this process, making it more accessible to new users. 

The BDA assumes the following directory structure:
- The placeholder `<project_name>` typically represents the name of the target material, optionally including its crystal structure.

```
    <project_name>
     │
     ├ bulk_supercell/ ── vasprun.xml
     │                 ├─ CONTCAR
     │                 ├─ LOCPOT
     │                  ...outputs
     │
     └ defects/ ── complete_def_en_ef_min.py
                ├─ run_complete.sh
                ├─ no_vatoms.py
                ├─ make_vAtoms_output.sh
                ├─ run_sxdefectalign_code.sh
                ├─ target_vertices_X_Rich.yaml
                ├─ target_vertices_Y_Rich.yaml
                ├─ Va_X_0/
                ├─ Va_X_1/
                ├─ Va_X_2/
                 ...
```
We recommend that users follow the same directory structure if possible.
The details of the workflow are explained step by step, using an example of GaN calculated with the PBE and HSE functional.

In BDA, there are five main scripts as shown in the directory structure:

`complete_def_en_ef_min.py`
`run_complete.sh`
`no_vatoms.py`
`make_vAtoms_output.sh`
`run_sxdefectalign_code.sh`



## 1. Finite-Size Correction ($E_{corr}$)

In charged-defect calculations, the defect artificially interacts with its periodic images due to the finite size of the simulation. These interactions distort the total energy, which then requires a correction. Christopher Freysoldt created a program to compute this correction from the periodic boundary conditions in charged-defect VASP calculations. It also accounts for the long-range electrostatic potentials assumed to follow a Gaussian distribution. The `run_sxdefectalign_code` script uses Freysolt's program called `sxdefectalign`, which must first be installed. 

**Download and setup**
Download at https://sxrepo.mpie.de/projects/sphinx-add-ons/files, then download sxdefectalign.bz2 and install with the following:
``` {.bash language="bash"}
bunzip2 sxdefectalign.bz2
chmod +x sxdefectalign
mv sxdefectalign ~/work/bin/
```
The sxdefectalign program is expected to be located in a directory named bin inside of your home work directory. If the bin folder does not exist, create it and move the Freysoldt `sxdefectalign` inside. The path is hard coded in the `run_sxdefectalign_code` script and must be updated if Freysold's `sxdefectalign` is installed elsewhere.

**Configuring the Script**
To use the `run_sxdefectalign_script` edits must first be made. The following describes what requires editing along with how to do so. 

-   Set the path to the bulk directory in line 13:

    ``` {.bash language="bash"}
    bulk="/path/to/BulkSupercell/"
    ```

-   Provide the dielectric tensor for your material in line 38:

    ``` {.bash language="bash"}
    --tensor 10.24,10.24,11.33
    ```

    Note: this should be the total dielectric tensor, both the electronic and ionic tensors added together. It follows Materials Project's convention.

-   Ensure POSCAR coordinates have defect centers in the first line:

    ``` {.bash language="bash"}
    0.083333 0.166666 0.499553 #defect center
    1.0
      12.8672158884    0.0000000000    0.0000000000
      -6.4336079442   11.1433358354    0.0000000000
       0.0000000000    0.0000000000   10.4781669866
    Ga N
    63 64
    direct
       0.0833333333    0.1666666667    0.9995534669 
       ...
    ```

    Note: when using PyDefect, the defect center does not appear in the first line of the POSCAR file. To locate it, refer to the defect_entry.json file, where it is specified as \"defect_center\".

-   Ensure directory names follow the convention of the script:

    ``` {.bash language="bash"}
    Va_Ga_-1    
    Va_Ga_-2    
    Va_Ga_-3    
    Va_Ga_0 
    Va_Ga_1
    ```
**Running the Script**

Make the script executable in the terminal and run:

``` {.bash language="bash"}
chmod +x run_sxdefectalign_code.sh
./run_sxdefectalign_code.sh
```

During execution, the script:

1.  Identifies all defect subdirectories in the working directory.

2.  Reads the bulk reference energy from the bulk OUTCAR file and writes it as the first entry in energies_correction.csv.

3.  For each defect directory that contains a LOCPOT file:

    1.  Determines the defect name and extracts the charge state from the directory name.

    2.  Locates the defect center from the first line of the corresponding POSCAR.

    3.  Converts the extracted charge into the sign convention required by sxdefectalign.

    4.  Executes sxdefectalign using the defect and bulk LOCPOT files and the user-specified dielectric tensor. This generates vAtoms.dat files for each directory, which will be used for the potential alignment correction. 

    5.  Extracts the correction energy from the sxdefectalign output.

	6. Extracts the correction energy from the sxdefectalign output.
	
    7.  Reads the defect total energy from the defect OUTCAR.

    8.  Appends the defect name, total energy, and correction energy to energies_correction.csv.

**Output**

Example of energies_correction.csv:

``` {.bash language="bash"}
Defect Name, Bulk Energy, Correction Energy
bulk, -779.26382452, 0
Va_Ga_0/, -769.12439871, 0
Va_Ga_-1/, -766.05547513, 0.1844
Va_Ga_1/, -771.58783550, 0.1844
Va_Ga_-2/, -762.69860429, 0.737601
Va_Ga_-3/, -759.18076603, 1.6596
...
```

## 2. Potential Alignment Correction ($\Delta V$)

This term is dependent on the long-range electrostatic potentials. 

The make_vAtoms_output script automatically searches all subdirectories in the working directory, reads the vAtoms.dat files if they exist, formats them as comma-separated values, and stores them in a file named vAtoms_output.csv. The vAtoms.dat files are created by running the `sxdefectalign` program. 


Each defect directory must contain a vAtoms.dat file generated by `run_sxdefectalign_code`. Directories without vAtoms.dat will be skipped automatically.
    



**Configuring the Script**

The script can be run without editing and does not require any modifications!


**Running the Script**

Make the script executable and run it in the terminal:
```
chmod +x make_vAtoms_output.sh  
./make_vAtoms_output.sh
```

-   The script prints nothing by default but will populate vAtoms_output.csv.
    
During execution, the script:

1.  Lists all subdirectories in the current working directory.
    
2.  Deletes any existing vAtoms_output.csv to prevent appending to old data.
    
3.  Writes a header row to the CSV.
    
4.  Loops through each defect directory:
    
    -   Checks for the presence of vAtoms.dat.
        
    -   Writes a directory marker into the CSV.
        
    -   Appends the contents of vAtoms.dat, formatting it as comma-separated values.
        
5.  Adds a final empty line and stop marker to the CSV for clarity.
   

**Output**

The final output is vAtoms_output.csv, which contains:

-   Header row: Column 1, Column 2, Column 3, Column 4, Column 5, ...

-   Directory markers (stop,<directory_name>) before each defect’s data.
    
-   Comma-separated vAtoms data for each defect.
    
-   Final stop marker at the end.
    

**Example snippet:**
```
Column 1, Column 2, Column 3, Column 4, Column 5  
stop,Va_Ga_0/  
9.87797,0,-0.0813866,-0.0813866,-0.00228509
6.12053,0,-0.0238853,-0.0238853,-3.09364
...
stop,Va_Ga_-1/  
9.87751,0.0402261,-0.0188545,-0.0590807,-0.00204068
6.08003,0.112015,0.0825346,-0.0294803,-3.05228
...  
stop


