# Tutorial of Behrend Defect Analyzer (BDA)

This page explains how to use the `BDA` code.

To follow this tutorial, defect calculations should already have been performed using [PyDefect](https://github.com/kumagai-group/pydefect) or another method [1]. This workflow also builds on scripts developed by [Zachery Willard](https://github.com/zacherywillard) [2] and [Evan Payne](https://github.com/EvanPayne22) [3], which handle defect energy extraction and formatting for BDA. If further information on formatting is required, refer to their repositories. Once these calculations are complete, charged defect formation energy analysis can begin.

The BDA is a tool used to calculate formation energies using the quantum simulation package known as [VASP](https://www.vasp.at/) [4]. It also uses [`sxdefectalign`](https://sxrepo.mpie.de/attachments/download/73/sxdefectalign-manual.pdf) to compute correction terms for charged defects [5], including finite-size and long-range potential energy corrections. The tool can then plot formation energy as a function of fermi energy.  These plots help identify the stability of defects, defect charge, and the effects of doping (i.e., changes in fermi energy). The BDA helps to automate this process. 


The formation energy of a defect is calculated as:

$$
E_\text{form}= E_\text{tot}^\text{defect} - E_\text{tot}^\text{bulk} - \sum_i n_i \mu_i + q(E_\text{vbm}+E_F+\Delta V) + E_\text{corr}
$$

Where:

- $E_\text{tot}^\text{defect}$ is the total energy of the defect supercell  
- $E_\text{tot}^\text{bulk}$ is the total energy of the pristine bulk  
- $n_i$ is the number of atoms added or removed  
- $\mu_i$ is the chemical potential of species $i$  
- $q$ is the defect charge  
- $E_F$ is the Fermi level  
- $E_\text{vbm}$ is the valence band maximum  
- $\Delta V$ is the potential alignment correction calculated by `sxdefectalign`  
- $E_\text{corr}$ is the finite-size and electrostatic correction obtained from `sxdefectalign`

For more details on formation energies and correction schemes in GaN, see Lyons and Van de Walle [6].

The BDA assumes the following directory structure:
- The placeholder `<project_name>` typically represents the name of the target material.

```
    <project_name>
     │
     ├ POSCAR #Bulk
     │
     ├ bulk_supercell/ ──
     │                 ├─ OUTCAR
     │                 ├─ LOCPOT
     │
     └ defects/ ── 
                ├─ run_complete.sh
                ├─ complete_def_en_ef_min.py
				├─ energies_final.py
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
We recommend that users follow the same directory structure if possible. The defect directories must also follow the naming convention shown. 

The BDA also assumes the following formatting for POSCAR files. The header must include the defect center. ie three real values representing the x,y,z coordinates of the defect in lattice units; the rest of the POSCAR should have the standard format see for exampl: [Materials Project](https://next-gen.materialsproject.org/materials) [7] and below:

    ```
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
   - Note: when using PyDefect, the defect center does not appear in the first line of the POSCAR file. To locate it, refer to the defect_entry.json file, where it is specified as \"defect_center\".

The details of the workflow are explained step by step, using an example of GaN calculated with the PBE and HSE functional.

# Step 1. Energy Corrections ($E_{corr}$, $\Delta V$)

In charged-defect VASP calculations, the defect artificially interacts with its periodic images due to the finite size of the simulation. A correction is required to gain the isolated defect formation energy. Christopher Freysoldt created a program to compute this correction ($E_{corr}$). It also generates data to determine the shift in the long-range electrostatic potential ($\Delta V$). The `run_sxdefectalign_code` script uses Freysoldt's program called `sxdefectalign`, which must first be installed. 

**Download and setup**
Download at https://sxrepo.mpie.de/projects/sphinx-add-ons/files, then download sxdefectalign.bz2 and install with the following:
``` {.bash language="bash"}
bunzip2 sxdefectalign.bz2
chmod +x sxdefectalign
mv sxdefectalign ~/work/bin/
```
The sxdefectalign program is expected to be located in a directory named bin inside of your work directory. If the bin folder does not exist, create it and move the Freysoldt `sxdefectalign` inside. The path is hard coded in the `run_sxdefectalign_code` script and must be updated if Freysoldt's `sxdefectalign` is installed elsewhere. The `make_vAtoms_output` script is also required. It automatically searches all subdirectories in the working directory, reads the `vAtoms.dat` files created by `run_sxdefectalign_code`, formats the vAtom data as comma-seperated values, and then stores them in a file named `vAtoms_output.csv`. 


**Configuring the Scripts**
To use the `run_sxdefectalign_script` edits must first be made as described below.

-   Set the path to the bulk directory in line 13:

    ``` {.bash language="bash"}
    bulk="/path/to/BulkSupercell/"
    ```

-   Provide the dielectric tensor for your material in line 38:

    ``` {.bash language="bash"}
    --tensor 10.24,10.24,11.33
    ```

    Note: this should be the total dielectric tensor, both the electronic and ionic tensors added together. It follows Materials Project's convention. Use the Materials Project's values or perform the calculation yourself.

The `make_vAtoms_output` script will not require any modification. 

## Running the Scripts

First make the `run_sxdefectalign_code` executable in the terminal and run:

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
	
    6.  Reads the defect total energy from the defect OUTCAR.

    7.  Appends the defect name, total energy, and correction energy to energies_correction.csv.


Then make `make_vAtoms_output` executable and run:
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
   

## Outputs
`run_sxdefectalign_code` produces `energies_correction.csv` which contains: 
  
- Header row: Defect Name, Bulk Energy, Correction Energy ($E_{Corr}$)  

-  Data rows: directory name, bulk energy, correction energy for each defect

**Example snippet:**
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
`make_vAtoms_output` produces `vAtoms_output.csv`, which contains:

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
```
## Energies Final

Once `energies_correction.csv` and `vAtoms_output.csv` are ready, use `energies_final.py` to combine these files with the reservoir energies and compute the potential alignment corrections (ΔV) for each defect. This script also calculates the standard deviation of ΔV based on the chosen set of atoms.  
  
### Script Arguments  
  
- `-poscar`: Path to the POSCAR file (default: `./POSCAR`)  
- `-vatoms`: Path to `vAtoms_output.csv` (default: `./vAtoms_output.csv`)  
- `-correction`: Path to `energies_correction.csv` (default: `./energies_correction.csv`)  
- `-percent`: Fraction of the furthest atoms used to compute ΔV (default: 0.8)  
- `-number`: Number of furthest atoms used for ΔV (default: -1; ignored if `-percent` is set)  
- `resen`: Per-atom bulk energies for each element in POSCAR order
  
**Note:** You can use either `-percent` or `-number` to select atoms for ΔV calculation. If both are provided, `-number` takes precedence.  
  
### Example Usage  
We recommend that users create a bash script to run the `energies_final.py` program. We will call it `run_energies_final.sh`. This makes keeping track of arguments easier. Here is an example:
```
#run energies_final.py	 Energy_per_atom Ga,N
python energies_final.py -2.91250895 -8.31707533 -percent 0.85 -poscar ./POSCAR -vatoms ./vAtoms_output.csv -correction ./energies_correction.csv
```
Take notice that the bulk energies per atom must come first, then any additional arguments. 
### Example Output

`energies_final.py` produces `energies_final.csv` with the following format:
```
Defect Name,Charge,Bulk Energy,Correction Energy,delta V,Std Deviation  
bulk,0.0,-779.26382452,0.0,0.0,0.0  
Va_Ga,0.0,-769.12439871,0.0,-0.1134629125,0.008629802275897968  
Va_Ga,-1.0,-766.05547513,0.1844,-0.10367685625,0.014714400975998342  
Va_Ga,1.0,-771.5878355,0.1844,-0.182173,0.012320662477115425  
Va_Ga,-2.0,-762.69860429,0.737601,-0.09537269999999999,0.02371790751483992  
Va_Ga,-3.0,-759.18076603,1.6596,-0.08617888124999998,0.03630872044557992
```
**Manual Creation**
This file can be created manually if these calculations have been done using another tool. If using PyDefect, the information can be found in the following files:
- **Defect Name**: Name of the defect directory (e.g., `Va_Ga_0/`).  
- **Charge**: Encoded in the defect directory name (e.g., `Va_Ga_-1/` → `-1`).  
- **Bulk Energy**: Extract from `OUTCAR` of the bulk calculation.  
- **Correction Energy ($E_\text{corr}$)**: Found in `defect_energy_info.yaml`.  
- **Potential Alignment (ΔV)**: Also from `defect_energy_info.yaml` (reported as alignment energy). Compute ΔV using:  $\Delta V =  \frac{E_\text{align}}{q}$ where $q$ is the defect charge.  
- **Standard Deviation of ΔV**: Not provided in PyDefect; can set to `0` if unknown.

# Step 2. Plotting


## References
[1] Yu Kumagai, Naoki Tsunoda, Akira Takahashi, and Fumiyasu Oba. Insights into oxygen vacancies from high-throughput first-principles calculations. *Phys. Rev. Materials*, 5:123803, 2021.  
[2] Zachery Willard. GitHub profile. https://github.com/zacherywillard, Accessed March 2026.  
[3] Evan Payne. GitHub profile. https://github.com/EvanPayne22, Accessed March 2026.  
[4] G. Kresse and J. Furthmüller. Efficient iterative schemes for ab initio total energy calculations using a plane-wave basis set. *Phys. Rev. B*, 54:11169–11186, 1996.  
[5] Christoph Freysoldt. Manual for sxdefectalign, version 3.0. Technical report, *MPI Fritz Haber Institute*, August 2022.  
[6] John L. Lyons and Chris G. Van de Walle. Computationally predicted energies and properties of defects in GaN. *npj Computational Materials*, 3:12, 2017.  
[7] Nubhav Jain, Shyue Ping Ong, Geoffroy Hautier, Wei Chen, William Davidson Richards, Stephen Dacek, Shreyas Cholia, Dan Gunter, David Skinner, Gerbrand Ceder, and Kristin A. Persson. The Materials Project: A materials genome approach to accelerating materials innovation. *APL Materials*, 1(1):011002, 2013.
