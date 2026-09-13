
# Tutorial of Behrend Defect Analyzer (BDA)

This page explains how to use the `BDA` code.

To follow this tutorial, defect calculations should already have been performed using [PyDefect](https://github.com/kumagai-group/pydefect) or another method [1]. This workflow also builds on scripts developed by [Zachery Willard](https://github.com/zacherywillard) [2] and [Evan Payne](https://github.com/EvanPayne22) [3], which handle defect energy extraction and formatting for BDA. If further information on formatting is required, refer to their repositories. Once these calculations are complete, charged defect formation energy analysis can begin.

The BDA is a tool used to calculate formation energies using the quantum simulation package known as [VASP](https://www.vasp.at/) [4]. It also uses [`sxdefectalign`](https://sxrepo.mpie.de/attachments/download/73/sxdefectalign-manual.pdf) to compute correction terms for charged defects [5], including finite-size and long-range potential energy corrections. The tool can then plot formation energy as a function of fermi energy.  These plots help identify the stability of defects, defect charge, and the effects of doping (i.e., changes in fermi energy). The BDA helps to automate this process. 


The formation energy of a defect is calculated as:

$$
E_\text{form}= E_\text{tot}^\text{defect} - E_\text{tot}^\text{bulk} - \sum_i \Delta n_i \mu_i^{eff} + q(E_\text{vbm}+E_F+\Delta V) + E_\text{corr}
$$

Where:

- $E_\text{tot}^\text{defect}$ is the total energy of the defect supercell  
- $E_\text{tot}^\text{bulk}$ is the total energy of the pristine bulk  
- $\Delta n_i$ is the number of atoms of species $i$ added or removed  
- $\mu_i^{eff}$ is the chemical potential of species $i$, defined as: $\mu_i^{eff}=\mu_i^{bulk}+\Delta \mu_i$
	- $\mu_i^{\mathrm{bulk}}$ is the bulk reservoir energy per atom of species $i$
	- $\Delta \mu_i$ is the relative chemical potential of species $i$ (e.g., Ga-rich)
- $q$ is the defect charge  
- $E_F$ is the Fermi level  
- $E_\text{vbm}$ is the valence band maximum  
- $\Delta V$ is the potential alignment correction calculated by `sxdefectalign`  
- $E_\text{corr}$ is the finite-size and electrostatic correction obtained from `sxdefectalign`

This workflow assumes users know all values except the energy correction terms. For more details on formation energies and correction schemes in GaN, see Lyons and Van de Walle [6].

The BDA assumes the following directory structure:
- The placeholder `<project_name>` typically represents the name of the target material.

```
	<project_name>/
	│
	├── bulk/
	│   ├── POSCAR
	│   ├── OUTCAR
	│   └── LOCPOT
	│
	└── defects/
	    ├── energies_final_vAtoms_plots.py
	    ├── formation_vs_fermi.py
	    ├── make_vAtoms_output.sh
	    ├── run_sxdefectalign.sh
	    ├── target_vertices_X_Rich.yaml
	    ├── target_vertices_Y_Rich.yaml
	    │ 
	    ├── <defect>_0/
	    │   └── POSCAR
	    ├── <defect>_1/
	    │   └── ...
	    ├── <defect>_-1/
	    │   └── ...
	    └── ...
```
We recommend that users follow the same directory structure. The defect directories must also follow the naming convention `<defect>_charge`!

The BDA also assumes the following formatting for POSCAR files. The header must include the defect center (i.e. three values representing the x,y,z coordinates of the defect in lattice units). The rest of the POSCAR should have the standard format, see for example: [Materials Project](https://next-gen.materialsproject.org/materials) [7] and below:

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
   
The plotting workflow requires a chemical potential input file (`target_vertices_<element>_Rich.yaml`) that defines the growth conditions used in the defect formation energy calculations (e.g., Ga-rich or N-rich). A detailed description of the file format and naming convention is provided in Step 2.

The rest of the details for the workflow are explained step by step, using GaN as an example, calculated with the PBE and HSE functional.
# Step 1. Energy Corrections ($E_{corr}$, $\Delta V$)

In charged-defect VASP calculations, the defect artificially interacts with its periodic images due to the finite size of the simulation. A correction is required to gain the isolated defect formation energy. Christopher Freysoldt created a program to compute this correction ($E_{corr}$). It also generates data to determine the shift in the long-range electrostatic potential ($\Delta V$). The `run_sxdefectalign` script uses Freysoldt's program called `sxdefectalign`, which must first be installed. 

**Download and setup**
Download at https://sxrepo.mpie.de/projects/sphinx-add-ons/files, then download sxdefectalign.bz2 and install with the following:
``` {.bash language="bash"}
bunzip2 sxdefectalign.bz2
chmod +x sxdefectalign
mv sxdefectalign ~/work/bin/
```
The `sxdefectalign` program is expected to be located in a directory named bin inside of your work directory. If the bin folder does not exist, create it and move the Freysoldt `sxdefectalign` inside. The path is hard coded in the `run_sxdefectalign` script and must be updated if Freysoldt's `sxdefectalign` is installed elsewhere. The `make_vAtoms_output` script is also required. It automatically searches all subdirectories in the working directory, reads the `vAtoms.dat` files created by `run_sxdefectalign`, formats the vAtom data as comma-seperated values, and then stores them in a file named `vAtoms_output.csv`. 


**Configuring the Scripts**
To use `run_sxdefectalign`, edits must first be made as described below.

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

First make the `run_sxdefectalign` executable in the terminal and run:

``` {.bash language="bash"}
chmod +x run_sxdefectalign.sh
./run_sxdefectalign.sh
```

During execution, the script prints the defect charge state, defect directory name, and the corresponding correction energy. An example output is shown below:
```
-0
Va_Ga_0/, -1072.37408904, 0
1
Va_Ga_-1/, -1069.98956324, 0.164549
-1
Va_Ga_1/, -1074.68454378, 0.164549
2
Va_Ga_-2/, -1067.30102614, 0.658195
...
```


Then make `make_vAtoms_output` executable and run:
```
chmod +x make_vAtoms_output.sh
./make_vAtoms_output.sh
```
 The script prints nothing by default but will populate vAtoms_output.csv.
  
## Outputs
`run_sxdefectalign` produces `energies_correction.csv` which contains: 
  
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
**Note:**

- Column 1: Distance from defect (Å), radial distance of each atom from the defect center. Used as the x-axis in ΔV plots.

- Column 2: Raw potential difference ($ΔV_{raw} = V_{\text{defect}} - V_{\text{bulk}}$).

- Column 3: Model long-range potential.  

- Column 4: Corrected potential ($\Delta V_{\text{aligned}} = \Delta V_{\text{raw}} - V_{\text{model}}$).

- Column 5: Weighting term, internal value from `sxdefectalign` (e.g., weighting or screening-related information).
## Energies Final and $\Delta V$ Plots

Once `energies_correction.csv` and `vAtoms_output.csv` are ready, use `energies_final_vAtoms_plots.py` to compute the potential alignment corrections ($\Delta V$) for each defect, create the `energies_final.csv` file, and create a `YAML` file containing relevant information on each neutral defect inside the `defects` directory. The `YAML` file is created by recursively parsing the defect directory for all `<defect>_0/POSCAR` files and storing information about each neutral defect, including its POSCAR path, elements, and atom counts. The script then determines all unique elements present across these neutral defect POSCARs. The chemical potentials must be provided for every unique element detected across the neutral defect POSCARs. Each chemical potential must be specified using the `Element=Energy` format. The order of the entries does not matter. For example, if the detected elements are `Ga` and `N`, both `-mu Ga=-2.91250895 N=-8.31707533` and `-mu N=-8.31707533 Ga=-2.91250895` are valid.
### Program Arguments  
- `-defectdirectory`: Directory containing neutral defects (default: `.`)
- `-defectyaml`: Output YAML for defect POSCAR info (default: `defect_poscars.yaml`)
- `-vatoms`: Path to `vAtoms_output.csv` (default: `./vAtoms_output.csv`)  
- `-correction`: Path to `energies_correction.csv` (default: `./energies_correction.csv`)  
- `-percent`: Fraction of the furthest atoms used to compute ΔV (default: 0.8)  
- `-number`: Number of furthest atoms used for ΔV (default: -1)
- `-mu`: One or more named chemical potentials in the form `Element=Energy`
- `-plotvatoms`: Boolean flag to generate vAtoms plots for all defects (default: `True`)
- `-vatomsxmin`: Minimum x-axis for vAtoms plots (default: -100, auto-scaled)
- `-vatomsxmax`: Maximum x-axis for vAtoms plots (default: -100, auto-scaled)
- `-vatomsymin`: Minimum y-axis for vAtoms plots (default: -100, auto-scaled)
- `-vatomsymax`: Maximum y-axis for vAtoms plots (default: -100, auto-scaled)
 
**Note:** You can use either `-percent` or `-number` to select atoms for $\Delta V$ calculation. If both are provided, `-number` takes precedence.  
  
### Example Usage  
We recommend that users create a small bash script to run the program. We will call it `run_energies_final_vAtoms_plots.sh`. This makes updating and keeping track of arguments easier. Here is an example:
```
#run energies_final_vAtoms_plots.py
python energies_final_vAtoms_plots.py -mu Ga=-2.91250895 N=-8.31707533 -percent 0.85
```
### Example Output

`energies_final_vAtoms_plots.py` produces `energies_final.csv` with the following format:
```
Defect Name, Charge, Bulk Energy, Correction Energy, Delta V, Std Deviation 
bulk,0,-779.26382452,0.0,0.0,0.0  
Va_Ga,0,-769.12439871,0.0,-0.1134629125,0.008629802275897968  
Va_Ga,-1,-766.05547513,0.1844,-0.10367685625,0.014714400975998342  
Va_Ga,1,-771.5878355,0.1844,-0.182173,0.012320662477115425  
Va_Ga,-2,-762.69860429,0.737601,-0.09537269999999999,0.02371790751483992  
Va_Ga,-3,-759.18076603,1.6596,-0.08617888124999998,0.03630872044557992
```
**Manual Creation of `energies_final.csv`**
This file can be created manually if these calculations have been done using another tool. If this is done, ensure that the order of the defects matches the order in `vAtoms.csv`. If using PyDefect, the information can be found in the following files:
- **Defect Name**: Name of the defect directory (e.g., `Va_Ga_0/`).  
- **Charge**: Encoded in the defect directory name (e.g., `Va_Ga_-1/` → `-1`).  
- **Bulk Energy**: Extract from `OUTCAR` of each calculation.  
- **Correction Energy ($E_\text{corr}$)**: Found in `defect_energy_info.yaml`.  
- **Potential Alignment (ΔV)**: Also from `defect_energy_info.yaml` (reported as alignment energy). Compute ΔV using:  $\Delta V =  \frac{E_\text{align}}{q}$ where $q$ is the defect charge.  
- **Standard Deviation of ΔV**: Not provided in PyDefect; can set to `0` if unknown.

`energies_final_vAtoms_plots.py` produces `defect_poscars.yaml` with the following format:
```
defects:
  Va_Ga:
    poscar: Va_Ga_0/POSCAR
    elements:
    - Ga
    - N
    counts:
    - 63
    - 64
  Va_N:
    poscar: Va_N_0/POSCAR
    elements:
    - Ga
    - N
    counts:
    - 64
    - 63
```
The program saves all ΔV plots in the `vAtomsImages` folder. These plots should be manually inspected to confirm their physical validity. An example is shown below:

<img src="images/vAtoms_for_Va_Ga_-3.png" alt="ΔV vs Radius for Va_Ga -3" width="600">

# Step 2. Plotting
Plotting the formation energy vs the fermi energy is a good way to qualitatively interpret which defects are most likely to be present when the material has a particular fermi energy. The `formation_vs_fermi.py` program will create these plots using files previously created in the tutorial. The program requires chemical potential input along with the `defect_poscars.yaml` file.
## Input Files

### Defect composition input

The program reads neutral defect compositions from `defect_poscars.yaml` by default. Each defect entry must provide `elements` and `counts` lists of equal length:
```
defects:
  Va_Ga:
    elements: [Ga, N]
    counts: [63, 64]
  Va_N:
    elements: [Ga, N]
    counts: [64, 63]
```
The element names and atom counts are compared with the bulk POSCAR to determine the composition change for each defect. The composition change is then used to calculate the chemical-potential contribution to the formation energy.
### Relative chemical potential input ($\Delta \mu$)
Relative chemical potentials are supplied through a YAML file. Each entry containing a `chem_pot` dictionary is treated as a chemical-potential phase or growth condition. For a single condition, use a file such as:

```
target: GaN
A:
  chem_pot:
    Ga: 0.0
    N: -1.31365
```

Multiple phase blocks may be included in one YAML file, provided that every block contains the same elements. The blocks are processed in file order. For systematic comparisons, separate files are still recommended, using names such as:
-   `target_vertices_Ga_Rich.yaml`
    
-   `target_vertices_N_Rich.yaml`

The `target` entry is retained as metadata; the program identifies chemical-potential conditions from the top-level blocks containing `chem_pot`.
### Bulk reservoir energies
The `-mu` values are the elemental bulk reservoir energies per atom. They must be supplied using the same element names as the `defect_poscars.yaml` file, in `Element=Energy` format. For example:
```
-mu Ga=-2.91250895 N=-8.31707533
```

The program validates that every element in the defect compositions has exactly one corresponding `-mu` value and rejects missing or extra elements, capitalization matters.
### Program Arguments  
-   `-plotsingledefect`: Generate individual plots for each defect (default: `False`)
-   `-defectyaml`: YAML containing defect compositions (default: `./defect_poscars.yaml`)
-  `-bulkposcar`: Path to the bulk POSCAR file (default: `../bulk/POSCAR`)  
-   `-correction`: Path to the final correction energies file (default: `./energies_final.csv`)
-   `-chempot`: Path to the chemical potential YAML file (default: `./target_vertices.yaml`)
-   `-ymax`: Maximum y-axis value for defect formation energy plot (default: `7`)
-   `-xmax`: Maximum x-axis value for defect plot (default: `-3`; replaced by the band gap)
-   `-ymin`: Minimum y-axis value for defect formation energy plot (default: `-7`)
-   `-xmin`: Minimum x-axis value for defect plot (default: `0`)
-   `-testfe`: Show defect charge-state information at a specified Fermi level (default: `-1`)
-   `-kT`: Thermal energy in eV used for occupation broadening (default: `0.05`)
-   `-printQ`: Print charge-state values at the intrinsic Fermi level (default: `False`)
-   `-colors`: List of colors for plotting (default: `["red", "blue", "orange", ...]`)
-   `-legloc`: Legend location identifier for plots (default: `8`)
-   `-hse`: `[band gap, VBM]` correction values for HSE calculations (default: `None`)
-   `--save_as`: Output filename prefix for generated plots (default: `combinedDefects`)
-   `-bg`: Band gap energy in eV (required)
-   `-vbm`: Valence band maximum offset in eV (required)
- `-mu`: One or more named chemical potentials in the form `Element=Energy`
## Example Usage  
We recommend that users create a small bash script to run the program. We will call it `run_formation_vs_fermi.sh`. This makes updating and keeping track of arguments easier. Here is an example:
```
#run formation_vs_fermi.py    Energy_per_atom Ga,N                                          						                    HSE_BG  HSE_VBM
python formation_vs_fermi.py -mu Ga=-2.91250895 N=-8.31707533 -bg 1.7378 -vbm 3.4099 -chempot target_vertices_Ga_Rich.yaml -ymin 0 -ymax 8 -hse 3.3212 2.3829 --save_as GaRich
python formation_vs_fermi.py -mu Ga=-2.91250895 N=-8.31707533 -bg 1.7378 -vbm 3.4099 -chempot target_vertices_N_Rich.yaml -ymin 0 -ymax 8 -hse 3.3212 2.3829 --save_as NRich
```
## Output

The program prints key defect information to the terminal, including:

- Number and names of elements in the defect system 
- Bulk and defect POSCAR compositions 
- Changes in elemental composition between the bulk and defect structures 
- Validation of defect energies in `energies_final.csv` 
- Effective chemical potentials (`μ_eff`) 
- Chemical potential terms and their contributions to the formation energy
- Defect formation energies at the VBM 
- Charge transition levels (Fermi energies) 
- Intrinsic Fermi defect level

An example output is shown below:

````
Chemical potential values: 2
energies_final.csv validation passed

Defects found in defect_poscars.yaml:
    Va_Ga
    Va_N

Phase: A
Defect: Va_Ga

Bulk POSCAR composition:
    Ga    64
    N     64

Defect POSCAR composition:
    Ga    63
    N     64

POSCAR Composition change:
    Ga    -1
    N     0

Effective chemical potentials:
    Ga    -2.912509 eV

Chemical potential terms:
    Ga: Delta_N = -1, mu = -2.912509, contribution = -2.912509 eV

Defect Formation Energies at VBM (1.999) in eV [phase: A]:
Va_Ga_0       7.670515 eV
Va_Ga_-1      9.069251 eV
Va_Ga_1       6.724370 eV
Va_Ga_-2     11.074738 eV
Va_Ga_-3     13.797562 eV

Transition from  1 to  0 at 0.94620 eV
Transition from  0 to -1 at 1.39880 eV
Transition from -1 to -2 at 2.00550 eV
Transition from -2 to -3 at 2.72290 eV

================================================================================

Phase: A
Defect: Va_N

Bulk POSCAR composition:
    Ga    64
    N     64

Defect POSCAR composition:
    Ga    64
    N     63

POSCAR Composition change:
    Ga    0
    N     -1

Effective chemical potentials:
    N     -9.630725 eV

Chemical potential terms:
    N: Delta_N = -1, mu = -9.630725, contribution = -9.630725 eV

Defect Formation Energies at VBM (1.999) in eV [phase: A]:
Va_N_0       1.573080 eV
Va_N_-1      4.582807 eV
Va_N_1      -1.173210 eV
Va_N_2      -1.654255 eV
Va_N_3      -1.907182 eV

Transition from  3 to  2 at 0.25300 eV
Transition from  2 to  1 at 0.48110 eV
Transition from  1 to  0 at 2.74630 eV
Transition from  0 to -1 at 3.00980 eV

Intrinsic Fermi Defect Level: 2.7463 eV

Formation energies saved to formation_energies.csv
Transition levels saved to transition_levels.csv
````

Two CSV files are created containing the calculated formation energies and transition levels.


`formation_energies.csv` contains the formation energy of every defect charge state at the selected VBM. Each row corresponds to one defect, charge state, and chemical potential phase.

For example:

```
Phase,Defect,Charge,Formation_Energy_eV
A,Va_Ga,0,7.670515
A,Va_Ga,-1,9.069251
A,Va_Ga,1,6.724370
...
```
 `transition_levels.csv` contains the Fermi energy at which the lowest-energy charge state changes from one charge state to another. Each row corresponds to one charge-state transition for a particular defect and chemical potential phase.

For example:
```
Phase,Defect,Charge_1,Charge_2,Transition_Fermi_Energy_eV
A,Va_Ga,1,0,0.94620
A,Va_Ga,0,-1,1.39880
A,Va_Ga,-1,-2,2.00550
...
```
The program will also store plots containing all defects in a directory named `chargeDefectPlots`. Here is an example of one such plot:
<img src="images/..." alt="" width="600">

If `-plotsingledefect True` is used, the program will also store individual defect plots in a `singleDefects` directory inside the `chargeDefectPlots` directory. Each individual plot is saved as a PNG file using the defect name. Examples of plots containing an individual defect and multiple defects are shown below:
<img src="images/..." alt="" width="600">

# References
[1] Yu Kumagai, Naoki Tsunoda, Akira Takahashi, and Fumiyasu Oba. Insights into oxygen vacancies from high-throughput first-principles calculations. *Phys. Rev. Materials*, 5:123803, 2021.  
[2] Zachery Willard. GitHub profile. https://github.com/zacherywillard, Accessed March 2026.  
[3] Evan Payne. GitHub profile. https://github.com/EvanPayne22, Accessed March 2026.  
[4] G. Kresse and J. Furthmüller. Efficient iterative schemes for ab initio total energy calculations using a plane-wave basis set. *Phys. Rev. B*, 54:11169–11186, 1996.  
[5] Christoph Freysoldt. Manual for sxdefectalign, version 3.0. Technical report, *MPI Fritz Haber Institute*, August 2022.  
[6] John L. Lyons and Chris G. Van de Walle. Computationally predicted energies and properties of defects in GaN. *npj Computational Materials*, 3:12, 2017.  
[7] Nubhav Jain, Shyue Ping Ong, Geoffroy Hautier, Wei Chen, William Davidson Richards, Stephen Dacek, Shreyas Cholia, Dan Gunter, David Skinner, Gerbrand Ceder, and Kristin A. Persson. The Materials Project: A materials genome approach to accelerating materials innovation. *APL Materials*, 1(1):011002, 2013
