# AutoDFT
A program for automated conformer searching and DFT calculations, designed
to be user-friendly for experimentalists with no computational background.

Conformer searching is performed using CREST (Pract, P.; Bohle, F.; Grimme, S.
Phys. Chem. Chem. Phys. 2020, 22, 7169). The lowest energy conformations within
a user-specified window are then optimized with DFT using Gaussian.
Scripts are provided to check the progress of Gaussian jobs and extract key
thermodynamic data with Goodvibes (Luchini, G.; Alegre-Requena, J. V.;
Funes-Ardoiz, I.; Paton, R. S. F1000 Research, 2020, 9, 291.).

The program is designed for use on a Linux cluster with a SLURM scheduler,
and has not been tested for use with Mac. The program is not compatible with
Windows.

## Installation (first time only)
Open a terminal window: Macs can use the built-in Terminal program, and PCs
can use the built-in Windows Terminal program. Connect to the computing cluster:
- ```ssh username@garibaldi.scripps.edu```

Install the micromamba package manager if not already installed:
- ```curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba```
- ```./bin/micromamba shell init -s bash -p ~/micromamba```
- ```source ~/.bashrc```

Create a micromamba environment and install packages (answer yes when prompted)
- ```micromamba create -n autodft```
- ```micromamba activate autodft```
- ```micromamba install -c conda-forge pandas goodvibes rdkit=2023.09.5 cclib xtb crest pyyaml```
- ```git clone git@github.com:shenvilab/autodft.git```
- ```cd autodft```
- ```pip install .```
- ```cd ..```

## Usage

### AutoDFT submission
This is performed using the ```run_autodft.py``` script, which takes molecules as 
SMILES strings. Alternatively, this can be performed using ```run_autodft_xyz.py```,
which takes an .xyz files of molecules as input, which is preferred for
transition metal complexes or organic structures for which automated generation
of a 3D structure fails (e.g. bicyclo[1.1.1]pentanes)

Connect to computing cluster and activate the autodft environment if not
done already:
```
micromamba activate autodft
```

Use of ```run_autodft.py``` is simply done as shown below, followed by answering
the provided prompts:
```
run_autodft.py
```

Use of ```run_autodft_xyz.py``` is analogously done with:
```
run_autodft_xyz.py
```

### Checking job status
The script ```gstatus.py``` is provided to help check the status of Gaussian
jobs submitted by AutoDFT or manually. Specify the name of all folders to check
for .log files, separated by spaces. Wildcards can be used as part of the folder
names. If no folders are specified, the script will look in the current
directory.

**Example:** See status of .log files in all folders (within the current one)
```
gstatus.py *
```
**Example:** See status of .log files in folders test1, test2, test3
```
gstatus.py test1 test2 test3
```
**Example:** See status of .log files in folders starting with test
```
gstatus.py test*
```
**Example:** See status of .log files in current directory
```
gstatus.py
```

### Compiling thermodynamic data
The script ```compile_results.py``` is provided to compile the key results for
all .log files in the specified directories. Specify the name of all folders to
check for .log files, separated by spaces. Wildcards can be used as part of the
folder names. If no folders are specified, the script will look in the current
directory.

**Example:** Compile results from .log files in all folders (within the current one)
```
gstatus.py *
```
**Example:** Compile results from .log files in folders test1, test2, test3
```
gstatus.py test1 test2 test3
```
**Example:** Compile results from .log files in folders starting with test
```
gstatus.py test*
```
**Example:** Compile results from all .log files in current directory
```
gstatus.py
```

Options:
- ```-a``` Use all conformations. (Default: False (use lowest only))
- ```-g``` Specify the name of the Goodvibes output file (Default: Goodvibes_output.csv)
- ```-o``` Specify the name of the summary file (Default: Goodvibes_output_summary.csv)

**Example:**
```
gstatus.py test* -a -g custom_filename.csv -o custom_summary_filename.csv
```

## Additional information

### Customization
If you would like, you can specify custom AutoDFT settings to use. To do this,
copy/paste the below text, modify as desired, and save it as a file named ```config.yaml```.

```
autodft_flow:         
  mem: 256mb                                  # Memory (for 'chaperone' job)
  processors: 1                               # Number of processors (for 'chaperone' job)
  time: '30:00:00'                            # Wall time (for 'chaperone' job). Should be at least the sum of
                                              # wall times for the CREST and Gaussian jobs

crest:
  mem: 12gb                                   # Memory
  method: gfn2                                # Method
  nprocshared: 16                             # Number of processors
  rm_extra_files: true                        # Remove extra CREST files
  script: run_crest.sh                        # Script name
  solvent: ''                                 # Solvent (do not include --gbsa flag)
  time: '6:00:00'                             # Wall time
  keywords: ''                                # Additional keywords


gaussian_input:
  version: '16'                               # Version of Gaussian (09 or 16)
  maxjobs: 10                                 # Maximum number of jobs to submit
  mem: 12gb                                   # Memory per job
  nprocshared: 12                             # Processors per job
  script: run_gaussian.sh                     # Script name
  time: '24:00:00'                            # Wall time per job
  write_chk: false                            # Whether to write and keep .chk files. For linked jobs, .chk files
                                              # will always be written, but this decides whether they are kept

gaussian_jobs:                                # Dashes indicate beginning of data for a job (can add or remove jobs as desired)
- basisset: GENECP                            # Basis set
  ecp_basisset_heavy: SDD                     # For split basis sets, the basis set for heavy atoms
  ecp_basisset_light: 6-31G(d)                # For split basis sets, the basis set for light atoms
  ecp_cutoff: 36                              # For split basis sets, the maximum atomic number for light atoms
  functional: B3LYP                           # Functional
  route: Opt Freq SCF=XQC                     # Route line (for linked jobs, do not add geom=allcheck guess=read)
- basisset: def2TZVP
  ecp_basisset_heavy: null
  ecp_basisset_light: null
  ecp_cutoff: null
  functional: M062X
  route: ''


goodvibes:
  conc: '1'                                   # Concentration (mol/L)
  f_cutoff: '100'                             # Frequency cutoff (cm-1)
  qs: truhlar                                 # Quasiharmonic oscillator approximation method
  keywords: ''                                # Additional keywords
```

Upload your ```config.yaml``` file to the computing cluster. When running
AutoDFT in the future, you can respond ```n``` to the prompt about using the
default configuration settings and specify the name of your ```config.yaml```
file (make sure you are in the folder containing that file).

### Getting SMILES strings
Draw the molecule in ChemDraw and highlight it. In the top menu:
**Edit > Copy as > SMILES**. Or use a keyboard shortcut
(**PC: Ctrl + Alt + C, Mac: Option + Command + C**).

#### Points of caution regarding SMILES strings

It's safer to draw the molecule in a 2-dimensional representation to help
ensure that stereochemistry is assigned properly. For fused rings, it's often
safest to draw the stereochemistry of hydrogens/substituents at the ring
junction, rather than to use endocyclic bonds.

Be especially careful about specifying stereochemistry when you have highly
bridged systems:

<p align="center">
![stereo_example1](images/stereo_example1.png)
</p>

Be careful about the directionality of wedges. Opposite directionality
corresponds to opposite stereochemistry:

<p align="center">
![stereo_example2](images/stereo_example2.png)
</p>

In certain cases, stereochemistry may be misintepreted as being at heteroatoms,
which can ultimately lead to invalid structures. Pay attention to any ChemDraw warnings:

<p align="center">
![stereo_example3](images/stereo_example3.png)
</p>