# AutoDFT
A program for automated conformer searching and DFT calculations.

Conformer searching is performed using CREST (Pract, P.; Bohle, F.; Grimme, S.
Phys. Chem. Chem. Phys. 2020, 22, 7169). The lowest energy conformations within
a user-specified window are then optimized with DFT using Gaussian.
Thermodynamic data is the extracted and processed with Goodvibes (Luchini, G.;
Alegre-Requena, J. V.; Funes-Ardoiz, I.; Paton, R. S. F1000 Research, 2020, 9,
291.).

The program is designed for use on a Linux cluster with a SLURM scheduler,
and has not been tested for use with Mac. The program is not compatible with
Windows.

## Installation (first time only)
Connect to the computing cluster as usual
- ```ssh username@garibaldi.scripps.edu```

Install the micromamba package manager if not already installed:
- ```curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba```
- ```./bin/micromamba shell init -s bash -p ~/micromamba```
- ```source ~/.bashrc```

Create a micromamba environment and install packages
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