import json
import os
import sys
from textwrap import dedent
from typing import Tuple

from rdkit import Chem

from autodft.config.config import Config
from autodft.utils.autodft_utils import charge_from_smiles, multiplicity_from_smiles


# CONSTANTS
COL_WIDTH1 = 55
COL_WIDTH2 = 15
DASHES = '-' * (COL_WIDTH1 + COL_WIDTH2)


class UserInputSMILES:
    """A class to read and store user inputs from the command line for
    an AutoDFT job using a SMILES string as input"""

    def __init__(self) -> None:

        # Info for all jobs
        print_settings_header()
        self.config = config_from_user()
        self.gaussian_jobs = gaussian_jobs(self.config)
        self.e_cutoff = e_cutoff_from_user()

        # Info for individual jobs
        self.print_mol_info_header()
        self.mols = mols_from_user()

        user_verification()
        self.export_json()
        self.write_sh_script()

    def print_mol_info_header(self) -> None:
        """Prints how to provide information about the molecule(s)."""

        print(dedent(f"""\
            {DASHES}
              Molecule information. Please provide the following:
            
                  Molecule name: The name that will become the base of all
                                 file names associated with this molecule.
                                 Example: mol1
                  SMILES string: Radical-bearing atoms require explicit
                                 hydrogens. Maximum 1 radical per molecule.
                               
              Enter a blank molecule name to stop adding structures.
            {DASHES}"""))

    def export_json(self) -> None:
        """Writes the user inputs to .json configuration files"""

        data = vars(self).copy()
        if data['config'] is not None:
            data['config'] = '../' + data['config']
        mols = data.pop('mols')
        for mol in mols:
            to_export = {}
            to_export.update(data)
            to_export.update(mol)
            write_json(to_export)

    def write_sh_script(self) -> None:
        """Writes shell scripts to run AutoDFT"""

        molnames = [m['molname'] for m in self.mols]
        for name in molnames:
            write_shell_script(molname=name, program='autodft_flow.py',
                               config_file=self.config)


class UserInputXYZ:
    """A class to read and store user inputs from the command line for
    an AutoDFT job using an XYZ file as input"""

    def __init__(self) -> None:

        # Info for all jobs
        print_settings_header()
        self.config = config_from_user()
        self.gaussian_jobs = gaussian_jobs(self.config)
        self.e_cutoff = e_cutoff_from_user()

        # Info for individual jobs
        self.print_mol_info_header()
        self.mols = xyzs_from_user()

        user_verification()
        self.export_json()
        self.write_sh_script()

    def print_mol_info_header(self) -> None:
        """Prints how to provide information about the molecule(s)."""

        print(dedent(f"""\
            {DASHES}
              Molecule information
              --------------------
                  Molecule name: The name that will become the base of all
                                 file names associated with this molecule.
                                 Example: mol1
                  XYZ file:      Include .xyz extension.
                  Charge:        Charge of molecule.
                  Multiplicity:  Spin multiplicity (number of unpaired
                                 electrons + 1).

              Enter a blank molecule name to stop adding structures.
            {DASHES}"""))

    def export_json(self) -> None:
        """Writes the user inputs to .json configuration files"""

        data = vars(self).copy()
        mols = data.pop('mols')
        for mol in mols:
            to_export = {}
            to_export.update(data)
            to_export.update(mol)
            write_json(to_export)

    def write_sh_script(self) -> None:
        """Writes shell scripts to run AutoDFT"""

        molnames = [m['molname'] for m in self.mols]
        for name in molnames:
            write_shell_script(molname=name, program='autodft_flow_xyz.py',
                               config_file=self.config)


# Helper functions

def print_settings_header() -> None:
    """Print a header for the information prompts"""

    template = f'{{:<{COL_WIDTH1}}}{{}}'
    print('\n' + DASHES)
    print(template.format('  Setting', 'Value'))
    print(DASHES)


def config_from_user() -> str:
    """Get configuration file to use for autodft job, or use default
    if not specified"""

    use_default = y_or_n('  o  Use default autodft settings? (y/n):')
    if use_default:
        return None

    msg = '  o  Name of .yaml configuration file:'
    file = prompt(msg)
    while (not file.endswith('.yaml')) or (file not in os.listdir()):
        if not file.endswith('.yaml'):
            print('     File name must end with .yaml:')
            file = prompt(msg)
        if file not in os.listdir():
            print('     File must be present in current directory:')
            file = prompt(msg)
    return file


def gaussian_jobs(config_file: str = None) -> list[dict]:
    """Returns a list of the Gaussian job(s) to use for constructing
    Gaussian input files"""

    config = Config.from_yaml(config_file)
    if use_default_l_o_t(config_file):
        return config.gaussian_jobs
    else:
        return alt_job_from_user()


def use_default_l_o_t(config_file: str = None) -> bool:
    """Get whether to use the default level of theory (default taken
    from .yaml config file if one has been specified)"""

    config = Config.from_yaml(config_file)
    default_jobs = config.gaussian_jobs
    total = len(default_jobs)

    print('  o  Default computational method:\n')
    for i, job in enumerate(default_jobs, start=1):
        print(f'       Job {i} of {total}:')
        print_l_o_t(job)
    print('')
    return y_or_n('     Use default method? (y/n):')


def print_l_o_t(job: dict[str]) -> str:
    """Takes a dictionary of information for a Gaussian job and 
    prints a string representation of the level of theory"""

    func, basis, route = job['functional'], job['basisset'], job['route']
    print(f'         {func}/{basis} {route}')

    if job['basisset'].lower() == 'genecp':
        light = job['ecp_basisset_light']
        heavy = job['ecp_basisset_heavy']
        cutoff = job['ecp_cutoff']
        print(f'         Basis set for atomic number <= {cutoff}: {light}')
        print(f'         Basis set for atomic number >  {cutoff}: {heavy}')


def individual_job_from_user() -> dict:
    """Get alternate specification of an individual Gaussian job from user
    (could be one of the linked jobs)"""

    job = {}
    tmp = f'      {{:<{COL_WIDTH1-6}}}'  # String template

    job['functional'] = input(tmp.format('o  Functional (Example: B3LYP):'))
    job['basisset'] = input(tmp.format('o  Basis set (Example: 6-31G(d)):'))
    print(tmp.format('o  Keywords (leave blank if none).'))
    job['route'] = input(tmp.format(
        '   Example: opt freq scrf=(smd,solvent=water)'))

    if job['basisset'].lower() != 'genecp':
        job['ecp_basisset_light'] = None
        job['ecp_basisset_heavy'] = None
        job['ecp_cutoff'] = None
    else:
        print(tmp.format('GENECP specified. Provide ECP information.'))
        job['ecp_basisset_light'] = input(
            tmp.format('o  Basis set for light atoms:'))
        job['ecp_basisset_heavy'] = input(
            tmp.format('o  Basis set for heavy atoms:'))
        job['ecp_cutoff'] = int(
            input(tmp.format('o  Maximum atomic number for light atoms:')))

    return job


def alt_job_from_user() -> dict:
    """Gets specifications for a Gaussian job specified by the user
    (which includes linked jobs, if the user desires)"""

    print('\n  Provide the information for the first job.')
    jobs = [individual_job_from_user()]
    while True:
        print('')
        msg = f'  {len(jobs)} job(s) provided. Add linked job? (y/n):'
        add_job = y_or_n(msg)
        if not add_job:
            break
        jobs.append(individual_job_from_user())

    return jobs


def mols_from_user() -> list[dict]:
    """Asks the user to specify molecule names and SMILES until a blank
    molecule name is entered. The molecule info is returned as a list
    of dictionaries (with each corresponding to one molecule)."""

    mols, molnames = [], []
    while True:
        print('  Molecule', len(mols) + 1)
        name = molname_from_user(not_allowed=molnames)
        if not name:
            break
        smiles, charge, multiplicity = mol_from_user_smiles()
        molnames.append(name)
        mols.append({'molname': name, 'smiles': smiles, 'charge': charge,
                     'multiplicity': multiplicity})

    if not mols:
        print('  No molecules have been provided. Exiting.')
        sys.exit(0)
    return mols


def molname_from_user(not_allowed: list[str]) -> str:
    """Get desired name to use as the root for all file names to create.
    Makes sure that the name does not match any of the specified names"""

    msg = f'  o  {"Name:":<14}'
    molname = input(msg)
    while os.path.isdir(molname) or molname in not_allowed:
        print('     This name has already been used.')
        molname = input(msg)
    return molname


def xyzs_from_user() -> list[dict]:
    """Asks the user to specify molecules (XYZ files with accompanying info)
    until a blank molecule name is entered. The molecule info is returned as a list
    of dictionaries (with each corresponding to one molecule)."""

    mols, molnames = [], []
    while True:
        print('  Molecule', len(mols) + 1)
        name = molname_from_user(not_allowed=molnames)
        if not name:
            break
        molnames.append(name)
        mols.append({'molname': name,
                     'xyzfile': xyzfile_from_user(),
                     'charge': charge_from_user(),
                     'multiplicity': multiplicity_from_user()})

    if not mols:
        print('  No molecules have been provided. Exiting.')
        sys.exit(0)
    return mols


def mol_from_user_smiles() -> Tuple[str, int, int]:
    """Get SMILES string, charge, multiplicity from user-provided SMILES"""

    msg = f'  o  {"SMILES:":<14}'
    smiles = input(msg)

    while not valid_smiles(smiles):
        if not valid_structure(smiles):
            print('     SMILES string does not correspond to a valid structure.')
        elif '*' in smiles:
            print('     SMILES string cannot contain wildcard (*)')
        else:
            print('     A maximum of 1 radical is allowed on the structure.')
            print('     Radical-bearing atoms require explicit hydrogens.')
        smiles = input(msg)

    smiles = Chem.CanonSmiles(smiles)
    charge = charge_from_smiles(smiles)
    multiplicity = multiplicity_from_smiles(smiles)
    return smiles, charge, multiplicity


def xyzfile_from_user() -> str:
    """Get the name of a .xyz file from the user"""

    msg = f'  o  {"XYZ file:":<14}'
    file = input(msg)
    while (not file.endswith('.xyz')) or (file not in os.listdir()):
        if not file.endswith('.xyz'):
            print('     File name must end with .xyz:')
            file = input(msg)
        if file not in os.listdir():
            print('     File must be present in current directory:')
            file = input(msg)
    return file


def e_cutoff_from_user() -> str:
    """Get energy cutoff in kcal/mol for conformer generation"""

    print('  o  Energy cutoff for conformers in kcal/mol')
    e_cutoff = prompt('       Recommended = 1 kcal/mol:')
    while not is_positive_float(e_cutoff):
        e_cutoff = prompt('     Energy cutoff must be a positive number:')
    return float(e_cutoff)


def charge_from_user() -> int:
    """Prompts the user for an integer input"""

    msg = f'  o  {"Charge:":<14}'
    num = input(msg)
    while not is_int(num):
        print('  Please input an integer: ')
        num = input(msg)
    return int(num)


def multiplicity_from_user() -> int:
    """Get spin multiplicity from user"""

    msg = f'  o  {"Multiplicity:":<14}'
    mult = input(msg)
    while not mult.isdigit() or int(mult) < 1:
        print('  Please input an integer 1 or greater: ')
        mult = input(msg)
    return int(mult)


def user_verification() -> None:
    """Ask user to verify provided inputs"""

    verified = y_or_n('  Confirm that the above information is correct (y/n):')
    if not verified:
        print('  Exiting. Please try again.\n')
        sys.exit(0)
    print(DASHES)


def valid_smiles(smiles: str) -> bool:
    """Returns whether a SMILES string is valid"""

    if not valid_structure(smiles):
        return False
    if '*' in smiles:
        return False
    if multiplicity_from_smiles(smiles) not in [1, 2]:
        return False
    return True


def valid_structure(smiles: str) -> bool:
    """Returns whether the SMILES string corresponds to a valid structure """

    return Chem.MolFromSmiles(smiles) is not None


def is_positive_float(string: str) -> bool:
    """Returns true if a string represents a positive float"""

    try:
        return float(string) > 0
    except ValueError:
        return False


def is_int(string: str) -> bool:
    """Returns true if a string represents an integer"""

    try:
        int(string)
        return True
    except ValueError:
        return False


def y_or_n(msg: str) -> bool:
    """Prompts a user for a yes or no answer"""

    response = prompt(msg).lower()
    while response not in ['y', 'n']:
        response = prompt('  Please input y or n: ').lower()
    return response == 'y'


def prompt(msg: str) -> str:
    """Prompts the user using the provided message, to which string
    formatting is applied"""

    return input(f'{msg:<{COL_WIDTH1}}')


def write_json(config_data: dict) -> None:
    """Write .json configuration file for autodft_flow.py"""

    molname = config_data['molname']
    os.mkdir(molname)
    json_path = f'{molname}/autodft_{molname}.json'
    with open(json_path, 'w') as f:
        json.dump(config_data, f, indent=4)


def write_shell_script(molname: str, program: str,
                       config_file: str = None) -> None:
    """Writes a shell script to execute autodft_flow, which will
    execute using the configuration data from the .json file.
    The program can be 'autodft_flow.py' or autodft_xyz_flow.py"""

    config = Config.from_yaml(config_file)
    sh_path = f'{molname}/autodft_{molname}.sh'
    
    with open(sh_path, 'w') as f:
        f.write(dedent(f"""\
            #!/bin/sh
            
            #SBATCH --nodes=1
            #SBATCH --ntasks=1
            #SBATCH --cpus-per-task={config.autodft_flow['processors']}
            #SBATCH --mem={config.autodft_flow['mem']}
            #SBATCH --time={config.autodft_flow['time']}
            #SBATCH --output=/dev/null
            #SBATCH --error=/dev/null
            #SBATCH --job-name={molname[-8:]}
            
            
            {program}
            """))
