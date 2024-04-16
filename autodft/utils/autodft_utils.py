from rdkit import Chem
from rdkit.Chem import Descriptors

import os
import json
import subprocess
import sys

from autodft.config.config import Config


def check_xtb() -> bool:
    """Checks whether the xTB package is installed"""

    try:
        subprocess.run(['xtb', '-h'], stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL)
    except FileNotFoundError:
        print('xTB is not installed. Do this with: conda install -c conda-forge xtb')
        print('    or: "micromamba install -c conda-forge xtb"')
        sys.exit()


def check_crest() -> bool:
    """Checks whether the CREST package is installed"""

    try:
        subprocess.run(['crest', '-h'], stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL)
        return True
    except FileNotFoundError:
        print('xTB is not installed. Do this with: "conda install -c conda-forge crest"')
        print('    or: "micromamba install -c conda-forge crest"')
        sys.exit()
        

def charge_from_smiles(smiles: str) -> int:
    """Returns the charge of a molecule based on its SMILES string"""

    mol = Chem.MolFromSmiles(smiles)
    return Chem.GetFormalCharge(mol)


def multiplicity_from_smiles(smiles: str) -> int:
    """Returns the spin multiplicity of a molecule based on its SMILES string"""

    mol = Chem.MolFromSmiles(smiles)
    return Descriptors.NumRadicalElectrons(mol) + 1


def input_data() -> dict:
    """Returns a dictionary of parameters specified in a
    .json configuration file"""

    json_files = [f for f in os.listdir() if f.endswith('.json')]
    if len(json_files) == 1:
        json_file = json_files[0]
        with open(json_file) as f:
            return json.load(f)
    else:
        raise Exception(
            'Invalid .json file(s) for autodft_flow. There must be one file named autodft_<molname>.json')


def write_autodft_parameters(log_file: str, inputs: dict, config: Config) -> None:
    """To a logging file, writes the autodft parameters that will be used,
    based on both user input and the supplied configuration information."""

    with open(log_file, 'a') as f:

        f.write(f'{" AutoDFT PARAMETERS ":*^75}\n')

        width = 40
        tmp = f'{{:<{width}}}{{}}\n'

        # Molecule info
        f.write('\nMOLECULE SPECIFICATIONS\n')
        f.write(tmp.format("o Name:", inputs['molname']))
        if 'smiles' in inputs:
            f.write(tmp.format("o SMILES string:", inputs['smiles']))
        if 'xyz_file' in inputs:
            f.write(tmp.format("o XYZ file:", inputs['xyz_file']))
        f.write(tmp.format("o Charge:", inputs['charge']))
        f.write(tmp.format("o Spin multiplicity:", inputs['multiplicity']))

        # CREST info
        f.write('\nCREST JOB\n')
        f.write(tmp.format("o Energy cutoff (kcal/mol):", inputs['e_cutoff']))
        f.write(tmp.format("o CREST method:", config.crest['method']))
        f.write(tmp.format("o CREST solvent:", config.crest['solvent']))
        f.write(tmp.format("o CREST keywords:", config.crest['keywords']))

        # Gaussian info
        f.write('\nGAUSSIAN JOBS\n')
        for i, job in enumerate(inputs['gaussian_jobs'], start=1):
            l_o_t = f'{job["functional"]}/{job["basisset"]} {job["route"]}'
            total = len(inputs['gaussian_jobs'])
            f.write(tmp.format(f'o Job {i} of {total}:', l_o_t))
            if job['basisset'].lower() == 'genecp':
                f.write(
                    '    ' + tmp.format('Basis set (light elements):', job['ecp_basisset_light']))
                f.write(
                    '    ' + tmp.format('Basis set (heavy elements):', job['ecp_basisset_heavy']))
                f.write(
                    '    ' + tmp.format('Max atomic # for light elements:', job['ecp_cutoff']))
        f.write(f'\n{" AutoDFT PROGRESS ":*^75}\n\n')


def cleanup(molname: str) -> None:
    """Delete .config and .sh files created by autodft to run autodft_flow"""

    try:
        os.remove(f'autodft_{molname}.json')
        os.remove(f'autodft_{molname}.sh')
    except FileNotFoundError:
        pass  # Possibly already deleted
