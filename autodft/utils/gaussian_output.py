import cclib
from rdkit import Chem
from rdkit.Chem import AllChem, rdDetermineBonds

import os
import logging
import itertools


logger = logging.getLogger(__name__)


# CONSTANTS
RMSD_THRESH_H_POLARONLY = 0.003
RMSD_THRESH_H_ALL = 0.015


def move_redundant(logfiles: list[str]) -> None:
    """Move Gaussian .log files for redundant structures into a new
    subdirectory."""

    subdir = 'duplicates'
    for logfile in find_redundant(logfiles):
        comfile = logfile.replace('.log', '.com')
        new_comfile = f'{subdir}/{comfile}'
        new_logfile = f'{subdir}/{logfile}'
        try:
            os.mkdir(subdir)
        except FileExistsError:
            pass
        try:
            os.rename(logfile, new_logfile)
            os.rename(comfile, new_comfile)
        except FileNotFoundError:
            pass


def find_redundant(logfiles: list[str]) -> list[str]:
    """From a list of .log files, find which ones are for redundant
    structures.

    Redundant structures are determined based on their RMS deviation,
    computed using structures without nonpolar hydrogens to speed
    up GetBestRMS"""

    optimized, polarH_only = optimized_mols(logfiles)
    files = optimized.keys()
    is_redundant = {f: False for f in files}
    
    for i, j in itertools.combinations(files, 2):
        if is_redundant[i]:
            continue
        try:
            rmsd = AllChem.GetBestRMS(optimized[i], optimized[j])
            thresh = RMSD_THRESH_H_POLARONLY if polarH_only else RMSD_THRESH_H_ALL
            if rmsd < thresh:
                is_redundant[j] = True
                logger.info(f'Duplicate structures found: {i} and {j}')
        except RuntimeError:
            # GetBestRMS fails when connectivity is different,
            # in which case structures are not redundant
            pass

    redundant_files = [f for f in is_redundant if is_redundant[f]]
    return redundant_files


def is_gaussian_logfile(filename):
    """Returns True if the file name is one of a Gaussian log file"""
    
    if not filename.lower().endswith('.log'):
        return False

    with open(filename) as f:
        return 'Entering Gaussian System' in f.readline()


def normal_termination(filename):
    """Returns True if the Gaussian job terminated normally"""

    if not is_gaussian_logfile:
        return False

    with open(filename) as f:
        return 'Normal termination' in f.readlines()[-1]


def error_termination(filename):
    """Returns True if the Gaussian job terminated with an error"""
    
    with open(filename) as f:
        return 'Error termination' in f.read()


def has_linked_jobs(filename):
    """Returns True if the .log file contains at least one linked job"""
    
    if not normal_termination(filename):
        raise InvalidGaussianLogFileError()
    
    with open(filename) as f:
        data = f.read()
    n_jobs = data.count('Entering Link 1')
    return n_jobs > 1


def optimized_mols(logfiles):
    """Read optimized structures from a list of Gaussian .log files.
    Returns a list of rdkit molecules for the optimized structures,
    and whether the molecule contains polar Hs only."""

    files = [f for f in logfiles if normal_termination(f)]
    opt_mols = {}

    for f in files:
        logdata = cclib.io.ccread(f)
        charge = logdata.charge
        multiplicity = logdata.mult
        xyz_writer = cclib.io.XYZWriter(logdata, lastgeom=True)
        xyz = xyz_writer.generate_repr()

        mol = Chem.MolFromXYZBlock(xyz)
        polarH_only = True
        try:
            if multiplicity == 2:
                # Workaround since DetermineBonds fails for radicals
                rdDetermineBonds.DetermineBonds(mol, charge=charge-1)
            else:
                rdDetermineBonds.DetermineBonds(mol, charge=charge)
        except ValueError:
            # Transition metal complexes fail with rdDetermineBonds
            opt_mols[f] = mol
            polarH_only = False
        else:
            opt_mols[f] = remove_nonpolar_Hs(mol)

    return opt_mols, polarH_only


def remove_atoms(mol, atoms_to_remove):
    """Removes the specified atom numbers from an rdkit molecule"""

    atoms_to_remove.sort(reverse=True)
    edit_mol = Chem.EditableMol(mol)
    for atom in atoms_to_remove:
        edit_mol.RemoveAtom(atom)
    return edit_mol.GetMol()


def remove_nonpolar_Hs(mol):
    """Removes all carbon-bound hydrogens from an rdkit molecule"""

    nonpolar_Hs = []
    for i, a in enumerate(mol.GetAtoms()):
        if a.GetSymbol() == 'H':
            bonded_to_C = any(
                [neighbor.GetSymbol() == 'C' for neighbor in a.GetNeighbors()])
            if bonded_to_C:
                nonpolar_Hs.append(i)
    return remove_atoms(mol, nonpolar_Hs)


class InvalidGaussianLogFileError(Exception):
    """Raised when a valid Gaussian .log file is required but
    is not provided"""
    pass
