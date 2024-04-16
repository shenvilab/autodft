#!/usr/bin/python3

import os
import subprocess

from autodft.helpers.user_input import UserInputXYZ
from autodft.utils.autodft_utils import check_xtb, check_crest


def main() -> None:

    check_xtb()
    check_crest()
    
    inputs = UserInputXYZ()
    molnames = [m['molname'] for m in inputs.mols]
    for name in molnames:
        xyzfile = next(m['xyzfile'] for m in inputs.mols if m['molname'] == name)
        os.rename(xyzfile, f'{name}/{xyzfile}')
        os.chdir(name)
        subprocess.run(['sbatch', f'autodft_{name}.sh'])
        os.chdir('..')


if __name__ == '__main__':
    main()
