#!/usr/bin/python3

import os
import subprocess

from autodft.helpers.user_input import UserInputSMILES
from autodft.utils.autodft_utils import check_xtb, check_crest


def main() -> None:

    check_xtb()
    check_crest()

    inputs = UserInputSMILES()
    molnames = [m['molname'] for m in inputs.mols]
    for name in molnames:
        os.chdir(name)
        subprocess.run(['sbatch', f'autodft_{name}.sh'])
        os.chdir('..')


if __name__ == '__main__':
    main()
