import pandas as pd

import subprocess
import os
import logging
from io import StringIO
from dataclasses import dataclass, field

from autodft.utils.files import files_in_dir


logger = logging.getLogger(__name__)


# Define constants

AU_TO_KCAL = 627.509541
G_kcal_name = 'qh-G(T)(_SPC) (kcal/mol)'
Grel_kcal_name = 'qh-G(T)(_SPC)_rel (kcal/mol)'


@dataclass
class GV_Executor:
    """Runs goodvibes with the specified options"""

    files: list[str] = field(default_factory=lambda: ['*.log'])
    qs: str = 'truhlar'
    f_cutoff: str = '100'
    conc: str = '1'
    spc: str = ''
    csv: bool = False
    logging_file: str = None

    def run(self, gv_output_name: str = '') -> None:
        """Runs goodvibes with the specified options"""

        # Collect arguments
        gv_args = self.files + ['--qs', self.qs, '-f', self.f_cutoff,
                                '-c', self.conc, '--imag']
        if self.spc:
            gv_args += ['--spc', self.spc]
        if self.csv:
            gv_args += ['--csv']

        # Run Goodvibes
        output = subprocess.run(['python', '-m', 'goodvibes'] + gv_args,
                                stdout=subprocess.PIPE, text=True)
        if self.logging_file is not None:
            with open(self.logging_file, 'a') as f:
                f.write(output.stdout + '\n')
        else:
            print(output.stdout)

        # Rename goodvibes output file, if desired
        if gv_output_name:
            original = 'Goodvibes_output.csv' if self.csv else 'Goodvibes_output.dat'
            os.replace(original, gv_output_name)


class GV_Results:
    """
    Processes a Goodvibes output file and the .log files used to create it.
    The goodvibes file must be a .csv file for the parsing to work properly
    """

    def __init__(self, datafile: str) -> None:
        """Creates a GV_Results object from a goodvibes output file (.csv)"""

        if not datafile.lower().endswith('.csv'):
            raise ValueError('The goodvibes output file must be a .csv type')

        self.datafile = datafile
        self.parsed = {'intro': '', 'stars': '',
                       'thermo_lines': '', 'error_lines': ''}
        self.df = None
        self._fix_column_names()
        self._parse_csv()

    def _fix_column_names(self) -> None:
        """Fixes a typo in the column names of the goodvibes output file"""

        filedata = ''
        with open(self.datafile) as f:
            filedata = f.read()
        with open(self.datafile, 'w') as f_out:
            f_out.write(filedata.replace(',im,freq', ',im_freq'))

    def _parse_csv(self) -> None:
        """Parses thermodynamic data (.csv format) and read into a dataframe"""

        with open(self.datafile) as f:
            is_thermo_data = False
            n_elements = None
            for line in f:
                if '*****' in line:
                    self.parsed['stars'] = line
                    continue
                if 'Warning! Couldn\'t find frequency information ...' in line:
                    self.parsed['error_lines'] += line
                    continue
                if line.startswith('   Structure,'):
                    is_thermo_data = True
                    n_elements = line.count(',') # equal because Goodvibes adds an extra comma
                if is_thermo_data:
                    if line.count(',') == n_elements:
                        line = line.replace(',\n', '\n')
                    self.parsed['thermo_lines'] += line[3:] # Remove spaces and/or bullet point
                else:
                    self.parsed['intro'] += line

        self.df = pd.read_csv(StringIO(self.parsed['thermo_lines']))
        self.g_name = [x for x in list(self.df) if 'qh-G(T)' in x][-1]

    # def g_min_hartree(self):
    #     """Returns the free energy of the lowest energy structure(s)
    #     in hartrees"""

    #     return min(list(self.df[self.g_name]))

    # def add_g_kcalmol(self) -> None:
    #     """Adds G_rel in kcal/mol to the dataframe of thermodynamic data"""

    #     self.df[G_kcal_name] = self.df[self.g_name] * AU_TO_KCAL

    # def add_grel_kcalmol(self) -> None:
    #     """Adds G_rel in kcal/mol to the dataframe of thermodynamic data"""

    #     g_kcal = self.df[self.g_name] * AU_TO_KCAL
    #     g_kcal_min = min(g_kcal)
    #     self.df[Grel_kcal_name] = self.df[self.g_name] * \
    #         AU_TO_KCAL - g_kcal_min

    # def mark_lowest(self) -> None:
    #     """Annotate the lowest energy structure as such"""

    #     self.df['Lowest'] = self.df[self.g_name].map(lambda x:
    #                                                  'yes' if x == self.g_min_hartree() else '')

    # def rename_lowest(self) -> None:
    #     """Adds '(lowest)' to the file name of the lowest free energy structure"""

    #     lowest_gv_entries = self.df.loc[self.df[self.g_name]
    #                                     == self.g_min_hartree()]
    #     lowest_file_roots = [file.removeprefix('o  ')
    #                          for file in lowest_gv_entries['Structure']]
    #     lowest_file_names = [file.removesuffix('(lowest)') + '.log'
    #                          for file in lowest_file_roots] # Remove any preexisting suffixes to avoid adding twice

    #     for filename in lowest_file_names:
    #         try:
    #             os.rename(filename, filename.replace('.log', '(lowest).log'))
    #         except FileNotFoundError:
    #             pass  # The file may have already been renamed

    # def write_csv(self) -> None:
    #     """Overwrite the original goodvibes output file (.csv),
    #     incorporating any added information"""

    #     df_data = self.df.to_csv(index=False, lineterminator='\n')
    #     df_lines = df_data.splitlines(keepends=True)
    #     table_header = df_lines[0]
    #     table_data = ''.join(df_lines[1:])

    #     with open(self.datafile, 'w') as f:
    #         f.write(self.parsed['intro'])
    #         f.write(table_header)
    #         f.write(self.parsed['stars'])
    #         f.write(table_data)
    #         f.write(self.parsed['error_lines'])
    #         f.write(self.parsed['stars'])


# def goodvibes_analysis(molname: str = None,
#                        output_file: str = None,
#                        qs: str = 'truhlar',
#                        conc: str = '1',
#                        f_cutoff: str = '100',
#                        linked: bool = False,
#                        add_g_kcalmol: bool = False,
#                        add_grel_kcalmol: bool = False,
#                        mark_lowest: bool = False,
#                        rename_lowest: bool = False                       
#                        ) -> None:
#     """Run and process Gaussian log files with Goodvibes"""

#     # Set Goodvibes output file name
#     if molname is None and output_file is None:
#         output_file = 'Goodvibes_output.csv'
#     if output_file is None:
#         output_file = f'{molname}_goodvibes_data.csv'

#     # Print status
#     logger.info('Extracting and correcting thermodynamic data with Goodvibes...')
#     logging_file = f'{molname}.out'
#     with open(logging_file, 'a') as f:
#         f.write(f'\n{" GOODVIBES OUTPUT ":*^75}\n\n')
    
#     # Run Goodvibes and process results
#     logfiles = files_in_dir('.log')
#     gv_executor = GV_Executor(files=logfiles, qs=qs, f_cutoff=f_cutoff,
#                               conc=conc, spc='link' if linked else '',
#                               csv=True, logging_file = logging_file)
#     gv_executor.run(output_file)
#     gv_results = GV_Results(output_file)

#     if add_g_kcalmol:
#         gv_results.add_g_kcalmol()

#     if add_grel_kcalmol:
#         gv_results.add_grel_kcalmol()

#     if mark_lowest:
#         gv_results.mark_lowest()

#     if rename_lowest:
#         gv_results.rename_lowest()

#     gv_results.write_csv()
#     logger.info('Goodvibes analysis complete.')

