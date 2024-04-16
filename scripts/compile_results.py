#!/usr/bin/python3

import pandas as pd
import numpy as np

import argparse
import subprocess
import os
from dataclasses import dataclass
from io import StringIO

from autodft.config.config import Config
from autodft.utils.files import glob_dirs, log_files, \
    molname_from_log, print_indent
from autodft.utils.gaussian_output import has_linked_jobs, \
    is_gaussian_logfile, normal_termination, error_termination


# CONSTANTS
AU_TO_KCAL = 627.509541


def cmdline_parser() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(dest='dirs', nargs='*', action='store',
                        help='Folder(s) to run script on.')
    parser.add_argument('-a', '--all', dest='use_all', action='store_true', default=0,
                        help='Optional: Compile all structures from each folder, not just the lowest energy ones')
    parser.add_argument('-o', '--output', dest='summary_file', action='store', default='Goodvibes_output_summary.csv',
                        help='Name of .csv file to write summarized Goodvibes data to')
    parser.add_argument('-g', '--gv_output', dest='gv_file', action='store', default='Goodvibes_output.csv',
                        help='Optional: Name of Goodvibes output file. Default: Goodvibes_output.csv')
    parser.add_argument('-c', '--config', dest='config_file', action='store', default=None,
                        help='Optional: Name of .config file contianing goodvibes settings')
    return parser.parse_args()


@dataclass(kw_only=True)
class GV_Summarizer:
    """Generates a .csv file containing summarized Goodvibes results for all
    .log files in the specified folders"""
    
    dirs: list[str] = None
    qs: str = 'truhlar'
    f_cutoff: str = '100'
    conc: str = '1'
    keywords: str = ''
    lowest_confs_only: bool = True
    gv_file: str = 'Goodvibes_output.csv'
    summary_file: str = 'Goodvibes_output_summary.csv'
    
    summary_df: pd.DataFrame = None
    
    @property
    def log_paths(self):
        """Returns a dictionary of completed, incomplete, and errored .log files"""
        
        logs = log_files(self.dirs)
        log_paths = [f'{folder}/{file}' for folder, files in logs.items()
                                        for file in files]
        log_paths = [f for f in log_paths if is_gaussian_logfile(f)]
        
        completed = [f for f in log_paths if normal_termination(f)]
        errored = [f for f in log_paths if error_termination(f)]
        incomplete = [f for f in log_paths if f not in completed + errored]

        return {'completed': completed, 'incomplete': incomplete, 'errored': errored}

    
    def run_goodvibes(self) -> None:
        """Runs goodvibes with the specified options"""
        
        # Check for consistent linked jobs
        completed_jobs = self.log_paths['completed']
        has_linked = [has_linked_jobs(x) for x in completed_jobs]
        consistent_linked = all([x == has_linked[0] for x in has_linked])
        if not consistent_linked:
            raise InconsistentLinkedJobsError()
        linked = has_linked[0]
        
        # Run Goodvibes
        command = ['python', '-m', 'goodvibes'] + completed_jobs \
                + ['--qs', self.qs, '-f', self.f_cutoff, '-c', self.conc,
                   '--imag', '--csv', '--output', 'temp']
        command += self.keywords.split()
        command += ['--spc', 'link'] if linked else []
        subprocess.run(command)
        os.replace('Goodvibes_temp.csv', self.gv_file)
    
    def summarize(self) -> None:
        """Creates a DataFrame containing only the key thermodynamic data:
        molecule name, imaginary frequencies, G (kcal/mol with all
        corrections applied)"""
        
        results = self.read_gv_csv()
        cols = list(results)
        g_name = [x for x in cols if 'qh-G(T)' in x][-1]
        results['G(kcal/mol)'] = results[g_name] * AU_TO_KCAL
        
        if self.lowest_confs_only:
            results['molname'] = results['Structure'].map(molname_from_log)
            lowest_confs_idx = results.groupby('molname')['G(kcal/mol)'].idxmin()
            results = results.loc[lowest_confs_idx]
        
        self.summary_df = results[['Structure', 'im_freq', 'G(kcal/mol)']]
    
    def read_gv_csv(self) -> None:
        """Parses .csv file of Goodvibes output and returns a dataframe"""
    
        self.fix_gv_columns()
        error_lines, thermo_lines = [], []
        
        with open(self.gv_file) as f:
            is_thermo_data = False
            elements_in_line = None
            for line in f:
                if '*****' in line:
                    continue
                if 'Warning! Couldn\'t find frequency information ...' in line:
                    error_lines += line
                    continue
                if line.startswith('   Structure,'):
                    is_thermo_data = True
                    elements_in_line = line.count(',') # equal because Goodvibes adds an extra comma
                if is_thermo_data:
                    if line.count(',') == elements_in_line:
                        line = line.replace(',\n', '\n')
                    thermo_lines += line[3:] # Remove spaces and/or bullet point
    
        csv_string = ''.join(thermo_lines)
        df = pd.read_csv(StringIO(csv_string))
        return df

    def report(self) -> None:
        """Print and write the summarized Goodvibes data"""
        
        if self.summary_df is None:
            raise NoResultsToReportError()
        
        # Generate messages about incomplete/failed jobs
        incomplete = self.log_paths['incomplete']
        errored = self.log_paths['errored']
        
        warning = ''
        if incomplete:
            warning += '\nThe following jobs did not complete and were not considered:\n'
            for log_path in incomplete:
                warning += f'o  {log_path.split("/")[-1]}\n'
        if errored:
            warning += '\nThe following jobs errored and were not considered:\n'
            for log_path in errored:
                warning += f'o  {log_path.split("/")[-1]}\n'
        warning += '\n' if warning else ''
        
        # Print messages
        dashes = '*' * 70
        template_header = '{:<45}{:>9}{:>16}'
        template_data = '{:<45}{:>9}{:>16.6f}'
        
        print_indent(dashes)
        print_indent(template_header.format('Structure', 'im_freq', 'G(kcal/mol)'))
        print_indent(dashes)
        for row in self.summary_df.itertuples(index=False):
            name, im_freq, g = row
            if np.isnan(im_freq):
                im_freq = ''
            print_indent(template_data.format(name, im_freq, g))
        print_indent(dashes)
        print_indent(warning)
                
        # Write results
        self.summary_df.to_csv(self.summary_file, index=False)
        with open(self.summary_file, 'a') as f:
            f.write(warning)
        
    def fix_gv_columns(self) -> None:
        """Fixes a typo in the column names of the goodvibes output file"""
        
        filedata = None
        with open(self.gv_file) as f:
            filedata = f.read()
        with open(self.gv_file, 'w') as f_out:
            f_out.write(filedata.replace(',im,freq', ',im_freq'))


class InconsistentLinkedJobsError(Exception):
    """Raised when some .log files contain linked Gaussian jobs but
    others do not"""
    pass


class NoResultsToReportError(Exception):
    """Raised when there are no summarized Goodvibes results to report"""
    pass


def main() -> None:
    args = cmdline_parser()
    dirs = glob_dirs(args.dirs)
    lowest_confs_only = not args.use_all
    config_file = args.config_file
    gv_file = args.gv_file
    summary_file = args.summary_file

    config = Config.from_yaml(config_file)
    gv_config = config.goodvibes    
    
    summarizer = GV_Summarizer(dirs=dirs,
                               qs=gv_config['qs'],
                               f_cutoff=gv_config['f_cutoff'],
                               conc=gv_config['conc'],
                               keywords=gv_config['keywords'],
                               gv_file=gv_file,
                               lowest_confs_only=lowest_confs_only,
                               summary_file=summary_file)    
    print('\nRunning Goodvibes...\n')
    summarizer.run_goodvibes()
    print('\nSummary of Goodvibes results:\n')
    summarizer.summarize()
    summarizer.report()
    
if __name__ == '__main__':
    main()
