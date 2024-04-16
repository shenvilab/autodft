#!/usr/bin/python3

import numpy as np

import argparse

from autodft.config.config import Config
from autodft.helpers.goodvibes_data import GV_Executor, GV_Results, AU_TO_KCAL
from autodft.utils.files import log_files, log_files_cwd, molname_from_log, glob_dirs
from autodft.utils.gaussian_output import has_linked_jobs


def cmdline_parser() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(dest='dirs', nargs='*', action='store',
                        help='Folder(s) to run script on.')
    parser.add_argument('-a', '--all', dest='use_all', action='store_true', default=0,
                        help='Optional: Compile all structures from each folder, not just the lowest energy ones')
    parser.add_argument('-o', '--output', dest='summary_file', action='store', default='Goodvibes_output_summary.csv',
                        help='Name of .csv file to write summarized Goodvibes data to')
    parser.add_argument('-g', '--gv_output', dest='gv_output', action='store', default='Goodvibes_output.csv',
                        help='Optional: Name of Goodvibes output file. Default: Goodvibes_output.csv')
    parser.add_argument('-c', '--config', dest='config_file', action='store', default=None,
                        help='Optional: Name of .config file contianing goodvibes settings')
    return parser.parse_args()


def iprint(msg: str) -> None:
    """Prints a message with an indent"""
    
    indent = ' ' * 3
    print(indent + msg)


class GV_Summary:
    """Looks in the specified folders for a goodvibes results files,
    the results of which are concatenated into a single .csv file.

    Note that this makes no attempt to check for consistent levels of theory
    or consistent Goodvibes specifications.
    """

    def __init__(self, dirs: list[str], summary_file: str,
                 gv_output: str='Goodvibes_output.csv',
                 config_file: str = 'Goodvibes_output_summary.csv') -> None:
        self.dirs = dirs
        self.gv_output = gv_output
        self.summary_file = summary_file
        self.config_file = config_file
        
        self.df = None
        self.selected_results = None

    def run(self) -> None:
        """Runs Goodvibes on all .log files in the specified directories
        and stores the results to a dataframe"""
        
        # Get log files
        logs = log_files(self.dirs) if self.dirs else log_files_cwd()        
        iprint('\nRunning Goodvibes...')
        log_paths = []
        for folder, files in logs.items():
            log_paths += [f'{folder}/{file}' for file in files]

        # Check for consistent linked jobs
        has_linked = [has_linked_jobs(x) for x in log_paths]
        consistent_linked = all([x == has_linked[0] for x in has_linked])
        if not consistent_linked:
            raise InconsistentLinkedJobsError()
        linked = has_linked[0]

        # Get configuration settings
        config = Config.from_yaml(self.config_file)
        gv_config = config.goodvibes
        qs = gv_config['qs']
        conc = gv_config['conc']
        f_cutoff = gv_config['f_cutoff']

        # Run Goodvibes
        gv_executor = GV_Executor(files=log_paths, qs=qs, f_cutoff=f_cutoff,
                                  conc=conc, csv=True,
                                  spc='link' if linked else '')
        gv_executor.run(self.gv_output)
        results = GV_Results(self.gv_output)
        self.df = results.df

    def process_results(self, use_all: bool=True) -> None:
        """Generates and stores a dataframe containing a summary of
        Goodvibes data (molecule name, imaginary freq., G (kcal/mol)"""
        
        results = self.df.copy()
        
        # Clean up and add requisite data
        cols = list(results)
        g_name = [x for x in cols if 'qh-G(T)' in x][-1]
        results['G(kcal/mol)'] = results[g_name] * AU_TO_KCAL
        
        # Keep lowest conformation of each molecule only
        if not use_all:
            results['molname'] = results['Structure'].map(molname_from_log)
            lowest_confs_idx = results.groupby('molname')['G(kcal/mol)'].idxmin()
            results = results.loc[lowest_confs_idx]
        
        # Get selected columns with qh-G(T) in kcal/mol
        results = results[['Structure', 'im_freq', 'G(kcal/mol)']]
        self.selected_results = results
        results.to_csv(self.summary_file, index=False)
        
    def print_gv_summary(self) -> None:
        """Prints a summary of the Goodvibes data"""
        
        dashes = '*' * 70
        template_header = '{:<45}{:>9}{:>16}'
        template_data = '{:<45}{:>9}{:>16.6f}'
        iprint(dashes)
        iprint(template_header.format('Structure', 'im_freq', 'G(kcal/mol)'))
        iprint(dashes)
        for row in self.selected_results.itertuples(index=False):
            name, im_freq, g = row
            if np.isnan(im_freq):
                im_freq = ''
            iprint(template_data.format(name, im_freq, g))
        iprint(dashes)
        
    def write_gv_summary(self) -> None:
        """Writes the summarized Goodvibes data to a .csv file"""
        
        self.selected_results.to_csv(self.summary_file, index=False)


class InconsistentLinkedJobsError(Exception):
    """Raised when some .log files contain linked Gaussian jobs but
    others do not"""
    pass


def main() -> None:
    args = cmdline_parser()
    dirs = glob_dirs(args.dirs)
    use_all = args.use_all
    config_file = args.config_file
    gv_output = args.gv_output
    summary_file = args.summary_file

    gv_summary = GV_Summary(dirs=dirs, gv_output=gv_output,
                            summary_file=summary_file,
                            config_file=config_file)
    gv_summary.run()
    
    print('Summary of Goodvibes results:\n')
    gv_summary.process_results(use_all=use_all)
    gv_summary.print_gv_summary()
    gv_summary.write_gv_summary()
    

if __name__ == '__main__':
    main()
