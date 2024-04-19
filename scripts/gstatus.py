#!/usr/bin/python3

import argparse

from autodft.utils.files import log_files, glob_dirs, print_indent


# CONSTANTS
INDENT = ' ' * 2
DASHES = '-' * 75


def cmdline_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(dest='dirs', nargs='*', action='store',
        help='Folders(s) to run script on, separated by spaces. Use * for wildcard. If none, will use current directory.')
    return parser.parse_args()


def status(logfile: str, folder: str = '.') -> str:
    """Returns a string of status information for a Gaussian .log file
    in the given directory. If no directory specified, will use the
    current directory"""

    # Initialize variables
    optimizations_done = 0
    opt_step_number = '-'
    jobs_completed = 0
    jobs_failed = '-'
    converged_maxforce = '-'
    converged_rmsforce = '-'
    converged_maxdisplacement = '-'
    converged_rmsdisplacement = '-'

    # Initialize variables for scans only
    is_scan = False
    max_found = False
    energies = []
    energy = None

    logfile_path = f'{folder}/{logfile}'
    with open(logfile_path) as f:
        for line in f:
            if 'Optimized Parameters' in line:
                optimizations_done += 1
                opt_step_number = '-'
                if is_scan:
                    energies.append(energy)
                    energy = None
            if 'Step number' in line:
                opt_step_number = line.split()[2]
            if 'Maximum Force' in line:
                converged_maxforce = line.split()[-1][0]
            if 'RMS     Force' in line:
                converged_rmsforce = line.split()[-1][0]
            if 'Maximum Displacement' in line:
                converged_maxdisplacement = line.split()[-1][0]
            if 'RMS     Displacement' in line:
                converged_rmsdisplacement = line.split()[-1][0]
            if 'Normal termination' in line:
                jobs_completed += 1
                opt_step_number = '-'
                converged_maxforce = converged_rmsforce = converged_maxdisplacement = converged_rmsdisplacement = '-'
            if 'Error termination' in line:
                jobs_failed = 'Y'
                opt_step_number = '-'
                converged_maxforce = converged_rmsforce = converged_maxdisplacement = converged_rmsdisplacement = '-'
            if ('#' in line) and any(x in line.lower() for x in ('opt', 'modredundant')):
                is_scan = True if not is_scan else is_scan
            if ('SCF Done:' in line) and is_scan:
                energy = float(line.split()[4])

    convergence_criteria = ''.join(
        [converged_maxforce, converged_rmsforce, converged_maxdisplacement, converged_rmsdisplacement])
    convergence_criteria = convergence_criteria.replace('N', '-')
    optimizations_done = '-' if not optimizations_done else optimizations_done
    jobs_completed = '-' if not jobs_completed else jobs_completed

    if is_scan:
        if energy is not None:  # Append energy of last structure to list, even of optimization of that structure has not yet completed
            energies.append(energy)
        if len(energies) >= 2:
            dE_signs = [e2 - e1 > 0 for e1, e2 in zip(energies, energies[1:])]
            max_found = any(
                [dE2 - dE1 == -1 for dE1, dE2 in zip(dE_signs, dE_signs[1:])])
        if max_found:
            optimizations_done = f'{optimizations_done}*'

    logfile_root = logfile.removesuffix('.log')
    status = '{:<35}{:>4}{:>8}{:>8}{:>8}{:>12}'.format(
        logfile_root[:36], optimizations_done, jobs_completed,
        jobs_failed, opt_step_number, convergence_criteria)
    return status


def print_header() -> None:
    """Prints the header"""

    print('')
    print_indent(DASHES)
    print_indent('{:<23}{:>16}{:>8}{:>8}{:>8}{:>12}'.format(
        'File', 'Opts', 'Jobs', 'Jobs', 'Opt', 'Converged?'))
    print_indent('{:<23}{:>16}{:>8}{:>8}{:>8}{:>12}'.format(
        '', 'done', 'done', 'failed', 'step', ''))
    print_indent(DASHES)


def print_footer() -> None:
    """Prints the footer"""

    print_indent(DASHES)
    print_indent(
        'Convergence info: max force, RMS force, max displacement, RMS displacement')
    print_indent(
        'For scan jobs: Star (*) indicates that an energy maximum has been found')
    print_indent(DASHES)
    print('')


def main() -> None:
    """Print status of all .log files in the user-specified folder"""

    args = cmdline_parser()
    dirs = glob_dirs(args.dirs)

    status_msgs = []
    logs = log_files(dirs)
    for folder, files in logs.items():
        status_msgs += [status(file, folder=folder) for file in files]

    if status_msgs:
        print_header()
        [print_indent(s) for s in status_msgs]
        print_footer()
    else:
        print('\n\n')


if __name__ == '__main__':
    main()
