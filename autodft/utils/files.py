import os
import sys
import glob


# CONSTANTS
INDENT = ' ' * 3


def print_indent(msg: str) -> None:
    """Prints an indented message"""

    for line in msg.splitlines():
        print(INDENT + line)


def log_files(dirs: list[str] = None) -> list[str]:
    """Takes a sequence of directory names and returns a dictionary of
    directories to .log files. If no directories are given, the current
    directory is used (denoted as '.')"""

    dirs = ['.'] if not dirs else [d for d in dirs if os.path.isdir(d)]
    logfile_dict = {}

    for d in dirs:
        logfiles_in_dir = [x for x in os.listdir(d) if x.endswith('.log')]
        logfiles_in_dir.sort(key=lambda x: (molname_from_log(x), filenum(x)))
        if logfiles_in_dir:
            logfile_dict[d] = logfiles_in_dir
    if not logfile_dict:
        print('\n  No .log files in the specified directory(s)\n')
        sys.exit()
    logfile_dict = dict(sorted(logfile_dict.items()))
    return logfile_dict


def filenum(filename: str) -> int:
    """For a file name that has a number before the extension,
    returns the number. Otherwise returns None.
    
    Example: filenum('Ni1a-conf2.log') # returns 2 """
    
    root = filename.split('.')[0]
    root = root.removesuffix('(lowest)')
    num = ''
    for char in root[::-1]:
        if not char.isdigit():
            break
        else:
            num = char + num
    try:
        return int(num)
    except ValueError:
        return None


def molname_from_log(filename: str) -> str:
    """"Returns the molecule name from a .log file name.
    If a label '-confxx' is present, it is removed."""

    root = filename.split('.')[0]
    num = filenum(filename)
    return root if num is None else root.removesuffix(f'-conf{num}')


def glob_dirs(dirs: list[str]) -> list[str]:
    """Returns a list of directories with all wildcards parsed"""
    
    matches = []
    for pattern in dirs:
        matches += glob.glob(pattern)
    match_dirs = [m for m in matches if os.path.isdir(m)]
    return match_dirs
    
    