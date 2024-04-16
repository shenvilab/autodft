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

    dirs = ['.'] if dirs is None else [d for d in dirs if os.path.isdir(d)]
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
    returns the number.
    
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
        raise FileWithoutNumberEndingError()


def molname_from_log(filename: str) -> str:
    """"Returns the molecule name from a .log file name formatted
    as molname-confxx.log."""

    root = filename.split('.')[0]
    num = str(filenum(filename))
    return root.removesuffix(num).removesuffix('-conf')


def glob_dirs(dirs: list[str]) -> list[str]:
    """Returns a list of directories with all wildcards parsed"""
    
    matches = []
    for pattern in dirs:
        matches += glob.glob(pattern)
    match_dirs = [m for m in matches if os.path.isdir(m)]
    return match_dirs
    
    
class FileWithoutNumberEndingError(Exception):
    """Raised when the file name (not including the extension) does not
    end with a number"""
    pass
