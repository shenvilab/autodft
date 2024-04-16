import os
import sys
import glob


def log_files(dirs: list[str]) -> list[str]:
    """Takes a sequence of directory names and returns a dictionary of
    directories to .log files"""

    dirs = [d for d in dirs if os.path.isdir(d)]
    logfile_dict = {}

    for d in dirs:
        logfiles = files_in_dir('.log', directory=d)
        if logfiles:
            logfile_dict[d] = logfiles
    if not logfile_dict:
        print('')
        print('  No .log files in the specified directory(s)')
        print('')
        sys.exit(0)
    logfile_dict = dict(sorted(logfile_dict.items()))
    return logfile_dict


def log_files_cwd() -> list[str]:
    """Returns a list of .log files in the current directory.
    
    For consistent output with log_files(), this list is given as the value
    of a dictionary whose key is the current directory (given as '.')"""

    logfiles = files_in_dir('.log')
    if not logfiles:
        print('')
        print('  No .log files in the current directory')
        print('')
        sys.exit(0)
    return {'.': logfiles}

    
def files_in_dir(extension: str, directory: str=None):
    """Returns an ordered list of numbered files with the given extension.
    For example: moleculename-conf1.log, moleculename-conf2.log, etc.
    Looks in current directory if none specified"""
    
    files = [x for x in os.listdir(directory) if x.endswith(extension)]
    files.sort(key=lambda x: (filebase(x), filenum(x)))
    return files


def filenum(filename: str) -> int:
    """For a file name that has a number before the extension,
    returns the number.
    
    Example: filenum('Ni1a-conf2.log') # returns 2 """
    
    root = filename.split('.')[0]
    root = root.split('(lowest)')[0]
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


def filebase(filename: str) -> str:
    """Returns the root of a file name with the following removed:
        -the extension if any
        -'(lowest)' at the end of the file name (excluding the extension)
        -a number suffix at the end of the file name
    
    Example: filebase('Ni1a-conf2.log') # returns 'Ni1a-conf'
    Example: filebase('Ni1a-conf5(lowest).log') # returns 'Ni1a-conf' """

    filename = filename.replace('(lowest).log', '.log')
    filename = filename.removesuffix('(lowest)')
    root = filename.split('.')[0]
    num = ''
    try:
        num = str(filenum(filename))
    except FileWithoutNumberEndingError:
        pass
    return root.removesuffix(num)


def molname_from_log(log: str) -> str:
    """Returns the molecule name from a .log file name formatted
    as molname-confxx.log"""
    
    base = filebase(log)
    return base.removesuffix('-conf')


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
