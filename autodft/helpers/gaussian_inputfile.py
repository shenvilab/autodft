import periodictable as pt

from dataclasses import dataclass
from typing import Tuple
from textwrap import dedent


@dataclass(kw_only=True)
class MolInfo:
    """A class with molecular information to write to a Gaussian input file.

    :param charge: Charge of the molecule
    :type charge: str
    :param multiplicity: Spin multiplicity of the molecule
    :type multiplicity: str
    :param cartesians:
        A string containing the elements and cartesian coordinates, formatted
        as in an .xyz file (but without the header used in .xyz files).
        The string can end with any number of blank lines (including 0)
    :type cartesians: str
    :param title:
        (Optional) Title to write to the input file.
        Default = 'Title card required'
    :type title: str
    :param elements:
        (Generated after initialization) Elements present in the molecule,
        sorted alphabetically
    :type elements: list[str]
    """

    charge: str
    multiplicity: str
    cartesians: str
    title: str = 'Title card required'

    elements: list[str] = None

    def __post_init__(self) -> None:
        self.elements = [line.split()[0]
                         for line in self.cartesians.splitlines()]
        self.elements = sorted(list(set(self.elements)))

    def __str__(self) -> str:
        """Returns a string of title, charge, multiplicity, cartesians
        formatted for a Gaussian input file"""

        return (f'{self.title}\n\n'
                f'{self.charge} {self.multiplicity}\n'
                f'{self.cartesians.rstrip()}\n\n')


@dataclass(kw_only=True)
class GaussianJob:
    """A class that provides the specifications for an individual Gaussian job

    :param functional:
        The DFT functional to use
    :type functional: str
    :param basisset:
        The basis set to use. Use 'genecp' for split basis sets
    :type basisset: str
    :param route:
        (Optional) Information to include in the route line separate from
        level of theory. Example: 'opt freq scrf=(smd,solvent=water)'
    :type route: str
    :param nprocshared:
        (Optional) Number of processors to request. Default = None
    :type nprocshared: str
    :param mem:
        (Optional) Amount of memory to request. Default = None
    :type mem: str
    :param write_chk:
        (Optional) Whether to write .chk files. Default = True
    :type write_chk: bool
    :param redundant_coords: 
        (Optional) Redundant coordinates for freezes/scans, with new lines
        denoted by '\n'. Example: 'B 24 45 F\nB 33 27 S 10 0.20000'.
        Default = None
    :type redundant_coords: str
    :param ecp_max_lowlevel_atomicnum:
        (Optional) Maximum atomic number to use with the low level of
        theory. Default = None
    :type ecp_max_lowlevel_atomicnum: int
    :param ecp_basisset_light:
        (Optional) Basis set to use for atomic numbers lower than and equal to
        the atomic number cutoff
    :type ecp_basisset_light: str
    :param ecp_basisset_heavy:
        (Optional) Basis set to use for atomic numbers greater than the
        atomic number cutoff
    :type ecp_basiset_high: str
    """

    functional: str
    basisset: str
    route: str = None
    nprocshared: str = None
    mem: str = None
    write_chk: bool = False
    redundant_coords: str = None
    ecp_max_lowlevel_atomicnum: int = None
    ecp_basisset_light: str = None
    ecp_basisset_heavy: str = None

    def redundant_coordblock(self) -> str:
        """Checks and formats the provided redundant coordinates"""

        if not self.redundant_coords:
            return ''

        lines = [line for line in self.redundant_coords.splitlines() if line]
        spaces = {'B': 3, 'A': 4, 'D': 5, 'F': 0, 'S': 2}
        for line in lines:
            n_spaces_actual = line.count(' ')
            n_spaces_correct = sum([spaces[char]
                                   for char in line if char in spaces.keys()])
            if n_spaces_actual != n_spaces_correct:
                raise InvalidRedundantCoordinatesError()

        return self.redundant_coords.rstrip() + '\n\n'

    def ecp_block_specified(self) -> bool:
        """Returns whether the ECP information is specified. Raises errors
        if the ECP information has inconsistencies"""

        n_ecp_specs = sum([bool(x) for x in
                           [self.ecp_max_lowlevel_atomicnum,
                            self.ecp_basisset_light,
                            self.ecp_basisset_heavy]])
        ecp_specified = n_ecp_specs == 3
        genecp_in_route = self.basisset.lower() == 'genecp'
        if n_ecp_specs not in [0, 3]:
            raise PartialECPSpecificationError()
        if ecp_specified != genecp_in_route:
            raise InconsistentECPPresenceError()
        return ecp_specified

    def light_heavy_elements(self, elements: list[str]) -> Tuple[str, str]:
        """Returns whitespace-delimited strings for the elements below and
        above the atomic number cutoff"""

        light = [e for e in elements if pt.elements.symbol(e).number
                                        <= self.ecp_max_lowlevel_atomicnum]
        heavy = [e for e in elements if pt.elements.symbol(e).number
                                        > self.ecp_max_lowlevel_atomicnum]
        light, heavy = ' '.join(light), ' '.join(heavy)
        return light, heavy

    def basisset_and_ecp(self, elements: list[str]) -> Tuple[str, str]:
        """Returns the basis set and ECP block to write to the com file.
        The appropriate basis set is provided if an ECP is specified but
        all elements fall above or below the atomic number cutoff."""

        if not self.ecp_block_specified():
            return self.basisset, ''

        light, heavy = self.light_heavy_elements(elements)
        if not light:
            return self.ecp_basisset_heavy, ''
        if not heavy:
            return self.ecp_basisset_light, ''

        ecp_block = dedent(f'''\
                           {light} 0
                           {self.ecp_basisset_light}
                           ****
                           {heavy} 0
                           {self.ecp_basisset_heavy}
                           ****
                           
                           {heavy} 0
                           {self.ecp_basisset_heavy}
                           
                           ''')
        return self.basisset, ecp_block

    def write(self, elements: list[str]) -> str:
        """Writes a string template for the GaussianJob.

        The string template has fields chk_line, route_extras,
        and mol_info to be filled in by a GaussianInputWriter object

        :param elements: A list of elements in the structure to write
        :type elements: list[str]
        :return: A string template for the GaussianJob
        :rtype: str
        """

        basisset, ecp_block = self.basisset_and_ecp(elements)
        redundant_coordblock = self.redundant_coordblock()

        output = ''
        output += f'%nprocshared={self.nprocshared}\n' if self.nprocshared else ''
        output += f'%mem={self.mem}\n' if self.mem else ''
        output += '{chk_line}'
        output += f'# {self.functional}/{basisset}'
        output += f' {self.route}{{route_extras}}\n\n' if self.route else '{route_extras}\n\n'
        output += '{mol_info}'
        output += f'{redundant_coordblock}{ecp_block}'

        return output


class GaussianInputWriter:
    """A class that stores GaussianJob objects and writes an input file where
    all of the individual jobs are linked together.

    :param jobs: A list of GaussianJob objects
    :type jobs: list[GaussianJob]

    Example usage:
    ::
        # Create test molecule
        mol1 = MolInfo(charge='0', multiplicity='1',
                       cartesians='C 0 0 0\nO 1 0 0\nH 0 1 0\n')

        # Create GaussianJobs
        job1 = GaussianJob(
            functional='b3lyp', basisset='genecp',
            route='opt freq scrf=(smd,solvent=water)',
            nprocshared='12', mem='2gb', write_chk=False,
            ecp_max_lowlevel_atomicnum=6,
            ecp_basisset_light='6-31g(d)', ecp_basisset_heavy='sdd')

        job2 = GaussianJob(
            functional='m062x', basisset='def2tzvp',
            route='scrf=(smd,solvent=water)', nprocshared='4', mem='1gb',
            write_chk=False, ecp_max_lowlevel_atomicnum=None,
            ecp_basisset_light=None, ecp_basisset_heavy=None)

        # Create GaussianInputWriter and write input file
        gwriter = GaussianInputWriter()
        gwriter.add(job1, job2)
        gwriter.write(mol1, 'testfile.com')
    """

    def __init__(self, *args: GaussianJob) -> None:
        """Constructor"""

        self.jobs = []
        [self.jobs.append(job) for job in args]

    @classmethod
    def from_dicts(cls, job_info_list: list[dict], nprocshared: int,
                   mem: str, write_chk: bool = False):
        """Creates a GaussianInputWriter from a dictionary of job 
        information supplied"""

        jobs = [GaussianJob(
            functional=job_info['functional'],
            basisset=job_info['basisset'],
            route=job_info['route'],
            ecp_basisset_light=job_info['ecp_basisset_light'],
            ecp_basisset_heavy=job_info['ecp_basisset_heavy'],
            ecp_max_lowlevel_atomicnum=job_info['ecp_cutoff'],
            nprocshared=nprocshared,
            mem=mem,
            write_chk=write_chk
        ) for job_info in job_info_list]
        return cls(*jobs)

    def add(self, *args: GaussianJob) -> None:
        """Adds any number of GaussianJob objects to the GaussianInputWriter

        :param *args: GaussianJob object(s)
        :type *args: GaussianJob
        :return: None
        :rtype: None
        """

        [self.jobs.append(job) for job in args]

    def write(self, mol_info: MolInfo, com_file: str):
        """Writes a Gaussian input file from the specified jobs and MolInfo

        :param mol_info:
            A MolInfo object containing information about the molecule
        :type mol_info: MolInfo
        :param com_file: The name of the .com file to write
        :type com_file: str
        :return: None
        :rtype: None
        """

        elements = mol_info.elements

        write_chk = any([job.write_chk for job in self.jobs]
                        ) or len(self.jobs) > 1
        chk_file = com_file.replace('.com', '.chk')
        chk_line = f'%chk={chk_file}\n' if write_chk else ''

        output = ''
        for i, job in enumerate(self.jobs):
            output += '\n--Link1--\n' if i > 0 else ''
            output += job.write(elements)
            output = output.replace(
                '{mol_info}', str(mol_info) if i == 0 else '')
            output = output.replace(
                '{route_extras}', ' geom=allcheck guess=read' if i > 0 else '')
        output = output.format(chk_line=chk_line) + '\n'

        with open(com_file, 'w') as f:
            f.write(output)


class PartialECPSpecificationError(Exception):
    """Raised when some but not all of the information to define the
    split basis set / effective core potential is provided"""
    pass


class InconsistentECPPresenceError(Exception):
    """Raised when the presence of ECP information differs between the
    basis set and the ECP-related attributes of the GaussianJob object"""
    pass


class InvalidRedundantCoordinatesError(Exception):
    """Raised when the syntax of redundant coordiantes (for freezes,
    scans) is invalid"""
    pass
