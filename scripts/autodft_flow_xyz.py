#!/usr/bin/python3

import logging

from autodft.config.config import Config
from autodft.helpers.crest_job import CrestJob
from autodft.helpers.gaussian_inputfile import GaussianInputWriter
from autodft.helpers.gaussian_manager import GaussianManager
from autodft.utils.autodft_utils import input_data, write_autodft_parameters, cleanup


def setup_logger(logging_file):
    logging.basicConfig(filename=logging_file, level=logging.INFO,
                        format='%(asctime)s   %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')
    logging.getLogger('cclib').propagate = False


def main():

    # Get user inputs and configuration settings
    inputs = input_data()
    molname = inputs['molname']
    
    # Get all other settings from .yaml config file
    config_file = inputs['config']
    config = Config.from_yaml(config_file)
    gaussian_config = config.gaussian_input

    logging_file = f'{molname}.out'
    setup_logger(logging_file)
    write_autodft_parameters(logging_file, inputs, config)

    try:

        # Initialization
        gwriter = GaussianInputWriter.from_dicts(
            job_info_list=inputs['gaussian_jobs'],
            nprocshared=gaussian_config['nprocshared'],
            mem=gaussian_config['mem'],
            write_chk=gaussian_config['write_chk'])
        logging.info('Initialized GaussianInputWriter.')
        gsubmitter = GaussianManager(
            molname=molname,
            charge=inputs['charge'],
            multiplicity=inputs['multiplicity'],
            gwriter=gwriter,
            **gaussian_config)
        logging.info('Initialized GaussianManager.')

        # Set up and run CREST job
        crest_job = CrestJob(
            molname=molname,
            xyzfile=inputs['xyzfile'],
            charge=inputs['charge'],
            multiplicity=inputs['multiplicity'],
            e_cutoff=inputs['e_cutoff'],
            **config.crest)

        # Set up and run Gaussian jobs
        gsubmitter.submit_from_crest_confs(crest_job=crest_job)

    except Exception:
        logging.exception('AutoDFT terminated with an error\n\n')
    
    else:
        logging.info('Normal termination of AutoDFT.')
            
    finally:
        cleanup(molname)


if __name__ == '__main__':
    main()
