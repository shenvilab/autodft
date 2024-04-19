import os
import yaml
from dataclasses import dataclass


# CONSTANTS
DEFAULT_CONFIG_FILE = 'default_parameters.yaml'


@dataclass
class Config:
    """A class that stores configuration information for autodft"""
    
    autodft_flow: dict
    crest: dict
    gaussian_input: dict
    gaussian_jobs: dict    
    goodvibes: dict


    @classmethod
    def default(cls):
        """Constructor that creates a Config object from
        the default .yaml configuration file"""
        
        script_dir = os.path.dirname(os.path.abspath(__file__))
        config_file_path = os.path.join(script_dir, DEFAULT_CONFIG_FILE)
        with open(config_file_path) as f:
            config_dict = yaml.safe_load(f)
            return cls(**config_dict)
    
    
    @classmethod
    def from_yaml(cls, yaml_file: str=None):
        """Constructor that creates a Config object from
        a custom .yaml configuration file"""
        
        if yaml_file is None:
            return cls.default()
            
        with open(yaml_file) as f:
            config_dict = yaml.safe_load(f)
            return cls(**config_dict)
        
        
