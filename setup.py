from setuptools import setup, find_packages

with open('README.md') as f:
    long_description=f.read()

setup(
    name='autodft',
    version='1.0.1',
    author='Shenvi Lab',
    author_email='sting@scripps.edu',
    url='https://github.com/shenvilab/autodft',
    description='Automated conformer searching and DFT calculations',
    long_description=long_description,
    long_description_content_type='text/markdown',
    packages=find_packages(),
    scripts=[
        'scripts/run_autodft.py',
        'scripts/run_autodft_xyz.py',
        'scripts/autodft_flow.py',
        'scripts/autodft_flow_xyz.py',
        'scripts/gstatus.py',
        'scripts/compile_results.py'
    ],
    package_data={'autodft': ['config/*.yaml']},
    data_files=[('', ['README.md', 'LICENSE.txt'])],
    include_package_data=True,
    zip_safe=False,
    install_requires=[
        'rdkit>=2022.03.1',
        'pandas',
        'goodvibes>=3.2',
        'cclib',
        'pyyaml>=6.0'],
    python_requires=">=3.10",
    license="GPL-v3",
)