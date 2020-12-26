import os
import glob
import setuptools
from distutils.core import setup

with open("README.md", 'r') as readme:
    long_description = readme.read()

setup(
    name='vivarium-cobra',
    version='0.0.3',
    packages=[
        'vivarium_cobra',
        'vivarium_cobra.processes',
        'vivarium_cobra.composites',
        'vivarium_cobra.experiments',
        'vivarium_cobra.library',
        'vivarium_cobra.models',
    ],
    author='Eran Agmon',
    author_email='eagmon@stanford.edu',
    url='https://github.com/vivarium-collective/vivarium-cobra',
    license='MIT',
    entry_points={
        'console_scripts': []},
    short_description='a vivarium wrapper for COBRApy',
    long_description=long_description,
    long_description_content_type='text/markdown',
    package_data={},
    include_package_data=True,
    install_requires=[
        'vivarium-core',
        'cobra',
        'Arpeggio',
    ],
)
