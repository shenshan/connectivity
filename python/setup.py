#!/usr/bin/env python

from distutils.core import setup

setup(
    name='cortical_map',
    version='0.1',
    author='Fabian Sinz',
    author_email='sinz@bcm.edu',
    description="Schemata for the analysis of the data recorded by Xiaolong and Shan. ",
    packages=['cortical_map'],
    requires=['numpy', 'pymysql', 'matplotlib','datajoint', 'commons'],
    license = "MIT",
)
