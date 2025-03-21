#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""
import os
import re
from setuptools import setup, find_packages


with open(os.path.join('cellmaps_imagedownloader', '__init__.py')) as ver_file:
    for line in ver_file:
        line = line.rstrip()
        if line.startswith('__version__'):
            version=re.sub("'", "", line[line.index("'"):])
        elif line.startswith('__description__'):
            desc = re.sub("'", "", line[line.index("'"):])
        elif line.startswith('__repo_url__'):
            repo_url = re.sub("'", "", line[line.index("'"):])
        elif line.startswith('__author__'):
            author = re.sub("'", "", line[line.index("'"):])
        elif line.startswith('__email__'):
            email = re.sub("'", "", line[line.index("'"):])

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = ['cellmaps_utils>=0.4.0,<1.0.0',
                'requests>=2.32.3,<3.0.0',
                'mygene>=3.2.2,<4.0.0',
                'pandas>=2.0.0,<3.0.0',
                'tqdm>=4.67.1,<5.0.0']

setup_requirements = [ ]

setup(
    author=author,
    author_email=email,
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
    ],
    description=desc,
    install_requires=requirements,
    license="MIT license",
    long_description=readme + '\n\n' + history,
    long_description_content_type='text/x-rst',
    include_package_data=True,
    keywords='cellmaps_imagedownloader',
    name='cellmaps_imagedownloader',
    packages=find_packages(include=['cellmaps_imagedownloader']),
    package_dir={'cellmaps_imagedownloader': 'cellmaps_imagedownloader'},
    package_data={'cellmaps_imagedownloader': ['readme_outputs.txt']},
    scripts=['cellmaps_imagedownloader/cellmaps_imagedownloadercmd.py'],
    setup_requires=setup_requirements,
    url=repo_url,
    version=version,
    zip_safe=False)
