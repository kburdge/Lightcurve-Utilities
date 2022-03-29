#!/usr/bin/env python


import sys
if sys.version < '2.6':
    raise ImportError("Python versions older than 2.6 are not supported.")
import glob
import os.path

from setuptools import (setup, find_packages)

# set basic metadata
PACKAGENAME = 'KBB_Utils'
DISTNAME = 'KBB_Utils'
AUTHOR = 'Kevin Burdge'
AUTHOR_EMAIL = 'kburdge@mit.edu'
LICENSE = 'GPLv3'



# -- dependencies -------------------------------------------------------------

setup_requires = [
    'setuptools',
]
install_requires = [
    'matplotlib',
    'numpy',
    'cuvarbase',
    'astropy',
]


# -- run setup ----------------------------------------------------------------

packagenames = find_packages()
scripts = glob.glob(os.path.join('bin', '*'))


setup(name=DISTNAME,
      provides=[PACKAGENAME],
      version='0.0.1',
      description=None,
      long_description=None,
      author=AUTHOR,
      author_email=AUTHOR_EMAIL,
      license=LICENSE,
      packages=packagenames,
      include_package_data=True,
      scripts=scripts,
      setup_requires=setup_requires,
      install_requires=install_requires,
      use_2to3=False,
      classifiers=[
          'Programming Language :: Python',
          'Development Status :: 3 - Alpha',
          'Intended Audience :: Science/Research',
          'Intended Audience :: End Users/Desktop',
          'Intended Audience :: Science/Research',
          'Natural Language :: English',
          'Topic :: Scientific/Engineering',
          'Topic :: Scientific/Engineering :: Astronomy',
          'Topic :: Scientific/Engineering :: Physics',
          'Operating System :: POSIX',
          'Operating System :: Unix',
          'Operating System :: MacOS',
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
      ],
)
print(packagenames)
