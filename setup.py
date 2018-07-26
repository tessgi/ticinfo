#!/usr/bin/env python
import sys
import os
from setuptools import setup

if "publish" in sys.argv[-1]:
    os.system("python setup.py sdist")
    os.system("twine upload dist/*")
    os.system("rm -rf dist/*")
    sys.exit()

# Load the __version__ variable without importing the package
exec(open('toco/version.py').read())

# Command-line tools
entry_points = {'console_scripts': [
    'toco = toco.toco:toco',
    'tocot = toco.toco:toco_simbad',
    'tococ = toco.toco:toco_coords',
    'tocon = toco.toco:toco_name',
]}

setup(name='ticinfo',
      version=__version__,
      description="Information for a TESS target",
      long_description="",
      author='Tom Barclay',
      author_email='tom@tombarclay.com',
      license='MIT',
      url='https://github.com/tessgi/toco',
      packages=['toco'],
      install_requires=[
                        'astroquery',
                        'astropy'
                        ],
      entry_points=entry_points,
      include_package_data=True,
      classifiers=[
          "Development Status :: 5 - Production/Stable",
          "License :: OSI Approved :: MIT License",
          "Operating System :: OS Independent",
          "Programming Language :: Python",
          "Intended Audience :: Science/Research",
          "Topic :: Scientific/Engineering :: Astronomy",
          ],
      )
