from setuptools import setup, find_packages
from numpy.distutils.core import setup, Extension

import glob
import subprocess

sys_includes = []
subprocess.call(['make', '-C', 'src'])
sys_includes += glob.glob('src/*.o')
sys_includes += glob.glob('src/*.so')
sys_includes += glob.glob('src/*.mod')

# interface for fortran code
ext1 = Extension(name    = 'fillit._fillBarnes',
                 sources = ['src/fillBarnes.pyf', 'src/fillBarnes.f90'],
                 extra_link_args = ['src/queues.o'])

ext2 = Extension(name    = 'fillit._fillWei',
                 sources = ['src/fillWei.pyf', 'src/fillWei.f90'],
                 extra_link_args = ['src/queues.o'])

ext3 = Extension(name    = 'fillit._fillZhou',
                 sources = ['src/fillZhou.pyf', 'src/fillZhou.f90'],
                 extra_link_args = ['src/queues.o'])

if __name__ == "__main__":
    setup(name = 'fillit',
          author            = "Tristan Salles",
          author_email      = "tristan.salles@sydney.edu.au",
          url               = "https://github.com/Geodels/fillit",
          version           = "0.0.1",
          description       = "A Python interface to compute pit filling using priority-flood on structured and unstructured meshes.",
          long_description = open('README.md').read(),
          long_description_content_type = "text/markdown",
          ext_modules       = [ext1,ext2,ext3],
          packages          = ['fillit'],
          package_data      = {'fillit': ['Notebooks/*ipynb',
                                          'Notebooks/Data/*'] ,'fillit':
                                          sys_includes},
          data_files=[('fillit',sys_includes)],
          include_package_data = True,
          classifiers       = ['Programming Language :: Python :: 2',
                               'Programming Language :: Python :: 2.6',
                               'Programming Language :: Python :: 2.7',
                               'Programming Language :: Python :: 3',
                               'Programming Language :: Python :: 3.3',
                               'Programming Language :: Python :: 3.4',
                               'Programming Language :: Python :: 3.5',
                               'Programming Language :: Python :: 3.6']
          )
