#!/usr/bin/env python3

import os
import shutil
import sys

from setuptools import setup
from setuptools.command.install import install


def readme():
    with open('README.md') as f:
        return f.read()


# Get the program version from another file.
__version__ = '0.0.0'
exec(open('badread/version.py').read())


class BadreadInstall(install):
    """
    The install process copies necessary files (like pre-built models) to the install location.
    """
    if __name__ == '__main__':
        def run(self):
            # Make sure we have permission to write the files.
            if os.path.isdir(self.install_lib) and not os.access(self.install_lib, os.W_OK):
                sys.exit('Error: no write permission for ' + self.install_lib + '  ' +
                         'Perhaps you need to use sudo?')
            if os.path.isdir(self.install_scripts) and not os.access(self.install_scripts, os.W_OK):
                sys.exit('Error: no write permission for ' + self.install_scripts + '  ' +
                         'Perhaps you need to use sudo?')

            install.run(self)

            # Copy error models to installation directory.
            error_models_source_dir = os.path.join('badread', 'error_models')
            error_models_dest_dir = os.path.join(self.install_lib, 'badread', 'error_models')
            if not os.path.exists(error_models_dest_dir):
                os.makedirs(error_models_dest_dir)
            shutil.copyfile(os.path.join(error_models_source_dir, 'nanopore.gz'),
                            os.path.join(error_models_dest_dir, 'nanopore.gz'))
            shutil.copyfile(os.path.join(error_models_source_dir, 'pacbio.gz'),
                            os.path.join(error_models_dest_dir, 'pacbio.gz'))

            # Copy qscore models to installation directory.
            qscore_models_source_dir = os.path.join('badread', 'qscore_models')
            qscore_models_dest_dir = os.path.join(self.install_lib, 'badread', 'qscore_models')
            if not os.path.exists(qscore_models_dest_dir):
                os.makedirs(qscore_models_dest_dir)
            shutil.copyfile(os.path.join(qscore_models_source_dir, 'nanopore.gz'),
                            os.path.join(qscore_models_dest_dir, 'nanopore.gz'))
            shutil.copyfile(os.path.join(qscore_models_source_dir, 'pacbio.gz'),
                            os.path.join(qscore_models_dest_dir, 'pacbio.gz'))


setup(name='Badread',
      version=__version__,
      description='Badread: a long read simulator that can mimic various kinds of read problems',
      long_description=readme(),
      long_description_content_type='text/markdown',
      url='https://github.com/rrwick/Badread',
      author='Ryan Wick',
      author_email='rrwick@gmail.com',
      license='GPLv3',
      packages=['badread'],
      install_requires=['edlib', 'numpy', 'scipy'],
      extras_require={'plot': ['matplotlib']},
      entry_points={"console_scripts": ['badread = badread.badread:main']},
      include_package_data=True,
      zip_safe=False,
      python_requires='>=3.6',
      cmdclass={'install': BadreadInstall})
