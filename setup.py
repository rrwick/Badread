#!/usr/bin/env python3

from setuptools import setup


def readme():
    with open('README.md') as f:
        return f.read()


# Get the program version from another file.
__version__ = '0.0.0'
exec(open('badread/version.py').read())


setup(name='Badread',
      version=__version__,
      description='Badread: a long read simulator that can mimic various kinds of read problems',
      long_description=readme(),
      url='https://github.com/rrwick/Badread',
      author='Ryan Wick',
      author_email='rrwick@gmail.com',
      license='GPLv3',
      packages=['badread'],
      entry_points={"console_scripts": ['badread = badread.badread:main']},
      include_package_data=True,
      zip_safe=False)
