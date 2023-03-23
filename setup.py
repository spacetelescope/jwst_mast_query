#!/usr/bin/env python
import os
import pkgutil
import subprocess
import sys
from setuptools import setup, find_packages, Extension, Command
from setuptools.command.test import test as TestCommand


# allows you to build sphinx docs from the package
# main directory with "python setup.py build_sphinx"

try:
    from sphinx.cmd.build import build_main
    from sphinx.setup_command import BuildDoc

    class BuildSphinx(BuildDoc):
        """Build Sphinx documentation after compiling C source files"""

        description = 'Build Sphinx documentation'

        def initialize_options(self):
            BuildDoc.initialize_options(self)

        def finalize_options(self):
            BuildDoc.finalize_options(self)

        def run(self):
            build_cmd = self.reinitialize_command('build_ext')
            build_cmd.inplace = 1
            self.run_command('build_ext')
            build_main(['-b', 'html', './docs', './docs/_build/html'])

except ImportError:
    class BuildSphinx(Command):
        user_options = []

        def initialize_options(self):
            pass

        def finalize_options(self):
            pass

        def run(self):
            print('!\n! Sphinx is not installed!\n!', file=sys.stderr)
            exit(1)

DOCS_REQUIRE = [
    'nbsphinx',
    'sphinx',
    'sphinx-automodapi',
    'sphinx-rtd-theme',
    'stsci-rtd-theme',
    'extension-helpers',
]
TESTS_REQUIRE = [
    'pytest',
]

setup(
    name='jwst_mast_query',
    description='Query MAST for JWST products, and download them',
    long_description=('Query MAST for JWST observations, select certain types of these products, provide a summary of these observations and products, and optionally download the selected products'),
    author='STScI (Rest)',
    author_email='arest@stsci.edu',
    url='https://github.com/spacetelescope/jwst_mast_query',
    keywords=['astronomy'],
    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
    packages = find_packages(),
    package_data={},
    use_scm_version=True,
    setup_requires=['setuptools_scm'],
    install_requires=[
        'astroquery>=0.4.2',
        'astropy',
        'ipython',
        'jinja2',
        'numpy',
        'pandas',
        'pyyaml',
        'scipy',
        'argparse'
    ],
    scripts=['jwst_mast_query/scripts/jwst_download.py'],
    include_package_data=True,
    extras_require={
        'docs': DOCS_REQUIRE,
        'test': TESTS_REQUIRE,
    },
    tests_require=TESTS_REQUIRE,
    cmdclass={
        'build_sphinx': BuildSphinx
    },)
