#! /usr/bin/env python

"""Tests for jwst_download.py script calls. These tests are closer to regression tests than unit tests

Author
------

B. Hilbert
"""
from glob import glob
import os
import pytest

import pandas as pd

import jwst_mast_query


# Get the directory where the config files for testing are located
CFG_DIR = os.path.join(os.path.dirname(jwst_mast_query.__file__), 'tests/testing_config_files')

# Determine if tests are being run on Github Actions
ON_GITHUB_ACTIONS = '/home/runner' in os.path.expanduser('~') or '/Users/runner' in os.path.expanduser('~')


commands = ['jwst_download.py -i nircam --propID 1068 --obsnums 3 4 --filetypes rate --Nobs_per_batch 3 --outrootdir ./ --date_select 2022-05-15 2022-05-16 --savetables 1068_obs34_query_cmd --skipdownload',
            'jwst_download.py -i niriss --propID 1063 --obsnums 108 --filetypes rate --Nobs_per_batch 3 --outrootdir ./ --date_select 2022-04-01 2022-04-02 --savetables 1063_obs108_query_cmd --skipdownload'
            ]
config_files = ['1068_obs34_query.cfg',
                '1063_obs108_query.cfg'
                ]
cmd_cfg = [(cmd, os.path.join(CFG_DIR, cfg)) for cmd, cfg in zip(commands, config_files)]


@pytest.mark.skipif(ON_GITHUB_ACTIONS, reason="No MAST token available when using Github Actions")
@pytest.mark.parametrize('command, config_file', cmd_cfg)
def test_jwst_download(command, config_file):
    """Test that a command line call and the equivalent call made using a config
    file return the same result.

    Parameters
    ----------
    command : str
        Command line call

    config_file : str
        Name of config file
    """
    # Run the command that uses no config file
    os.system(command)

    # Run the command that does use the config file
    os.system(f'jwst_download.py -c {config_file}')

    # Compare the resulting tables
    base = os.path.basename(config_file).split('.')[0]
    cmd_table = pd.read_fwf(f'{base}_cmd.selprod.txt', header=0)
    cfg_table = pd.read_fwf(f'{base}_cfg.selprod.txt', header=0)
    assert sorted(cmd_table['obsID'].values)  == sorted(cfg_table['obsID'].values)

    # Remove table files
    for filename in glob(f'{base}*txt'):
        os.remove(filename)
