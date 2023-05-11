#! /usr/bin/env python

"""Tests for jwst_download.py script calls. These tests are closer to regression tests than unit tests

Author
------

B. Hilbert
"""
from glob import glob
import os
import pytest

import numpy as np
import pandas as pd

import jwst_mast_query


# Get the directory where the config files for testing are located
CFG_DIR = os.path.join(os.path.dirname(jwst_mast_query.__file__), 'tests/testing_config_files')

# Determine if tests are being run on Github Actions
ON_GITHUB_ACTIONS = '/home/runner' in os.path.expanduser('~') or '/Users/runner' in os.path.expanduser('~')


commands = ['jwst_download.py -i nircam --propID 1068 --obsnums 3 4 --filetypes rate --Nobs_per_batch 3 --outrootdir ./ --date_select 2022-05-15 2022-05-16 --savetables 1068_obs34_query_cmd --skipdownload',
            'jwst_download.py -i niriss --propID 1063 --obsnums 108 --filetypes rate --Nobs_per_batch 3 --outrootdir ./ --date_select 2022-04-01 2022-04-02 --savetables 1063_obs108_query_cmd --skipdownload',
            'jwst_download.py -i nircam --propID 1068 --obsnums 3 --Nobs_per_batch 3 --outrootdir ./ --date_select 2022-05-15 2022-05-16 --savetables 1068_obs3_query_with_guidestars_cmd --skipdownload --guidestars',
            'jwst_download.py -i nircam --propID 1068 --obsnums 3 --Nobs_per_batch 3 --outrootdir ./ --date_select 2022-05-15 2022-05-16 --savetables 1068_obs3_query_guidestars_only_cmd --skipdownload --guidestar_data_only']

config_files = ['1068_obs34_query.cfg',
                '1063_obs108_query.cfg',
                '1068_obs3_query_with_guidestars.cfg',
                '1068_obs3_query_guidestars_only.cfg'
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


@pytest.mark.skipif(ON_GITHUB_ACTIONS, reason="No MAST token available when using Github Actions")
def test_guidestar_download():
    """Test that guidestar data can be downloaded if requested

    Parameters
    ----------
    command : str
        Command line call
    """
    sci_plus_gs_cmd = 'jwst_download.py -i nircam --propID 1068 --obsnums 3 --Nobs_per_batch 3 --outrootdir ./ --date_select 2022-05-15 2022-05-16 --savetables sci_plus_gs --skipdownload --guidestars'
    gs_only_cmd = 'jwst_download.py -i nircam --propID 1068 --obsnums 3 --Nobs_per_batch 3 --outrootdir ./ --date_select 2022-05-15 2022-05-16 --savetables gs_only --skipdownload --guidestar_data_only'
    sci_only_cmd = 'jwst_download.py -i nircam --propID 1068 --obsnums 3 --Nobs_per_batch 3 --outrootdir ./ --date_select 2022-05-15 2022-05-16 --savetables sci_only --skipdownload'

    # Run the commands
    os.system(sci_plus_gs_cmd)
    os.system(gs_only_cmd)
    os.system(sci_only_cmd)

    # Check for the existence of guidestar files in each
    both_table = pd.read_fwf('sci_plus_gs.selprod.txt', header=0)
    gs_table = pd.read_fwf('gs_only.selprod.txt', header=0)
    sci_table = pd.read_fwf('sci_only.selprod.txt', header=0)

    both_file = np.array(['gs' in e for e in both_table['productFilename']])
    gs_file = np.array(['gs' in e for e in gs_table['productFilename']])
    sci_file = np.array(['gs' not in e for e in sci_table['productFilename']])

    assert np.all(gs_file)
    assert np.all(sci_file)
    assert np.sum(both_file) < len(both_table['productFilename'].values)

    # Remove table files
    bases = ['sci_plus_gs', 'gs_only', 'sci_only']
    for base in bases:
        for filename in glob(f'{base}*txt'):
            os.remove(filename)
