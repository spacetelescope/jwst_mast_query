#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed July 28 11:01:23 2021

@author: arest
"""
from astroquery.mast import Mast
from astropy.time import Time
import numpy as np
import pandas as pd
import astropy.io.fits as fits
import astroquery
import argparse,os,sys,re,io
import yaml

from jwst_query import query_mast


# MAST API documentation:
# https://mast.stsci.edu/api/v0/pyex.html

# MAST CAOM Field Descriptions:
# https://mast.stsci.edu/api/v0/_c_a_o_mfields.html

# File naming schemes:
# https://jwst-pipeline.readthedocs.io/en/latest/jwst/data_products/file_naming.html
# https://jira.stsci.edu/browse/JSDP-1778


if astroquery.__version__<'0.4.2':
    raise RuntimeError("astroquery version=%s, at least 0.4.2 required!" % astroquery.__version__)

class download_mast(query_mast):
    def __init__(self):
        '''
        Set server and web service endpoints for astroquery
        '''
        query_mast.__init__(self)

    def define_options(self,**args):
        parser = query_mast.define_options(self,**args)

        # options for output directory file
        if 'JWSTDOWNLOAD_OUTDIR' in os.environ and os.environ['JWSTDOWNLOAD_OUTDIR'] != '':
            outrootdir = os.environ['JWSTDOWNLOAD_OUTDIR']
        else:
            outrootdir = '.'
        #print('Default output rootdir: %s' % outrootdir)
    
        parser.add_argument('-d','--download', action='store_true', default=False, help='download the data.  (default=%(default)s)')
        parser.add_argument('-o','--outrootdir', default=outrootdir, help=('output root directory.''(default=%(default)s)'))
        parser.add_argument('--outsubdir', default=None, help=('subdir added to the output root directory. If None, then propID is used (default=%(default)s)'))
        parser.add_argument('--clobber', action='store_true', default=False, help='existing files are overwritten, otherwise they are skipped (default=%(default)s)')

        return(parser)

    def fetch_products(self, ix_selected_products=None):
        '''
        Retrieve specified data products from the MAST web service
        '''
        
        if ix_selected_products is None:
            ix_selected_products = self.ix_selected_products
        
        uri_list_str = '&'.join(["uri=" + x for x in self.productTable.t.loc[ix_selected_products,'dataURI']])
        print(uri_list_str)      

#        uri_list_str = '&'.join(["uri=" + x for x in product_table["dataURI"]])
        now_str = ''.join(str(x).zfill(2) for x in Time.now().ymdhms).split('.')[0]
        
        shellscript_filename = self.outdir+'/mastDownload_' + now_str + '.sh'
        makepath4file(shellscript_filename)
        
        print('Saving shell script:',shellscript_filename)        
        #print(uri_list_str)
        download_url = self.JwstObs._portal_api_connection.MAST_BUNDLE_URL + ".sh?" + uri_list_str
        #print('vvv',download_url)
        self.JwstObs._download_file(download_url, shellscript_filename)

if __name__ == '__main__':

    download = download_mast()
    parser = download.define_options()
    args = parser.parse_args()

    # the arguments are saved in download.params
    download.get_arguments(args)
    if download.verbose>2:
        print('params:', download.params)
        
    # set the outdir, based on arguments --outrootdir, --outsubdir. If no outsubdir is specified, the propID is used.
    download.set_outdir()

    # get the MJD limits, based on --mjdlimits --lockbacktime --mjdmin --mjdmax --lre3 --lre4
    mjd_min, mjd_max = download.get_mjd_limits()

        

        

