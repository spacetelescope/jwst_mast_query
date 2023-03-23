#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed July 28 11:01:23 2021

@author: arest
"""
#from astroquery.mast import Mast
from astropy.time import Time
#import numpy as np
#import pandas as pd
import astroquery
import os,sys,shutil,re
#import yaml
from astropy.utils.data import download_file

from jwst_mast_query.jwst_query import query_mast
from jwst_mast_query.pdastro import makepath4file,AnotB
#from pdastro import makepath4file,AnotB


# MAST API documentation:
# https://mast.stsci.edu/api/v0/pyex.html

# MAST CAOM Field Descriptions:
# https://mast.stsci.edu/api/v0/_c_a_o_mfields.html

# File naming schemes:
# https://jwst-pipeline.readthedocs.io/en/latest/jwst/data_products/file_naming.html
# https://jira.stsci.edu/browse/JSDP-1778


if astroquery.__version__<'0.4.2':
    raise RuntimeError("astroquery version=%s, at least 0.4.2 required!" % astroquery.__version__)

def rmfile(filename,raiseError=True):
    " if file exists, remove it "
    if os.path.lexists(filename):
        os.remove(filename)
        if os.path.isfile(filename):
            if raiseError == 1:
                raise RuntimeError('ERROR: Cannot remove %s' % filename)
            else:
                return(1)
    return(0)


class download_mast(query_mast):
    def __init__(self):
        '''
        Set server and web service endpoints for astroquery
        '''
        query_mast.__init__(self)

    def define_options(self,**args):
        parser = query_mast.define_options(self,**args)

        # options for output directory file
        #if 'JWSTDOWNLOAD_OUTDIR' in os.environ and os.environ['JWSTDOWNLOAD_OUTDIR'] != '':
        #    outrootdir = os.environ['JWSTDOWNLOAD_OUTDIR']
        #else:
        #    outrootdir = '.'
        #print('Default output rootdir: %s' % outrootdir)

        parser.add_argument('-n','--skipdownload', action='store_true', default=False, help='skip downloading the data.  (default=%(default)s)')
#        parser.add_argument('-o','--outrootdir', default=outrootdir, help=('output root directory. (default=%(default)s)'))
#        parser.add_argument('--outsubdir', default=None, help=('subdir added to the output root directory. If None, then propID is used (default=%(default)s)'))
        parser.add_argument('--clobber', action='store_true', default=None, help='existing files are overwritten, otherwise they are skipped (default=%(default)s)')

        return(parser)

    def download_product(self, url, outfilename, clobber=False):
        """
        Download the product url and save it to outfilename.

        Parameters
        ----------
        clobber : boolean, optional
            Re-download even if exists.

        Returns
        -------
        dlcode,dlstr:
            0,None: file does not exist
            1,'downloaded': file was succesfully downloaded
            2,'exists':     file already exists, skipped
            3,'ERRORNone':  outfilemame is None
            4,'ERRORempty': outfilemame is empty string
            6,'ERRORdel':   file could not be removed before download
            7,'ERRORdownload': Error during download
            9,'ERRORnotexist': After teh file was downloaded and saved, it doesn't exist
        """
        if outfilename is None:
            print('ERROR: outfilename is None!')
            return(3,'ERRORNone')
        if outfilename=='':
            print('ERROR: outfilename is empty string!')
            return(4,'ERRORempty')

#        if os.path.exists(outfilename):
#            if not clobber and :
#            print(f'WARNING: {outfilename} exists and clobber=False, thus skipping re-downloading it!')
#            return(2,'exists')

        makepath4file(outfilename)

        if rmfile(outfilename,raiseError=False):
            print('ERROR: could not remove old file {outfilename}')
            return(6,'ERRORdel')

        tstart = Time.now()
        if self.verbose: print(f'Dowloading file {url} to {outfilename}')

        try:
            self.JwstObs._download_file(url, outfilename)
        except Exception as e:
            msg = str(e)
            print(f'ERROR: for file {outfilename}: {msg}')
            return(7,'ERRORdownload')

        if not os.path.isfile(outfilename):
            print(f'ERROR: {outfilename} does not exist!!')
            return(9,'ERRORnotexist')

        tend = Time.now()
        if self.verbose: print('time passed for download process: {:.2f} seconds'.format((tend-tstart).to_value('sec')))

        return(1,'downloaded')



    def download_products(self, productTable=None, ix_selected_products=None, clobber=None, ask_confirm_download=True):

        if productTable is None:
           productTable=self.productTable

        if ix_selected_products is None:
            ix_selected_products = self.ix_selected_products

        if clobber is None:
            clobber=self.params['clobber']

#        if 'dl_code' not in self.imoutcols: self.imoutcols.append('dl_code')
#        if 'dl_str' not in self.imoutcols: self.imoutcols.append('dl_str')

        ixs_exist = productTable.ix_equal('dl_code',2,indices = ix_selected_products)

        skip_string = '### None of these files are currently present in the output directory. Downloading all.'
        if not self.params['clobber']:
            ixs_download = AnotB(ix_selected_products,ixs_exist)
            #print(f'\n###############################\n### Downloading {len(ixs_download)} files')
            if len(ixs_exist)>0:
                skip_string = f'### skipping {len(ixs_exist)} files since they already exist'
        else:
            ixs_download = ix_selected_products
            #print(f'\n###############################\n### Downloading {len(ixs_download)} files')
            if len(ixs_exist)>0:
                skip_string = f'### clobbering {len(ixs_exist)} files that already exist'

        if len(ixs_download) == 0:
            if len(ix_selected_products)>0:
                print('###############################\n### Nothing to download!!!')
                print(skip_string)
            else:
                print('###############################\n### No Files Matching Query. Nothing to download.')
            sys.exit(0)
        else:
            print(f'\n###############################\n### Downloading {len(ixs_download)} files')
            print(skip_string)

        print(f'Outdir: {self.outrootdir}/<proposal_id> where <proposal_id> is the value in the proposal_id column')

        if ask_confirm_download and len(ixs_download)>0:
            do_it = input('Do you want to continue and download these products [y/n]?  ')
            if do_it.lower() in ['y','yes']:
                pass
            elif do_it.lower() in ['n','no']:
                print('OK, stopping....')
                sys.exit(0)
            else:
                print(f'Hmm, \'{do_it}\' is neither yes or no. Don\'t know what to do, so stopping ....')
                sys.exit(0)
        else:
            if len(ixs_download)<=0:
                print('############## Nothing to download!!!!')
                return(0)

        counter = 1
        successcounter=0
        failedcounter=0

        mastURL = self.JwstObs._portal_api_connection.MAST_DOWNLOAD_URL + "?uri="

        ixs_download = productTable.ix_sort_by_cols(['proposal_id','obsnum'],indices=ixs_download)

        for ix in ixs_download:
            print(f'\n### Downloading #{counter} out of {len(ixs_download)} files (status: {successcounter} successful, {failedcounter} failed): {os.path.basename(productTable.t.loc[ix,"outfilename"])}')

            url = mastURL + productTable.t.loc[ix,'dataURI']
            outfilename = productTable.t.loc[ix,'outfilename']

            (dl_code,dl_str) = self.download_product(url, outfilename, clobber=clobber)

            productTable.t.loc[ix,'dl_code']=dl_code
            productTable.t.loc[ix,'dl_str']=dl_str
            if dl_code<=2:
                successcounter+=1
            else:
                failedcounter+=1
            counter+=1

        print(f'\n### Download complete: {successcounter} successful, {failedcounter} failed\n###############################')
        return(0)
