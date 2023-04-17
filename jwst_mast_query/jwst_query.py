#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed July 28 11:01:23 2021

@author: arest, bhilbert
"""
from astroquery.mast import Mast, Observations
from astropy.time import Time
from copy import deepcopy
import numpy as np
import pandas as pd
#import astropy.io.fits as fits
#from astropy.nddata import bitmask
import astroquery
import argparse,os,sys,re,copy,shutil
import yaml
from PIL import Image

# pdastroclass is wrapper around pandas.
from jwst_mast_query.pdastro import pdastroclass,unique,AnotB,AorB,AandB,split_commonpath

from jwst_mast_query.utils.constants import MAST_OBS_MODES, PARAM_DEFAULTS

# MAST API documentation:
# https://mast.stsci.edu/api/v0/pyex.html

# MAST CAOM Field Descriptions:
# https://mast.stsci.edu/api/v0/_c_a_o_mfields.html

# File naming schemes:
# https://jwst-pipeline.readthedocs.io/en/latest/jwst/data_products/file_naming.html
# https://jira.stsci.edu/browse/JSDP-1778


if astroquery.__version__<'0.4.2':
    raise RuntimeError("astroquery version=%s, at least 0.4.2 required!" % astroquery.__version__)

# environment variables:

# $MAST_API_TOKEN: needed to access the MAST database. If the token is not set
# with $MAST_API_TOKEN, then it needs to be specified with --token.
# it looks like that for publically available data (like the LRE data), no token is needed...

# $JWST_QUERY_CFGFILE: the config file is optional. It an be specified with
# $JWST_QUERY_CFGFILE or with --configfile

# List of allowed sca (detector) names
ALL_SCAS = ['a1', 'a2', 'a3', 'a4', 'a5', 'along',
            'b1', 'b2', 'b3', 'b4', 'b5', 'blong',
            'guider1','guider2',
            'nrs1','nrs2',
            'mirimage', 'mirifulong', 'mirifushort',
            'nis']


def getlimits(lims):
    if lims is None or len(lims)==0:
        return(None)
    if len(lims)==1:
        if re.search('\+$',lims[0]):
            lims[0]=re.sub('\+$','',lims[0])
            return([lims[0],None])
        elif  re.search('\-$',lims[0]):
            lims[0]=re.sub('\-$','',lims[0])
            return([None,lims[0]])
        else:
            return([lims[0],lims[0]])
    elif len(lims)==2:
        return([lims[0],lims[1]])
    else:
        raise RuntimeError(f'limits can have only 2 entries, more given! {lims}')

def addlink2string(s,link):
    return('<a href="%s">%s</a>' % (link,s))

def imagestring4web(imagename, width=None, height=None):
    imstring = '<img src="%s"' % imagename
    if height != None:
        imstring += f' height={height}'
        #if isinstance(height,int): height = str(height)
        #imstring += 'height=%s' % height
    if width != None:
        imstring += f' width={width}'
        #if isinstance(width,int): width = str(width)
        #imstring += 'width=%s' % width
    imstring +='>'
    return(imstring)

def image4table(imagename, TNimagename=None, width=None, height=None):
    if TNimagename is None:
        TNimagename=imagename
    return(addlink2string(imagestring4web(TNimagename,width=width,height=height),imagename))


def rmfile(filename,raiseError=1,gzip=False):
    " if file exists, remove it "
    if os.path.lexists(filename):
        os.remove(filename)
        if os.path.isfile(filename):
            if raiseError == 1:
                raise RuntimeError('ERROR: Cannot remove %s' % filename)
            else:
                return(1)
    if gzip and os.path.lexists(filename+'.gz'):
        os.remove(filename+'.gz')
        if os.path.isfile(filename+'.gz'):
            if raiseError == 1:
                raise RuntimeError('ERROR: Cannot remove %s' % filename+'.gz')
            else:
                return(2)
    return(0)


class query_mast:
    def __init__(self):
        '''
        Set server and web service endpoints for astroquery
        '''
        self.JwstObs = Mast()
        server = 'https://mast.stsci.edu'
        self.JwstObs._portal_api_connection.MAST_REQUEST_URL = server + "/portal_jwst/Mashup/Mashup.asmx/invoke"
        self.JwstObs._portal_api_connection.MAST_DOWNLOAD_URL = server + "/jwst/api/v0.1/download/file"
        self.JwstObs._portal_api_connection.COLUMNS_CONFIG_URL = server + "/portal_jwst/Mashup/Mashup.asmx/columnsconfig"
        self.JwstObs._portal_api_connection.MAST_BUNDLE_URL = server + "/jwst/api/v0.1/download/bundle"


        #self.RET_COLUMNS = ['proposal_id','dataURL','obsid','obs_id','t_min','t_exptime']
        self.SERVICES = {
                'SI_search': 'Mast.Jwst.Filtered.',
                'Caom_search':'Mast.Caom.Filtered',
                'Product_search':'Mast.Caom.Products'
                }

        self.verbose = 0
        self.debug = 0

        # delete all tables and outdir, set the default parameters
        self.reset()

    def reset(self):
        """
        delete all tables and outdir, set the default parameters into self.params

        Returns
        -------
        None.

        """
        self.outrootdir = None

        # observation and product table
        self.obsTable = pdastroclass()
        self.productTable = pdastroclass()
        # summary table: statistics/info for each propID/obsnum pair
        self.summary = pdastroclass()

        # indices of the selected products
        # They are sorted by self.params['sortcols_productTable'] (can also be set in config file)
        self.ix_selected_products = None

        # indices for obsTable, sorted by self.params['sortcols_obsTable'] (can also be set in config file)
        self.ix_obs_sorted = None
        # indices for summary table, sorted by self.params['sortcols_summaryTable'] (can also be set in config file)
        self.ix_summary_sorted = None

        # self.params will be populated from define_options()
        self.params = deepcopy(PARAM_DEFAULTS)

    def define_options(self, parser=None, usage=None, conflict_handler='resolve'):
        if parser is None:
            parser = argparse.ArgumentParser(usage=usage, conflict_handler=conflict_handler)

        parser.add_argument('-i', '--instrument', type=str, default=PARAM_DEFAULTS['instrument'], choices=['niriss','nircam','nirspec','miri','fgs'],
                            help='Specify the instrument (nircam, nirspec, niriss, miri, fgs)  (default=%(default)s)')
        parser.add_argument('--obsmode', type=str, nargs='+', default=PARAM_DEFAULTS['obsmode'], help=('Optionally specify the observation modes to query for. '
                                                                                                       'If provided, MAST will be queried only for these types of '
                                                                                                       'observations. e.g. assuming instrument is "nircam", '
                                                                                                       'obsmode: ["image", "wfss"] will query MAST for NIRCAM/IMAGE '
                                                                                                       'and NIRCAM/WFSS data. If left empty, the query will include '
                                                                                                       'all modes. For a list of valid modes, see the example config '
                                                                                                       'file in the repository. Note that as of 24 March 2023, MAST '
                                                                                                       'is in the process of reprocessing all data to add the mode '
                                                                                                       'values to the instrument names. Until this process is '
                                                                                                       'complete, some older data may not include the mode name, '
                                                                                                       'and therefore will not be found via a query that includes '
                                                                                                       'the mode name. (default=%(default)s)'))
        parser.add_argument('-v','--verbose', default=PARAM_DEFAULTS['verbose'], action='count')
        parser.add_argument('--propID', type=int, nargs="+", default=PARAM_DEFAULTS['propID'], help=('Specify the JWST proposal number to query for. This can be '
                                                                                                     'a 5 digit integer, or for a smaller number, an integer '
                                                                                                     'with or without prepended zeros. e.g. 1409 or 01409.'))
        parser.add_argument('--obsnums', type=int, nargs='+', default=PARAM_DEFAULTS['obsnums'], help=('Specify the observation numbers within the propID. This '
                                                                                                       'can be a single number, or a bracketed list of numbers. '
                                                                                                       'e.g. 3 or [3, 103] (default=%(default)s)'))
        parser.add_argument('--guidestars', action='store_true', default=PARAM_DEFAULTS['guidestars'], help=('If this is set to True, guidestar products are '
                                                                                                             'also included. Note: there are a **lot** of guide '
                                                                                                             'star products. We recommend you set to True only if '
                                                                                                             'really needed! (default=%(default)s)'))
        parser.add_argument('-f','--filetypes',  type=str, nargs="+", default=PARAM_DEFAULTS['filetypes'], help=('List of file types to select in the product table, '
                                                                                                                 'e.g., _uncal.fits or _uncal.jpg. If no suffix is '
                                                                                                                 'given, .fits is appended. If only letters, then _ '
                                                                                                                 'and .fits are added. For example, "uncal" gets '
                                                                                                                 'expanded to _uncal.fits. Typical image filetypes '
                                                                                                                 'are uncal, rate, rateints, cal. For downloading a '
                                                                                                                 'single file type, the brackets must still surround '
                                                                                                                 'the file suffix, as the script expects a list. A '
                                                                                                                 'relatively complete list of options includes: '
                                                                                                                 '["_segm.fits", "_asn.json", "_pool.csv", "_i2d.jpg", '
                                                                                                                 '"_thumb.jpg", "_cat.ecsv", "_i2d.fits", "_uncal.fits", '
                                                                                                                 '"_uncal.jpg", "_cal.fits", "_trapsfilled.fits", '
                                                                                                                 '"_cal.jpg", "_rate.jpg", "_rateints.jpg", '
                                                                                                                 '"_trapsfilled.jpg", "_rate.fits", "_rateints.fits"] '
                                                                                                                 'See the JWST calibration pipeline documentation for '
                                                                                                                 'a complete list. (default=%(default)s)'))
        parser.add_argument('--calib_levels',  type=int, nargs="+", default=PARAM_DEFAULTS['calib_levels'], help=('Only select products with the specified '
                                                                                                                   'calibration levels (calib_level column in '
                                                                                                                   'productTable) (default=%(default)s)'))
        parser.add_argument('--Nobs_per_batch', type=int, default=PARAM_DEFAULTS['Nobs_per_batch'], help=('For large programs, the query for the products can time '
                                                                                                          'out, and in these cases it is better to split up the query '
                                                                                                          'into Nobs_per_batch observations per batch. For example, '
                                                                                                          'if there are 22 observations, and Nobs_per_batch=4, then '
                                                                                                          'there will be 6 batches, 5 batches with 4 observations, '
                                                                                                          'and the last batch with 2. (default=%(default)s)'))
        parser.add_argument('-c','--configfile', type=str, default=PARAM_DEFAULTS['configfile'], help=('optional config file. default is set to $JWST_QUERY_CFGFILE. '
                                                                                                       'Use -vvv to see the full list of all set parameters.'))
        parser.add_argument('--login', default=PARAM_DEFAULTS['login'], nargs=2, help='username and password for login')
        parser.add_argument('--token', type=str, default=PARAM_DEFAULTS['token'], help=('MAST API token. You can also set the token in the environment '
                                                                                        'variable \$MAST_API_TOKEN'))
        parser.add_argument('--outrootdir', default=PARAM_DEFAULTS['outrootdir'], help=('The base directory for the downloaded products. This can be the name of '
                                                                                        'an environment variable (such as JWSTDOWNLOAD_OUTDIR), or '
                                                                                        'a path. Note the preceding "$" in the case where an environment variable is '
                                                                                        'given. If you include the --outrootdir command line argument when calling '
                                                                                        'jwst_query.py or jwst_download.py, that value will override any value '
                                                                                        'provided in the config file. (default=%(default)s)'))
        parser.add_argument('--outsubdir', default=PARAM_DEFAULTS['outsubdir'], help=('Any additional directory to add to the base directory. This can be used to '
                                                                                      'customize the organization of the downloaded products. (default=%(default)s)'))
        parser.add_argument('--obsnum2outsubdir', default=PARAM_DEFAULTS['obsnum2outsubdir'], help=('If True, the observation number will be used to create a '
                                                                                                    'subdirectory into which the appropriate files will be placed. '
                                                                                                    'e.g. $JWSTDOWNLOAD_OUTDIR/01410/obsnum23/'))
        parser.add_argument('--propIDs_obsnum2outsubdir', default=PARAM_DEFAULTS['propIDs_obsnum2outsubdir'], help=('Specify list of propID for which obsnum2outsubdir '
                                                                                                                    'is True. Only for the listed proposal numbers will '
                                                                                                                    'observation numbers be used to create a '
                                                                                                                    'subdirectory into which the appropriate files will '
                                                                                                                    'be placed. e.g. $JWSTDOWNLOAD_OUTDIR/01410/obsnum23/'))
        parser.add_argument('--skip_propID2outsubdir', action='store_true', default=PARAM_DEFAULTS['skip_propID2outsubdir'], help=('By default, the APT proposal ID '
                                                                                                                                   'is added as a subdir to the output '
                                                                                                                                   'directory. You can skip this with '
                                                                                                                                   'this option (default=%(default)s)'))
        parser.add_argument('--skip_check_if_outfile_exists', action='store_true', default=PARAM_DEFAULTS['skip_check_if_outfile_exists'],
                            help=('By default, the script queries MAST for products, and then checks if each file already exists in the output directory or not '
                                  '("dl_code" and "dl_str" columns). For large numbers of products, this can take time. By setting this option to True, this '
                                  'check can be skipped. (default=%(default)s)'))
        parser.add_argument('--skip_check_filesize', action='store_true', default=PARAM_DEFAULTS['skip_check_filesize'], help=('Don\'t check if output files have the '
                                                                                                                               'correct filesize. This makes it faster '
                                                                                                                               'for large lists, but files might be '
                                                                                                                               'corrupted.'))
        parser.add_argument('-e','--obsid_select', nargs="+", default=PARAM_DEFAULTS['obsid_select'], help=('Specify obsid range applied to "obsID" and "parent_obsid" '
                                                                                                            'columns in the PRODUCT table. If single value, then exact '
                                                                                                            'match. If single value has "+" or "-" at the end, then it '
                                                                                                            'is a lower and upper limit, respectively. If two values, '
                                                                                                            'then range. Examples: 385539+, 385539-, 385539 385600 '
                                                                                                            '(default=%(default)s)'))
        parser.add_argument('-l','--obsid_list', nargs="+", default=PARAM_DEFAULTS['obsid_list'], help=('Specify list of obsid applied to "obsID" and "parent_obsid" '
                                                                                                        'columns in the PRODUCT table. examples: 385539 385600 385530 '
                                                                                                        '(default=%(default)s)'))
        parser.add_argument('--sca', nargs="+", default=PARAM_DEFAULTS['sca'],
                            choices=['a1','a2','a3','a4','a5','along','b1','b2','b3','b4','b5','blong',
                                     'guider1','guider2','nrs1','nrs2','mirimage','mirifulong','mirifushort','nis'],
                                     help='Specify list of sca\'s to select. a5=along, b5=blong')
        parser.add_argument('--mastcolumns_obsTable', nargs="+", default=PARAM_DEFAULTS['mastcolumns_obsTable'], help=('Core columns returned from MAST to the obsTable. '
                                                                                                                       '(default=%(default)s)'))
        parser.add_argument('--outcolumns_obsTable', nargs="+", default=PARAM_DEFAULTS['outcolumns_obsTable'], help=('Output columns for the obsTable. '
                                                                                                                     '(default=%(default)s)'))
        parser.add_argument('--sortcols_obsTable', nargs="+", default=PARAM_DEFAULTS['sortcols_obsTable'], help=('The obsTable is sorted based on these columns. '
                                                                                                                       'The defaults sort the table in the order the '
                                                                                                                       'observations were observed. (default=%(default)s)'))
        parser.add_argument('--outcolumns_productTable', nargs="+", default=PARAM_DEFAULTS['outcolumns_productTable'], help=('List of columns to be shown in product '
                                                                                                                             'table (default=%(default)s)'))
        parser.add_argument('--sortcols_productTable', nargs="+", default=PARAM_DEFAULTS['sortcols_productTable'], help=('The productTable is sorted based on these '
                                                                                                                         'columns. The default sorts the table based '
                                                                                                                         'on calibration level. (default=%(default)s)'))
        parser.add_argument('--sortcols_summaryTable', nargs="+", default=PARAM_DEFAULTS['sortcols_summaryTable'], help=('The summary table is sorted based on these '
                                                                                                                         'columns. The default sorts the table in the '
                                                                                                                         'order the observations were observed. '
                                                                                                                         '(default=%(default)s)'))
        parser.add_argument('--jpg_separate_subdir', default=PARAM_DEFAULTS['jpg_separate_subdir'], help=('If True, downloaded jpgs are saved in separate "jpg" '
                                                                                                          'subdirectories, along side the fits files. (default=%(default)s)'))

        time_group = parser.add_argument_group("Time constraints for the observation/product search")
        time_group.add_argument('-l', '--lookbacktime', type=float, default=PARAM_DEFAULTS['lookbacktime'], help=('Lookback time in days. The script will query MAST '
                                                                                                                  'over a time from the lookback time to the present '
                                                                                                                  'moment. Note that all other time parameters '
                                                                                                                  '(date_select, etc) override the lookback time. '
                                                                                                                  '(default=%(default)s)'))
        time_group.add_argument('-d','--date_select', nargs="+", default=PARAM_DEFAULTS['date_select'], help=('Specify date range (MJD or isot format) applied to '
                                                                                                              '"dateobs_center" column. If single value, then exact '
                                                                                                              'match. If single value has "+" or "-" at the end, then '
                                                                                                              'it is a lower and upper limit, respectively. date_select '
                                                                                                              'will override the lookbacktime. Examples: '
                                                                                                              '58400+, 58400-,2020-11-23+, 2020-11-23 2020-11-25  '
                                                                                                              '(default=%(default)s)'))
        parser.add_argument('-s', '--savetables', type=str, default=PARAM_DEFAULTS['savetables'], help=('Save the tables (selected products, obsTable, summary with '
                                                                                                        'suffix selprod.txt, obs.txt, summary.txt, respectively) with '
                                                                                                        'the specified string as basename. Tables are saved in the '
                                                                                                        'same output directory as the data.  If no string is '
                                                                                                        'provided, the tables are not saved. (default=%(default)s)'))
        parser.add_argument('-w', '--makewebpages', action='store_true', default=PARAM_DEFAULTS['makewebpages'],
                            help=('Make webpages for the products for each propID containing info and images of the retrieved data. Note that this option is not '
                                  'currently in the config file, and must be specified at the command line. (default=%(default)s)'))
        parser.add_argument('--webpage_tablefigsize_width', default=PARAM_DEFAULTS['webpage_tablefigsize_width'], help=("Width in pixels of the resized jpg images "
                                                                                                                        "to be inserted into the index.html summary file."
                                                                                                                        "Recommended: 100-150, or don't specify if "
                                                                                                                        "webpage_mkthumbnails, since then the size of "
                                                                                                                        "the thumbnails are used by default. "
                                                                                                                        "(default=%(default)s)"))
        parser.add_argument('--webpage_tablefigsize_height', default=PARAM_DEFAULTS['webpage_tablefigsize_height'], help=("Height in pixels of the resized jpg images "
                                                                                                                          "to be inserted into the index.html summary "
                                                                                                                          "file. If left undefined, the height will be "
                                                                                                                          "determined from webpage_thumbnail_width and "
                                                                                                                          "the aspect ratio of the original image.  "
                                                                                                                          "(default=%(default)s)"))
        parser.add_argument('--webpage_level12_jpgs', nargs="+", default=PARAM_DEFAULTS['webpage_level12_jpgs'], help=('List of filetypes whose thumbnails will be shown in '
                                                                                                            'index.html. (default=%(default)s)'))
        parser.add_argument('--webpage_fitskeys2table', nargs="+", default=PARAM_DEFAULTS['webpage_fitskeys2table'], help=('List of fits header keywords that should '
                                                                                                                           'be copied to the table. (default=%(default)s)') )
        parser.add_argument('--webpage_cols4table', nargs="+", default=PARAM_DEFAULTS['webpage_cols4table'], help=('Columns to be shown in index.html. (default=%(default)s)'))
        parser.add_argument('--webpage_sortcols', nargs="+", default=PARAM_DEFAULTS['webpage_sortcols'], help=('Columns to sort index.html by. (default=%(default)s)'))
        parser.add_argument('--webpage_mkthumbnails', default=PARAM_DEFAULTS['webpage_mkthumbnails'], help=('If True, a thumbnail jpg is created for each of the '
                                                                                                            'jpg products listed in webpage_level12_jpgs. '
                                                                                                            '(default=%(default)s)'))
        parser.add_argument('--webpage_thumbnails_overwrite', default=PARAM_DEFAULTS['webpage_thumbnails_overwrite'], help=('If True, remake thumbnails even if they '
                                                                                                                            'already exist (default=%(default)s)'))
        parser.add_argument('--webpage_thumbnails_width', default=PARAM_DEFAULTS['webpage_thumbnails_width'], help=('Width in pixels of the resized jpg images '
                                                                                                                    'to be inserted into the index.html summary '
                                                                                                                    'file. (default=%(default)s)'))
        parser.add_argument('--webpage_thumbnails_height', default=PARAM_DEFAULTS['webpage_thumbnails_height'], help=('Height in pixels of the resized jpg images '
                                                                                                                      'to be inserted into the index.html summary '
                                                                                                                      'file. (default=%(default)s)'))
        parser.add_argument('--skipdownload', action='store_true', default=PARAM_DEFAULTS['skipdownload'], help='If set, no files will be downloaded. (default=%(default)s)')

        return(parser)

    def get_arguments(self, args, configfile=None):
        '''

        Parameters
        ----------
        args : list
            pass the command line arguments to self.params.
            makes sure that the propID has the correct format (5 digit string)
        configfile : string, optional
            Config filename. The default is None. If None, then
            $JWST_QUERY_CONFIGFILE is used if exists.

        Returns
        -------
        None.

        '''

        def subenvvarplaceholder(paramsdict):
            """ Loop through all string parameters and substitute environment variables. environment variables have the form $XYZ """
            envvarpattern=re.compile('\$(\w+)')

            for param in paramsdict:
                if isinstance(paramsdict[param], str):
                    envvarnames=envvarpattern.findall(paramsdict[param])
                    if envvarnames:
                        for name in envvarnames:
                            if not (name in os.environ):
                                if param=='outrootdir':
                                    # this is a special hack: if outrootdir is not defined, it is set to '.'
                                    print(f'WARNING! environment variable JWSTDOWNLOAD_OUTDIR in {param} not defined! Thus setting {param} to None...')
                                    paramsdict[param] = None
                                else:
                                    raise RuntimeError("environment variable %s used in config file, but not set!" % name)
                            else:
                                envval=os.environ[name]
                                subpattern='\$%s' % (name)
                                paramsdict[param] = re.sub(subpattern,envval,paramsdict[param])
                elif isinstance(paramsdict[param], dict):
                #elif type(dict[param]) is types.DictType:
                    # recursive: sub environment variables down the dictionary
                    subenvvarplaceholder(paramsdict[param])
            return(0)

        # get the parameters from the config file
        if args.configfile is not None:
            if not os.path.isfile(args.configfile):
                raise RuntimeError('config file %s does not exist!' % (args.configfile))
            print(f'Loading config file {args.configfile}')
            cfgparams = yaml.load(open(args.configfile,'r'), Loader=yaml.FullLoader)
            self.params.update(cfgparams)
            self.params['configfile'] = args.configfile

            if args.verbose>2:
                print('\n### CONFIG FILE PARAMETERS:')
                for p in cfgparams:
                    print('config file args: setting %s to' % (p),cfgparams[p])

            # Make sure propID and obsnums, if they have been entered in the config file,
            # are lists, rather than integers. When entered via command line, they are
            # guaranteed to be lists.
            list_keys = ['propID', 'obsnums']
            for key in list_keys:
                if key in self.params.keys():
                    if not isinstance(self.params[key], list) and self.params[key] is not None:
                        self.params[key] = [self.params[key]]

            # Ensure that the propID is a 5-digit string
            if 'propID' in self.params.keys():
                if self.params['propID'] is not None:
                    self.params['propID'] = ['%05d' % (int(self.params['propID'][0]))]

        subenvvarplaceholder(self.params)

        # Make sure obsmode is a list, even if it's just [None]
        list_keys = ['obsmode']
        for key in list_keys:
            if key in self.params.keys():
                if not isinstance(self.params[key], list):
                    self.params[key] = [self.params[key]]

        # Go through optional parameters.
        # 'None' does not overwrite previously set parameters (either default or config file)
        if args.verbose>2:
            print('\n### OPTIONAL COMMAND LINE PARAMETERS:')
        argsdict = vars(args)
        for arg in argsdict:

            # skip config file
            if arg=='configfile': continue

            if argsdict[arg] is not None and argsdict[arg] != []:
                if args.verbose>2:
                    print('optional args: setting %s to %s' % (arg,argsdict[arg]))
                self.params[arg]=argsdict[arg]
            else:
                if arg not in  self.params:
                    self.params[arg]=None

        if 'obsnums' not in self.params.keys():
            self.params['obsnums']=None

        # propID
        if self.params['propID'] is not None:
            if len(self.params['propID'])>1:
                self.params['obsnums'] = self.params['propID'][1:]
            self.params['propID'] = '%05d' % (int(self.params['propID'][0]))
        self.verbose = self.params['verbose']

        if self.verbose>2:
            print('\n### FINAL PARAMETERS:')
            for p in self.params:
                print(f'{p}: {self.params[p]}')

        print('INSTRUMENT: ', self.params['instrument'])
        print('obsmode: ', self.params['obsmode'])
        print('propID: ',self.params['propID'])
        print('obsnums: ',self.params['obsnums'])

        return(0)

    def login(self,token=None,raiseErrorFlag=False):
        if token is None:
            token = self.params['token']

        if token is None:
            if raiseErrorFlag:
                raise RuntimeError("No token!! Cannot login! you can set the token with \$API_MAST_TOKEN, as 'token' in the config file, or with --token ")
            else:
                print("Warning: no token given, could not login, only public requests possible")
                return(1)
        else:
            self.JwstObs.login(token=token, store_token=False)
        return(0)

    def obsid_select(self,obsid_select,
                     productTable=None,
                     ix_selected_products=None):

        if productTable is None:
           productTable=self.productTable

        if ix_selected_products is None:
            ix_selected_products = self.ix_selected_products

        # parse trailing '+' and '-', and get limits
        limits = getlimits(obsid_select)
        if limits is None:
            return(ix_selected_products)

        for i in range(len(limits)):
            if limits[i] is not None: limits[i]=int(limits[i])

        ixs = self.productTable.ix_not_null('obsID',indices=ix_selected_products)
        ixs_keep1 = self.productTable.ix_inrange('obsID',limits[0],limits[1],indices=ixs)
        ixs_keep2 = self.productTable.ix_inrange('parent_obsid',limits[0],limits[1],indices=ixs)
        ixs_keep = unique(AorB(ixs_keep1,ixs_keep2))
        print(f'obsid cut {limits[0]} - {limits[1]}: keeping {len(ixs_keep)} from {len(ixs)}')
        return(ixs_keep)

    def obsid_list(self,obsid_list,
                     productTable=None,
                     ix_selected_products=None):

        if productTable is None:
           productTable=self.productTable

        if ix_selected_products is None:
            ix_selected_products = self.ix_selected_products

        if obsid_list is None or len(obsid_list)==0:
            return(ix_selected_products)

        ixs_keep = []
        for obsid in obsid_list:
            obsid = int(obsid)
            ixs2add = self.productTable.ix_equal('obsID',obsid,indices=ix_selected_products)
            ixs_keep.extend(ixs2add)
            ixs2add = self.productTable.ix_equal('parent_obsid',obsid,indices=ix_selected_products)
            ixs_keep.extend(ixs2add)
        ixs_keep = unique(ixs_keep)
        print(f'obsid list cut: keeping {len(ixs_keep)} from {len(ix_selected_products)}')
        return(ixs_keep)


    def select_sca(self, sca_list,
                   productTable=None,
                   ix_selected_products=None):

        if productTable is None:
           productTable=self.productTable

        if ix_selected_products is None:
            ix_selected_products = self.ix_selected_products

        if sca_list is None or len(sca_list)==0:
            return(ix_selected_products)

        if isinstance(sca_list,str):
            sca_list=[sca_list]

        ixs_keep = []
        for sca in sca_list:
            if sca=='a5': sca='along'
            if sca=='b5': sca='blong'
            ixs2add = self.productTable.ix_equal('sca',sca,indices=ix_selected_products)
            ixs_keep.extend(ixs2add)
        ixs_keep = unique(ixs_keep)
        print(f'select sca {sca_list}: keeping {len(ixs_keep)} from {len(ix_selected_products)}')
        return(ixs_keep)

    def select_obsnums(self, obsnum_list,
                       productTable=None,
                       ix_selected_products=None):

        if productTable is None:
           productTable=self.productTable

        if ix_selected_products is None:
            ix_selected_products = self.ix_selected_products

        if obsnum_list is None or len(obsnum_list)==0:
            return(ix_selected_products)

        if isinstance(obsnum_list,str):
            obsnum_list=[obsnum_list]

        ixs_keep = []
        for obsnum in obsnum_list:
            ixs2add = self.productTable.ix_equal('obsnum',obsnum,indices=ix_selected_products)
            ixs_keep.extend(ixs2add)
        ixs_keep = unique(ixs_keep)
        print(f'select obsnum {obsnum_list}: keeping {len(ixs_keep)} from {len(ix_selected_products)}')
        return(ixs_keep)


    def get_mjd_limits(self, lookbacktime=None, date_select=None):
        """Determine the earliest and latest dates to use for the query

        Parameters
        ----------
        lookbacktime : int
            Number of days before the current date to search

        date_select : str
            String of starting and stopping dates to use

        Returns
        -------
        mjd_min : float
            MJD of the earliest time to search
        mjd_max : float
            MJD of the latest time to search
        """
        if lookbacktime is None: lookbacktime = self.params['lookbacktime']

        if date_select is None: date_select = self.params['date_select']

        # If date_select is a string, split it into a list.
        # When provided at the command line, date_select is a list.
        # When provided in cfg file, date_select is a string.
        if isinstance(date_select, str):
            date_select = date_select.split(' ')

        mjd_min = mjd_max = None
        # parse trailing '+' and '-', and get limits
        limits = getlimits(date_select)
        if limits is not None:
            # Convert dates into MJD if necessary
            for i in (0,1):
                if limits[i] is not None:
                    try:
                        limits[i] = float(limits[i])
                    except:
                        limits[i]= Time(limits[i], format='isot').to_value('mjd')
            mjd_min = limits[0]
            mjd_max = limits[1]
        else:
            if (mjd_min is None):
                if self.params['lookbacktime']>0.0:
                    if self.verbose>1: print('Note: Looking back %.1f days' % (self.params['lookbacktime']))
                    mjd_min = Time.now().mjd - self.params['lookbacktime']

        if mjd_max is None:
            mjd_max = Time.now().mjd+0.1

        return(mjd_min, mjd_max)

    def observation_query(self, propID=None, instrument=None, mjd_min=None, mjd_max=None, token=None):
        '''
        Perform query for observations matching JWST instrument and program ID
        that began after a particular date.

        Parameters
        ----------
        propID : int
            Proposal ID number
        instrument : str
            Instrument name, e.g. 'nircam'
        mjd_min : float
            MJD of the earliest date to search
        mjd_max : float
            MJD of the latest date to search

        Returns
        -------
        self.obsTable
        '''

        if instrument is None:
            instrument = self.params['instrument']

        print('MJD range:',mjd_min, mjd_max)

        # Work-around for MAST update, where instrument name being changed to
        # a string of the format "<instrument>/<mode>" rather than "<instrument>"
        # e.g "MIRI/IMAGE" rather than "MIRI". Wildcards are not yet supported
        # by astroquery when working with the JwstFiltered data, so we need to
        # get a list of all possible modes associated with the proposal, and then
        # query for each separately

        # To get the list of modes, we can query using Observations, which does support
        # wildcard characters, but returnes only limited information. In this case, we
        # care only about the instrument_name values in the results.

        # Check to see if the user entered a list of observing modes. If not, we
        # query all modes for the instrument via wildcard
        if self.params['obsmode'][0] is None or len(self.params['obsmode']) == 0:
            print(f'No obsmode given. Querying for all files for {instrument}.')
            inst_list = list(set(Observations.query_criteria(instrument_name=f'{instrument}*',
                                                             t_min=[mjd_min, mjd_max])['instrument_name']))
        else:
            initial_inst_list = [f'{instrument.upper()}/{mode.upper()}' for mode in self.params['obsmode']]
            inst_list = deepcopy(initial_inst_list)

            # Check against the set of allowed modes, and alert the user to any invalid entries
            # The MAST search will run without crashing if there is a bad entry. It will simply
            # ignore the entry.
            valid_modes = [mode.split('/')[1].lower() for mode in MAST_OBS_MODES[instrument]]
            for mode in initial_inst_list:
                if mode not in MAST_OBS_MODES[instrument]:
                    mode_only = mode.split('/')[1].lower()
                    print(f'\n\nWARNING: "{mode_only}" is not a valid instrument mode for the MAST search.')
                    print('This entry will be ignored in the query.')
                    print(f'Valid {instrument} modes are: {valid_modes}\n\n')
                    inst_list.remove(mode)

            # If all of the user-input obsmodes are invalid, then we fall back to using the wildcard
            if len(inst_list) == 0:
                print('All entered modes are invalid. Falling back to query across all possible obsmodes.')
                inst_list = list(set(Observations.query_criteria(instrument_name=f'{instrument}*',
                                                             t_min=[mjd_min, mjd_max])['instrument_name']))

        # If inst_list is empty, then there are no observations to query for
        if len(inst_list) == 0:
            return(0)

        columns = ','.join(self.params['mastcolumns_obsTable'])
        service = self.SERVICES['Caom_search']
        params = {"columns":columns,
                  "filters":[
                  {"paramName":"obs_collection","values":["JWST"]},
                  {"paramName":"instrument_name","values":inst_list},
                  {"paramName":"t_min",
                    "values":[{"min":mjd_min,"max":mjd_max}]},
                  ]}

        # Only add propID entry if not None. If it is set to None it doesn't work!
        if propID is not None:
            params["filters"].append({"paramName":"proposal_id","values":[propID]})

        if self.verbose:
            print('\n#### Querying Obstable....')
        if self.verbose>1:
            print('query params:',params)

        if token is not None:
            self.obsTable.t = self.JwstObs.service_request(service, params, mast_token=token).to_pandas()
        else:
            self.obsTable.t = self.JwstObs.service_request(service, params).to_pandas()
        if self.verbose:
            print('Obstable obtained: length:',len(self.obsTable.t))

        if len(self.obsTable.t)==0:
            print('WARNING!! No observations found!!')
            return(0)

        # Find the obsnum # from the filename if possible.
        ixs = self.obsTable.getindices()
        obsnumsearch = re.compile('^jw\d{5}(\d{3})\d{3}\_')
        self.obsTable.t['obsnum']=pd.NA
        for ix in ixs:
            m = obsnumsearch.search(self.obsTable.t.loc[ix,'obs_id'])
            if m is not None:
                self.obsTable.t.loc[ix,'obsnum']=int(m.groups()[0])
            else:
                self.obsTable.t.loc[ix,'obsnum']=pd.NA

        self.obsTable.mjd2dateobs('t_min','date_min')

        self.ix_obs_sorted = self.obsTable.ix_sort_by_cols(self.params['sortcols_obsTable'])

        # make proposal_id integer
        if 'proposal_id' in self.obsTable.t.columns:
            self.obsTable.t['proposal_id']=self.obsTable.t['proposal_id'].astype('int')
        if 'obsid' in self.obsTable.t.columns:
            self.obsTable.t['obsid']=self.obsTable.t['obsid'].astype('int')

        self.obsTable.t.style.set_properties(subset=["obs_id"], **{'text-align': 'left'})
        #self.obsTable.default_formatters['obs_id']='{:<}'.format


        if self.verbose>1: print('Obstable columns:',self.obsTable.t.columns)

        return self.obsTable

    def fix_obsnum(self, obsTable=None, productTable=None):
        """
        fix obsnum=NA in productTable and obsTable

        Parameters
        ----------
        obsTable : pdastroclass, optional
            Table with observations. The default is None. If None, then self.obsTable is used
        productTable : pdastroclass, optional
            Table with products. The default is None. If None, then self.productTable is used
        Returns
        -------
        None.

        """
        if obsTable is None:
            obsTable=self.obsTable

        if productTable is None:
           productTable=self.productTable

        self.obsTable.t['obsnum']=pd.Series(self.obsTable.t['obsnum'],dtype=pd.Int32Dtype())

        ixs_obsTable = obsTable.getindices()
        for ix_obsTable in ixs_obsTable:
            # find all products for a given observation
            obsid = obsTable.t.loc[ix_obsTable,'obsid']
            if self.verbose>2:
                print('### obsid',obsid)
            ixs_prodTable = productTable.ix_equal('parent_obsid',obsid)
            if self.verbose>2:
                productTable.write(indices=ixs_prodTable)

            # If the obsnum cannot be determined from the obsTable 'obs_id', then get it from the products!
            if obsTable.t.loc[ix_obsTable,'obsnum'] is pd.NA:
                # remove None entries
                ixs = productTable.ix_not_null('obsnum',indices=ixs_prodTable)
                if len(ixs)>0:
                    obsnum = unique(productTable.t.loc[ixs,'obsnum'])
                    if len(obsnum)==0:
                        obsTable.write(indices=ixs_obsTable)
                        productTable.write(indices=ixs_prodTable)
                        raise RuntimeError('No obsnum for these observations and products!')
                    """
                    if len(obsnum)>1:
                        print('\n BUUUGGGGG!!!!')
                        obsTable.write(indices=[ix_obsTable])
                        productTable.write(indices=ixs)
                        print('More than one obsnum!',obsnum)
                        raise RuntimeError('BUGGG!!!!!')
                    """
                    if len(obsnum)==1:
                        obsTable.t.loc[ix_obsTable,'obsnum'] = obsnum
                    else:
                        if self.debug:
                            obsTable.write(indices=[ix_obsTable])
                            productTable.write(indices=ixs)
                            print('More than one obsnum!',obsnum)
                            raise RuntimeError('BUGGG!!!!!')

                else:
                    print('Warning: could not determine obsnum for observation %s' % (obsTable.t.loc[ix_obsTable,'obsid']))

            # If there are entries in productTable for which 'obsnum' is None (happens to some calib_level=3 products) fill them up
            ixs_null = productTable.ix_is_null('obsnum',indices=ixs_prodTable)
            if len(ixs_null)>0:
                productTable.t.loc[ixs_null,'obsnum'] = obsTable.t.loc[ix_obsTable,'obsnum']


    def product_query_bkp(self, obsTable=None, guidestars=None):
        '''
        Perform query for data products based on obs_id's in observation table
        '''

        if obsTable is None:
            obsTable=self.obsTable

        if guidestars is None:
            guidestars = self.params['guidestars']
            if guidestars is None:
                guidestars = False


        # query MAST for all products for the obsid's
        obsids = ','.join(obsTable.t['obsid'].astype('str'))
        service = self.SERVICES['Product_search']
        params = {"obsid":obsids,
                  "columns":['type','productType'],
                  "format":"json"
                  }
        if self.verbose:
            print('\n#### Querying ProductTable....')
        self.productTable.t = self.JwstObs.service_request(service, params).to_pandas()
        if self.verbose:
            print('productTable obtained: length:',len(self.productTable.t))
            #print('productTable columns:',self.productTable.t.columns)

        # remove guide stars if wanted...
        if not guidestars:
            ixs_products = self.productTable.getindices()
            gs_text = '_gs-'
            ixs_gs = self.productTable.ix_matchregex('productFilename',gs_text)
            ixs_keep_products = AnotB(ixs_products,ixs_gs)
            self.productTable.t=self.productTable.t.loc[ixs_keep_products].reset_index()
            if self.verbose:
                print('Removing %d guide star products from a total of %d products, %d left' % (len(ixs_gs),len(ixs_products),len(ixs_keep_products)))

        # Fill the suffix column with the suffix of the form _bla1.bla2, e.g. _uncal.fits
        # This will later be used to figure out
        self.productTable.t['filetype'] = self.productTable.t['productFilename'].str.extract(r'(\_[a-zA-Z0-9]+\.[a-zA-Z0-9]+)$')

        scas = list(np.repeat(np.nan, len(self.productTable.t['obs_id'])))
        for i, oid in enumerate(self.productTable.t['obs_id']):
            match = [word for word in ALL_SCAS if word in oid]
            if len(match) > 0:
                scas[i] = match[0]
        self.productTable.t['sca'] = scas

        # Find the obsnum # from the filename if possible.
        ixs = self.productTable.getindices()
        obsnumsearch = re.compile('^jw\d{5}(\d{3})\d{3}\_')
        self.productTable.t['obsnum']=pd.NA
        for ix in ixs:
            m = obsnumsearch.search(self.productTable.t.loc[ix,'obs_id'])
            if m is not None:
                self.productTable.t.loc[ix,'obsnum']=int(m.groups()[0])
            else:
                self.productTable.t.loc[ix,'obsnum']=pd.NA

        # make proposal_id integer
        self.productTable.t['proposal_id']=self.productTable.t['proposal_id'].astype('int')
        # make obsid integer
        if 'obsID' in self.productTable.t.columns:
            self.productTable.t['obsID']=self.productTable.t['obsID'].astype('int')
        if 'parent_obsid' in self.productTable.t.columns:
            self.productTable.t['parent_obsid']=self.productTable.t['parent_obsid'].astype('int')

        self.fix_obsnum(obsTable=obsTable, productTable=self.productTable)

        if self.verbose: print('productTable columns:',self.productTable.t.columns)
        if self.verbose:
            allfiletypes = unique(self.productTable.t['filetype'])
            print('List of all filetypes of obtained products:',allfiletypes)

        return self.productTable

    def product_query(self, obsTable=None, guidestars=None, Nobs_per_batch=2):
        '''
        Perform query for data products based on obs_id's in observation table
        '''

        if obsTable is None:
            obsTable=self.obsTable

        if guidestars is None:
            guidestars = self.params['guidestars']
            if guidestars is None:
                guidestars = False


        # query MAST for all products for the obsid's
        obsids = ','.join(obsTable.t['obsid'].astype('str'))
        service = self.SERVICES['Product_search']
        #params = {"obsid":obsids,
        #          "columns":['type','productType'],
        #          "format":"json"
        #          }
        if Nobs_per_batch is not None:
            if self.verbose:
                print(f'\n#### Querying ProductTable in Nobs={Nobs_per_batch} per batch ....')

            obsid_list = list(obsTable.t['obsid'].astype('str'))
            obsids_chunks = [obsid_list[i:i+Nobs_per_batch] for i in range(0, len(obsid_list), Nobs_per_batch)]
            tmptables = []
            if self.verbose:
                print(f'{len(obsids_chunks)} batches')
            for obsids_chunk in obsids_chunks:
                obsids = ','.join(obsids_chunk)
                if self.verbose:
                    print(f'Querying products for obsids={obsids}')
                params = {"obsid":obsids,
                          "columns":['type','productType'],
                          "format":"json"
                          }
                tmptable = self.JwstObs.service_request(service, params).to_pandas()

                # Remove duplicate rows based on the productFilename value
                tmptable.drop_duplicates(subset='productFilename', keep='first', inplace=True, ignore_index=True)

                if self.verbose:
                    print(f'After filtering duplicate entries, {len(tmptable)} products found for this chunk of observations')
                tmptables.append(tmptable)

            self.productTable.t = pd.concat(tmptables, axis=0, ignore_index=True)

            # Remove duplicate rows that come from different queries
            self.productTable.t.drop_duplicates(subset='productFilename', keep='first', inplace=True, ignore_index=True)
        else:
            params = {"obsid":obsids,
                      "columns":['type','productType'],
                      "format":"json"
                      }
            if self.verbose:
                print('\n#### Querying ProductTable....')
            self.productTable.t = self.JwstObs.service_request(service, params).to_pandas()
        #self.productTable.write()
        if self.verbose:
            print(f'productTable obtained with {len(self.productTable.t)} entries')

        # remove guide stars if wanted...
        if not guidestars:
            ixs_products = self.productTable.getindices()
            gs_text = '_gs-'
            ixs_gs = self.productTable.ix_matchregex('productFilename',gs_text)
            ixs_keep_products = AnotB(ixs_products,ixs_gs)
            self.productTable.t=self.productTable.t.loc[ixs_keep_products].reset_index()
            if self.verbose:
                print('Removing %d guide star products from a total of %d products, %d left' % (len(ixs_gs),len(ixs_products),len(ixs_keep_products)))

        # Fill the suffix column with the suffix of the form _bla1.bla2, e.g. _uncal.fits
        # This will later be used to figure out
        self.productTable.t['filetype'] = self.productTable.t['productFilename'].str.extract(r'(\_[a-zA-Z0-9]+\.[a-zA-Z0-9]+)$')

        scas = list(np.repeat(np.nan, len(self.productTable.t['obs_id'])))
        for i, oid in enumerate(self.productTable.t['obs_id']):
            match = [word for word in ALL_SCAS if word in oid]
            if len(match) > 0:
                scas[i] = match[0]
        self.productTable.t['sca'] = scas

        #bla = self.productTable.t['obs_id'].str.extract(r'^jw\d{8}(\d{3})_(\d{2})')
        #print(bla,bla.columns,bla[0])
        # get the visit number
        self.productTable.t['visit']= self.productTable.t['obs_id'].str.extract(r'^jw\d{8}(\d{3})_\d{2}\d{1}\d{2}_')
        #self.productTable.t['vstgrp']= self.productTable.t['obs_id'].str.extract(r'^jw\d{8}\d{3}_(\d{2})\d{1}\d{2}_')

        # mark all products that are fits or jpg in separate columns
        self.productTable.t['fits']=self.productTable.t['jpg']=False
        ixs_fits = self.productTable.ix_matchregex('filetype',r'\.fits$')
        self.productTable.t.loc[ixs_fits,'fits']=True
        ixs_fits = self.productTable.ix_matchregex('filetype',r'\.jpg$')
        self.productTable.t.loc[ixs_fits,'jpg']=True

        # Find the obsnum # from the filename if possible.
        ixs = self.productTable.getindices()
        obsnumsearch = re.compile('^jw\d{5}(\d{3})\d{3}\_')
        self.productTable.t['obsnum']=pd.NA
        for ix in ixs:
            m = obsnumsearch.search(self.productTable.t.loc[ix,'obs_id'])
            if m is not None:
                self.productTable.t.loc[ix,'obsnum']=int(m.groups()[0])
            else:
                self.productTable.t.loc[ix,'obsnum']=pd.NA

        # make proposal_id integer
        self.productTable.t['proposal_id']=self.productTable.t['proposal_id'].astype('int')
        # make obsid integer
        if 'obsID' in self.productTable.t.columns:
            self.productTable.t['obsID']=self.productTable.t['obsID'].astype('int')
        if 'parent_obsid' in self.productTable.t.columns:
            self.productTable.t['parent_obsid']=self.productTable.t['parent_obsid'].astype('int')

        self.fix_obsnum(obsTable=obsTable, productTable=self.productTable)

        if self.verbose: print('productTable columns:',self.productTable.t.columns)
        if self.verbose:
            allfiletypes = unique(self.productTable.t['filetype'])
            print('List of all filetypes of obtained products:',allfiletypes)

        return self.productTable


    def product_filter(self, productTable=None, filetypes=None, calib_levels=None, gs_omit=True):
        '''
        Filter the list of products based on the filetype (or better suffix)
        '''
        if productTable is None:
           productTable=self.productTable

        if filetypes is None:
            filetypes=self.params['filetypes']

        if calib_levels is None:
            calib_levels=self.params['calib_levels']

        if self.verbose:
            print('\n#### Selecting products....')

        # if necessary, add leading '_' and suffix '.fits'
        if filetypes is not None:
            for i in range(len(filetypes)):
                if filetypes[i] in ['fits','jpg']:
                    filetypes[i]+='$'
                    continue
                if re.search('^_',filetypes[i]) is None:
                    filetypes[i] = '_'+filetypes[i]
                if re.search('\.',filetypes[i]) is None:
                    filetypes[i] += '.fits'

            print('allowed filetype list:',filetypes)

        ix_products = self.productTable.getindices()
        # remove guide stars if wanted...
        if gs_omit:
            gs_text = '_gs-'
            ix_gs = self.productTable.ix_matchregex('productFilename',gs_text)
            if self.verbose>1:
                print('Removing %d guide star products from a total of %d products, %d left' % (len(ix_gs),len(ix_products),len(ix_products)-len(ix_gs)))
            ix_products = AnotB(ix_products,ix_gs)


        if calib_levels is not None:
            ix_calib_level = []
            for calib_level in calib_levels:
                ixs = self.productTable.ix_equal('calib_level',calib_level,indices=ix_products)
                if self.verbose>1:
                    print('Keeping %d products with calib_level %d' % (len(ixs),calib_level))
                ix_calib_level.extend(ixs)
            if self.verbose>1:
                print('Keeping %d out of %d products with the correct calib_levels' % (len(ix_calib_level),len(ix_products)))
            ix_products = ix_calib_level

        # Loop trough the filetypes and get all entries from ix_products list
        self.ix_selected_products = []
        if filetypes is not None:
            for filetype in filetypes:

                # make sure the '.' in the regular expression is literal, and also add '$' to the end
                regex=re.sub('\.','\.',filetype)+'$'

                # get the indices for matching filetype ...
                ix_matching_filetype = self.productTable.ix_matchregex('filetype',regex,indices=ix_products)
                # ... and add them to the list of good indices
                self.ix_selected_products.extend(ix_matching_filetype)

                if self.verbose>1:
                    print('%d products with filetype regex matching %s' % (len(ix_matching_filetype),regex))
            if self.verbose:
                print('%d products with correct filetypes left' % (len(self.ix_selected_products)))
        else:
            self.ix_selected_products = ix_products

        # get the list of valid filetypes that exist
        self.params['filetypes'] = unique(self.productTable.t.loc[self.ix_selected_products,'filetype'])

        self.ix_selected_products = self.productTable.ix_sort_by_cols(self.params['sortcols_productTable'],indices=self.ix_selected_products)


        return(self.ix_selected_products)

    def update_obstable_indices(self, ixs_prod, ix_obs_sorted=None, obsTable=None, productTable=None):
        if obsTable is None:
            obsTable=self.obsTable

        if productTable is None:
           productTable=self.productTable

        if ix_obs_sorted is None:
            ix_obs_sorted = self.ix_obs_sorted

        # get all parent_obsids of the selected products
        obsids = unique(productTable.t.loc[ixs_prod,'parent_obsid'])

        # get all keys in observation table that are associated with these obsids
        new_ix_obs_sorted = []
        for obsid in obsids:
            new_ix_obs_sorted.extend(self.obsTable.ix_equal('obsid',obsid,indices=ix_obs_sorted))

        # Keep the order from ix_obs_sorted!
        new_ix_obs_sorted = AandB(ix_obs_sorted,new_ix_obs_sorted)

        return(new_ix_obs_sorted)

    def get_outdir(self,productTable, ix,
                   outdir=None,
                   skip_propID2outsubdir=False,
                   obsnum2outsubdir=False,
                   jpg_separate_subdir=False):

        if outdir is None:
            outdir = self.outrootdir
        if outdir is None:
            raise RuntimeError("outdir is not defined!")
        outdir = outdir.rstrip("/")


        if not skip_propID2outsubdir:
            outdir += f'/{productTable.t.loc[ix,"proposal_id"]:05d}'
            #outdir += '/{:05d}'.format(productTable.t.loc[ix,"proposal_id"])

        # shameless hack: for NIRCam
        if obsnum2outsubdir or (int(productTable.t.loc[ix,"proposal_id"]) in self.params["propIDs_obsnum2outsubdir"]):
            if productTable.t.loc[ix,"obsnum"] is pd.NA:
                outdir += '/NA'
            else:
                outdir += f'/obsnum{productTable.t.loc[ix,"obsnum"]:02d}'

        if jpg_separate_subdir and productTable.t.loc[ix,"jpg"]:
            outdir += '/jpg'

        return(outdir)


    def mk_outfilename(self, productTable, ix,
                       outdir=None,
                       skip_propID2outsubdir=False,
                       obsnum2outsubdir=False,
                       jpg_separate_subdir=False,
                       skip_check_if_outfile_exists=False,
                       skip_check_filesize=False):

        # get the output directory
        outdir = self.get_outdir(productTable, ix, outdir=outdir,
                                 skip_propID2outsubdir=skip_propID2outsubdir,
                                 obsnum2outsubdir=obsnum2outsubdir,
                                 jpg_separate_subdir=jpg_separate_subdir)


        #if outdir is None:
        #    outdir = self.outrootdir
        #if outdir is None:
        #    raise RuntimeError("outdir is not defined!")
        #outdir.rstrip("/")


        #if not skip_propID2outsubdir:
        #    outdir += '/{:05d}'.format(productTable.t.loc[ix,"proposal_id"])

        outfilename = productTable.t.loc[ix,'productFilename']

        fullfilename = os.path.abspath(f'{outdir}/{outfilename}')
        productTable.t.loc[ix,'outfilename']=fullfilename

        if not skip_check_if_outfile_exists:
            if os.path.exists(fullfilename):
                if not skip_check_filesize:
#                    print('VVV',productTable.t.loc[ix,'outfilename'],productTable.t.loc[ix,'size'],os.path.getsize(productTable.t.loc[ix,'outfilename']))
                    outfile_filesize = os.path.getsize(productTable.t.loc[ix,'outfilename'])
                    if outfile_filesize != int(productTable.t.loc[ix,'size']):
                        productTable.t.loc[ix,'dl_code'] = 4
                        productTable.t.loc[ix,'dl_str'] = 'ERRfilesize'
                    else:
                        productTable.t.loc[ix,'dl_code'] = 2
                        productTable.t.loc[ix,'dl_str'] = 'exists'
                else:
                    productTable.t.loc[ix,'dl_code'] = 2
                    productTable.t.loc[ix,'dl_str'] = 'exists'

            else:
                productTable.t.loc[ix,'dl_code'] = 0
                productTable.t.loc[ix,'dl_str'] = None


    def set_outrootdir(self,outrootdir=None,outsubdir=None):
        self.outrootdir = outrootdir
        if self.outrootdir is None:
            self.outrootdir = self.params['outrootdir']
        # if outdir is not defined, use '.'
        if self.outrootdir is None or self.outrootdir == '':
            self.outrootdir = '.'

        if outsubdir is not None and outsubdir!='':
            self.outrootdir  += f'/{outsubdir}'
        elif self.params['outsubdir'] is not None and self.params['outsubdir']!='':
            self.outrootdir  += f'/{self.params["outsubdir"]}'

        self.outrootdir = os.path.abspath(self.outrootdir)
        return(self.outrootdir)


    def mk_outfilenames(self,
                        productTable=None,
                        ix_selected_products=None,
                        outdir=None,
                        skip_propID2outsubdir=False,
                        obsnum2outsubdir=False,
                        jpg_separate_subdir=False,
                        skip_check_if_outfile_exists=False,
                        skip_check_filesize=False):

        if productTable is None:
           productTable=self.productTable

        if ix_selected_products is None:
            ix_selected_products = self.ix_selected_products



        productTable.t.loc[ix_selected_products,'outfilename'] = None
        productTable.t.loc[ix_selected_products,'dl_code'] = pd.NA
        productTable.t.loc[ix_selected_products,'dl_str'] = None
        productTable.t['dl_code'] = productTable.t['dl_code'].astype(pd.Int32Dtype())
#        if 'outfilename' not in self.imoutcols: self.imoutcols.append('outfilename')
#        if 'dl_code' not in self.imoutcols: self.imoutcols.append('dl_code')
#        if 'dl_str' not in self.imoutcols: self.imoutcols.append('dl_str')
        for ix in ix_selected_products:
            self.mk_outfilename(productTable,
                                ix,
                                outdir=outdir,
                                skip_propID2outsubdir=skip_propID2outsubdir,
                                obsnum2outsubdir=obsnum2outsubdir,
                                jpg_separate_subdir=jpg_separate_subdir,
                                skip_check_if_outfile_exists=skip_check_if_outfile_exists,
                                skip_check_filesize=skip_check_filesize)
        return(self.ix_selected_products)


    def update_obsTable_with_selectedProducts(self, obsTable=None, productTable=None, ix_selected_products=None, filetypes=None):
        """
        links product and observation Table. Counts how many products of which type for each observation.

        Parameters
        ----------
        obsTable : pdastroclass, optional
            Table with observations. The default is None. If None, then self.obsTable is used
        productTable : pdastroclass, optional
            Table with products. The default is None. If None, then self.productTable is used
        ix_selected_products : indeces, optional
            indices of selected products for productTable. The default is None. If None, then self.ix_selected_products is used

        Returns
        -------
        None.

        """
        if obsTable is None:
            obsTable=self.obsTable

        if productTable is None:
           productTable=self.productTable

        if ix_selected_products is None:
            ix_selected_products = self.ix_selected_products

        if filetypes is None:
            filetypes=self.params['filetypes']
        # add the filetypes to the output columns of obsTable
        for filetype in filetypes:
            if filetype not in self.params['outcolumns_obsTable']:
                self.params['outcolumns_obsTable'].append(filetype)

        for filetype in filetypes:
            obsTable.t[filetype]=None

        ixs_obsTable = obsTable.getindices()
        for ix_obsTable in ixs_obsTable:
            # find all products for a given observation
            obsid = obsTable.t.loc[ix_obsTable,'obsid']
            if self.verbose>2:
                print('### obsid',obsid)
            ixs_prodTable = productTable.ix_equal('parent_obsid',obsid,indices=ix_selected_products)
            if self.verbose>2:
                productTable.write(indices=ixs_prodTable)

            # count the product filetypes for each observations
            for ix_prodTable in ixs_prodTable:
                if obsTable.t.loc[ix_obsTable,productTable.t.loc[ix_prodTable,'filetype']] is None:
                    obsTable.t.loc[ix_obsTable,productTable.t.loc[ix_prodTable,'filetype']] = 0
                obsTable.t.loc[ix_obsTable,productTable.t.loc[ix_prodTable,'filetype']] += 1

        return(0)

    def mk_summary_tables(self, obsTable=None, ix_obs_sorted=None, filetypes=None):
        if obsTable is None:
            obsTable=self.obsTable

        if ix_obs_sorted is None:
            ix_obs_sorted = self.ix_obs_sorted


        if filetypes is None:
            filetypes=self.params['filetypes']

        # Add the filetypes to the summary table column
        columns_summary = ['proposal_id','obsnum']
        columns_summary.extend(filetypes)
        columns_summary.extend(['date_start'])

        self.summary = pdastroclass(columns=columns_summary)
        self.summary.t['obsnum']=pd.Series(dtype=pd.Int32Dtype())
        self.summary.t['proposal_id']=pd.Series(dtype=pd.Int64Dtype())

        #Loop through each propID
        propIDs = unique(obsTable.t.loc[ix_obs_sorted,'proposal_id'])
        for propID in propIDs:
            ixs_propID =  obsTable.ix_equal('proposal_id',propID,indices=ix_obs_sorted)

            ixs_propID_notnull = obsTable.ix_not_null('obsnum', indices=ixs_propID)
            obsnums = unique(obsTable.t.loc[ixs_propID_notnull,'obsnum'])
            if len(ixs_propID)!=len(ixs_propID_notnull):
                obsnums.append(pd.NA)
            if self.verbose>2:
                print('\n#### propID:',propID,' obsnums:',obsnums)
            # Loop through echo obsnum for the given propID, and make a summary entry
            for obsnum in obsnums:
                newline_obsnum={'proposal_id':propID}
                if obsnum is pd.NA:
                    ixs_obsnum = obsTable.ix_is_null('obsnum',indices=ixs_propID)
                else:
                    ixs_obsnum = obsTable.ix_equal('obsnum',obsnum,indices=ixs_propID_notnull)
                if self.verbose>2:
                    obsTable.write(indices=ixs_obsnum)
                newline_obsnum['obsnum']=obsnum

                #if obsnum is None:
                #    newline_obsnum['obsnum']=None
                #else:
                #    newline_obsnum['obsnum']=obsnum

                for filetype in filetypes:
                    newline_obsnum[filetype]=obsTable.t.loc[ixs_obsnum,filetype].sum()

                # get the start and end date for an observation
                date_start = obsTable.t.loc[ixs_obsnum,'date_min'].min()
                newline_obsnum['date_start']=date_start

                self.summary.newrow(newline_obsnum)

        #make the type integer for the filetype, proposal_id, and obsnum columns
        for filetype in filetypes:
            self.summary.t[filetype] = self.summary.t[filetype].astype('int')
        self.summary.t['proposal_id'] = self.summary.t['proposal_id'].astype('int')


        #self.summary.t['obsnum'] = self.summary.t['obsnum'].astype('int')

        # sort summary table
        self.ix_summary_sorted = self.summary.ix_sort_by_cols(self.params['sortcols_summaryTable'])

        return(0)

    def TNjpgname(self,jpgname):
        TNname = re.sub('jpg$','TN.jpg',jpgname)
        if TNname==jpgname:
            raise RuntimeError(f'thumbnail name {TNname} is equal to the original name {jpgname}!')
        return(TNname)

    def mk_thumbnail(self, jpgname, TNjpgname, tn_width, tn_height):
        #basewidth = 300
        img = Image.open(jpgname)
        if tn_width is None and tn_height is None:
            raise RuntimeError('Cannot make thumbnail, at least one of width and height has to be specified')
        if tn_height is None:
            wpercent = (tn_width / float(img.size[0]))
            tn_height = int((float(img.size[1]) * float(wpercent)))
        elif tn_width is None:
            hpercent = (tn_height / float(img.size[1]))
            tn_width = int((float(img.size[0]) * float(hpercent)))

        img = img.resize((tn_width, tn_height), Image.ANTIALIAS)
        img.save(TNjpgname)


    def mk_thumbnails(self, productTable, ixs,
                      jpg_separate_subdir=None,figcols=None,
                      tn_width=None, tn_height=None, overwrite=None):

        if tn_width is None:
            tn_width=self.params['webpage_thumbnails_width']

        if tn_height is None:
            tn_height=self.params['webpage_thumbnails_height']

        if jpg_separate_subdir is None:
            jpg_separate_subdir=self.params['jpg_separate_subdir']

        if figcols is None:
           figcols =self.params['webpage_level12_jpgs']

        if overwrite is None:
           overwrite =self.params['webpage_thumbnails_overwrite']

        ixs_uncal  = productTable.ix_equal('filetype','_uncal.fits',indices=ixs)

        if len(ixs_uncal)==0:
            print('WARNING: no entries to do jpg thumbnails for! returning ')
            return(0)

        suffixes = []
        for figcol in figcols:
            suffix = figcol
            if not re.search('^_',suffix): suffix=f'_{suffix}'
            if not re.search('\.jpg$',suffix): suffix+='.jpg'
            suffixes.append(suffix)

        #make the thumbnails
        print(f'Making thumbnails for {suffixes}')
        for ix in ixs_uncal:
            outdir = os.path.dirname(productTable.t.loc[ix,'outfilename'])
            if jpg_separate_subdir: outdir += '/jpg'
            for suffix in suffixes:
                jpgname = f"{outdir}/{productTable.t.loc[ix,'obs_id']}{suffix}"
                TNjpgname = self.TNjpgname(jpgname)
                if overwrite or not os.path.isfile(TNjpgname):
                    if os.path.isfile(jpgname):
                        if self.verbose>2: print(f'Making thumbnail {TNjpgname} for {jpgname}')
                        self.mk_thumbnail(jpgname, TNjpgname, tn_width, tn_height)
        print('Making thumbnails done!')
        return(0)

    def mk_webpage4ixs(self, ixs, productTable,
                       htmlname,asciiname,description,
                       jpg_separate_subdir=None,
                       width=None, height=None,
                       tn_width=None, tn_height=None,
                       figcols=None, fitskeys2table=None,
                       sortcols=None, mkthumbnails=None):

        if width is None:
            width=self.params['webpage_tablefigsize_width']

        if height is None:
            height=self.params['webpage_tablefigsize_height']

        if jpg_separate_subdir is None:
            jpg_separate_subdir=self.params['jpg_separate_subdir']

        if figcols is None:
           figcols =self.params['webpage_level12_jpgs']

        if fitskeys2table is None:
           fitskeys2table = self.params['webpage_fitskeys2table']

        if sortcols is None:
           sortcols =self.params['webpage_sortcols']

        if mkthumbnails is None:
           mkthumbnails =self.params['webpage_mkthumbnails']

        ixs_uncal  = productTable.ix_equal('filetype','_uncal.fits',indices=ixs)
        ixs_uncal =  productTable.ix_sort_by_cols(sortcols,indices=ixs_uncal)

        if len(ixs_uncal)==0:
            print(f'WARNING: no entries for {description}! Removing {htmlname}, {asciiname}, and returning...')
            rmfile(htmlname)
            rmfile(asciiname)
            return(0)

        suffixes = []
        for figcol in figcols:
            suffix = figcol
            if not re.search('^_',suffix): suffix=f'_{suffix}'
            if not re.search('\.jpg$',suffix): suffix+='.jpg'
            suffixes.append(suffix)


        # copy the fitsheader keywords to the table
        productTable.fitsheader2table('outfilename',indices=ixs_uncal,
                                      optionalfitskey=fitskeys2table)

        if mkthumbnails:
            self.mk_thumbnails(productTable, ixs_uncal,
                               jpg_separate_subdir=jpg_separate_subdir,figcols=figcols,
                               tn_width=tn_width, tn_height=tn_height)

        #add the figures to the table
        for ix in ixs_uncal:
            outdir = os.path.dirname(productTable.t.loc[ix,'outfilename'])
            if jpg_separate_subdir: outdir += '/jpg'

            # get the subdir of outdir that is not common with the html dir
            (commonpath,subdir)=split_commonpath(outdir, os.path.dirname(htmlname))
            if commonpath == '': raise RuntimeError(f'Somethings is wrong, htmldir={os.path.dirname(htmlname)} has no common path with outdir={outdir}')

            for suffix,figcol in zip(suffixes,figcols):
                if  subdir=='':
                    jpgname = f"{productTable.t.loc[ix,'obs_id']}{suffix}"
                else:
                    jpgname = f"{subdir}/{productTable.t.loc[ix,'obs_id']}{suffix}"

                if mkthumbnails:
                    TNjpgname = self.TNjpgname(jpgname)
                else:
                    TNjpgname=None

                # make a thumbnail that links to the full size image
                productTable.t.loc[ix,figcol]=image4table(jpgname,TNimagename=TNjpgname,width=width,height=height)
                #productTable.t.loc[ix,figcol]=addlink2string(imagestring4web(jpgname,width=None,height=p),jpgname)

        outcols4html = self.params['webpage_cols4table']

        # Make the ascci table(without the jpg cols), and save it.
        # Use a try/except to insulate against the case where requested suffixes for download in webpage_level12_jpgs
        # are not included in webpage_cols4table
        outcols4ascci = copy.deepcopy(outcols4html)
        for figcol in figcols:
            try:
                outcols4ascci.remove(figcol)
            except ValueError:
                pass
        print(f'Writing {description} pandas ascii to {asciiname}')
        print(outcols4ascci)
        productTable.write(filename=asciiname, indices=ixs_uncal, columns=outcols4ascci)

        # write the table to index.html
        print(f'writing {description} html to {htmlname}')
        f=open(htmlname,'w')
        s_asciilink = addlink2string('ascii-table',os.path.basename(asciiname))
        f.writelines([f'Level 1+2 Products for {description} ({s_asciilink} here)'])
        (errorflag,lines)=productTable.write(return_lines=True, indices=ixs_uncal, columns=outcols4html, htmlflag=True, htmlsortedtable=True, escape=False)
        if errorflag: raise RuntimeError(f'Soemthing went wrong when doing the webpage table for {description}')
        f.writelines(lines)
        f.close()

        # Make sure sortable.js is in the html directory
        htmldir = os.path.dirname(htmlname)
        jsfilename = os.path.join(os.path.realpath(os.path.dirname(os.path.join('../', __file__))), 'sortable.js')
        dest_jsfilename = f'{htmldir}/sortable.js'
        if not os.path.isfile(dest_jsfilename):
            if not os.path.isfile(jsfilename):
                raise RuntimeError(f'java script {jsfilename} for sortable tables does not exist!')
            shutil.copy(jsfilename, dest_jsfilename)

        return(0)

    def mk_webpages4propIDs(self, productTable=None, ix_selected_products=None,
                            jpg_separate_subdir=None,
                            width=None, height=None,
                            figcols=None, fitskeys2table=None,
                            sortcols=None, mkthumbnails=None):
        if productTable is None:
            # make a deep copy so that the thumbnail columns are not copied into the self.productTable
            productTable= copy.deepcopy(self.productTable)

        if ix_selected_products is None:
            ix_selected_products = self.ix_selected_products

        propIDs = unique(productTable.t.loc[ix_selected_products,'proposal_id'])
        for propID in propIDs:
            if self.verbose:
                print(f'Making webpage for propID={propID}')
            propID = int(propID)
            htmldir = f'{self.outrootdir}/{propID:05d}'
            htmlname = f'{htmldir}/index.html'
            asciiname = f'{htmldir}/{propID:05d}.products.txt'
            # description is used for the title of the webpage etc
            description=f'propID={propID}'
            # get all indices for propID
            ixs_propID = productTable.ix_equal('proposal_id',propID,indices=ix_selected_products)

            self.mk_webpage4ixs(ixs_propID, productTable,
                                htmlname,asciiname,description,
                                jpg_separate_subdir=jpg_separate_subdir,
                                width=width, height=height,
                                figcols=figcols, fitskeys2table=fitskeys2table,
                                sortcols=sortcols,mkthumbnails=mkthumbnails)


    def mk_all_tables(self, filetypes=None, showtables=True):

        if filetypes is not None:
            self.params['filetypes'] = filetypes

        # get the MJD limits, based on --mjdlimits --lockbacktime --mjdmin --mjdmax
        mjd_min, mjd_max = self.get_mjd_limits()

        # get the observations: stored in self.obsTable
        self.observation_query(propID=self.params['propID'], mjd_min = mjd_min, mjd_max = mjd_max)

        if len(self.obsTable.t)==0:
            print('\n################################\nNO OBSERVATIONS FOUND! exiting....\n################################')
            self.ix_selected_products = []
            return(0)

        if self.verbose>1:
            print(self.obsTable.t)

        # get the products:  stored in self.productTable
        # this list contains in general several entries for each observations.
        # If not otherwise specified, uses self.obsTable as starting point
        self.product_query(Nobs_per_batch=self.params['Nobs_per_batch'])
        if self.verbose>1:
            print(self.productTable.t)

        # select the products to download
        # If not otherwise specified, uses self.productTable as starting point
        # uses self.filetypes for filtering, which is set by optional parameters.
        # the selected products indices are put into self.ix_selected_products
        self.product_filter()

        # get selected obsnums if specified
        self.ix_selected_products = self.select_obsnums(self.params['obsnums'])

        # make some obsid cuts!
        self.ix_selected_products = self.obsid_select(self.params['obsid_select'])
        self.ix_selected_products = self.obsid_list(self.params['obsid_list'])

        # get selected SCAs if specified
        self.ix_selected_products = self.select_sca(self.params['sca'])

        # only keep entries in the obstable that are parent_obsid in product table
        self.ix_obs_sorted = self.update_obstable_indices(self.ix_selected_products)

        # define the output filenames, and check if they exist.
        self.mk_outfilenames(skip_propID2outsubdir=self.params['skip_propID2outsubdir'],
                             obsnum2outsubdir=self.params['obsnum2outsubdir'],
                             jpg_separate_subdir=self.params['jpg_separate_subdir'],
                             skip_check_if_outfile_exists=self.params['skip_check_if_outfile_exists'],
                             skip_check_filesize=self.params['skip_check_filesize'])

        # update obsTable with selected products: count the different filetypes, and also update the obsnum if None in obsTable
        self.update_obsTable_with_selectedProducts()

        # print the table to screen if either verbose or not downloading
        if showtables:
            print('\n######################\n### Selected Products:\n######################')
            self.productTable.write(indices=self.ix_selected_products, columns=self.params['outcolumns_productTable'])

            print('\n#################\n### Observations:\n#################')
            self.obsTable.write(columns=self.params['outcolumns_obsTable'],indices=self.ix_obs_sorted)

        # make summary table
        self.mk_summary_tables()
        if showtables:
            print('\n##########################\n### Summary propID/obsnum:\n##########################')
            self.summary.write(indices=self.ix_summary_sorted)

        # save the tables if wanted
        if self.params['savetables'] is not None:
            prodfile = os.path.join(self.outrootdir, self.params['savetables']+'.selprod.txt')
            obsfile = os.path.join(self.outrootdir, self.params['savetables']+'.obs.txt')
            summaryfile = os.path.join(self.outrootdir, self.params['savetables']+'.summary.txt')

            self.productTable.write(filename=prodfile, indices=self.ix_selected_products, columns=self.params['outcolumns_productTable'], verbose=2)
            self.obsTable.write(filename=obsfile, columns=self.params['outcolumns_obsTable'], indices=self.ix_obs_sorted, verbose=2)
            self.summary.write(filename=summaryfile, indices=self.ix_summary_sorted, verbose=2)

        return(0)



if __name__ == '__main__':

    query = query_mast()
    parser = query.define_options()
    args = parser.parse_args()

    # the arguments are saved in query.params
    query.get_arguments(args)
    if query.verbose>2:
        print('params:', query.params)

    # use arguments or $API_MAST_TOKEN to login
    query.login(raiseErrorFlag=False)

    # self.outrootdir is set depending on outrootdir and outsubdir in cfg file or through the options --outrootdir and --outsubdir
    query.set_outrootdir()
    if query.verbose: print(f'Outdir: {query.outdir}')

    # make the tables
    query.mk_all_tables()

    # make the webpages
    if args.makewebpages:
        query.mk_webpages4propIDs()



