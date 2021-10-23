#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed July 28 11:01:23 2021

@author: arest, bhilbert
"""
from astroquery.mast import Mast
from astropy.time import Time
import numpy as np
import pandas as pd
#import astropy.io.fits as fits
#from astropy.nddata import bitmask
import astroquery
import argparse,os,sys,re,types
import yaml

# pdastroclass is wrapper around pandas.
from pdastro import pdastroclass,unique,AnotB


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
                'Caom_search':'Mast.Caom.Filtered.JwstOps',
                'Product_search':'Mast.Caom.Products.JwstOps'
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
        self.outdir = None

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

        # columns returned from MAST to the obsTable
        #self.mastcolumns_obsTable = ['proposal_id','dataURL','obsid','obs_id','t_min','t_exptime']

        # output columns for the tables. Note that the columns for the individual filetypes
        # are automatically added to the obsTable
        #self.outcolumns_productTable = ['proposal_id','obsnum','obsID','parent_obsid','obs_id','dataproduct_type','productFilename','filetype','calib_level','size','description']
        #self.outcolumns_obsTable = ['proposal_id','obsnum','obsid','obs_id','t_min','t_exptime','date_min']

        # self.params will be populated with the arguments
        self.params = {}

        # columns returned from MAST to the obsTable
        # These are the default values, they can be changed in the config file
        self.params['mastcolumns_obsTable']=['proposal_id','dataURL','obsid','obs_id','t_min','t_exptime']

        # output columns for the tables. Note that the columns for the individual filetypes
        # are automatically added to the obsTable
        # These are the default values, they can be changed in the config file
        self.params['outcolumns_productTable']=['proposal_id','obsnum','obsID','parent_obsid','obs_id','dataproduct_type','productFilename','filetype','calib_level','size','description']
        self.params['outcolumns_obsTable']=['proposal_id','obsnum','obsid','obs_id','t_min','t_exptime','date_min']

        # The productTable is sorted based on these columns  (can also be set in config file)
        self.params['sortcols_productTable']=['calib_level','filetype','obsID']
        # The obsTable is sorted based on these columns  (can also be set in config file)
        self.params['sortcols_obsTable']=['date_min','proposal_id','obsnum']
        # The summary table is sorted based on these columns (can also be set in config file)
        self.params['sortcols_summaryTable']=['date_start','proposal_id','obsnum']

        self.params['instrument']='nircam'
        self.params['filetypes']=['uncal']
        self.params['guidestars']=False
        self.params['lookbacktime']=1.0

    def define_options(self,parser=None,usage=None,conflict_handler='resolve'):
        if parser is None:
            parser = argparse.ArgumentParser(usage=usage,conflict_handler=conflict_handler)

        # default for config file, if available
        if 'JWST_QUERY_CFGFILE' in os.environ and os.environ['JWST_QUERY_CFGFILE'] != '':
            cfgfilename = os.environ['JWST_QUERY_CFGFILE']
        else:
            cfgfilename = None

        # default for token is $MAST_API_TOKEN
        if 'MAST_API_TOKEN' in os.environ:
            defaulttoken = os.environ['MAST_API_TOKEN']
        else:
            defaulttoken = None


        parser.add_argument('-i', '--instrument', type=str, default=None, choices=['niriss','nircam','nirspec','miri','fgs'], help='Instrument.  (default=%(default)s)')

        parser.add_argument('-v','--verbose', default=0, action='count')
        parser.add_argument('--propID', type=int, default=None, help='Search for data for this proposal ID(=APT #) only.')
        parser.add_argument('--guidestars', action='store_true', default=None, help='Don\'t skip guidestars. By default, they are skipped')
        parser.add_argument('-f','--filetypes',  type=str, nargs="+", default=None, help=('List of product filetypes to get, e.g., _uncal.fits or _uncal.jpg. If only letters, then _ and .fits are added, for example uncal gets expanded to _uncal.fits. Typical image filetypes are uncal, rate, rateints, cal (default=%(default)s)'))
        parser.add_argument('--calib_levels',  type=int, nargs="+", default=None, help=('Only select products with the specified calibration levels (calib_level column in productTable) (default=%(default)s)'))

        parser.add_argument('-c','--configfile', type=str, default=cfgfilename, help='optional config file. default is set to $JWST_QUERY_CFGFILE. Use -vvv to see the full list of all set parameters.')
        parser.add_argument('--login', default=None, nargs=2, help='username and password for login')
        parser.add_argument('--token', type=str, default=defaulttoken, help='MAST API token. You can also set the token in the environment variable \$MAST_API_TOKEN')

        parser.add_argument('--outrootdir', default=None, help='output root directory (default=%(default)s)')
        parser.add_argument('--outsubdir', default=None, help='osubdir added to utput root directory (default=%(default)s)')
        parser.add_argument('--skip_propID2outsubdir', action='store_true', default=None, help='By default, the APT proposal ID is added as a subdir to the output directory. You can skip this with this option (default=%(default)s)')
        parser.add_argument('--skip_check_if_outfile_exists', action='store_true', default=None,help='Don\'t check if output files exists. This makes it faster for large lists')

        parser.add_argument('-e','--obsid_select', nargs="+", default=[], help='Specify obsid range applied to "obsID" column in the PRODUCT table. If single value, then exact match. If single value has "+" or "-" at the end, then it is a lower and upper limit, respectively. If two values, then range. Examples: 385539+, 385539-, 385539 385600 (default=%(default)s)')
        parser.add_argument('-l','--obsid_list', nargs="+", default=[], help='Specify list of obsid applied to "obsID" column in the PRODUCT table. examples: 385539 385600 385530 (default=%(default)s)')

        time_group = parser.add_argument_group("Time constraints for the observation/product search")

        time_group.add_argument('-l', '--lookbacktime', type=float, default=None, help='lookback time in days.')
#        time_group.add_argument('--mjd_min', type=float, default=None, help='minimum MJD. overrides lookback time.')
#        time_group.add_argument('--mjd_max', type=float, default=None, help='maximum MJD. overrides lookback time.')
#        time_group.add_argument('-m', '--mjd_limits', default=None, type=float, nargs=2, help='specify the MJD limits. overrides lookback time and mjd_min/max optional arguments.')
#        time_group.add_argument('-d', '--date_limits', default=None, type=str, nargs=2, help='specify the date limits (ISOT format). overrides lookback time and mjd* optional arguments.')
        time_group.add_argument('-d','--date_select', nargs="+", default=[], help='Specify date range (MJD or isot format) applied to "dateobs_center" column. If single value, then exact match. If single value has "+" or "-" at the end, then it is a lower and upper limit, respectively. Examples: 58400+, 58400-,2020-11-23+, 2020-11-23 2020-11-25  (default=%(default)s)')

        time_group.add_argument('--lre3', action='store_true', default=None, help='Use the LRE-3 date limits. Overrides lookback and mjd* options.')
        time_group.add_argument('--lre4', action='store_true', default=None, help='Use the LRE-4 date limits. Overrides lookback and mjd* options.')
        time_group.add_argument('--lre5', action='store_true', default=None, help='Use the LRE-5 date limits. Overrides lookback and mjd* options.')
        time_group.add_argument('--lre6', action='store_true', default=None, help='Use the LRE-6 date limits. Overrides lookback and mjd* options.')

        parser.add_argument('-s', '--savetables', type=str, default=None, help='save the tables (selected products, obsTable, summary with suffix selprod.txt, obs.txt, summary.txt, respectively) with the specified string as basename (default=%(default)s)')


        return(parser)

    def get_arguments(self, args, configfile = None):
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
                                raise RuntimeError("environment variable %s used in config file, but not set!" % name)
                            envval=os.environ[name]
                            subpattern='\$%s' % (name)
                            paramsdict[param] = re.sub(subpattern,envval,paramsdict[param])
                elif isinstance(paramsdict[param], dict):
                #elif type(dict[param]) is types.DictType:
                    # recursive: sub environment variables down the dictiionary
                    subenvvarplaceholder(paramsdict[param])
            return(0)


        # get the parameters from the config file
        if args.configfile is not None:
            #cfgparams = yaml.load_file(args.configfile)
            if not os.path.isfile(args.configfile):
                raise RuntimeError('config file %s does not exist!' % (args.configfile))
            print(f'Loading config file {args.configfile}')
            cfgparams = yaml.load(open(args.configfile,'r'), Loader=yaml.FullLoader)
            self.params.update(cfgparams)
            subenvvarplaceholder(self.params)

            if args.verbose>2:
                print('\n### CONFIG FILE PARAMETERS:')
                for p in cfgparams:
                    print('config file args: setting %s to' % (p),cfgparams[p])

        # Go through optional parameters.
        # 'None' does not overwrite previously set parameters (either default or config file)
        if args.verbose>2:
            print('\n### OPTIONAL COMMAND LINE PARAMETERS:')
        argsdict = vars(args)
        for arg in argsdict:

            # skip config file
            if arg=='configfile': continue

            if argsdict[arg] is not None:
                if args.verbose>2:
                    print('optional args: setting %s to %s' % (arg,argsdict[arg]))
                self.params[arg]=argsdict[arg]
            else:
                if arg not in  self.params:
                    self.params[arg]=None

        # propID
        if self.params['propID'] is not None:
            self.params['propID'] = '%05d' % (int(self.params['propID']))
        self.verbose = self.params['verbose']

        if self.verbose>2:
            print('\n### FINAL PARAMETERS:')
            for p in self.params:
                print(p,self.params[p])

        print('INSTRUMENT:',self.params['instrument'])
        #sys.exit(0)

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
            self.JwstObs.login(token=token, store_token=True)
        return(0)

    def obsid_select(self,obsid_select,
                     productTable=None,
                     ix_selected_products=None):

        if productTable==None:
           productTable=self.productTable

        if ix_selected_products is None:
            ix_selected_products = self.ix_selected_products

        # parse trailing '+' and '-', and get limits
        limits = getlimits(obsid_select)
        if limits is None:
            return(ix_selected_products)

        for i in range(len(limits)):
            if limits[i] is not None: limits[i]=int(limits[i])

        ixs = self.productTable.ix_remove_null('obsID',indices=ix_selected_products)
        ixs_keep = self.productTable.ix_inrange('obsID',limits[0],limits[1],indices=ixs)
        print(f'obsid cut {limits[0]} - {limits[1]}: keeping {len(ixs_keep)} from {len(ixs)}')
        return(ixs_keep)

    def obsid_list(self,obsid_list,
                     productTable=None,
                     ix_selected_products=None):

        if productTable==None:
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
        print(f'obsid list cut: keeping {len(ixs_keep)} from {len(ix_selected_products)}')
        return(ixs_keep)


    def get_mjd_limits(self, lookbacktime=None, date_select=None, lre3=False, lre4=False, lre5=False, lre6=False):
        if lookbacktime is None: lookbacktime = self.params['lookbacktime']

        if date_select is None: date_select = self.params['date_select']
        if not lre3: lre3 = self.params['lre3']
        if not lre4: lre4 = self.params['lre4']
        if not lre5: lre5 = self.params['lre5']
        if not lre6: lre6 = self.params['lre6']

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
        elif lre3:
            if self.verbose: print('setting mjd limits to LRE-3!')
            mjd_min = Time('2021-05-19', format='iso').mjd
            mjd_max = mjd_min+7
        elif lre4:
            if self.verbose: print('setting mjd limits to LRE-4!')
            mjd_min = Time('2021-06-14', format='iso').mjd
            mjd_max = mjd_min+7
        elif lre5:
            if self.verbose: print('setting mjd limits to LRE-5!')
            mjd_min = Time('2021-08-08', format='iso').mjd
            mjd_max = mjd_min+7
        elif lre6:
            if self.verbose: print('setting mjd limits to LRE-6!')
            mjd_min = Time('2021-10-17', format='iso').mjd
            mjd_max = mjd_min+10
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
        '''

        if instrument is None:
            instrument = self.params['instrument']

        print('MJD range:',mjd_min, mjd_max)

        columns = ','.join(self.params['mastcolumns_obsTable'])
        service = self.SERVICES['Caom_search']
        params = {"columns":columns,
                  "filters":[
                  {"paramName":"obs_collection","values":["JWST"]},
                  {"paramName":"instrument_name","values":[instrument]},
#                  {"paramName":"proposal_id","values":[propID]},
                  {"paramName":"t_min",
                    "values":[{"min":mjd_min,"max":mjd_max}]},
                  ]}

        # Only add propID entry if not None. If it is set to None it doesn't work!
        if propID is not None:
            params["filters"].append({"paramName":"proposal_id","values":[propID]})

        if self.verbose:
            print('\n#### Querying Obstable....')
        if token is not None:
            self.obsTable.t = self.JwstObs.service_request(service, params,mast_token=token).to_pandas()
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
        self.obsTable.t['obsnum']=None
        for ix in ixs:
            m = obsnumsearch.search(self.obsTable.t.loc[ix,'obs_id'])
            if m is not None:
                self.obsTable.t.loc[ix,'obsnum']=int(m.groups()[0])
            else:
                self.obsTable.t.loc[ix,'obsnum']=None

        self.obsTable.mjd2dateobs('t_min','date_min')

        self.ix_obs_sorted = self.obsTable.ix_sort_by_cols(self.params['sortcols_obsTable'])

        # make proposal_id integer
        if 'proposal_id' in self.obsTable.t.columns:
            self.obsTable.t['proposal_id']=self.obsTable.t['proposal_id'].astype('int')
        if 'obsid' in self.obsTable.t.columns:
            self.obsTable.t['obsid']=self.obsTable.t['obsid'].astype('int')

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
        if obsTable==None:
            obsTable=self.obsTable

        if productTable==None:
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
                ixs = productTable.ix_remove_null('obsnum',indices=ixs_prodTable)
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


    def product_query(self, obsTable=None, guidestars=None):
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

        # fill the suffix column with the suffix of the form _bla1.bla2, e.g. _uncal.fits
        # This will later be used to figure out
        self.productTable.t['filetype'] = self.productTable.t['productFilename'].str.extract(r'(\_[a-zA-Z0-9]+\.[a-zA-Z0-9]+)$')

        # Find the obsnum # from the filename if possible.
        ixs = self.productTable.getindices()
        obsnumsearch = re.compile('^jw\d{5}(\d{3})\d{3}\_')
        self.productTable.t['obsnum']=None
        for ix in ixs:
            m = obsnumsearch.search(self.productTable.t.loc[ix,'obs_id'])
            if m is not None:
                self.productTable.t.loc[ix,'obsnum']=int(m.groups()[0])
            else:
                self.productTable.t.loc[ix,'obsnum']=None

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
        if productTable==None:
           productTable=self.productTable

        if filetypes==None:
            filetypes=self.params['filetypes']

        if calib_levels==None:
            calib_levels=self.params['calib_levels']

        if self.verbose:
            print('\n#### Selecting products....')

        # if necessary, add leading '_' and suffix '.fits'
        for i in range(len(filetypes)):
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
            print('VVVVVVVVVVv',calib_levels)
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
        for filetype in filetypes:

            # make sure the '.' in the regular expression is literal, and also add '$' to the end
            regex=re.sub('\.','\.',filetype)+'$'

            # get the indices for matching filetype ...
            ix_matching_filetype = self.productTable.ix_matchregex('filetype',regex,indices=ix_products)
            # ... and add them to the list of good indices
            self.ix_selected_products.extend(ix_matching_filetype)

            if self.verbose>1:
                print('%d products with filetype regex matching %s' % (len(ix_matching_filetype),regex))

        self.ix_selected_products = self.productTable.ix_sort_by_cols(self.params['sortcols_productTable'],indices=self.ix_selected_products)


        if self.verbose:
            print('%d products with correct filetypes left' % (len(self.ix_selected_products)))

        return(self.ix_selected_products)

    def mk_outfilename(self, productTable, ix,
                       outdir=None,
                       skip_propID2outsubdir=False,
                       obsnum2outsubdir=False,
                       #info2filename=False,
                       skip_check_if_outfile_exists=False):


        if outdir is None:
            outdir = self.outdir
        if outdir is None:
            raise RuntimeError("outdir is not defined!")
        outdir.rstrip("/")


        if not skip_propID2outsubdir:
            #outdir += f'/{productTable.t.loc[ix,"proposal_id"]}'
            outdir += '/{:05d}'.format(productTable.t.loc[ix,"proposal_id"])

        if obsnum2outsubdir:
            if productTable.t.loc[ix,"obsnum"] is pd.NA:
                outdir += '/NA'
            else:
                outdir += f'/{productTable.t.loc[ix,"obsnum"]}'


        outfilename = productTable.t.loc[ix,'productFilename']

        fullfilename = os.path.abspath(f'{outdir}/{outfilename}')
        productTable.t.loc[ix,'outfilename']=fullfilename

        if not skip_check_if_outfile_exists:
            if os.path.exists(fullfilename):
                productTable.t.loc[ix,'dl_code'] = 2
                productTable.t.loc[ix,'dl_str'] = 'exists'
            else:
                productTable.t.loc[ix,'dl_code'] = 0
                productTable.t.loc[ix,'dl_str'] = None


    def set_outdir(self,outrootdir=None,outsubdir=None):
        self.outdir = outrootdir
        if self.outdir is None:
            self.outdir = self.params['outrootdir']
        # if outdir is not defined, use '.'
        if self.outdir is None or self.outdir == '':
            self.outdir = '.'


        if outsubdir is not None and outsubdir!='':
            self.outdir  += f'/{outsubdir}'
        elif self.params['outsubdir'] is not None and self.params['outsubdir']!='':
            self.outdir  += f'/{self.params["outsubdir"]}'

        self.outdir = os.path.abspath(self.outdir)
        return(self.outdir)

    def mk_outfilenames(self,
                        productTable=None,
                        ix_selected_products=None,
                        outdir=None,
                        skip_propID2outsubdir=False,
#                        info2filename=False,
                        skip_check_if_outfile_exists=False):

        if productTable==None:
           productTable=self.productTable

        if ix_selected_products is None:
            ix_selected_products = self.ix_selected_products



        productTable.t.loc[ix_selected_products,'outfilename'] = None
        productTable.t.loc[ix_selected_products,'dl_code'] = np.nan
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
                                skip_check_if_outfile_exists=skip_check_if_outfile_exists)
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
        if obsTable==None:
            obsTable=self.obsTable

        if productTable==None:
           productTable=self.productTable

        if ix_selected_products is None:
            ix_selected_products = self.ix_selected_products

        if filetypes==None:
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

    def mk_summary_tables(self, obsTable=None, filetypes=None):
        if obsTable==None:
            obsTable=self.obsTable

        if filetypes==None:
            filetypes=self.params['filetypes']

        # Add the filetypes to the summary table column
        columns_summary = ['proposal_id','obsnum']
        columns_summary.extend(filetypes)
        columns_summary.extend(['date_start'])

        self.summary = pdastroclass(columns=columns_summary)
        self.summary.t['obsnum']=pd.Series(dtype=pd.Int32Dtype())
        self.summary.t['proposal_id']=pd.Series(dtype=pd.Int64Dtype())

        #Loop through each propID
        propIDs = unique(obsTable.t['proposal_id'])
        for propID in propIDs:
            ixs_propID =  obsTable.ix_equal('proposal_id',propID)

            ixs_propID_notnull = obsTable.ix_remove_null('obsnum', indices=ixs_propID)
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

    def mk_all_tables(self, filetypes=None, showtables=True):

        if filetypes is not None:
            self.params['filetypes'] = filetypes

        # get the MJD limits, based on --mjdlimits --lockbacktime --mjdmin --mjdmax --lre3 --lre4 --lre5 --lre6
        mjd_min, mjd_max = self.get_mjd_limits()

        # get the observations: stored in self.obsTable
        self.observation_query(self.params['propID'], mjd_min = mjd_min, mjd_max = mjd_max)

        if len(self.obsTable.t)==0:
            print('\n################################\nNO OBSERVATIONS FOUND! exiting....\n################################')
            return(0)

        if self.verbose>1:
            print(self.obsTable.t)


        # get the products:  stored in self.productTable
        # this list contains in general several entries for each observations.
        # If not otherwise specified, uses self.obsTable as starting point
        self.product_query()
        if self.verbose>1:
            print(self.productTable.t)
            #print(self.productTable.columns)

        # select the products to download
        # If not otherwise specified, uses self.productTable as starting point
        # uses self.filetypes for filtering, which is set by optional parameters.
        # the selected products indices are put into self.ix_selected_products
        self.product_filter()

        # make some obsid cuts!
        self.ix_selected_products = self.obsid_select(self.params['obsid_select'])
        self.ix_selected_products = self.obsid_list(self.params['obsid_list'])


        # definte the output filenames, and check if they exist.
        self.mk_outfilenames(skip_propID2outsubdir=self.params['skip_propID2outsubdir'],
                             skip_check_if_outfile_exists=self.params['skip_check_if_outfile_exists'])

        # update obsTable with selected products: count the different filetypes, and also update the obsnum if None in obsTable
        self.update_obsTable_with_selectedProducts()

        # print the table to screen if either verbose or not downloading
        if showtables:
            print('\n######################\n### Selected Products:\n######################')
            self.productTable.write(indices=self.ix_selected_products,columns=self.params['outcolumns_productTable'])

            print('\n#################\n### Observations:\n#################')
            self.obsTable.write(columns=self.params['outcolumns_obsTable'],indices=self.ix_obs_sorted)

        # make summary table
        self.mk_summary_tables()
        if showtables:
            print('\n##########################\n### Summary propID/obsnum:\n##########################')
            self.summary.write(indices=self.ix_summary_sorted)

        # save the tables if wanted
        if self.params['savetables'] is not None:
            self.productTable.write(filename=self.params['savetables']+'.selprod.txt',indices=self.ix_selected_products,columns=self.params['outcolumns_productTable'],verbose=2)
            self.obsTable.write(filename=self.params['savetables']+'.obs.txt',columns=self.params['outcolumns_obsTable'],indices=self.ix_obs_sorted,verbose=2)
            self.summary.write(filename=self.params['savetables']+'.summary.txt',indices=self.ix_summary_sorted,verbose=2)

        return(0)

if __name__ == '__main__':

    query = query_mast()
    parser = query.define_options()
    args = parser.parse_args()

    # the arguments are saved in query.params
    query.get_arguments(args)
    if query.verbose>2:
        print('params:', query.params)

    query.login(raiseErrorFlag=False)

    # make the tables
    query.mk_all_tables()




