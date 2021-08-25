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
#import astropy.io.fits as fits
from astropy.nddata import bitmask
import astroquery
import argparse,os,sys,re,io
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

        parser.add_argument('-i', '--instrument', type=str, default=None, choices=['niriss','nircam','nirspec','miri','fgs'], help='Instrument.  (default=%(default)s)')

        parser.add_argument('-v','--verbose', default=0, action='count')
        parser.add_argument('--propID', type=int, default=None, help='Download data for this program ID.')
        parser.add_argument('--guidestars', action='store_true', default=None, help='Don\'t skip guidestsars  (default=%(default)s)')
        parser.add_argument('-f','--filetypes',  type=str, nargs="+", default=None, help=('List of product filetypes to get, e.g., _uncal.fits or _uncal.jpg. If only letters, then _ and .fits are added, for example uncal gets expanded to _uncal.fits. Typical image filetypes are uncal, rate, rateints, cal (default=%(default)s)'))

        parser.add_argument('-c','--configfile', type=str, default=cfgfilename, help='optional config file. default is set to $JWST_QUERY_CFGFILE. Use -vvv to see the full list of all set parameters.')
        parser.add_argument('--token', type=str, default=None, help='MAST API token. You can also set the token in the environment variable $MAST_API_TOKEN')


        parser.add_argument('-d','--download', action='store_true', default=None, help='download the data.  (default=%(default)s)')
        #parser.add_argument('-o','--outrootdir', default=outrootdir, help=('output root directory.''(default=%(default)s)'))
        #parser.add_argument('--outsubdir', default=None, help=('subdir added to the output root directory. If None, then propID is used (default=%(default)s)'))
        #parser.add_argument('--clobber', action='store_true', default=False, help='existing files are overwritten, otherwise they are skipped (default=%(default)s)')

        parser.add_argument('-m', '--mjd_limits', default=None, type=float, nargs=2, help='specify the MJD limits. overrides lookback time and mjd_min/max optional arguments(default=%(default)s)')
        parser.add_argument('-l', '--lookbacktime', type=float, default=None, help='lookback time in days (default=%(default)s)')
        parser.add_argument('--mjd_min', type=float, default=None, help='minimum MJD. overrides lookback time  (default=%(default)s)')
        parser.add_argument('--mjd_max', type=float, default=None, help='maximum MJD. (default=%(default)s)')
        
        parser.add_argument('--lre3', action='store_true', default=None, help='Use the LRE-3 mjd limits (default=%(default)s)')
        parser.add_argument('--lre4', action='store_true', default=None, help='Use the LRE-4 mjd limits (default=%(default)s)')
        parser.add_argument('--lre5', action='store_true', default=None, help='Use the LRE-5 mjd limits (default=%(default)s)')
        parser.add_argument('--lre6', action='store_true', default=None, help='Use the LRE-6 mjd limits (default=%(default)s)')
        
        parser.add_argument('-s', '--savetables', type=str, default=None, help='save the tables (selected products, obsTable, summary with suffix selprod.txt, obs.txt, summary.txt, respectively) with the specified string as basename (default=%(default)s)')

        
        return(parser)

    def get_arguments(self, args, configfile = None): 
        '''

        Parameters
        ----------
        args : list
            pass the command line arguments to self.params.
            make sure propID has the correct format (5 digit string)
        configfile : string, optional
            Config filename. The default is None. If None, then 
            $JWST_QUERY_CONFIGFILE is used if exists.

        Returns
        -------
        None.

        '''

        # get the parameters from the config file
        if args.configfile is not None:
            #cfgparams = yaml.load_file(args.configfile)
            if not os.path.isfile(args.configfile):
                raise RuntimeError('config file %s does not exist!' % (args.configfile))
            cfgparams = yaml.load(open(args.configfile,'r'), Loader=yaml.FullLoader) 
            self.params.update(cfgparams)
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

    def set_outdir(self, outrootdir=None, outsubdir=None):
        if outrootdir is None: outrootdir = self.params['outrootdir']
        if outsubdir is None: outsubdir = self.params['outsubdir']
        if outsubdir is None: outsubdir = self.params['propID']
        
        self.outdir = '%s/%s' % (outrootdir,outsubdir)


    def get_mjd_limits(self, mjd_limits=None, lookbacktime=None, mjd_min=None, mjd_max=None, lre3=False, lre4=False, lre5=False):
        if mjd_limits is None: mjd_limits = self.params['mjd_limits']
        if lookbacktime is None: lookbacktime = self.params['lookbacktime']
        if mjd_min is None: mjd_min = self.params['mjd_min']
        if mjd_max is None: mjd_max = self.params['mjd_max']
        if not lre3: lre3 = self.params['lre3']
        if not lre4: lre4 = self.params['lre4']
        if not lre5: lre5 = self.params['lre5']
        
        if self.params['mjd_limits'] is not None:
            mjd_min =  self.params['mjd_limits'][0]
            mjd_max =  self.params['mjd_limits'][1]

            
        if lre3:
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
        else:
            if (mjd_min is None):
                if self.params['lookbacktime']>0.0:
                    if self.verbose>1: print('Note: Looking back %.1f days' % (self.params['lookbacktime']))
                    mjd_min = Time.now().mjd - self.params['lookbacktime']
            
        if mjd_max is None:
            mjd_max = Time.now().mjd+0.1

        return(mjd_min, mjd_max)
            
    def observation_query(self, propID=None, instrument=None, mjd_min=None, mjd_max=None):
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
        if self.params['token'] is not None:
            self.obsTable.t = self.JwstObs.service_request(service, params,mast_token=self.params['token']).to_pandas()
        else:
            self.obsTable.t = self.JwstObs.service_request(service, params).to_pandas()
        if self.verbose: 
            print('Obstable obtained: length:',len(self.obsTable.t))
            
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
        self.obsTable.t['proposal_id']=self.obsTable.t['proposal_id'].astype('int')

        if self.verbose>1: print('Obstable columns:',self.obsTable.t.columns)

        return self.obsTable    

    def product_query(self, obsTable=None):
        '''
        Perform query for data products based on obs_id's in observation table
        '''
        
        if obsTable==None:
            obsTable=self.obsTable
        
        # query MAST for all products for the obsid's
        obsids = ','.join(obsTable.t['obsid'])
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

        if self.verbose: print('productTable columns:',self.productTable.t.columns)
        if self.verbose: 
            allfiltypes = unique(self.productTable.t['filetype'])
            print('List of all filetypes of obtained products:',allfiltypes)
        
        return self.productTable
    
    def product_filter(self, productTable=None, filetypes=None, gs_omit=True):
        '''
        Filter the list of products based on the filetype (or better suffix)
        '''
        if productTable==None:
           productTable=self.productTable
        
        if filetypes==None:
            filetypes=self.params['filetypes']

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
                
            # If the obsnum cannot be determined from the obsTable 'obs_id', then get it from the products!
            if obsTable.t.loc[ix_obsTable,'obsnum'] is None:
                # remove None entries
                ixs = productTable.ix_remove_null('obsnum',indices=ixs_prodTable)
                obsnum = unique(productTable.t.loc[ixs,'obsnum'])
                if len(obsnum)==0:
                    obsTable.write(indices=ix_obsTable)
                    productTable.write(indices=ixs_prodTable)
                    raise RuntimeError('No obsnum for these observations and products!')
                if len(obsnum)>1:
                    print('More than one obsnum!',obsnum)
                    raise RuntimeError('BUGGG!!!!!')
                obsTable.t.loc[ix_obsTable,'obsnum'] = obsnum
                
            # If there are entries in productTable for which 'obsnum' is None (happens to some calib_level=3 products) fill them up
            ixs_null = productTable.ix_is_null('obsnum',indices=ixs_prodTable)
            if len(ixs_null)>0:
                productTable.t.loc[ixs_null,'obsnum'] = obsTable.t.loc[ix_obsTable,'obsnum'] 
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
            
        #Loop thhrough each propID
        propIDs = unique(obsTable.t['proposal_id'])
        for propID in propIDs:
            ixs_propID =  obsTable.ix_equal('proposal_id',propID)
            obsnums = unique(obsTable.t.loc[ixs_propID,'obsnum'])
            if self.verbose>2: 
                print('\n#### propID:',propID,' obsnums:',obsnums)
            # Loop through echo obsnum for the given propID, and make a summary entry
            for obsnum in obsnums:
                ixs_obsnum = obsTable.ix_equal('obsnum',obsnum,indices=ixs_propID)
                if self.verbose>2: 
                    obsTable.write(indices=ixs_obsnum)
                newline_obsnum={'proposal_id':propID}
                if obsnum is None:
                    newline_obsnum['obsnum']=None
                else:
                    newline_obsnum['obsnum']=obsnum
                    
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
        self.summary.t['obsnum'] = self.summary.t['obsnum'].astype('int')

        # sort summary table
        self.ix_summary_sorted = self.summary.ix_sort_by_cols(self.params['sortcols_summaryTable'])
   
        return(0)

if __name__ == '__main__':

    query = query_mast()
    parser = query.define_options()
    args = parser.parse_args()

    # the arguments are saved in query.params
    query.get_arguments(args)
    if query.verbose>2:
        print('params:', query.params)
        
    # get the MJD limits, based on --mjdlimits --lockbacktime --mjdmin --mjdmax --lre3 --lre4
    mjd_min, mjd_max = query.get_mjd_limits()

    # get the observations: query.obsTable
    query.observation_query(query.params['propID'], mjd_min = mjd_min, mjd_max = mjd_max)
    if query.verbose>0:
        print(query.obsTable.t)
        
    if len(query.obsTable.t)>0:
        # get the products: query.productTable
        # this list contains in general several entries for each observations.
        # If not otherwise specified, uses query.obsTable as starting point
        query.product_query()
        if query.verbose>1:
            print(query.productTable.t)
            #print(query.productTable.columns)
        
        # select the products to download
        # If not otherwise specified, uses query.productTable as starting point
        # uses query.filetypes for filtering, which is set by optional parameters.
        query.product_filter()
        
        # update obsTable with selected products: count the different filetypes, and also update the obsnum if None in obsTable
        query.update_obsTable_with_selectedProducts()

        # print the table to screen if either verbose or not downloading 
        if query.verbose>0 or (not args.download):
            print('\n#####################\n### Selected Products:\n#####################')
            query.productTable.write(indices=query.ix_selected_products,columns=query.params['outcolumns_productTable'])

            print('\n####################\n### Observations:\n####################')
            query.obsTable.write(columns=query.params['outcolumns_obsTable'],indices=query.ix_obs_sorted)

        # make summary table
        query.mk_summary_tables()
        if query.verbose>0 or (not args.download):
            print('\n##########################\n### Summary propID/obsnum:\n##########################')
            query.summary.write(indices=query.ix_summary_sorted)

        # save teh tables if wanted
        if query.params['savetables'] is not None:
            query.productTable.write(filename=query.params['savetables']+'.selprod.txt',indices=query.ix_selected_products,columns=query.params['outcolumns_productTable'],verbose=2)
            query.obsTable.write(filename=query.params['savetables']+'.obs.txt',columns=query.params['outcolumns_obsTable'],indices=query.ix_obs_sorted,verbose=2)
            query.summary.write(filename=query.params['savetables']+'.summary.txt',indices=query.ix_summary_sorted,verbose=2)
            

        sys.exit(0)

        if args.download:
            query.fetch_products()
    else:
        print('\n################################\nNO OBSERVATIONS FOUND! exiting....\n################################')
        

        

