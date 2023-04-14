"""Globally defined and used variables for the ``jwst_mast_query`` project.

"""
import os


MAST_OBS_MODES = {'miri': ['MIRI/CORON', 'MIRI/IFU', 'MIRI/IMAGE', 'MIRI/SLIT', 'MIRI/SLITLESS', 'MIRI/TARGACQ'],
                  'nircam': ['NIRCAM/CORON', 'NIRCAM/GRISM', 'NIRCAM/IMAGE', 'NIRCAM/TARGACQ', 'NIRCAM/WFSS'],
                  'niriss': ['NIRISS/AMI', 'NIRISS/IMAGE', 'NIRISS/SOSS', 'NIRISS/WFSS'],
                  'nirspec': ['NIRSPEC/BOTS', 'NIRSPEC/IFU', 'NIRSPEC/MSA', 'NIRSPEC/SLIT']
                  }

# Default for config file, if available
if 'JWST_QUERY_CFGFILE' in os.environ and os.environ['JWST_QUERY_CFGFILE'] != '':
    cfgfilename = os.environ['JWST_QUERY_CFGFILE']
else:
    cfgfilename = None

# Default for token is $MAST_API_TOKEN
if 'MAST_API_TOKEN' in os.environ:
    defaulttoken = os.environ['MAST_API_TOKEN']
else:
    defaulttoken = None

PARAM_DEFAULTS = {'mastcolumns_obsTable': ['proposal_id','dataURL','obsid','obs_id','t_min','t_exptime','instrument_name'],
                  'outcolumns_productTable': ['proposal_id','obsnum','obsID','parent_obsid','obs_id','sca','visit','dataproduct_type','productFilename','filetype','calib_level','size','description'],
                  'outcolumns_obsTable': ['proposal_id','obsnum','obsid','obs_id','t_min','t_exptime','date_min','instrument_name'],
                  'sortcols_productTable': ['calib_level','filetype','obsID'],
                  'sortcols_obsTable': ['date_min','proposal_id','obsnum'],
                  'sortcols_summaryTable': ['date_start','proposal_id','obsnum'],
                  'configfile': cfgfilename,
                  'token': defaulttoken,
                  'instrument': 'nircam',
                  'obsmode': None,
                  'propID': None,
                  'obsnums': None,
                  'obsid_select': [],
                  'obsid_list': [],
                  'sca': None,
                  'verbose': 0,
                  'filetypes': ['uncal'],
                  'guidestars': False,
                  'lookbacktime': 1,
                  'calib_levels': None,
                  'Nobs_per_batch': 2,
                  'obsnum2outsubdir': None,
                  'propIDs_obsnum2outsubdir': [],
                  'date_select': [],
                  'savetables': None,
                  'makewebpages': False,
                  'outrootdir': None,
                  'outsubdir': None,
                  'skip_propID2outsubdir': False,
                  'skip_check_if_outfile_exists': False,
                  'skip_check_filesize': False,
                  'jpg_separate_subdir': False,
                  'webpage_tablefigsize_width': None,
                  'webpage_tablefigsize_height': None,
                  'webpage_level12_jpgs': ['_uncal.jpg', '_dark.jpg', '_rate.jpg', '_rateints.jpg', '_trapsfilled.jpg', '_cal.jpg', '_crf.jpg'],
                  'webpage_fitskeys2table': ['TARG_RA', 'TARG_DEC', 'FILTER', 'PUPIL', 'READPATT', 'NINTS', 'NGROUPS', 'NFRAMES', 'DATE-BEG',
                                             'DATE-END', 'EFFINTTM', 'EFFEXPTM'],
                  'webpage_cols4table': ['proposal_id', 'obsnum', 'visit', 'obsID', 'parent_obsid', 'sca', 'FILTER', 'PUPIL', 'READPATT', 'uncal',
                                         'dark', 'rate', 'rateints', 'cal', 'TARG_RA', 'TARG_DEC', 'NINTS', 'NGROUPS', 'NFRAMES', 'DATE-BEG', 'DATE-END',
                                         'EFFINTTM', 'EFFEXPTM', 'size', 'obs_id', 'outfilename'],
                  'webpage_sortcols': ['proposal_id', 'obsnum', 'visit', 'sca'],
                  'webpage_mkthumbnails': True,
                  'webpage_thumbnails_overwrite': False,
                  'webpage_thumbnails_width': 120,
                  'webpage_thumbnails_height': None,
                  'login': None,
                  'skipdownload': False
                  }