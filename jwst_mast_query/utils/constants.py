"""Globally defined and used variables for the ``jwst_mast_query`` project.

"""

MAST_OBS_MODES = {'miri': ['MIRI/CORON', 'MIRI/IFU', 'MIRI/IMAGE', 'MIRI/SLIT', 'MIRI/SLITLESS', 'MIRI/TARGACQ'],
                  'nircam': ['NIRCAM/CORON', 'NIRCAM/GRISM', 'NIRCAM/IMAGE', 'NIRCAM/TARGACQ', 'NIRCAM/WFSS'],
                  'niriss': ['NIRISS/AMI', 'NIRISS/IMAGE', 'NIRISS/SOSS', 'NIRISS/WFSS'],
                  'nirspec': ['NIRSPEC/BOTS', 'NIRSPEC/IFU', 'NIRSPEC/MSA', 'NIRSPEC/SLIT']
                  }
