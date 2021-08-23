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
from astropy.nddata import bitmask
import astroquery
import argparse,os,sys,re,io
from scipy.interpolate import interp1d


# MAST API documentation:
# https://mast.stsci.edu/api/v0/pyex.html

# MAST CAOM Field Descriptions:
# https://mast.stsci.edu/api/v0/_c_a_o_mfields.html

# File naming schemes in this page
# https://jwst-pipeline.readthedocs.io/en/latest/jwst/data_products/file_naming.html
# https://jira.stsci.edu/browse/JSDP-1778


if astroquery.__version__<'0.4.2':
    raise RuntimeError("astroquery version=%s, at least 0.4.2 required!" % astroquery.__version__)


def makepath(path,raiseError=True):
    if path == '':
        return(0)
    if not os.path.isdir(path):
        os.makedirs(path)
        if not os.path.isdir(path):
            if raiseError:
                raise RuntimeError('ERROR: Cannot create directory %s' % path)
            else:
                return(1)
    return(0)

def makepath4file(filename,raiseError=True):
    path = os.path.dirname(filename)
    if not os.path.isdir(path):
        return(makepath(path,raiseError=raiseError))
    else:
        return(0)
    
#https://numpy.org/doc/stable/reference/routines.set.html
def AorB(A,B):
    if len(A) == 0:
        return(B)
    if len(B) == 0:
        return(A)
    return(np.union1d(A,B))

def AandB(A,B,assume_unique=False):
    return(np.intersect1d(A,B,assume_unique=assume_unique))

def AnotB(A,B,keeporder=False):
    if keeporder:
        # This is slower, but keeps order
        out=[]
        for i in A: 
            if not(i in B):
                out.append(i)
        return(out)
    return(np.setdiff1d(A,B))

def not_AandB(A,B):
    return(np.setxor1d(A,B))

def unique(A):
    unique = []
    for a in A:
        if a not in unique:
            unique.append(a)
    return unique


class pdastroclass:
    def __init__(self,hexcols=[],hexcol_formatters={},**kwargs):
        self.t = pd.DataFrame(**kwargs)
   
        self.verbose = 0
        
        # if self.auto_convert_dtypes==True, then before .to_string() is called in self.write, self.t.convert_dtypes() is run 
        #self.auto_convert_dtypes = True

        # special cases: columns in hexadecimal format. 
        self.hexcols=hexcols
        self.hexcol_formatters=hexcol_formatters


        # if self.default_dtypeMapping != None, then the dtype mapping is applied to the table .to_string() is called in self.write
        # This makes sure that formatters don't throw errors if the type of a column got changed to float or object during
        # one of the table operations
        self.default_dtypeMapping = None
        # example:
        # self.default_dtypeMapping = {'counter':np.int64}
        
        # default_formatters are passed to to_string() in self.write
        self.default_formatters = {}
        # example:
        # self.default_formatters = {'MJD':'{:.6f}'.format,'counter':'{:05d}'.format}
       
        # dictionary for the splines. arguments are the y columns of the spline
        self.spline={}


    def load_lines(self,lines,sep='\s+',**kwargs):
        errorflag = self.load_spacesep(io.StringIO('\n'.join(lines)),sep=sep,**kwargs)
        return(errorflag)

    def load_cmpfile(self,filename,**kwargs):
        """
        load old frankenstein format of cmp file 

        Parameters
        ----------
        filename : 
            cmp filename.

        Returns
        -------
        errorflag and fits header

        """
        cmphdr = fits.getheader(filename)
        s = ''
        for i in range(1,int(cmphdr['NCOLTBL'])+1):
            s+=' %s' % cmphdr['COLTBL%d' % i]
        s = re.sub('Xpos','X',s)
        s = re.sub('Ypos','Y',s)
        print(s)
        lines = open(filename,'r').readlines()
        lines[0]=s
        errorflag = self.load_spacesep(io.StringIO('\n'.join(lines)))
        return(errorflag,cmphdr)
        
    def load_spacesep(self,filename,test4commentedheader=True,namesMapping=None,roundingMapping=None,
                      hexcols=None,auto_find_hexcols=True,
                      na_values=['None','-','--'],verbose=False,**kwargs):
        
        kwargs['delim_whitespace']=True

        #also test for commented header to make it compatible to old format.
        self.load(filename,na_values=na_values,test4commentedheader=test4commentedheader,
                  namesMapping=namesMapping,roundingMapping=roundingMapping,
                  hexcols=hexcols,auto_find_hexcols=auto_find_hexcols,verbose=verbose,**kwargs)

        return(0)

    def load(self,filename,raiseError=True,test4commentedheader=False,namesMapping=None,roundingMapping=None,
             hexcols=None,auto_find_hexcols=True,verbose=False,**kwargs):
        #self.t = ascii.read(filename,format='commented_header',delimiter='\s',fill_values=[('-',0),('--',0)])

        try:
            if verbose: print('Loading %s' % filename)
            self.t = pd.read_table(filename,**kwargs)
        except Exception as e:
            print('ERROR: could not read %s!' % filename)
            if raiseError:
                raise RuntimeError(str(e))
            return(1)
        
        if test4commentedheader:
            # This is to make it compatible to my old-style commented header files!
            # commented header: make sure it doesn't count the '#' as a column!
            if self.t.columns[0]=='#':
                renamemapping = {}
                for i in range(len(self.t.columns)-1):
                    renamemapping[self.t.columns[i]]=self.t.columns[i+1]
                renamemapping[self.t.columns[-1]]='__delme'
                self.t = self.t.rename(columns=renamemapping)
                self.t = self.t.drop(columns=['__delme'])

            # make sure the # is not kept in column name!
            if self.t.columns[0][0]=='#':
                self.t = self.t.rename(columns={self.t.columns[0]:self.t.columns[0][1:]})
            
        if hexcols is None: hexcols=[]
        hexcols.extend(self.hexcols)
        self.formattable(namesMapping=namesMapping,roundingMapping=roundingMapping,hexcols=hexcols,auto_find_hexcols=auto_find_hexcols)
        
        
        return(0)

    def write(self,filename=None,indices=None,columns=None,formatters=None,raiseError=True,overwrite=True,verbose=False, 
              index=False, makepathFlag=True,convert_dtypes=False,hexcols=None,**kwargs):

        # make sure indices are converted into a valid list
        indices=self.getindices(indices)
        
        # make sure columns are converted into a valid list
        columns=self.getcolnames(columns)

        # make the path to the file if necessary
        if not (filename is None):
            if makepathFlag:
                if makepath4file(filename,raiseError=False):
                    errorstring='ERROR: could not make directory for %s' % filename
                    if raiseError:
                        raise RuntimeError(errorstring)
                    #print(errorstring)
                    return(1)

            # if overwrite, then remove the old file first
            if os.path.isfile(filename):
                if overwrite:
                    os.remove(filename)
                    if os.path.isfile(filename):
                        errorstring='ERROR: could not clobber %s' % filename
                        if raiseError:
                            raise RuntimeError(errorstring)
                        #print(errorstring)
                        return(2)
                else:
                    print('Warning: file exists, not deleting it, skipping! if you want to overwrite, use overwrite option!')
                    return(0)
        
        # Fix the dtypes if wanted. THIS IS DANGEROUS!!!
        if convert_dtypes:
            # Using pandas convert_dtypes converts the columns into panda types, e.g. Int64. These types
            # are different than the numpy types (e.g., int64). The numpy types work with '{}'.format, the panda
            # int types do not!!! Therefore I set convert_dtype to False by default.
            self.t=self.t.convert_dtypes()
        if not(self.default_dtypeMapping is None):
            self.formattable(dtypeMapping=self.default_dtypeMapping)            
        
        # if both formatters and self.defaultformatters are None, then no formatting. formatters supersedes self.defaultformatters
        if formatters is None:
            formatters = self.default_formatters

        # hexcols
        if hexcols is None: hexcols=[]
        hexcols.extend(self.hexcols)
        if len(hexcols)>0:
            if formatters is None: formatters={}
            for hexcol in self.hexcols:
                if not(hexcol in self.t.columns):
                    #Nothing to do yet!
                    continue
                if not (hexcol in formatters):
                    if isinstance(self.t[hexcol].dtype,np.int64):
                        formatters[hexcol]='0x{:08x}'.format
                    else:
                        formatters[hexcol]='0x{:04x}'.format
            
        if verbose>1 and not(filename is None): print('Saving %d rows into %s' % (len(indices),filename))
        if len(indices)==0:
            # just save the header
            if filename is None:
                print(' '.join(columns)+'\n')
            else:
                if columns is None:
                    columns = []
                open(filename,'w').writelines(' '.join(columns)+'\n')
        else:
            if filename is None:
                print(self.t.loc[indices].to_string(index=index, columns=columns, formatters=formatters, **kwargs))
            else:
                self.t.loc[indices].to_string(filename, index=index, columns=columns, formatters=formatters, **kwargs)

        if not (filename is None):
            # some extra error checking...
            if not os.path.isfile(filename):
                errorstring='ERROR: could not save %s' % filename
                if raiseError:
                    raise RuntimeError(errorstring)
                #print(errorstring)
                return(3)
                
        return(0)
        
    def formattable(self,namesMapping=None,roundingMapping=None,dtypeMapping=None,hexcols=None,auto_find_hexcols=False):
        
        if not(namesMapping is None):
            self.t = self.t.rename(columns=namesMapping)
            
        if not(roundingMapping is None):
            self.t = self.t.round(roundingMapping)

        if not(dtypeMapping is None):
            for col in dtypeMapping:
                self.t[col] = self.t[col].astype(dtypeMapping[col])
        

        # hexadecimal columns
        hexpattern = re.compile('^0x[a-fA-F0-9]+$')
        # if wanted, find columns that are in hexadecimal string format, and convert it interger. 
        # These columns are added to self.hexcols
        if auto_find_hexcols and len(self.t)>0:
            for col in self.t.columns:
                if self.t[col].dtype == object and isinstance(self.t.at[0,col],str):
                    if not(hexpattern.search(self.t.at[0,col]) is None):
                        # add column to self.hexcols, so that it stays a hexcol when writing and saving 
                        if not(col in self.hexcols):
                            self.hexcols.append(col)
                        # convert it to int
                        self.t[col] = self.t[col].apply(int, base=16)
        # Go through hexcols, and check if they need conversion.
        # These columns are also added to self.hexcols
        if not(hexcols is None):
            for hexcol in hexcols:
                # add column to self.hexcols, so that it stays a hexcol when writing and saving 
                if not(hexcol in self.hexcols):
                    self.hexcols.append(hexcol)
                if not(hexcol in self.t.columns):
                    # nothing to do (yet)!
                    continue
                # if the column is still a string starting with '0x', convert it to int
                if self.t[hexcol].dtype == object and isinstance(self.t.at[0,hexcol],str):
                    if not(hexpattern.search(self.t.at[0,hexcol]) is None):
                         self.t[hexcol] = self.t[hexcol].apply(int, base=16)
               
                
            
        return(0)
    
    def getindices(self,indices=None):
        """make indices conform (input can be None,([list],), int, str, or list). The output is a list """
        
        #If indices is None, return all values
        if indices is None: 
            return(self.t.index.values)
        
        # If indices=([indiceslist],), then it needs to be changed to indices=indiceslist
        if isinstance(indices,tuple):
            if len(indices)==0:
                #return an empty list instead of empty tuple
                return([])
            # teh first entry is a list, then these are the relevant indices!
            if isinstance(indices[0],list) or isinstance(indices[0],np.ndarray):
                return(indices[0])
            else:
                return(list(indices))
        
        # if the passed value is an integer or str, make it a list!
        if isinstance(indices,int) or isinstance(indices,str) or isinstance(indices,float):
            return([indices])           
        
        indices=np.array(indices)
        
        #if not (isinstance(indices,list) or isinstance(indices,np.ndarray)):
        #    raise RuntimeError("Can't convert this to an indices list!",type(indices),indices)
            
        return(indices)            

    def getcolnames(self,colnames=None):
        """Return a list of all colnames of colnames=None. If colnames=string, return a list"""
        if (colnames is None):
             colnames = self.t.columns[:] 
        elif isinstance(colnames,str):
            if colnames.lower()=='all':
                colnames = self.t.columns[:]
            else:
                colnames=[colnames]
        return(colnames)
            

    def ix_remove_null(self,colnames=None,indices=None):
        # get the indices based on input.
        indices=self.getindices(indices)
        
        # get the column names over which to iterate
        colnames=self.getcolnames(colnames)
        
        for colname in colnames:
            #print('XXX',indices)
            (notnull,) = np.where(pd.notnull(self.t.loc[indices,colname]))
            indices = indices[notnull]
            #print('YYY',notnull)
        return(indices)

    def ix_is_null(self,colnames=None,indices=None):
        # get the indices based on input.
        indices=self.getindices(indices)
        
        # get the column names over which to iterate
        colnames=self.getcolnames(colnames)
        
        for colname in colnames:
            #print('XXX',indices)
            (null,) = np.where(pd.isnull(self.t.loc[indices,colname]))
            indices = indices[null]
            #print('YYY',null)
        return(indices)

    def ix_equal(self,colnames,val,indices=None):
        # get the indices based on input.
        
        # use isnull() if val is None
        if val is None:
            indices = self.ix_is_null(colnames,indices=indices)
            return(indices)

        indices=self.getindices(indices)        
        
        # get the column names over which to iterate
        colnames=self.getcolnames(colnames)
        for colname in colnames:
            (keep,) = np.where(self.t.loc[indices,colname].eq(val))
            indices = indices[keep]
            
        return(indices)
        
    def ix_inrange(self,colnames=None,lowlim=None,uplim=None,indices=None,
                   exclude_lowlim=False,exclude_uplim=False):

        # get the indices based on input.
        indices=self.getindices(indices)
        
        # get the column names over which to iterate
        colnames=self.getcolnames(colnames)
        #print(colnames)
        for colname in colnames:
            if not(lowlim is None):
                if exclude_lowlim:
                    (keep,) = np.where(self.t.loc[indices,colname].gt(lowlim))
                else:
                    (keep,) = np.where(self.t.loc[indices,colname].ge(lowlim))
                indices = indices[keep]
                #print('lowlim cut:',keep)

            if not(uplim is None):
                if exclude_uplim:
                    (keep,) = np.where(self.t.loc[indices,colname].lt(uplim))
                else:
                    (keep,) = np.where(self.t.loc[indices,colname].le(uplim))
                indices = indices[keep]
                #print('uplim cut:',keep)
        return(indices)
    
    def ix_outrange(self,colnames=None,lowlim=None,uplim=None,indices=None,
                    exclude_lowlim=False,exclude_uplim=False):

        # get the indices based on input.
        indices=self.getindices(indices)        
        
        # get the column names over which to iterate
        colnames=self.getcolnames(colnames)
        
        #print('BBB',indices)
        for colname in colnames:
            if not(lowlim is None):
                if exclude_lowlim:
                    (keeplow,) = np.where(self.t.loc[indices,colname].lt(lowlim))
                else:
                    (keeplow,) = np.where(self.t.loc[indices,colname].le(lowlim))
                #print('lowlim cut:',keeplow)
            else:
                keeplow=[]

            if not(uplim is None):
                if exclude_uplim:
                    (keepup,) = np.where(self.t.loc[indices,colname].gt(uplim))
                else:
                    (keepup,) = np.where(self.t.loc[indices,colname].ge(uplim))
                #print('uplim cut:',keepup)
            else:
                keepup=[]
                
            indices = indices[AorB(keeplow,keepup)]
            
        return(indices)
    
    def ix_unmasked(self,maskcol,maskval=None,indices=None):

        # get the indices based on input.
        indices=self.getindices(indices)  
        
        if maskval is None:
            (keep,) = np.where(self.t.loc[indices,maskcol].eq(0))
        else:
            (keep,) = np.where(bitmask.bitfield_to_boolean_mask(self.t.loc[indices,maskcol].astype('int'),ignore_flags=~maskval,good_mask_value=True))
        indices = indices[keep]
        return(indices)           
    
    def ix_masked(self,maskcol,maskval=None,indices=None):

        # get the indices based on input.
        indices=self.getindices(indices)  
        
        if maskval is None:
            (keep,) = np.where(self.t.loc[indices,maskcol].ne(0))
        else:
            (keep,) = np.where(bitmask.bitfield_to_boolean_mask(self.t.loc[indices,maskcol].astype('int'),ignore_flags=~maskval))
        indices = indices[keep]
        return(indices)    
    
    def ix_matchregex(self,col,regex,indices=None):
        # get the indices based on input.
        indices=self.getindices(indices)  
        
        (keep,) = np.where((self.t.loc[indices,col].str.contains(regex)==True))
        #bla = self.t[col].str.contains(regex)==True

        indices = indices[keep]
        return(indices)    

    def ix_sort_by_cols(self,cols,indices=None):

        # get the indices based on input.
        indices=self.getindices(indices)  
        
        # get the column names (makes sure that it is a list)
        cols=self.getcolnames(cols)

        ix_sorted = self.t.loc[indices].sort_values(cols).index.values

        return(ix_sorted)

    def newrow(self,dicti=None):
        self.t = self.t.append(dicti,ignore_index=True)
        return(self.t.index.values[-1])
        
    def add2row(self,index,dicti):
        self.t.loc[index,list(dicti.keys())]=list(dicti.values())
        return(index)

    def fitsheader2table(self,fitsfilecolname,indices=None,requiredfitskeys=None,optionalfitskey=None,raiseError=True,skipcolname=None,headercol=None,ext=None,extname=None):

        indices = self.getindices(indices)        

        # initialize columns if necessary
        if requiredfitskeys!=None:
            for fitskey in requiredfitskeys:
                if not (fitskey in self.t.columns):
                    self.t[fitskey]=None
        if optionalfitskey!=None:
            for fitskey in optionalfitskey:
                if not (fitskey in self.t.columns):
                    self.t[fitskey]=None

        if headercol!=None and (not (headercol in self.t.columns)):
            self.t[headercol]=None

        # loop through the images
        for index in indices:
            header = fits.getheader(self.t.loc[index,fitsfilecolname],ext=ext,extname=extname)
            if headercol!=None:
                self.t[headercol]=header
                
            if skipcolname!=None:
                self.t.loc[index,skipcolname]=False
            if requiredfitskeys!=None:
                for fitskey in requiredfitskeys:
                    if fitskey in header:
                        self.t.loc[index,fitskey]=header[fitskey]
                    else:
                        if raiseError:
                            raise RuntimeError("fits key %s does not exist in file %s" % (fitskey,self.t[fitsfilecolname][index]))
                        else:
                            self.t.loc[index,fitskey]=None
                            if skipcolname!=None:
                                 self.t.loc[index,skipcolname]=True
                                 
            if optionalfitskey!=None:
                for fitskey in optionalfitskey:
                    if fitskey in header:
                        self.t.loc[index,fitskey]=header[fitskey]
                    else:
                        self.t.loc[index,fitskey]=None

    def dateobs2mjd(self,dateobscol,mjdcol,timeobscol=None,indices=None,tformat='isot'):
        indices = self.getindices(indices)  
        if len(indices)==0:
            return(0)

        if not (mjdcol in self.t.columns):
            self.t.loc[indices,mjdcol]=None

        if timeobscol!=None:
            dateobslist = list(self.t.loc[indices,dateobscol]+'T'+self.t.loc[indices,timeobscol])
        else:
            dateobslist = list(self.t.loc[indices,dateobscol])

        dateobjects = Time(dateobslist,  format=tformat, scale='utc')
        mjds = dateobjects.mjd
        self.t.loc[indices,mjdcol]=mjds
        
    def mjd2dateobs(self,mjdcol,dateobscol,indices=None,tformat='isot'):
        indices = self.getindices(indices)  
        if len(indices)==0:
            return(0)

        if not (dateobscol in self.t.columns):
            self.t.loc[indices,dateobscol]=None


        mjdlist = list(self.t.loc[indices,mjdcol])

        dateobjects = Time(mjdlist, format='mjd')
        dateobs = dateobjects.to_value('isot')
        self.t.loc[indices,dateobscol]=dateobs
        
    def calc_color(self,f1,df1,f2,df2,outcolor=None,outcolor_err_nameformat='e(%s)',indices=None,color_formatter='{:.3f}'.format):
        indices = self.getindices(indices)  
        if len(indices)==0:
            return(0)
        
        # get color and uncertainty column names
        if outcolor is None: outcolor = '%s-%s' % (f1,f2)
        outcolor_err = outcolor_err_nameformat % outcolor
        # delete all previous results
        self.t.loc[indices,outcolor]= np.nan
        self.t.loc[indices,outcolor_err]= np.nan

        # Only use indices for which both filters have uncertainties, i.e. they are not upper limits
        ix_good = self.ix_remove_null(df1,indices=indices)
        ix_good = self.ix_remove_null(df2,indices=ix_good)

        self.t[outcolor] = self.t.loc[ix_good,f1] - self.t.loc[ix_good,f2]
        self.t[outcolor_err] = np.sqrt(np.square(self.t.loc[ix_good,df1]) + np.square(self.t.loc[ix_good,df2]))
        
        if not(color_formatter is None):
            if self.default_formatters is None:
                self.default_formatters = {}
            self.default_formatters[outcolor]=color_formatter
            self.default_formatters[outcolor_err]=color_formatter
        
        return(0)

    
    def flux2mag(self,fluxcol,dfluxcol,magcol,dmagcol,indices=None,
                 zpt=None,zptcol=None, upperlim_Nsigma=None):

        indices = self.getindices(indices)  
        if len(indices)==0:
            return(0)
        
        # delete all previous results
        self.t.loc[indices,magcol]= np.nan
        self.t.loc[indices,dmagcol]= np.nan

        # if upperlim_Nsigma is not None, then upper limits are calculated. 
        if upperlim_Nsigma is None:
            indices_mag = indices
            # no upper limits
            indices_ul = indices_ul_negative = []
        else:
           # calculate the S/N
            self.t.loc[indices,'__tmp_SN']=self.t.loc[indices,fluxcol]/self.t.loc[indices,dfluxcol]
            # Get the indices with valid S/N
            indices_validSN = self.ix_remove_null('__tmp_SN',indices=indices)
 
            # get the indices for which the S/N>=upperlim_Nsigma
            indices_mag = self.ix_inrange(['__tmp_SN'],upperlim_Nsigma,indices=indices_validSN)
            # all the other indices can only be used for upper limits
            indices_ul = AnotB(indices_validSN,indices_mag) 
            # get the indices for which the S/N is negative
            indices_ul_negative = self.ix_inrange(['__tmp_SN'],None,0.0,indices=indices_ul)
            #self.t = self.t.drop(columns=['__tmp_SN'])
            
        # Calculate the mags and uncertainties
        self.t.loc[indices_mag,magcol]= -2.5*np.log10(self.t.loc[indices_mag,fluxcol])
        self.t.loc[indices_mag,dmagcol] = 2.5 / np.log(10.0) * self.t.loc[indices_mag,dfluxcol]/self.t.loc[indices_mag,fluxcol]

        # Are there upper limits to be calculated? 
        if len(indices_ul)>0:
            # just to be sure, set dmags to nan
            self.t.loc[indices_ul,dmagcol] = np.nan
            
            indices_ul_positive = AnotB(indices_ul,indices_ul_negative) 
            if len(indices_ul_positive)>0:
                # If flux is positive: upper limit = flux + Nsigma * dflux
                self.t.loc[indices_ul_positive,magcol] = -2.5*np.log10(self.t.loc[indices_ul_positive,fluxcol]+upperlim_Nsigma*self.t.loc[indices_ul_positive,dfluxcol])

            if len(indices_ul_negative)>0:
                # If flux is negative: upper limit = Nsigma * dflux
                self.t.loc[indices_ul_negative,magcol] = -2.5*np.log10(upperlim_Nsigma*self.t.loc[indices_ul_negative,dfluxcol])
        
        if not(zpt is None):        
            self.t.loc[indices,magcol]+=zpt
        if not(zptcol is None):        
            self.t.loc[indices,magcol]+=self.t.loc[indices,zptcol]


    def initspline(self,xcol,ycol,indices=None,
                   kind='cubic',bounds_error=False,fill_value='extrapolate', 
                   **kwargs):
        if not (xcol in self.t.columns):
            raise RuntimeError("spline: x column %s does not exist in table!" % xcol)
        if not (ycol in self.t.columns):
            raise RuntimeError("spline: y column %s does not exist in table!" % ycol)

        # make sure there are no nan values
        indices = self.ix_remove_null(colnames=[xcol,ycol])        

        # initialize the spline and save it in self.spline with the key ycol
        self.spline[ycol]= interp1d(self.t.loc[indices,xcol],self.t.loc[indices,ycol],
                                    kind=kind,bounds_error=bounds_error,fill_value=fill_value,**kwargs)
        
    def getspline(self,xval,ycol):
        if not(ycol in self.spline):
            raise RuntimeError('Spline for column %s is not defined!' % ycol)
        return(self.spline[ycol](xval))
    
    def plot(self,indices=None,*args,**kwargs):
        indices = self.getindices(indices)
        ax = self.t.loc[indices].plot(*args,**kwargs)        
        return(ax)



class jwstMASTclass:
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
        
        
        self.RET_COLUMNS = ['proposal_id','dataURL','obsid','obs_id','t_min','t_exptime']
        self.SERVICES = {
                'SI_search': 'Mast.Jwst.Filtered.',
                'Caom_search':'Mast.Caom.Filtered.JwstOps',
                'Product_search':'Mast.Caom.Products.JwstOps'
                }
        
        # self.params will be populated with the arguments
        self.params = {}
        
        # delete all tables and outdir
        self.reset()
        
    def reset(self):
        self.outdir = None    
        
        # observation and product table
        self.obsTable = pdastroclass()
        self.productTable = pdastroclass()
        # summary table: statistics/info for each propID/obsnum pair
        self.summary = pdastroclass()

        # indices of the selected products
        # They are sorted by self.sortcols_productTable
        self.ix_selected_products = None
        self.sortcols_productTable = ['calib_level','filetype','obsID']
    
        # indices for obsTable, sorted by date_min
        self.ix_obs_sorted = None
        # indices for summary table, sorted by date_start
        self.ix_summary_sorted = None

        # output columns for the tables. Note that the columns for the individual filetypes
        # are automatically added to the obsTable
        self.outcolumns_productTable = ['proposal_id','obsnum','obsID','parent_obsid','obs_id','dataproduct_type','productFilename','filetype','calib_level','size','description']
        self.outcolumns_obsTable = ['proposal_id','obsnum','obsid','obs_id','t_min','t_exptime','date_min']




    def define_options(self,parser=None,usage=None,conflict_handler='resolve'):
        if parser is None:
            parser = argparse.ArgumentParser(usage=usage,conflict_handler=conflict_handler)

        # options for config file
        if 'JWSTDOWNLOAD_OUTDIR' in os.environ and os.environ['JWSTDOWNLOAD_OUTDIR'] != '':
            outrootdir = os.environ['JWSTDOWNLOAD_OUTDIR']
        else:
            outrootdir = '.'
        #print('Default output rootdir: %s' % outrootdir)
    
        parser.add_argument('-v','--verbose', default=0, action='count')
        parser.add_argument('--propID', type=str, default=None, help='Download data for this program ID.')

        parser.add_argument('-o','--outrootdir', default=outrootdir, help=('output root directory.''(default=%(default)s)'))
        parser.add_argument('--outsubdir', default=None, help=('subdir added to the output root directory. If None, then propID is used (default=%(default)s)'))

        parser.add_argument('-f','--filetypes',  type=str, nargs="+", default=['uncal'], help=('List of product filetypes to get, e.g., _uncal.fits or _uncal.jpg. If only letters, then _ and .fits are added, for example uncal gets expanded to _uncal.fits. Typical image filetypes are uncal, rate, rateints, cal (default=%(default)s)'))

        parser.add_argument('--clobber', action='store_true', default=False, help='existing files are overwritten, otherwise they are skipped (default=%(default)s)')
        parser.add_argument('-d','--download', action='store_true', default=False, help='download the data.  (default=%(default)s)')
        parser.add_argument('--guidestars', action='store_true', default=False, help='Don\'t skip guidestsars  (default=%(default)s)')

        parser.add_argument('-m', '--mjd_limits', type=float, nargs=2, help='specify the MJD limits. overrides lookback time and mjd_min/max optional arguments(default=%(default)s)')
        parser.add_argument('-l', '--lookbacktime', type=float, default=15, help='lookback time in days (default=%(default)s)')
        parser.add_argument('--mjd_min', type=float, default=None, help='minimum MJD. overrides lookback time  (default=%(default)s)')
        parser.add_argument('--mjd_max', type=float, default=None, help='maximum MJD. (default=%(default)s)')
        
        parser.add_argument('--lre3', action='store_true', default=False, help='Use the LRE-3 mjd limits (default=%(default)s)')
        parser.add_argument('--lre4', action='store_true', default=False, help='Use the LRE-4 mjd limits (default=%(default)s)')
        parser.add_argument('--lre5', action='store_true', default=False, help='Use the LRE-4 mjd limits (default=%(default)s)')
        
        parser.add_argument('-i', '--instrument', type=str, default='nircam', choices=['niriss','nircam','nirspec','miri','fgs'], help='Instrument.  (default=%(default)s)')

        parser.add_argument('-s', '--savetables', type=str, default=None, help='save the tables (selected products, obsTable, summary with suffix selprod.txt, obs.txt, summary.txt, respectively) with the specified string as basename (default=%(default)s)')

        
        return(parser)

    def get_arguments(self, args): 
        '''

        Parameters
        ----------
        args : list
            pass the command line arguments to self.params.
            make sure propID has the correct format (5 digit string)

        Returns
        -------
        None.

        '''
        self.params.update(vars(args))
        if self.params['propID'] is not None:
            self.params['propID'] = '%05d' % (int(self.params['propID']))
        self.verbose = self.params['verbose']
        
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
         
        columns = ','.join(self.RET_COLUMNS)
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

        self.ix_obs_sorted = self.obsTable.ix_sort_by_cols(['date_min'])


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
        
        self.ix_selected_products = self.productTable.ix_sort_by_cols(self.sortcols_productTable,indices=self.ix_selected_products)

        
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
            if filetype not in self.outcolumns_obsTable:
                self.outcolumns_obsTable.append(filetype)

        for filetype in filetypes:
            obsTable.t[filetype]=None

        ixs_obsTable = obsTable.getindices()
        for ix_obsTable in ixs_obsTable:
            obsid = obsTable.t.loc[ix_obsTable,'obsid']
            if self.verbose>2: 
                print('### obsid',obsid)
            ixs_prodTable = productTable.ix_equal('parent_obsid',obsid,indices=ix_selected_products)
            if self.verbose>2: 
                productTable.write(indices=ixs_prodTable)
                
            for ix_prodTable in ixs_prodTable:
                if obsTable.t.loc[ix_obsTable,productTable.t.loc[ix_prodTable,'filetype']] is None:
                    obsTable.t.loc[ix_obsTable,productTable.t.loc[ix_prodTable,'filetype']] = 0
                obsTable.t.loc[ix_obsTable,productTable.t.loc[ix_prodTable,'filetype']] += 1
                
            if obsTable.t.loc[ix_obsTable,'obsnum'] is None:
                obsnum = unique(productTable.t.loc[ixs_prodTable,'obsnum'])
                if len(obsnum)==0:
                    obsTable.write(indices=ix_obsTable)
                    productTable.write(indices=ixs_prodTable)
                    raise RuntimeError('No obsnum for these observations and products!')
                if len(obsnum)>1:
                    print('More than one obsnum!',obsnum)
                    raise RuntimeError('BUGGG!!!!!')
                obsTable.t.loc[ix_obsTable,'obsnum'] = obsnum
            
    def mk_summary_tables(self, obsTable=None, filetypes=None):
        if obsTable==None:
            obsTable=self.obsTable

        if filetypes==None:
            filetypes=self.params['filetypes']
            
        columns_summary = ['proposal_id','obsnum']
        columns_summary.extend(filetypes)
        columns_summary.extend(['date_start'])
        self.summary = pdastroclass(columns=columns_summary)
            
        propIDs = unique(obsTable.t['proposal_id'])
        for propID in propIDs:
            ixs_propID =  obsTable.ix_equal('proposal_id',propID)
            obsnums = unique(obsTable.t.loc[ixs_propID,'obsnum'])
            if self.verbose>2: 
                print('\n#### propID:',propID,' obsnums:',obsnums)
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
                
                
        
        #self.summary.write()

        #make the type integer for the filetype columns
        for filetype in filetypes:
            self.summary.t[filetype] = self.summary.t[filetype].astype('int')
        #ixs = self.summary.ix_is_null('obsnum')
        #self.summary.t.loc[ixs,'obsnum'] = None
        #self.summary.t['obsnum'] = self.summary.t['obsnum'].astype('int')
        self.summary.t['proposal_id'] = self.summary.t['proposal_id'].astype('int')
        self.summary.t['obsnum'] = self.summary.t['obsnum'].astype('int')

        self.ix_summary_sorted = self.summary.ix_sort_by_cols(['date_start'])
   
        return(0)

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

    jwstMAST = jwstMASTclass()
    parser = jwstMAST.define_options()
    args = parser.parse_args()

    # the arguments are saved in jwstMAST.params
    jwstMAST.get_arguments(args)
    if jwstMAST.verbose>2:
        print('params:', jwstMAST.params)
        
    # set the outdir, based on arguments --outrootdir, --outsubdir. If no outsubdir is specified, the propID is used.
    jwstMAST.set_outdir()

    # get the MJD limits, based on --mjdlimits --lockbacktime --mjdmin --mjdmax --lre3 --lre4
    mjd_min, mjd_max = jwstMAST.get_mjd_limits()

    # get the observations: jwstMAST.obsTable
    jwstMAST.observation_query(jwstMAST.params['propID'], mjd_min = mjd_min, mjd_max = mjd_max)
    if jwstMAST.verbose>0:
        print(jwstMAST.obsTable.t)
        
    if len(jwstMAST.obsTable.t)>0:
        # get the products: jwstMAST.productTable
        # this list contains in general several entries for each observations.
        # If not otherwise specified, uses jwstMAST.obsTable as starting point
        jwstMAST.product_query()
        if jwstMAST.verbose>1:
            print(jwstMAST.productTable.t)
            #print(jwstMAST.productTable.columns)
        
        # select the products to download
        # If not otherwise specified, uses jwstMAST.productTable as starting point
        # uses jwstMAST.filetypes for filtering, which is set by optional parameters.
        jwstMAST.product_filter()
        
        # update obsTable with selected products: count the different filetypes, and also update the obsnum if None in obsTable
        jwstMAST.update_obsTable_with_selectedProducts()

        # print the table to screen if either verbose or not downloading 
        if jwstMAST.verbose>0 or (not args.download):
            print('\n#####################\n### Selected Products:\n#####################')
            jwstMAST.productTable.write(indices=jwstMAST.ix_selected_products,columns=jwstMAST.outcolumns_productTable)

            print('\n####################\n### Observations:\n####################')
            jwstMAST.obsTable.write(columns=jwstMAST.outcolumns_obsTable,indices=jwstMAST.ix_obs_sorted)

            
            
        # make summary table
        jwstMAST.mk_summary_tables()
        if jwstMAST.verbose>0 or (not args.download):
            print('\n##########################\n### Summary propID/obsnum:\n##########################')
            jwstMAST.summary.write(indices=jwstMAST.ix_summary_sorted)

        # save teh tables if wanted
        if jwstMAST.params['savetables'] is not None:
            jwstMAST.productTable.write(filename=jwstMAST.params['savetables']+'.selprod.txt',indices=jwstMAST.ix_selected_products,columns=jwstMAST.outcolumns_productTable,verbose=2)
            jwstMAST.obsTable.write(filename=jwstMAST.params['savetables']+'.obs.txt',columns=jwstMAST.outcolumns_obsTable,indices=jwstMAST.ix_obs_sorted,verbose=2)
            jwstMAST.summary.write(filename=jwstMAST.params['savetables']+'.summary.txt',indices=jwstMAST.ix_summary_sorted,verbose=2)
            

        sys.exit(0)

        if args.download:
            jwstMAST.fetch_products()
    else:
        print('\n################################\nNO OBSERVATIONS FOUND! exiting....\n################################')
        

        

