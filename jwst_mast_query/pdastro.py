#!/usr/bin/env python
'''
wrapper around pandas with convenience functions to ease handling of tables
A. Rest
'''
import sys,os,re,types,copy,io
import numpy as np
from astropy.time import Time
import astropy.io.fits as fits
import pandas as pd
from pandas.core.dtypes.common import is_object_dtype,is_float_dtype,is_string_dtype,is_integer_dtype

from astropy.nddata import bitmask
from astropy import units as u
from astropy.coordinates import SkyCoord

from scipy.interpolate import interp1d

def split_commonpath(path1,path2, strip_basenames=False, raiseError=True):
    """
    compare two paths, and return the common path, and whatever is left from
    path1 that is not common

    Parameters
    ----------
    path1 : string
    path2 : string
    strip_basenames: if True, then the basenames are removed before the split
    raiseError: if True
    Returns
    -------
    (commonpath,not_common)

    """
    if strip_basenames:
        path1 = os.path.basename(path1)
        path2 = os.path.basename(path2)

    commonpath = os.path.commonpath([path1,path2])
    if len(commonpath)==0:
        return('',path1)
    not_common = path1[len(commonpath):]
    not_common = not_common.lstrip("/")
    commonpath = commonpath.rstrip("/")
    return(commonpath,not_common)

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

def rmfiles(filenames,raiseError=1,gzip=False):
    if not (type(filenames) is list):
        raise RuntimeError("List type expected as input to rmfiles")
    errorflag = 0
    for filename in filenames:
        errorflag |= rmfile(filename,raiseError=raiseError,gzip=gzip)
    return(errorflag)


#https://numpy.org/doc/stable/reference/routines.set.html
def AorB(A,B):
    if len(A) == 0:
        return(B)
    if len(B) == 0:
        return(A)
    return(np.union1d(A,B))

def AandB(A,B,assume_unique=False,keeporder=False):
    if keeporder:
        # This is slower, but keeps order
        out=[]
        for i in A:
            if i in B:
                out.append(i)
        return(out)
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

def order_A_like_B(A,B):
    """
    Parameters
    ----------
    A : list-like object
    B : list-like object

    Returns
    -------
    returns A but ordered like the following:
        all entries that are also in B come first, in the order of B
        all entries that are not in B come next, in the original order

    """
    tmp = AandB(A,B,keeporder=True)
    notB = AnotB(A,tmp,keeporder=True)
    tmp.extend(notB)
    return(tmp)

def unique(A):
    unique = []
    for a in A:
        if a not in unique:
            unique.append(a)
    return unique

def radec2coord(ra, dec):
    unit = [u.deg, u.deg]
    if ':' in str(ra):
        # Assume input RA/DEC are hour/degree
        unit[0] = u.hourangle

    try:
        coord = SkyCoord(ra, dec, frame='icrs', unit=unit)
        return(coord)
    except ValueError:
        print(f'ERROR: cannot interpret: RA={ra} DEC={dec}')
        return(None)



class pdastroclass:
    def __init__(self,hexcols=[],hexcol_formatters={},**kwargs):
        self.t = pd.DataFrame(**kwargs)

        self.verbose = 0

        # if self.auto_convert_dtypes==True, then before .to_string() is called in self.write, self.t.convert_dtypes() is run
        #self.auto_convert_dtypes = True

        # special cases: columns in hexadecimal format.
        self.hexcols=hexcols
        self.hexcol_formatters=hexcol_formatters

        # if a file is successfully loaded, the filename is saved in this varviable
        self.filename = None


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

        # add list of columns to be skipped when using write function
        self.skipcols = []

        # dictionary for the splines. arguments are the y columns of the spline
        self.spline={}


    def load_lines(self,lines,sep='\s+',**kwargs):
        #errorflag = self.load_spacesep(io.StringIO('\n'.join(lines)),sep=sep,**kwargs)
        errorflag = self.load_spacesep(io.StringIO('\n'.join(lines)),**kwargs)
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
                      hexcols=None,auto_find_hexcols=True, delim_whitespace=True,
                      na_values=['None','-','--'],verbose=False,**kwargs):

        #kwargs['delim_whitespace']=True

        #also test for commented header to make it compatible to old format.
        self.load(filename,na_values=na_values,test4commentedheader=test4commentedheader,
                  namesMapping=namesMapping,roundingMapping=roundingMapping,delim_whitespace=delim_whitespace,
                  hexcols=hexcols,auto_find_hexcols=auto_find_hexcols,verbose=verbose,**kwargs)

        return(0)

    def load(self,filename,raiseError=True,test4commentedheader=False,namesMapping=None,roundingMapping=None,
             hexcols=None,auto_find_hexcols=True,verbose=False,delim_whitespace=True,**kwargs):
        #self.t = ascii.read(filename,format='commented_header',delimiter='\s',fill_values=[('-',0),('--',0)])

        try:
            if verbose: print('Loading %s' % filename)
            self.t = pd.read_table(filename,delim_whitespace=delim_whitespace,**kwargs)
            self.filename = filename
        except Exception as e:
            print('ERROR: could not read %s!' % filename)
            if raiseError:
                raise RuntimeError(str(e))
            self.filename = None
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

    def write(self,filename=None,indices=None,columns=None,formatters=None,
              raiseError=True,overwrite=True,verbose=False,
              commentedheader=False,
              index=False, makepathFlag=True,convert_dtypes=False,
              hexcols=None,skipcols=None,
              return_lines=False,
              htmlflag=False, htmlsortedtable=False, **kwargs):

        # make sure indices are converted into a valid list
        indices=self.getindices(indices)

        # make sure columns are converted into a valid list
        columns=self.getcolnames(columns)
        if skipcols is None: skipcols = self.skipcols
        if len(skipcols)>0:
            columns=AnotB(columns,skipcols,keeporder=True)

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
                    if raiseError:
                        raise RuntimeError(f'file {filename} already exists! use overwrite option...')

                    print(f'Warning: file {filename} already exists, not deleting it, skipping! if you want to overwrite, use overwrite option!')
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

        # ugly hack: backwards compatibility to old commented header texttable format:
        # rename first column temporarily to have a '#'
        if commentedheader:
            columns=list(columns)
            self.t.rename(columns={columns[0]:f'#{columns[0]}'},inplace=True)
            columns[0]=f'#{columns[0]}'

        if len(indices)==0:
            if columns is None: columns = []
            # just save the header
            lines = [' '.join(columns)+'\n']
        else:
            if htmlflag:
                lines = self.t.loc[indices].to_html(index=index, columns=columns, formatters=formatters, **kwargs)
                if htmlsortedtable:
                    lines = re.sub('^\<table','<script type="text/javascript" src="sortable.js"></script>\n<table class="sortable" id="anyid" ',lines)
            else:
                lines = self.t.loc[indices].to_string(index=index, columns=columns, formatters=formatters, **kwargs)
        if filename is None:
            if not return_lines:
                print(lines)
        else:
            open(filename,'w').writelines(lines)


            """
            if filename is None:
                if return_lines:
                    return(0,self.t.loc[indices].to_string(index=index, columns=columns, formatters=formatters, **kwargs))
                else:
                    print(self.t.loc[indices].to_string(index=index, columns=columns, formatters=formatters, **kwargs))
            else:
                if not htmlflag:
                    self.t.loc[indices].to_string(filename, index=index, columns=columns, formatters=formatters, **kwargs)
                else:
                    if htmlsortedtable:
                        lines = self.t.loc[indices].to_html(index=index, columns=columns, formatters=formatters, **kwargs)
                        lines = re.sub('^\<table','<script type="text/javascript" src="sortable.js"></script>\n<table class="sortable" id="anyid" ',lines)
                        f = open(filename,"w")
                        f.writelines(lines)
                        f.close()
                    else:
                        self.t.loc[indices].to_html(filename, index=index, columns=columns, formatters=formatters, **kwargs)
                    #self.t.loc[indices].to_html(filename, index=index, columns=columns, formatters=formatters, **kwargs)
            """

        # reverse ugly hack
        if commentedheader:
            columns[0]=re.sub('^\#','',columns[0])
            self.t.rename(columns={"#"+columns[0]:columns[0]},inplace=True)

        if filename is not None:
            # some extra error checking...
            if not os.path.isfile(filename):
                errorstring='ERROR: could not save %s' % filename
                if raiseError:
                    raise RuntimeError(errorstring)
                #print(errorstring)
                return(3)
        if return_lines:
            return(0,lines)
        else:
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
        print('ix_remove_null deprecated, replace with ix_not_null')
        return(self.ix_not_null(colnames=colnames,indices=indices))

    def ix_not_null(self,colnames=None,indices=None):
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
        row_df = pd.DataFrame(dicti, index=[len(self.t)])
        self.t = pd.concat([self.t, row_df])
        return(self.t.index.values[-1])

    def add2row(self,index,dicti):
        self.t.loc[index,list(dicti.keys())]=list(dicti.values())
        return(index)

    def fitsheader2table(self,fitsfilecolname,indices=None,requiredfitskeys=None,optionalfitskey=None,
                         raiseError=True,verify='silentfix',skipcolname=None,headercol=None,ext=None,extname=None,
                         prefix=None,suffix=None):
        def fitskey2col(fitskey,prefix=None,suffix=None):
            col = fitskey
            if prefix is not None: col = prefix+col
            if suffix is not None: col += suffix
            return(col)


        indices = self.getindices(indices)

        # initialize columns if necessary
        if requiredfitskeys!=None:
            for fitskey in requiredfitskeys:
                col = fitskey2col(fitskey,prefix=prefix,suffix=suffix)
                if not (col in self.t.columns):
                    self.t[col]=None
        if optionalfitskey!=None:
            for fitskey in optionalfitskey:
                col = fitskey2col(fitskey,prefix=prefix,suffix=suffix)
                if not (col in self.t.columns):
                    self.t[col]=None

        if headercol!=None and (not (headercol in self.t.columns)):
            self.t[headercol]=None

        # loop through the images
        for index in indices:
            #header = fits.getheader(self.t.loc[index,fitsfilecolname],ext=ext,extname=extname)
            # It was impossible to use 'verify' with getheader...
            hdu = fits.open(self.t.loc[index,fitsfilecolname],ext=ext,extname=extname,output_verify="silentfix")
            if verify is not None: hdu.verify(verify)
            header = hdu[0].header
            if headercol!=None:
                self.t[headercol]=header

            if skipcolname!=None:
                self.t.loc[index,skipcolname]=False
            if requiredfitskeys!=None:
                for fitskey in requiredfitskeys:
                    col = fitskey2col(fitskey,prefix=prefix,suffix=suffix)
                    if fitskey in header:
                        self.t.loc[index,col]=header[fitskey]
                    else:
                        if raiseError:
                            raise RuntimeError("fits key %s does not exist in file %s" % (fitskey,self.t[fitsfilecolname][index]))
                        else:
                            self.t.loc[index,col]=None
                            if skipcolname!=None:
                                 self.t.loc[index,skipcolname]=True

            if optionalfitskey!=None:
                for fitskey in optionalfitskey:
                    col = fitskey2col(fitskey,prefix=prefix,suffix=suffix)
                    if fitskey in header:
                        self.t.loc[index,col]=header[fitskey]
                    else:
                        self.t.loc[index,col]=None

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

    def radeccols_to_SkyCoord(self,racol=None,deccol=None,indices=None, frame='icrs',
                              assert_0_360_limits=True,assert_pm90_limits=True):
        if (racol is None) and (deccol is None):
            raise RuntimeError('You need to specify at least one of racol or deccol')

        indices = self.getindices(indices)
        for col in [racol,deccol]:
            if col is None: continue
            ixs_null = self.ix_is_null(col,indices)
            if len(ixs_null)>0:
                self.write(indices=ixs_null)
                raise RuntimeError(f'Null values for column {col} in above row(s)')

        if len(indices)==0:
            print('Warning, trying to assert ra/dec columns are decimal for 0 rows')
            return([],[])

        # fill ra and dec
        ra = np.full(len(indices),np.nan)
        dec = np.full(len(indices),np.nan)
        if racol is not None: ra = np.array(self.t.loc[indices,racol])
        if deccol is not None: dec = np.array(self.t.loc[indices,deccol])

        # check if all cols are already in numerical format. If yes, all good! This can speed things up!!
        numflag=True
        for col in [racol,deccol]:
            if col is None: continue
            if not(col in self.t.columns):
                raise RuntimeError(f'Column {col} is not in columns {self.t.columns}')
            if not(is_float_dtype(self.t[col]) or is_integer_dtype(self.t[col])):
                numflag=False

        if numflag:
            # no conversion needed!
            # check ra dec limits if wanted
            if (ra is not None) and assert_0_360_limits:
                ra = np.where(ra<0.0,ra+360.0,ra)
                ra = np.where(ra>=360.0,ra-360.0,ra)
            if (dec is not None) and assert_pm90_limits:
                dec = np.where(dec<-90.0,dec+180.0,dec)
                dec = np.where(dec> 90.0,dec-180.0,dec)

            coord = SkyCoord(ra, dec, frame=frame, unit=(u.deg,u.deg))
        else:
            # convert the strings into decimal degrees
            unit = [u.deg,u.deg]

            check_line_by_line = False
            # check format of racol if necessary
            if racol is not None:
                hexagesimal = unique([(re.search(':',x) is not None) for x in ra])
                if len(hexagesimal)==1:
                    # if hexagesimal true, set unit of RA to hourangle!!
                    if hexagesimal[0]:
                        unit[0] = u.hourangle
                elif len(hexagesimal)==2:
                    # both formats? then check line by line!
                    print('Warning: it looks like there are inconsistent formats in RA column {racol}! checking line by line for sexagesimal')
                    check_line_by_line = True
                else:
                    raise RuntimeError('Something is wrong here when trying to determine if RA col {racol} is sexagesimal! ')


            if check_line_by_line:
                coord = np.full(len(indices),np.nan)
                hexagesimal = [(re.search(':',x) is not None) for x in ra]
                raunits=np.where(hexagesimal,u.hourangle,u.deg)
                for i in range(len(ra)):
                    coord[i] = SkyCoord(ra[i], dec[i], frame=frame, unit=(raunits[i],u.deg))
            else:
                coord = SkyCoord(ra, dec, frame=frame, unit=unit)
            # no need to check ra dec limits, already done in SkyCoord

        return(indices,coord)


    def assert_radec_cols_decimal_degrees(self,racol=None,deccol=None,
                                          outracol=None,outdeccol=None,
                                          indices=None,coordcol=None,
                                          assert_0_360_limits=True,
                                          assert_pm90_limits=True):

        (indices,coord) = self.radeccols_to_SkyCoord(racol=racol,deccol=deccol,indices=indices,
                                                     assert_0_360_limits=assert_0_360_limits,assert_pm90_limits=assert_pm90_limits)

        if outracol is None: outracol = racol
        if outdeccol is None: outdeccol = deccol
        if outracol is not None: self.t.loc[indices,outracol]=coord.ra.degree
        if outdeccol is not None: self.t.loc[indices,outdeccol]=coord.dec.degree
        if coordcol is not None:
            self.t.loc[indices,coordcol]=coord
            # Don't write coordcol
            if not (coordcol in self.skipcols): self.skipcols.append(coordcol)

        return(0)

    def assert_radec_cols_sexagesimal(self,racol=None,deccol=None,
                                      outracol=None,outdeccol=None,
                                      indices=None,coordcol=None,
                                      precision=3,
                                      assert_0_360_limits=True,
                                      assert_pm90_limits=True):

        (indices,coord) = self.radeccols_to_SkyCoord(racol=racol,deccol=deccol,indices=indices,
                                                     assert_0_360_limits=assert_0_360_limits,
                                                     assert_pm90_limits=assert_pm90_limits)

        if outracol is None: outracol = racol
        if outdeccol is None: outdeccol = deccol

        # list of ra/dec pairs
        hmsdms = map(lambda x: x.split(' '), coord.to_string(sep=':', style='hmsdms', precision=precision))
        # unzip pairs into a list of ra and a list of dec
        radeclist = list(zip(*hmsdms))

        if outracol is not None: self.t.loc[indices,outracol]=radeclist[0]
        if outdeccol is not None: self.t.loc[indices,outdeccol]=radeclist[1]
        if coordcol is not None:
            self.t.loc[indices,coordcol]=coord
            # Don't write coordcol
            if not (coordcol in self.skipcols): self.skipcols.append(coordcol)

        return(0)


    def radeccols_to_coord(self,racol,deccol,coordcol,indices=None,
                           assert_0_360_limits=True,assert_pm90_limits=True):

        (indices,coord) = self.radeccols_to_SkyCoord(racol=racol,deccol=deccol,indices=indices,
                                                     assert_0_360_limits=assert_0_360_limits,assert_pm90_limits=assert_pm90_limits)

        self.t.loc[indices,coordcol]=coord
        # Don't write coordcol
        if not (coordcol in self.skipcols): self.skipcols.append(coordcol)

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
            indices_validSN = self.ix_not_null('__tmp_SN',indices=indices)
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


class pdastrostatsclass(pdastroclass):
    def __init__(self,**kwargs):
        pdastroclass.__init__(self,**kwargs)
        self.reset()
        self.set_statstring_format()
        self.c4_smalln = [0.0, 0.0, 0.7978845608028654, 0.8862269254527579, 0.9213177319235613, 0.9399856029866251, 0.9515328619481445]

    def c4(self,n):
        #http://en.wikipedia.org/wiki/Unbiased_estimation_of_standard_deviation
        if n<=6:
            return(self.c4_smalln[n])
        else:
            return(1.0 - 1.0/(4.0*n) - 7.0/(32.0*n*n) - 19.0/(128.0*n*n*n))

    def reset(self):
        self.statparams = {}
        for k  in ['mean','mean_err','stdev','stdev_err','X2norm','ix_good','ix_clip']:
            self.statparams[k]=None
        for k  in ['Ngood','Nclip','Nchanged','Nmask','Nnan','converged','i']:
            self.statparams[k]=0
        self.calc_stdev_X2_flag = True

    def set_statstring_format(self,format_floats='{:f}',format_ints='{:d}',format_none='{}', format_X2norm='{:.2f}'):
        self.str_format1 = "mean:%s(%s) stdev:%s(%s) X2norm:%s Nchanged:%s Ngood:%s Nclip:%s" % (format_floats,format_floats,format_floats,format_floats,format_X2norm,format_ints,format_ints,format_ints)
        self.str_format_none = "mean:%s(%s) stdev:%s(%s) X2norm:%s Nchanged:%s Ngood:%s Nclip:%s" % (format_none,format_none,format_none,format_none,format_none,format_none,format_none,format_none)
        #self.str_format2 = "mean:%s(%s) stdev:%s X2norm:%s Nchanged:%s Ngood:%s Nclip:%s" % (format_floats,format_floats,format_floats,format_floats,format_ints,format_ints,format_ints)

    def statstring(self):
        if self.statparams['mean'] is None or self.statparams['stdev'] is None:
            formatstring = "WARNING! i:{:02d} "+ self.str_format_none
        else:
            formatstring = "i:{:02d} "+ self.str_format1

        s = formatstring.format(self.statparams['i'],self.statparams['mean'],self.statparams['mean_err'],
                                self.statparams['stdev'],self.statparams['stdev_err'],self.statparams['X2norm'],
                                self.statparams['Nchanged'],self.statparams['Ngood'],self.statparams['Nclip'])
        return(s)



    def calcaverage_errorcut(self,datacol, noisecol, indices=None,
                             mean=None,Nsigma=None,medianflag=False,
                             return_ix=False, verbose=0):

        # get the indices based on input.
        indices=self.getindices(indices)

        # If N-sigma cut and second iteration (i.e. we have a stdev from the first iteration), skip bad measurements.
        if not(Nsigma is None) and not(mean is None):
            ix_good_bkp = copy.deepcopy(self.statparams['ix_good'])
            (keep,) = np.where(np.absolute(self.t.loc[indices,datacol]-mean)<=Nsigma*self.t.loc[indices,noisecol])
            ix_good  = indices[keep]
        else:
            ix_good_bkp = None
            ix_good  = indices

        if verbose>3 and not(ix_good_bkp is None) and not(ix_good is None):
            print('{} good data after sigma clipping, {} clipped'.format(len(ix_good),len(indices)-len(ix_good)))
            self.write(indices=ix_good)
            #self.write(columns=['objID','filter','psfFlux','psfFluxErr','psfMag','psfMagErr','psfMagErr_tot'],indices=ix_good)

        Ngood = len(ix_good)

        if Ngood>1:
            if medianflag:
                mean = self.t.loc[ix_good,datacol].median()
                if verbose>1: print('median: {:f}'.format(mean))
                stdev =  np.sqrt(1.0/(Ngood-1.0)*np.sum(np.square(self.t.loc[ix_good,datacol] - mean)))/self.c4(Ngood)
                mean_err = self.t.loc[ix_good,noisecol].median()/np.sqrt(Ngood-1)
                #mean_err = stdev/np.sqrt(Ngood-1)

            else:
                c1 = np.sum(1.0*self.t.loc[ix_good,datacol]/np.square(self.t.loc[ix_good,noisecol]))
                c2 = np.sum(1.0/np.square(self.t.loc[ix_good,noisecol]))
                mean = c1/c2
                mean_err = np.sqrt(1.0/c2)
                stdev = self.t.loc[ix_good,datacol].std()

            stdev_err = 1.0*stdev/np.sqrt(2.0*Ngood)
            X2norm = 1.0/(Ngood-1.0)*np.sum(np.square((self.t.loc[ix_good,datacol] - mean)/self.t.loc[ix_good,noisecol]))

        else:
            if Ngood==1:
                mean = self.t.loc[ix_good[0],datacol]*1.0
                mean_err = self.t.loc[ix_good[0],noisecol]*1.0
            else:
                mean = None
                mean_err = None

            X2norm   = None
            stdev     = None
            stdev_err = None

        self.statparams['ix_good']=ix_good
        self.statparams['Ngood']=Ngood
        self.statparams['ix_clip']=AnotB(indices,ix_good)
        self.statparams['Nclip']=len(indices) - Ngood
        if not(ix_good_bkp is None):
            self.statparams['Nchanged'] = len(not_AandB(ix_good_bkp,ix_good))
        else:
            self.statparams['Nchanged'] = 0

        self.statparams['mean']      = mean
        self.statparams['stdev']     = stdev
        self.statparams['mean_err']  = mean_err
        self.statparams['stdev_err'] = stdev_err
        self.statparams['X2norm']    = X2norm

        if Ngood<1:
            return(1)
        return(0)

    def calcaverage_sigmacut(self,datacol, noisecol=None, indices=None,
                             mean=None,stdev=None,Nsigma=None,
                             percentile_cut=None,percentile_Nmin=3,
                             medianflag=False,
                             return_ix=False, verbose=0,rescol='__tmp_residuals'):

        # get the indices based on input.
        indices=self.getindices(indices)
        if len(indices)==0:
            print('Warning!! no data passed for sigma cut!')
            self.reset()
            return(2)

        ix_good_bkp = None
        if (percentile_cut is None) or (len(indices)<=percentile_Nmin):
            # If N-sigma cut and second iteration (i.e. we have a stdev from the first iteration), skip bad measurements.
            if not(Nsigma is None) and not(stdev is None) and not(mean is None):
                ix_good_bkp = copy.deepcopy(self.statparams['ix_good'])
                (keep,) = np.where(np.absolute(self.t.loc[indices,datacol]-mean)<=Nsigma*stdev)
                ix_good  = indices[keep]
                if verbose>3:
                    print('good data after sigma clipping:')
                    self.write(indices=ix_good)
                    #self.write(columns=['objID','filter','psfFlux','psfFluxErr','psfMag','psfMagErr'],indices=ix_good)
            else:
                if verbose>3:
                    print('No sigma clipping yet...')
                ix_good  = indices
        else:
            # percentile clipping!!!
            if mean is None:
                if medianflag:
                    median = self.t.loc[indices,datacol].median()
                    if verbose>1: print('median: {:f}'.format(median))
                    mean = median
                else:
                    mean = self.t.loc[indices,datacol].mean()
                    if verbose>1: print('mean: {:f}'.format(mean))

            self.t[rescol]=np.absolute(self.t.loc[indices,datacol]-mean)
            max_residual = np.percentile(np.absolute(self.t.loc[indices,datacol]-mean),percentile_cut)
            if verbose:
                print('%f percentile cut: max residual for cut: %f' % (percentile_cut,max_residual))
            ix_good  = self.ix_inrange(rescol,None,max_residual,exclude_uplim=True)


            #print('all')
            #self.write(columns=['objID','psfMag','psfMagErr',rescol],indices=indices)
            #print('good')
            #self.write(columns=['objID','psfMag','psfMagErr',rescol],indices=ix_good)

            # make sure the percentile clipping does not cut too much away for small numbers: always keep the best percentile_Nmin
            if len(ix_good)<percentile_Nmin:
                print('Warning: %d<%d made it through the percentile cut, too few, taking the %d with the least residuals' % (len(ix_good),percentile_Nmin,percentile_Nmin))
                residuals = np.sort(self.t.loc[indices,rescol])
                max_residual = residuals[percentile_Nmin-1]
                ix_good  = self.ix_inrange(rescol,None,max_residual,exclude_uplim=False)
                if len(ix_good)<percentile_Nmin:
                    raise RuntimeError('%d<%d in percentile cut!' % (len(ix_good),percentile_Nmin))

            if verbose>3:
                print('good data after percentile clipping:')
                self.write(indices=ix_good)
                #self.write(columns=['objID','filter','psfFlux','psfFluxErr','psfMag','psfMagErr',rescol],indices=ix_good)


            #sys.exit(0)

            #residuals = np.sort(np.absolute(self.t.loc[indices,datacol]-mean))
            #print(residuals)
            #max_residual = np.percentile(np.absolute(self.t.loc[indices,datacol]-mean),percentile_cut)

            #sys.exit(0)


        Ngood = len(ix_good)
        if Ngood>1:
            if medianflag:
                median = self.t.loc[ix_good,datacol].median()
                #mean = scipy.median(self.t.loc[ix_good,datacol])
                if verbose>1: print('median: {:f}'.format(median))
                stdev =  np.sqrt(1.0/(Ngood-1.0)*np.sum(np.square(self.t.loc[ix_good,datacol] - median)))/self.c4(Ngood)
                mean = median
            else:
                mean = self.t.loc[ix_good,datacol].mean()
                if verbose>1: print('mean: {:f}'.format(mean))
                stdev = self.t.loc[ix_good,datacol].std()

            mean_err = stdev/np.sqrt(Ngood-1)
            stdev_err = 1.0*stdev/np.sqrt(2.0*Ngood)
            if noisecol is None:
                X2norm = 1.0/(Ngood-1.0)*np.sum(np.square((self.t.loc[ix_good,datacol] - mean)/stdev))
            else:
                X2norm = 1.0/(Ngood-1.0)*np.sum(np.square((self.t.loc[ix_good,datacol] - mean)/self.t.loc[ix_good,noisecol]))
        else:
            if Ngood==1:
                mean = self.t.loc[ix_good[0],datacol]*1.0
                if noisecol is None:
                    mean_err = None
                else:
                    mean_err = self.t.loc[ix_good[0],noisecol]*1.0
            else:
                mean = None
                mean_err = None

            X2norm   = None
            stdev     = None
            stdev_err = None

        self.statparams['ix_good']=ix_good
        self.statparams['Ngood']=Ngood
        self.statparams['ix_clip']=AnotB(indices,ix_good)
        self.statparams['Nclip']=len(indices) - Ngood
        if not(ix_good_bkp is None):
            self.statparams['Nchanged'] = len(not_AandB(ix_good_bkp,ix_good))
        else:
            self.statparams['Nchanged'] = 0

        self.statparams['mean']      = mean
        self.statparams['stdev']     = stdev
        self.statparams['mean_err']  = mean_err
        self.statparams['stdev_err'] = stdev_err
        self.statparams['X2norm']    = X2norm

        if Ngood<1:
            return(1)
        return(0)

    def calcaverage_sigmacutloop(self,datacol, indices=None, noisecol=None, sigmacutFlag=False,
                                 maskcol=None, maskval=None,
                                 removeNaNs = True,
                                 Nsigma=3.0,Nitmax=10,verbose=0,
                                 percentile_cut_firstiteration=None,
                                 median_firstiteration=True):
        """
         median_firstiteration: in the first iteration, use the median instead the mean. This is more robust if there is a population of bad measurements
        """

        if noisecol is None:
            sigmacutFlag = True

        # get the indices based on input.
        indices=self.getindices(indices)

        self.reset()

        # exclude data if wanted
        if maskcol!=None:
            Ntot = len(indices)
            indices = self.ix_unmasked(maskcol,maskval=maskval,indices=indices)
            self.statparams['Nmask']= Ntot-len(indices)
            if verbose>1: print('Keeping {:d} out of {:d}, skippin {:d} because of masking in column {} (maskval={})'.format(len(indices),Ntot,Ntot-len(indices),maskcol,maskval))
        else:
            self.statparams['Nmask']= 0

        # remove null values if wanted
        if removeNaNs:
            colnames = [datacol]
            if not(noisecol is None): colnames.append(noisecol)
            if not(maskcol is None): colnames.append(maskcol)
            Ntot = len(indices)
            indices = self.ix_remove_null(colnames,indices=indices)
            self.statparams['Nnan']= Ntot-len(indices)
            if verbose>1: print('Keeping {:d} out of {:d}, skippin {:d} because of null values in columns {:s}'.format(len(indices),Ntot,Ntot-len(indices),",".join(colnames)))
        else:
            self.statparams['Nnan']= 0


        while ((self.statparams['i']<Nitmax) or (Nitmax==0)) and (not self.statparams['converged']):
            # median only in first iteration and if wanted
            medianflag = median_firstiteration and (self.statparams['i']==0) and (Nsigma!=None)
            # percentile_cut only in first interation and if wanted
            percentile_cut = None
            if (self.statparams['i']==0):
                percentile_cut = percentile_cut_firstiteration

            if sigmacutFlag:
                errorflag = self.calcaverage_sigmacut(datacol, indices=indices, noisecol=noisecol,
                                                      mean = self.statparams['mean'], stdev = self.statparams['stdev'],
                                                      Nsigma=Nsigma,
                                                      medianflag=medianflag, percentile_cut=percentile_cut,
                                                      verbose=verbose)
            else:
                errorflag = self.calcaverage_errorcut(datacol, noisecol, indices=indices,
                                                      mean = self.statparams['mean'],
                                                      Nsigma=Nsigma, medianflag=medianflag, verbose=verbose)

            #if self.statparams['i']==0:
            #    self.statparams['stdev']=0.05

            if verbose>2:
                print(self.statstring())

            # Not converged???
#            if errorflag or self.statparams['stdev']==None or self.statparams['stdev']==0.0 or self.statparams['mean']==None:
            if errorflag or self.statparams['stdev']==None or (self.statparams['stdev']==0.0 and sigmacutFlag) or self.statparams['mean']==None:
                self.statparams['converged']=False
                break
            # Only do a sigma cut if wanted
            if Nsigma == None or Nsigma == 0.0:
                self.statparams['converged']=True
                break
            # No changes anymore? If yes converged!!!
            if (self.statparams['i']>0) and (self.statparams['Nchanged']==0) and (not medianflag):
                self.statparams['converged']=True
                break
            self.statparams['i']+=1
            print()

        if not(self.statparams['converged']):
            if self.verbose>1:
                print('WARNING! no convergence!')

        return(not self.statparams['converged'])

    """
    def colnames4params(self,columns=None,colmapping={},prefix='',suffix='',skipcols=[]):
        cols=[]
        keys=[]
        for k in ['mean','mean_err','stdev','stdev_err','X2norm','Ngood','Nclip','Nmask','Nnan','converged','i']:
            if k in skipcols:
                continue
            outcol=k
            if k in colmapping:
                outcol=colmapping[k]
            keys.append(k)
            cols.append('{}{}{}'.format(prefix,outcol,suffix))
        return(keys,cols)
"""


    def intializecols4statparams(self,params=None,prefix='',suffix='',parammapping={},skipparams=[],
                                 addformat2defaultformatter=True,format4outvals='{:.3f}',
                                 setcol2None=True):
        outcols=[]
        outparams=[]

        if params is None:
            params=['mean','mean_err','stdev','stdev_err','X2norm','Ngood','Nclip','Nmask','Nnan','converged','i']

        if addformat2defaultformatter and (self.default_formatters is None):
                    self.default_formatters = {}

        for param in params:
            if param in skipparams:
                continue

            outparams.append(param)

            colbase=param
            if param in parammapping:
                colbase=parammapping[param]

            outcol='{}{}{}'.format(prefix,colbase,suffix)
            outcols.append(outcol)

            if addformat2defaultformatter:
                if not(re.search('^N',param) is None) or param =='i':
                    self.default_formatters[outcol]='{:d}'.format
                elif param=='X2norm':
                    self.default_formatters[outcol]='{:.2f}'.format
                elif param in ['converged']:
                    pass
                else:
                    self.default_formatters[outcol]=format4outvals.format
            if setcol2None:
                if param in ['mean','mean_err','stdev','stdev_err','X2norm']:
                    self.t[outcol]=np.nan
                else:
                    self.t[outcol]=None

        return(set(zip(outparams,outcols)))

    def statresults2table(self,statparams,param2columnmapping,destindex=None):
        resultdict={}
        #statparams = copy.deepcopy(statparams)
        # loop through keys,outcols, and assign values
        for (param,outcol) in param2columnmapping:
            if destindex is None:
                resultdict[outcol]:statparams[param]
            else:
                self.t.loc[destindex,outcol]=statparams[param]

        # either put the data into a new row or into an existing one
        if destindex is None:
            outindex = self.t.newrow(resultdict)
        else:
            outindex = destindex
            #self.t.loc[destindex,outcols]=vals
        return(outindex)

