# jwst_mast_query

This repository contains software for querying MAST for JWST observations and products, displaying the results in ascii tables, summarizing the files in a local html file, and downloading the selected products to a local directory.

`jwst_mast_query` uses the MAST API to query MAST. It first queries over the input time frame and retrieves all of the matching "observations". In this case the script uses a definition of "observation" that matches that in APT. An observation can be composed of one or more exposures, from one or more detectors of a single instrument. For each observation, the script identifies all of the individual "products" or files, to download. There can be many products for a given observation, including science data files in various stages of calibration, guide star observations, jpg files of the observations, etc. This means if the user asks for too many observations, the number of products can be so large that the MAST query times out. For example: The pre-launch LRE5 rehearsal has 264 NIRCam observations, which have in turn >17,000 products! About 12,000 of these are guide star products, about 5,000 are science observation products, and 500 are raw, uncalibrated science exposures. `jwst_mast_query` is able to download all NIRCam LRE5 products, but times out when downloading NIRSpec LRE5 products, as NIRSpec had more observations during the rehearsal. The easiest way to avoid time-outs is to remove unwanted products when defining the query.


## Installation

### Install the latest released version
`pip install jwst_mast_query`

### Install the development version
`pip install git+https://github.com/spacetelescope/jwst_mast_query.git`

### Server Installation
If you encounter the error below when running, try installing `pip install keyrings.alt`
and repeating the installation of jwst_mast_query.

> keyring.errors.NoKeyringError: No recommended backend was available. Install a recommended 3rd party backend package; or, install the keyrings.alt package if you want to use the non-recommended backends. See https://pypi.org/project/keyring for details.


## Environment Variables
These environment variables are not required, but they make life easier:

- **MAST_API_TOKEN**: Set this variable to your MAST token, and you get automatically logged in. You can also pass your token with --token. NOTE: if you don't use your MAST token for 10 days, it will become invalid, and you have to get a new one from here: https://auth.mast.stsci.edu/token
- **JWST_QUERY_CFGFILE**: Set this variable to your config file (yaml), and it gets loaded automatically. There is a default config file jwst_query.cfg. The config file can also be supplied at run time using --config and a path to a local config file.

- **JWSTDOWNLOAD_OUTDIR**: The yaml config file accepts environment variable in the format $XYZ.

- In the default config file, "outrootdir: $JWSTDOWNLOAD_OUTDIR". This allows different people to use the same config file, but store the images in different locations.


## Use

The jwst_mast_query package contains two tools, both designed to be called from the command line. The first is **jwst_download.py**. This script will query MAST for files matching the input parameters, and then download those files into a local output directory specified by the user. There is also an option to download associated jpgs from MAST and create a simple index.html file with a table containing these preview images and basic observation information. If requested, the script will also save ascii tables with information on the files identified by the query.

The other script is **jwst_query.py**. This script will query MAST for files matching the input parameters. **jwst_download.py** wraps around **jwst_query.py** and provides the abiltiy to download the files identified by **jwst_query.py**.  

> **NOTE: jwst_download.py** will check to see if the if the files returned by the MAST query already exist in the specified output directory. If they do, it will not download the files. It checks on a file-by-file basis, and also checks if the file is complete. Only missing or incomplete files will be downloaded, meaning that you can safely run **jwst_download.py** multiple times with the same query, and you will end up only downloading missing/new data, saving time.

## Inputs

Inputs to the command line call of **jwst_download.py** are optional, but we recommend at a minimum using a configuration file to define input parameters. To get started, download/copy the [configuration file](https://github.com/spacetelescope/jwst_mast_query/blob/main/jwst_mast_query/jwst_query.cfg) from the repository. For a detailed description of the options in the configuration file, see the **Config File** section below.

Below we show an example of a typical call to **jwsst_download.py**, with a few options specified. In this case, we specify verbose mode (`-v`) in order to get more details printed to the screen as the command runs. In addition, we want to locate all data from JWST proposal 1410 (`--propID`). We specify the name of the config file using `--config`. Since there is no path given, it is assumed that the config file is in the current working directory. We set the `--lookbacktime` to 3 days. This means that **jwst_download.py** will only search for observations taken within the last 3 days. Finally, we request data only from NIRCam's NRCA1 and NRCA2 detectors using the `--sca` option.

`jwst_download.py -v --propID 1410 --config jwst_query.cfg --lookbacktime 3 --sca a1 a2`

### Common command line options

- --lookbacktime : The number of days before the present to use as the beginning of the query

- --propID : The JWST proposal number. e.g. `--propID 1409` or `--propID 01409`

- --obsnums : Optional observation number or numbers within propID to retrieve. e.g. `--obsnums 3 103` would retrieve files only from observations 3 and 103.

- --makewebpages : If set, an index.html file is generated, containing info and images of the retrieved data. Note that this option is not currently in the config file, and must be specified at the command line.


These are just a few of the options that can be set. The rest are specified in the config file provided in the call. Config file details are given below.
For more example calls, see the **Examples** section below.


### Config File and Input Options


#### Hierarchy
Most config parameters get their values from one of three places: defaults, the config file, and command line arguments. 

When **jwst_download.py** or **jwst_query.py** are run, first the config file is read in, and all parameters in the config file are saved in self.params, overwriting the existing default values. Then, any command line arguments that are provided are added to self.params, overwriting the config file parameters.

#### Config file

The config file is a convenient way to set and keep track of parameter values. The table below lists the contents of the config file and defines each parameter.


| Parameter Name: default value                                                                       | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
|-----------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| instrument: nircam                                                                                  | Specify the instrument (nircam, nirspec, niriss, miri, fgs)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
| propID                                                                                  | Specify the JWST proposal number to query for. This can be a 5 digit integer, or for a smaller number, an integer with or without prepended zeros. e.g. 1409 or 01409.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            |
| obsnums                                                                                  | Specify the observation numbers within the propID. This can be a single number, or a bracketed list of numbers. e.g. 3 or [3, 103]                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            |
| outrootdir: $JWSTDOWNLOAD_OUTDIR                                                                    | The base directory for the downloaded products. This can be the name of an environment variable (such as JWSTDOWNLOAD_OUTDIR in this example), or a path. Note the preceding "$" in the case where an environment variable is given. If you include the --outrootdir command line argument when calling jwst_query.py or jwst_download.py, that value will override the value provided here.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            |
| outsubdir:                                                                                          | Any additional directory to add to the base directory. This can be used to customize the organization of the downloaded products.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| skip_propID2outsubdir: False                                                                        | Don't use the propID in the output directory. This can be used for custom output directory, e.g., to put the products from several APT files into one directory                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |
| obsnum2outsubdir: True                                                                              | If True, the observation number will be used to create a subdirectory into which the appropriate files will be placed. e.g. $JWSTDOWNLOAD_OUTDIR/01410/obsnum23/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        |
| propIDs_obsnum2outsubdir: [1409]                                                                    | If True, only for the listed proposal numbers will observation numbers will be used to create a subdirectory into which the appropriate files will be placed. e.g. $JWSTDOWNLOAD_OUTDIR/01410/obsnum23/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| skip_check_if_outfile_exists: False                                                                 | The script queries MAST for products, and then checks if each file already exists in the output directory or not ("dl_code" and "dl_str" columns). For large numbers of products, this can take time. By setting this option to True, this check can be skipped.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        |
| Nobs_per_batch: 2                                                                                   | For large programs, the query for the products can time out, and in these cases it is better to split up the query into Nobs_per_batch observations per batch. For example, if there are 22 observations, and Nobs_per_batch=4, then there will be 6 batches, 5 batches with 4 observations, and the last batch with 2.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        |
| obsmode: ['image', 'wfss']                                                                          | Optionally specify the observation modes to query for. If provided, MAST will be queried only for these types of observations. e.g. assuming instrument is 'nircam', obsmode: ['image', 'wfss'] will query MAST for NIRCAM/IMAGE and NIRCAM/WFSS data. If left empty, the query will include all modes. For a list of valid modes, see the example config file in the repository. Note that as of 24 March 2023, MAST is in the process of reprocessing all data to add the mode values to the instrument names. Until this process is complete, some older data may not include the mode name, and therefore will not be found via a query that includes the mode name.                                                                                                                                                                                                                                                                                                                                                                      |
| filetypes: ['uncal']                                                                                | List of file types to select in the product table, e.g., \_uncal.fits or \_uncal.jpg. If no suffix is given, .fits is appended. If only letters, then \_ and .fits are added. For example, 'uncal' gets expanded to \_uncal.fits. Typical image filetypes are uncal, rate, rateints, cal. For downloading a single file type, the brackets must still surround the file suffix, as the script expects a list. A relatively complete list of options includes: ['\_segm.fits', '\_asn.json', '\_pool.csv', '\_i2d.jpg', '\_thumb.jpg', '\_cat.ecsv', '\_i2d.fits', '\_uncal.fits', '\_uncal.jpg', '\_cal.fits', '\_trapsfilled.fits', '\_cal.jpg', '\_rate.jpg', '\_rateints.jpg', '\_trapsfilled.jpg', '\_rate.fits', '\_rateints.fits'] See the [JWST calibration pipeline documentation](https://jwst-pipeline.readthedocs.io/en/latest/jwst/introduction.html#pipeline-step-suffix-definitions) for a complete list. |
| guidestars: False                                                                                   | If guidestars is set to True, guidestar products are also included. Note: there are a **lot** of guide star products. We recommend you set to True only if really needed!                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| lookbacktime: 1.0                                                                                   | Lookback time in days. The script will query MAST over a time from the lookback time to the present moment. Note that all other time parameters (date_select, etc) override the lookback time.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| mastcolumns_obsTable: ['proposal_id', 'dataURL', 'obs_id', 't_min','t_exptime']                     | Core columns returned from MAST to the obsTable                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |
| outcolumns_productTable                                                                             | List of columns to be shown in product table, e.g., ['proposal_id', 'obsnum', 'obsID', 'parent_obsid', 'obs_id', 'dataproduct_type', 'productFilename', 'filetype', 'calib_level', 'size', 'outfilename', 'dl_code', 'dl_str']                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| outcolumns_obsTable: ['proposal_id', 'obsnum', 'obsid', 'obs_id', 't_min', 't_exptime', 'date_min'] | Output columns for the obsTable.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        |
| sortcols_productTable: ['calib_level','filetype','obsID']                                           | The productTable is sorted based on these columns. The default sorts the table based on calibration level.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
| sortcols_obsTable: ['date_min','proposal_id','obsnum']                                              | The obsTable is sorted based on these columns. The defaults sort the table in the order the observations were observed                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |
| sortcols_summaryTable: ['date_start','proposal_id','obsnum']                                        | The summary table is sorted based on these columns. The default sorts the table in the order the observations were observed                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
| date_select: []                                                                                     | Specify date range (MJD or isot format) applied to "dateobs_center" column. If single value, then only exact matches will be returned. If a single value has "+" or "-" at the end, then it is a lower and upper limit, respectively. Examples: 58400+, 58400-, 2020-11-23+,2020-11-23 2020-11-25                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| savetables:                                                                                         | Save the tables (selected products, obsTable, summary with suffix selprod.txt, obs.txt, summary.txt, respectively) with the specified string as basename. Filenames will then be *xxxx.summary.txt*, *xxxx.obs.txt*, and *xxxx.selprod.txt*, where xxxx is the savetables value. Tables are saved in the same output directory as the data.  If no string is provided, the tables are not saved.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     |
| webpage_thumbnail_width: 100                                                                        | Width in pixels of the resized jpg images to be inserted into the index.html summary file. |
| webpage_thumbnail_height:                                                                           | Height in pixels of the resized jpg images to be inserted into the index.html summary file. If left undefined, the height will be determined from webpage_thumbnail_width and the aspect ratio of the original image.              |
| webpage_level12_jpgs: ['\_uncal.jpg','\_dark.jpg','\_rate.jpg','\_rateints.jpg','\_trapsfilled.jpg','\_cal.jpg','\_crf.jpg'] | Filetypes whose thumbnails will be shown in index.html. |


## Outputs

The primary output of **jwst_download.py** are the downloaded files themselves. By default, the downloaded files are saved into the directory:

\<outrootdir\>\<outsubdir\>\<proposal number\>obsnum\<XX\>

where:

proposal number is the 5-digit, zero-padded APT number. i.e. the proposal ID.
outsubdir is an optional user-specified additional subdirectory.
XX is the observation number, as specified in the proposal.


In addition to querying MAST and downloading the selected files, `jwst_mast_query` saves several ASCII files containing tables with details of the files. These are the Summary table, the obsTable, and the productTable. 

By default, the summary table will be saved in the output directory with the name *.summary.txt* This table contains a high-level summary of the data identified in the query. The example below shows that the table contains the propposal ID, observation number, number of uncal.fits files found for each observation, and the date of the beginning of the observation.

    proposal_id  obsnum  \_uncal.fits         date_start
        1410       1            1 2022-02-13T22:58:32.378
        1410       3            1 2022-02-14T20:22:58.543
        1410      14            1 2022-02-14T20:34:26.670


The obsTable contains data about each observation found to have data matching the query. In the example below we see that this table contains the individual observation IDs for the files, as well as the exposure time for each.

    proposal_id  obsnum  obsid                obs_id                  t_min      t_exptime       date_min              \_uncal.fits
        1410       1    71672220 jw01410001001_02101_00001_guider1 59623.957319    161.052   2022-02-13T22:58:32.378         1
        1410      14    71673297 jw01410014001_02101_00001_guider1 59624.857253    161.052   2022-02-14T20:34:26.670         1
        1410       3    71673298 jw01410003001_02101_00001_guider1 59624.849289    161.052   2022-02-14T20:22:58.543         1


The table of selected products gives even more details about each downloaded product. In the example below we see there is information on the type of each data product, the calibration pipelne level through which the file has been run, and the location to which the file has been saved. 

    proposal_id obsnum   obsID  parent_obsid                   obs_id             sca     dataproduct_type    filetype  calib_level size                       outfilename                                  dl_code  dl_str
       1410       1    71672220  71672220     jw01410001001_02101_00001_guider1 guider1      image          \_uncal.fits     1    100713600 /jwst_data/01410/jw01410001001_02101_00001_guider1_uncal.fits        0     NaN
       1410       3    71673298  71673298     jw01410003001_02101_00001_guider1 guider1      image          \_uncal.fits     1    100713600 /jwst_data/01410/jw01410003001_02101_00001_guider1_uncal.fits        0     NaN
       1410      14    71673297  71673297     jw01410014001_02101_00001_guider1 guider1      image          \_uncal.fits     1    100713600 /jwst_data/01410/jw01410014001_02101_00001_guider1_uncal.fits        0     NaN



## Examples

Get all NIRCam NRCA1 and NRCA2 files for proposal 1410 taken in the last 3 days

    jwst_download.py -v --propID 1410 --config jwst_query.cfg --lookbacktime 3 --sca a1 a2


### Specify a proposal ID and specific observation numbers

Download the fits and jpg files for JWST proposal 1138, taken in the last 1 day and create an index.html summary file.

    jwst_download.py -v -c jwst_query.cfg --outrootdir /jwst_data -l 1 --propID 01138 --makewebpages --filetypes jpg fits

Download the fits files for JWST proposal 1409, observations 3 and 103 only, taken in the last 2 days.

    jwst_download.py -v --config jwst_query.cfg --lookbacktime 2 --propID 1409 --obsnums 3 103

Download the fits and jpg files for JWST proposal 1138, observation 4 only, taken in the last 1 day and create an index.html summary file.

    jwst_download.py -v -c jwst_query.cfg --outrootdir /jwst_data -l 1 --propID 01138 --obsnums 4 --makewebpages --filetypes jpg fits

Download the fits and jpg files for JWST proposal 1138, observations 4, 5, and 7 only, taken in the last 1 day and create an index.html summary file.

    jwst_download.py -v -c jwst_query.cfg --outrootdir /jwst_data -l 1 --propID 01138 --obsnums 4 5 7 --makewebpages --filetypes jpg fits


### Specify dates

Get all files for proposal 743 with an observation date between Aug 11, 2021 16:49:49 and Aug 12, 2021 16:49:49

    jwst_download.py -v  -c jwst_query.cfg --propID 743 --date_select 2021-08-11T16:49:49 2021-08-12T16:49:49

Get all files for proposal 743 with an observation date of Aug 11, 2021 16:49:49 or later

    jwst_download.py -v  -c jwst_query.cfg --propID 743 --date_select 2021-08-11T16:49:49+

Get all files for proposal 743 with an observation date of Aug 11, 2021 16:49:49 or earlier

    jwst_download.py -v  -c jwst_query.cfg --propID 743 --date_select 2021-08-11T16:49:49-

Get all files for proposal 743 with an observation date of MJD 59430.0 or later

    jwst_download.py -v  -c jwst_query.cfg --propID 743 --date_select 59430.0+
 
Download only the jpg (not the fits) files from the last 5 days, and create a table of results, saved into index.html

    jwst_download.py -v -c jwst_query.cfg --outrootdir /my_jwst_data -lookbacktime 5 --makewebpages --filetype jpg



[//]: # (jwst_download.py --lre5 -v  -c jwst_query.cfg --propID 743 --obsid_select 55783016+   But no one will know this obsid for their project)

[//]: # (jwst_download.py --lre5 -v  -c jwst_query.cfg --propID 743 --obsid_select 55783016 55783018  But no one will know this obsid for their project)

[//]: # (jwst_download.py --lre5 -v  -c jwst_query.cfg --propID 743 --obsid_list 55783016 55783018  But no one will know this obsid for their project)

