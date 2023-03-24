#! /usr/bin/env python

"""Script that calls the jwst_download.py module. This is the top-level script that users should call.
It will parse arguments, make tables, query MAST, download files, and create webpages.
"""
import sys
from jwst_mast_query import jwst_download as jdl


def main():
    download = jdl.download_mast()
    parser = download.define_options()
    args = parser.parse_args()

    # The config file is loaded into self.params, and overwritten with the arguments that are not None
    download.get_arguments(args)
    if download.verbose>2:
        print('params:', download.params)

    # Use arguments or $API_MAST_TOKEN to login
    download.login(raiseErrorFlag=True)

    # self.outrootdir is set depending on outrootdir and outsubdir in cfg file or through the options --outrootdir and --outsubdir
    download.set_outrootdir()
    if download.verbose: print(f'Outdir: {download.outrootdir}')

    # Make the tables, but don't show them yet, since the output files need to be updated first
    if download.mk_all_tables(showtables=True):
        sys.exit(0)

    if not args.skipdownload and len(download.ix_selected_products)>0:
        download.download_products()
        print('\n######################\n### Downloaded Selected Products:\n######################')
        download.productTable.write(indices=download.ix_selected_products,columns=download.params['outcolumns_productTable'])
    else:
        if len(download.ix_selected_products)<=0:
            print('############## Nothing selected!!!! exiting...')
            sys.exit(0)

    # Make the webpages
    if args.makewebpages:
        download.mk_webpages4propIDs()


if __name__ == '__main__':
    main()