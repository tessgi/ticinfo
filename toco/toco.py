from __future__ import absolute_import
from . import logger
import argparse

from astroquery.mast import Catalogs
from astroquery.simbad import Simbad

import sys
import warnings
warnings.filterwarnings('ignore', category=UserWarning, append=True)


MAXTIC = 10000003844


class Target(object):

    def __init__(self, tic):
        self.validate_tic(tic)
        self.tic = int(tic)

    def validate_tic(self, tic):
        try:
            tic = int(tic)
        except:
            logger.error("Not a valid TIC number")
            sys.exit(1)

        if tic > MAXTIC:
            logger.error("Not a valid TIC number")
            sys.exit(1)

        if tic < 0:
            logger.error("Not a valid TIC number")
            sys.exit(1)

    def query(self):
        catalogData = Catalogs.query_criteria(catalog="Tic", ID=self.tic)
        return catalogData


def print_results(tic=12350):
    target = Target(tic)
    catalogData = target.query()

    print(catalogData[['ID', 'ra', 'dec', 'pmRA', 'pmDEC',
                       'eclong', 'eclat', 'Tmag', 'Vmag', 'Kmag', 'Teff',
                       'rad', 'mass', 'd', ]])
    print()


def toco(args=None):
    """
    exposes toco to the command line
    """
    if args is None:
        parser = argparse.ArgumentParser(
            description="Information for a TESS target")
        parser.add_argument('tic', type=float,
                            help="TICID of target")
        args = parser.parse_args(args)
        args = vars(args)
    tic = args['tic']

    _output = print_results(tic)
