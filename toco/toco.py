from __future__ import absolute_import
from . import logger
import argparse

from astropy.coordinates import SkyCoord, get_constellation
from astroquery.mast import Catalogs
from astroquery.simbad import Simbad
from astropy import units as u

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
        catalogData = Catalogs.query_object('TIC' + str(self.tic), catalog="TIC")
        return catalogData


def print_results(tic=12350, simbad_search=False):
    target = Target(tic)
    catalogData = target.query()

    print(catalogData[['ID', 'ra', 'dec', 'pmRA', 'pmDEC',
                       'eclong', 'eclat', 'Tmag', 'Vmag', 'Kmag', 'Teff',
                       'rad', 'mass', 'd', ]][0])

    if simbad_search:
        skobj = SkyCoord(ra=catalogData['ra'] * u.degree,
                         dec=catalogData['dec'] * u.degree,
                         frame='icrs')
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            result_table = Simbad.query_region(skobj, radius=20 * u.arcsec)

        if result_table is None:
            logger.warning("No Simbad target resolved")
        else:
            print()
            print('Target name: {}'.format(
                result_table['MAIN_ID'][0].decode('utf-8')))
            print("The target is in constellation {}".format(get_constellation(
                skobj)[0]))

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


def toco_simbad(args=None):
    """
    like toco but prints some more info
    """
    if args is None:
        parser = argparse.ArgumentParser(
            description="Information for a TESS target")
        parser.add_argument('tic', type=float,
                            help="TICID of target")
        args = parser.parse_args(args)
        args = vars(args)
    tic = args['tic']
    _output = print_results(tic, simbad_search=True)


def toco_name(args=None):
    """
    like toco but prints some more info
    """
    if args is None:
        parser = argparse.ArgumentParser(
            description="Information for a TESS target")
        parser.add_argument('name', nargs='+',
                            help="Name of the target")
        args = parser.parse_args(args)
        args.name = ' '.join(args.name)
        args = vars(args)
    name = args['name']
    tic = get_tic_name(name)
    _output = print_results(tic, simbad_search=True)



def toco_coords(args=None):
    """
    like tocot but starts with coords
    """
    if args is None:
        parser = argparse.ArgumentParser(
            description="Information for a TESS target")
        parser.add_argument('ra', type=float,
                            help="Right Ascension in decimal degrees (J2000).")
        parser.add_argument('dec', type=float,
                            help="Declination in decimal degrees (J2000).")
        args = parser.parse_args(args)
        args = vars(args)
    ra = args['ra']
    dec = args['dec']
    tic = get_tic_radec(ra, dec)
    _output = print_results(tic, simbad_search=True)


def get_tic_radec(ra, dec):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        catalogData = Catalogs.query_region('{} {}'.format(
            ra, dec),
            catalog='Tic', radius=0.006)

    try:
        return catalogData['ID'][0]
    except IndexError:
        logger.error("No TIC target at those coordiantes")
        sys.exit(1)

def get_tic_name(name):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        result_table = Simbad.query_object(name)
    if result_table is None:
        logger.error("Target name failed to resolve, please check")
        sys.exit(1)

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        ra_sex = result_table['RA'][0]
        dec_sex = result_table['DEC'][0]
        catalogData = Catalogs.query_region(SkyCoord
            (ra_sex, dec_sex, unit=(u.hour, u.deg)),
            catalog='Tic', radius=0.006)

    try:
        return catalogData['ID'][0]
    except IndexError:
        logger.error("No TIC target at those coordiantes")
        sys.exit(1)
