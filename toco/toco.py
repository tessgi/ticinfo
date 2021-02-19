from __future__ import absolute_import
from . import logger
import argparse

from numpy import ones, array
from astropy.coordinates import SkyCoord, get_constellation
from astroquery.mast import Catalogs
from astroquery.simbad import Simbad
from astroquery.mast import Observations
from astroquery.exceptions import ResolverError, NoResultsWarning
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
        catalogData = Catalogs.query_object(
            f'TIC {self.tic}', catalog="TIC", radius=3 * u.arcsec)
        return catalogData

    def get_obs(self):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=NoResultsWarning)

            products = Observations.query_criteria(
                objectname=f'TIC {self.tic}',
                t_exptime=(19, 1801),
                project='TESS', obs_collection='TESS',
                radius=.0001 * u.arcsec)

            if len(products) == 0:
                return [[], [], []]
            # get the 2min data
            mask_2 = ones(len(products), dtype=bool)
            mask_2 &= array([(texp > 119) for texp in products['t_exptime']])
            mask_2 &= array(
                [prov.lower() == 'spoc' for prov in products['provenance_name']])
            mask_2 &= array([filename.endswith('_lc.fits') if isinstance(
                filename, str) else False for filename in products['dataURL']])

            obs_sectors_2 = sorted(products[mask_2]['sequence_number'])

            # get the FFI data
            mask_ffi = ones(len(products), dtype=bool)
            mask_ffi &= array([(texp > 121) for texp in products['t_exptime']])
            mask_ffi &= array(
                [(dataproduct == 'image') for
                 dataproduct in products['dataproduct_type']])
            mask_ffi &= array(
                [prov.lower() == 'spoc' for prov in products['provenance_name']])

            obs_sectors_ffi = sorted(products[mask_ffi]['sequence_number'])

            # get the 20 s data
            mask_20 = ones(len(products), dtype=bool)
            mask_20 &= array([(texp < 119) for texp in products['t_exptime']])
            mask_20 &= array(
                [prov.lower() == 'spoc' for prov in products['provenance_name']])
            mask_20 &= array([filename.endswith('_lc.fits') if isinstance(
                filename, str) else False for filename in products['dataURL']])

            obs_sectors_20 = sorted(products[mask_20]['sequence_number'])
        obs_sectors = [obs_sectors_2, obs_sectors_ffi, obs_sectors_20]
        return obs_sectors


def print_results(tic=12350, simbad_search=False, data_search=False):
    target = Target(tic)
    catalogData = target.query()[0]

    catalogData['ra'] = catalogData['ra'].round(5)
    catalogData['dec'] = catalogData['dec'].round(5)
    catalogData['eclong'] = catalogData['eclong'].round(5)
    catalogData['eclat'] = catalogData['eclat'].round(5)
    catalogData['pmRA'] = catalogData['pmRA'].round(2)
    catalogData['pmDEC'] = catalogData['pmDEC'].round(2)
    catalogData['Tmag'] = catalogData['Tmag'].round(2)
    catalogData['Vmag'] = catalogData['Vmag'].round(2)
    catalogData['Kmag'] = catalogData['Kmag'].round(2)

    print(catalogData[['ID', 'ra', 'dec', 'pmRA', 'pmDEC',
                       'eclong', 'eclat', 'Tmag', 'Vmag', 'Kmag', 'Teff',
                       'rad', 'mass', 'd', ]])

    if simbad_search:
        skobj = SkyCoord(ra=catalogData['ra'] * u.degree,
                         dec=catalogData['dec'] * u.degree,
                         frame='icrs')
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')

            customSimbad = Simbad()
            customSimbad.add_votable_fields(
                'ra(2;A;ICRS;J2000;2000)', 'dec(2;D;ICRS;J2000;2000)')
            customSimbad.remove_votable_fields('coordinates')

            # try different search radii, be fast if possible
            for i in [5, 10, 20]:
                result_table = customSimbad.query_region(
                    skobj, radius=i * u.arcsec)
                if result_table is None:
                    continue
                else:
                    break

        if result_table is None:
            logger.warning("No Simbad target resolved")
        else:
            print()
            print('Target name: {}'.format(
                result_table['MAIN_ID'][0]))
        print("The target is in constellation {}".format(get_constellation(
            skobj)))

    if data_search:
        obs_sectors = target.get_obs()
        obs2, obsffi, obs20 = obs_sectors

        print(f'FFI data at MAST for sectors:   {sorted(list(set(obsffi)))}')
        print(f'2-min data at MAST for sectors: {sorted(list(set(obs2)))}')
        print(f'20-s data at MAST for sectors:  {sorted(list(set(obs20)))}')

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
        parser.add_argument('-s', '--showdata', action='store_true',
                            help='List sectors with TESS data')
        args = parser.parse_args(args)
        args = vars(args)
    tic = args['tic']
    data_search = args['showdata']
    _output = print_results(tic, data_search=data_search)


def toco_simbad(args=None):
    """
    like toco but prints some more info
    """
    if args is None:
        parser = argparse.ArgumentParser(
            description="Information for a TESS target")
        parser.add_argument('tic', type=float,
                            help="TICID of target")
        parser.add_argument('-s', '--showdata', action='store_true',
                            help='List sectors with TESS data')
        args = parser.parse_args(args)
        args = vars(args)
    tic = args['tic']
    data_search = args['showdata']
    _output = print_results(tic, simbad_search=True, data_search=data_search)


def toco_name(args=None):
    """
    like toco but prints some more info
    """
    if args is None:
        parser = argparse.ArgumentParser(
            description="Information for a TESS target")
        parser.add_argument('name', nargs='+',
                            help="Name of the target")
        parser.add_argument('-s', '--showdata', action='store_true',
                            help='List sectors with TESS data')
        args = parser.parse_args(args)
        args.name = ' '.join(args.name)
        args = vars(args)
    name = args['name']
    tic = get_tic_name(name)
    data_search = args['showdata']
    _output = print_results(tic, simbad_search=True, data_search=data_search)


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
        parser.add_argument('-s', '--showdata', action='store_true',
                            help='List sectors with TESS data')
        args = parser.parse_args(args)
        args = vars(args)
    ra = args['ra']
    dec = args['dec']
    data_search = args['showdata']
    tic = get_tic_radec(ra, dec)
    _output = print_results(tic, simbad_search=True, data_search=data_search)


def get_tic_radec(ra, dec):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        catalogData = Catalogs.query_region('{} {}'.format(
            ra, dec),
            catalog='Tic', radius=1 * u.arcsec)

    try:
        return catalogData['ID'][0]
    except IndexError:
        logger.error("No TIC target at those coordiantes")
        sys.exit(1)


def get_tic_name(name):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        customSimbad = Simbad()
        customSimbad.add_votable_fields(
            'ra(2;A;ICRS;J2000;2000)', 'dec(2;D;ICRS;J2000;2000)')
        customSimbad.remove_votable_fields('coordinates')
        result_table = customSimbad.query_object(name)
    if result_table is None:
        logger.error("Target name failed to resolve, please check")
        sys.exit(1)

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        ra_sex = result_table['RA_2_A_ICRS_J2000_2000'][0]
        dec_sex = result_table['DEC_2_D_ICRS_J2000_2000'][0]
        catalogData = Catalogs.query_region(SkyCoord
                                            (ra_sex, dec_sex, unit=(u.hour, u.deg)),
                                            catalog='Tic', radius=0.006)

    try:
        return catalogData['ID'][0]
    except IndexError:
        logger.error("No TIC target at those coordiantes")
        sys.exit(1)
