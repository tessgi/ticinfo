from __future__ import absolute_import

import argparse
import sys
import warnings

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord, get_constellation
from astroquery.exceptions import NoResultsWarning, ResolverError
from astroquery.mast import Catalogs, Observations
from astroquery.simbad import Simbad
from astropy.time import Time
from numpy import array, ones

from tesswcs import pointings, WCS



_next_cycle = pointings['Cycle'].max()
future_mask = pointings['Cycle'] == _next_cycle
past_mask = np.asarray(pointings['End']) < Time.now().jd
_ra = np.asarray(pointings['RA'])
_dec = np.asarray(pointings['Dec'])
_roll = np.asarray(pointings['Roll'])
_sectors = np.asarray(pointings['Sector'])

from . import get_logger

# warnings.filterwarnings('ignore', category=UserWarning, append=True)

logger = get_logger()

MAXTIC = 10005000540


class Target(object):
    def __init__(self, tic):
        self.validate_tic(tic)
        self.tic = int(tic)

    def validate_tic(self, tic):
        try:
            tic = int(tic)
        except:
            logger.error("[red bold]Not a valid TIC number[/]", extra={"markup":True})
            sys.exit(1)

        if tic > MAXTIC:
            logger.error("[red bold]Not a valid TIC number[/]", extra={"markup":True})
            sys.exit(1)

        if tic < 0:
            logger.error("[red bold]Not a valid TIC number[/]", extra={"markup":True})
            sys.exit(1)

    def query(self):
        try:
            catalogData = Catalogs.query_object(
                f"TIC {self.tic}", catalog="TIC", radius=3 * u.arcsec
            )
            return catalogData
        except ResolverError:
            logger.error(
                "[bold red]Target name failed to resolve, please check the input.[/bold red]",
                extra={"markup": True},
            )
            sys.exit(1)

    def get_obs(self):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=NoResultsWarning)

            products = Observations.query_criteria(
                objectname=f"TIC {self.tic}",
                t_exptime=(19, 1801),
                project="TESS",
                obs_collection="TESS",
                dataproduct_type="timeseries",  # discard images
                radius=0.0001 * u.arcsec,
            )

            if len(products) == 0:
                return [[], [], []]
            # get the 2min data
            mask_2 = ones(len(products), dtype=bool)
            mask_2 &= array([(texp > 119) for texp in products["t_exptime"]])
            mask_2 &= array(
                [prov.lower() == "spoc" for prov in products["provenance_name"]]
            )
            mask_2 &= array(
                [
                    filename.endswith("_lc.fits")
                    if isinstance(filename, str)
                    else False
                    for filename in products["dataURL"]
                ]
            )

            obs_sectors_2 = sorted(products[mask_2]["sequence_number"])

            # get the FFI data
            mask_ffi = ones(len(products), dtype=bool)
            mask_ffi &= array([(texp > 121) for texp in products["t_exptime"]])
            mask_ffi &= array(
                [
                    (dataproduct == "image")
                    for dataproduct in products["dataproduct_type"]
                ]
            )
            mask_ffi &= array(
                [prov.lower() == "spoc" for prov in products["provenance_name"]]
            )

            obs_sectors_ffi = sorted(products[mask_ffi]["sequence_number"])

            # get the 20 s data
            mask_20 = ones(len(products), dtype=bool)
            mask_20 &= array([(texp < 119) for texp in products["t_exptime"]])
            mask_20 &= array(
                [prov.lower() == "spoc" for prov in products["provenance_name"]]
            )
            mask_20 &= array(
                [
                    filename.endswith("lc.fits") if isinstance(filename, str) else False
                    for filename in products["dataURL"]
                ]
            )

            obs_sectors_20 = sorted(products[mask_20]["sequence_number"])
        obs_sectors = [obs_sectors_2, obs_sectors_ffi, obs_sectors_20]
        return obs_sectors


def show_results(tic=12350, simbad_search=False, data_search=False, future_search=False, long_table=False):
    target = Target(tic)
    catalogData = target.query()[0]
    skobj = SkyCoord(
        ra=catalogData["ra"] * u.degree,
        dec=catalogData["dec"] * u.degree,
        frame="icrs",
    )
    if simbad_search:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            customSimbad = Simbad()
            customSimbad.add_votable_fields(
                "ra(2;A;ICRS;J2000;2000)", "dec(2;D;ICRS;J2000;2000)"
            )
            customSimbad.remove_votable_fields("coordinates")

            # try different search radii, be fast if possible
            for i in [5, 10, 20]:
                result_table = customSimbad.query_region(skobj, radius=i * u.arcsec)
                if result_table is None:
                    continue
                else:
                    break

        if result_table is None:
            logger.warning("No Simbad target resolved")
        else:
            logger.info(
                f"Target name: [green bold]{result_table['MAIN_ID'][0]}[/green bold]",
                extra={"markup": True},
            )
        logger.info(
            f"[italic]The target is in constellation {get_constellation(skobj)}[/italic]\n",
            extra={"markup": True},
        )

    catalogData["ra"] = catalogData["ra"].round(5)
    catalogData["dec"] = catalogData["dec"].round(5)
    catalogData["eclong"] = catalogData["eclong"].round(5)
    catalogData["eclat"] = catalogData["eclat"].round(5)
    catalogData["pmRA"] = catalogData["pmRA"].round(2)
    catalogData["pmDEC"] = catalogData["pmDEC"].round(2)
    catalogData["Tmag"] = catalogData["Tmag"].round(2)
    catalogData["Vmag"] = catalogData["Vmag"].round(2)
    catalogData["Kmag"] = catalogData["Kmag"].round(2)
    catalogData["plx"] = catalogData["plx"].round(2)

    #RAH - modified code to print output expected from the tutorial
    logger.info(
        f"""{catalogData[['ID', 'ra', 'dec', 'plx', 'Tmag', 'Vmag', 'Kmag', 'Teff',
                       'rad', 'mass']]}\n"""
    )

    if data_search:
        logger.start_spinner("Searching for TESS Data Product Availability...")

        obs_sectors = target.get_obs()
        obs2, obsffi, obs20 = obs_sectors
        past_sectors = get_observability(skobj, _sectors[past_mask], _ra[past_mask], _dec[past_mask], _roll[past_mask])
        logger.stop_spinner()
        logger.info(f"[bold]Stored[/bold]: Object observable in archived TESS FFIs during sectors:   {past_sectors}", extra={'markup':True})
        logger.info(f"[bold]Stored[/bold]: FFI light curve data at MAST for sectors:   {sorted(list(set(obsffi)))}", extra={'markup':True})
        logger.info(f"[bold]Stored[/bold]: 2-min light curve data at MAST for sectors: {sorted(list(set(obs2)))}", extra={'markup':True})
        logger.info(f"[bold]Stored[/bold]: 20-s light curve data at MAST for sectors:  {sorted(list(set(obs20)))}", extra={'markup':True})
    
    if future_search:
        logger.start_spinner(f"Determining observability in Cycle {_next_cycle}...")
        future_sectors = get_observability(skobj, _sectors[future_mask], _ra[future_mask], _dec[future_mask], _roll[future_mask])
        logger.stop_spinner()
        logger.info(f"[bold]Future[/bold]: Will be observable in Cycle {_next_cycle} sectors:  {future_sectors}", extra={'markup':True})

    if long_table:
        logger.info(f"\n[bold]Additional Object Information:[/bold]", extra={'markup':True})
        logger.info(f"""{catalogData[['pmRA', 'pmDEC', 'eclong', 'eclat']]}\n""", extra={'markup':True})

def toco(args=None):
    """
    exposes toco to the command line
    """
    if args is None:
        parser = argparse.ArgumentParser(description="Information for a TESS target")
        parser.add_argument(
            "input",
            nargs="+",
            help="Input target, can be TIC ID, name, or RA/Dec coords in decimal.",
        )
        parser.add_argument(
            "-s",
            "--stored",
            action="store_true",
            help="Will query MAST to find what data is available. Will list sectors with TESS data for target in FFIs, 2 minute TPFs, and 20s TPFs.",
        )
        parser.add_argument(
            "-f",
            "--future",
            action="store_true",
            help=f"Will determine which sectors the target is obserable in during the next observing cycle (Cycle {_next_cycle})",
        )
        parser.add_argument(
            "-l",
            "--long",
            action="store_true",
            help=f"Will output an additional table that contains ecliptic lat/long and pm RA/Dec",
        )
        args = parser.parse_args(args)
        args = vars(args)
    tic = args["input"]
    logger.start_spinner("Searching for target...")
    if len(tic) == 1:
        tic = tic[0]
        if isinstance(tic, str):
            if tic.startswith("TIC") | tic.startswith("tic"):
                tic = tic[3:].strip()

    elif len(tic) == 2:
        if tic[0] in ["TIC", "tic"]:
            tic = tic[1]
        #RAH: The below strips numbers of non number inputs like +-., for RA and Dec
        #This now allows the program to process input RA and DECs    
        elif np.all([s.lstrip('+-.,').isnumeric for s in tic]):
            tic = get_tic_radec(*tic)
        else:
            tic = get_tic_name(" ".join(tic))
    else:
        tic = get_tic_name("".join(tic))

    if not tic.isnumeric():
        tic = get_tic_name(tic)
    logger.stop_spinner()
    show_results(tic, simbad_search=True, data_search=args["stored"], future_search=args["future"], long_table=args["long"])


def toco_name(args=None):
    logger.info(
        "[bold red]`tocon` is now deprecated, you can have the same functionality from calling `toco`[/]",
        extra={"markup": True},
    )
    toco(args)


def toco_simbad(args=None):
    logger.info(
        "[bold red]`tocot` is now deprecated, you can have the same functionality from calling `toco`[/]",
        extra={"markup": True},
    )
    toco(args)


def toco_coords(args=None):
    logger.info(
        "[bold red]`tococ` is now deprecated, you can have the same functionality from calling `toco`[/]",
        extra={"markup": True},
    )
    toco(args)


def get_tic_radec(ra, dec):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        catalogData = Catalogs.query_region(
            "{} {}".format(ra, dec), catalog="Tic", radius=1 * u.arcsec
        )

    try:
        return catalogData["ID"][0]
    except IndexError:
        logger.error("No TIC target at those coordiantes")
        logger.stop_spinner()
        sys.exit(1)


def get_tic_name(name):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        customSimbad = Simbad()
        customSimbad.add_votable_fields(
            "ra(2;A;ICRS;J2000;2000)", "dec(2;D;ICRS;J2000;2000)"
        )
        customSimbad.remove_votable_fields("coordinates")
        result_table = customSimbad.query_object(name)
    if result_table is None:
        logger.error(
            "[bold red]Target name failed to resolve, please check the input.[/bold red]",
            extra={"markup": True},
        )
        logger.stop_spinner()
        sys.exit(1)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ra_sex = result_table["RA_2_A_ICRS_J2000_2000"][0]
        dec_sex = result_table["DEC_2_D_ICRS_J2000_2000"][0]
        catalogData = Catalogs.query_region(
            SkyCoord(ra_sex, dec_sex, unit=(u.hour, u.deg)), catalog="Tic", radius=0.006
        )

    try:
        return catalogData["ID"][0]
    except IndexError:
        logger.error("No TIC target at those coordiantes")
        sys.exit(1)

def get_observability(coord, sectors, ra, dec, roll):
    """Find out if target will be observable in this Cycle"""
    observable_sectors = []
    for sdx in range(len(ra)):
        observable = False
        for camera in np.arange(1, 5):
            for ccd in np.arange(1, 5):
                wcs = WCS.predict(ra[sdx], dec[sdx], roll[sdx], camera=camera, ccd=ccd)
                observable |= wcs.footprint_contains(coord)
        if observable:
            observable_sectors.append([sectors[sdx]])
    if observable_sectors == []:
        return []
    else:
        return np.hstack(observable_sectors)

