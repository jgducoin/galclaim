#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import errno
from astropy.table import Column
import copy
import healpy as hp
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import ascii
from astropy.table import Table
import numpy as np
import utils
import plotting
import argparse

"""
Created on Tue Apr 21 22:33:58 2020

@author: Ducoin Jean-Grégoire ducoin@iap.fr
         David Corre
"""
"""
This code provide tools to calculate the association false alarm rate between a
given transient error box and a given galaxy.
"""


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


class Association:

    """
    Dedicated class to calculate the p value of a given association 
    """

    def __init__(
        self,
        RA,
        Dec,
        error_radius,
        catalog,
        transient,
        transient_error_radius,
        verbose,
    ):
        """
        Initialisation of the class

        Input parameters
        ----------------
        RA : float
            Right ascension of the center of the transient error box
        Dec : float
            Declination of the center of the transient error box        
        error_radius : float
            Radius of the circle of the error box.
            To be provided in deg!!!
        catalog : string
            catalog you want to consider, can be "Pan-STARRS", "HSC" or "AllWISE"
        """

        self.transient_RA = RA
        self.transient_Dec = Dec
        self.transient_RA_err = error_radius
        self.transient = transient
        self.error_radius = transient_error_radius
        self.verbose = verbose
        if catalog not in ["Pan-STARRS", "AllWISE", "HSC"]:
            raise ValueError(
                "the catalog you want to consider is not in the list of \
possible catalog, it can be 'Pan-STARRS', 'HSC'  or 'AllWISE' "
            )

        self.catalog = catalog

        # first_query to get all compatible objects with the transient error box
        self.query_transient_box = self.query_box(
            RA=RA, Dec=Dec, radius=error_radius, catalog=catalog
        )

    def query_box(self, RA, Dec, radius, catalog):
        """Query the provided position with the given radius and catalog"""

        coords = SkyCoord(
                    RA,
                    Dec,
                    unit=(u.degree,u.degree),
                    frame='icrs')
        
        if catalog == "Pan-STARRS" or catalog == 'AllWISE':
            query = utils.query_catalog(
                coords,
                radius * u.degree,
                catalog,
                radius_min = 1 *u.arcsec
            )

        elif catalog == "HSC":
            query = utils.hsccone_full(
                RA,
                Dec,
                radius,
                release="v3",
                format="csv",
                magtype="magauto",
                columns=None,
                verbose=False,
            )
            
            if len(query) == 0:
                
                query = None

        # Check that the query returns something
        return query
    
    
    def get_sigma_localy_catalog(self):
        """
        Calcuate the sigma parameter according to Chrimes18 (https://arxiv.org/pdf/1804.08971.pdf) 

        In th case of HSC the catalog is not all sky,
        in order to avoid border image issue instead of onsidering a big circle
        aroud the transient we take all the sources in the same image (202*202 arcsec^2)

        Input parameters
        ----------------
        """

        # Security if no object are found compatible with the transient error box
        
        if self.query_transient_box is None:
            return

        if self.catalog == "HSC":
            Ncounts, dummy_radius, self.col_find_list = utils.get_N_in_image(
                self, self.transient_RA, self.transient_Dec
            )

            # for HSC this is the image size (see get_N_in_image function)
            self.total_area = (202 / 3600) ** 2 - (np.pi * 30. / 3600 ** 2)


        elif self.catalog == 'Pan-STARRS' or self.catalog == "AllWISE":
            big_radius = 0.05  # 3 arcmin

            # add a security, if the big radius is smaller than the transient error radius
            # (often true with BAT detection only), take the transient error radius in those cases
            if big_radius < self.transient_RA_err:
                big_radius = self.transient_RA_err

            query = self.query_box(
                RA=self.transient_RA,
                Dec=self.transient_Dec,
                radius=big_radius,
                catalog=self.catalog,
            )

            Ncounts = utils.count_galaxies_in_query(
                self.query_transient_box, query, self.transient_RA, self.transient_Dec, self.catalog
            )

            # the second terme is the center of the filed supressed (donuts)
            self.total_area = (np.pi * big_radius ** 2) - (np.pi * 30. / 3600 ** 2)


        self.sigma_list = copy.deepcopy(Ncounts)
        for filt, N_list in self.sigma_list.items():
            self.sigma_list[filt] = N_list / self.total_area


    def compute_pval(self):
        """
        This function compute the final p value of the association

        Input parameters
        ----------------

        """
        # Security if no object are found compatible with the transient error box
        if self.query_transient_box is None:
            if self.verbose:
                print("No compatible object found in the transient error box")
            self.pval_list = None
            return

        compute_distance_list(self)

        nu = copy.deepcopy(self.sigma_list)


        for filt, nu_list in nu.items():
            nu[filt] = (
                np.pi
                * (np.array(self.query_transient_box["dist_to_center"] ** 2))
                * self.sigma_list[filt]
            )

        self.pval_list = copy.deepcopy(nu)

        for filt, pval in self.pval_list.items():

            exp_term = np.exp(-nu[filt])
            self.pval_list[filt] = (1 - (exp_term)).tolist()
            self.pval_list[filt] = np.round(self.pval_list[filt],7)


    def parse_result(self, catalog):
        """
        This function parse all the computation result to an readable table

        Input parameters
        ----------------

        """

        # security if no object compatible with the transient error box
        if self.pval_list is None:
            self.table = None
            return

        if catalog == 'Pan-STARRS':

            self.table = self.query_transient_box[
                "objName",
                "objID",
                "raMean",
                "decMean",
                "gKronMag",
                "gKronMagErr",
                "rKronMag",
                "rKronMagErr",
                "iKronMag",
                "iKronMagErr",
                "zKronMag",
                "zKronMagErr",
                "yKronMag",
                "yKronMagErr",
                "objType",
                "dist_to_center",
            ]
            
            #rounding all the columns that we need to round
            self.table["gKronMag"] = np.round(self.table["gKronMag"],3)
            self.table["gKronMagErr"] = np.round(self.table["gKronMagErr"],3)
            self.table["rKronMag"] = np.round(self.table["gKronMag"],3)
            self.table["rKronMagErr"] = np.round(self.table["gKronMagErr"],3)
            self.table["iKronMag"] = np.round(self.table["gKronMag"],3)
            self.table["iKronMagErr"] = np.round(self.table["gKronMagErr"],3)
            self.table["zKronMag"] = np.round(self.table["gKronMag"],3)
            self.table["zKronMagErr"] = np.round(self.table["gKronMagErr"],3)
            self.table["yKronMag"] = np.round(self.table["gKronMag"],3)
            self.table["yKronMagErr"] = np.round(self.table["gKronMagErr"],3)
            
        elif catalog == 'AllWISE':
            self.table = self.query_transient_box[
                'designation',
                'ra',
                'dec',
                'w1mpro',
                'w1sigmpro',
                'w2mpro',
                'w2sigmpro',
                'w3mpro',
                'w3sigmpro',
                'w4mpro',
                'w4sigmpro',
                'objType',
                'dist_to_center'
            ]
            
            
        #type_list = self.table['objType']
        self.table.rename_column('objType', 'type')
        # go back to arcsec for the distance_to_center
        self.table["dist_to_center"] = self.table["dist_to_center"] * 3600
        self.table["dist_to_center"] = np.round(self.table["dist_to_center"],3)
        #self.table["type"] = type_list

        for filt, pval in self.pval_list.items():
            self.table[filt + "_pval"] = pval

        self.table.sort("dist_to_center")
        
        


    def parse_result_HSC(self):
        """
        This function parse all the computation result to an readable table

        Input parameters
        ----------------

        """
        # security if no object compatible with the transient error box
        if self.pval_list == None:
            self.table = None
            return

        type_list = utils.object_type_HSC(self.query_transient_box)

        for i in range(len(type_list)):
            if type_list[i] == "undefined":
                self.pval_list[i] = "None"

        self.table = self.query_transient_box[
            "MatchID", "MatchRA", "MatchDec", "Filter", "dist_to_center"
        ]

        # change the dtype of the Filter columns, unknown error while saving to ecsv otherwise (when there is only one filter)
        new_filt_col = []
        # print('self.table["Filter"] =', self.table['Filter'])
        for i in range(len(self.table["Filter"])):

            new_filt_col += [self.table["Filter"][i][0]]

        # print("new_filt_col =", new_filt_col)

        self.table["Filter"] = np.array(new_filt_col)

        # go back to arcsec for the distance_to_center
        self.table["dist_to_center"] = self.table["dist_to_center"] * 3600
        self.table["dist_to_center"] = np.round(self.table["dist_to_center"],3)

        mag_col = []

        for i in range(len(self.query_transient_box)):

            for col_name in self.query_transient_box.columns:
                if col_name.endswith(self.table["Filter"][i]):
                    col_find = col_name


            mag_col += [self.query_transient_box[col_find][i]]

        self.table["mag"] = np.round(mag_col,3)
        self.table["type"] = type_list

        for filt, pval in self.pval_list.items():
            self.table[filt + "_pval"] = pval

        self.table.sort("dist_to_center")


def compute_distance_list(self):
    """
    We have to compute the distance between each objects and the transient error box center
    """

    col_dist = []

    transient_center = SkyCoord(
        self.transient_RA,
        self.transient_Dec,
        unit=(u.degree,u.degree),
        frame='icrs')

    if self.catalog == "Pan-STARRS":

        for i in range(len(self.query_transient_box)):

            RA_object = self.query_transient_box["raMean"][i]
            Dec_object = self.query_transient_box["decMean"][i]

            dist_object = hp.rotator.angdist(
                [float(RA_object), float(Dec_object)],
                [float(self.transient_RA), float(self.transient_Dec)],
                lonlat=True,
            )
            # go back from rad to deg
            dist_object = (180 / np.pi) * dist_object[0]

            col_dist += [dist_object]

    if self.catalog == "HSC":

        for i in range(len(self.query_transient_box)):

            RA_object = self.query_transient_box["MatchRA"][i]
            Dec_object = self.query_transient_box["MatchDec"][i]

            dist_object = hp.rotator.angdist(
                [float(RA_object), float(Dec_object)],
                [float(self.transient_RA), float(self.transient_Dec)],
                lonlat=True,
            )
            # go back from rad to deg
            dist_object = (180 / np.pi) * dist_object[0]

            col_dist += [dist_object]

    if self.catalog == "AllWISE":
        coords = SkyCoord(
                    self.query_transient_box["ra"],
                    self.query_transient_box["dec"],
                    unit=(u.degree,u.degree),
                    frame='icrs')
        sep = coords.separation(transient_center)

        col_dist = sep.value 

    self.query_transient_box["dist_to_center"] = col_dist

    return


def do_association(
    transients,
    outputdir,
    catalogs,
    do_plot_dist_pval=True,
    verbose=False,
):
    """
    Final function call by ../crosmatch/HostHunt.py to compute pval for each transients


    """

    # Create output dir if not exists
    mkdir_p(outputdir)
    mkdir_p(outputdir+'/plots')
    
    print("computing pval for each transients ...")

    for i in range(len(transients)):

        transient = transients["ID"][i]

        RA = transients['RA'][i]
        DEC = transients['Dec'][i]
        transient_error_radius = float(transients['pos_err'][i])

        # we consider candidate host inside an 30 arcsec radius circle
        transient_candidate_radius = 30 / 3600

        for catalog in catalogs:

            # Initialyse the class
            Association_transient = Association(
                RA,
                DEC,
                transient_candidate_radius,
                catalog,
                transient,
                transient_error_radius,
                verbose,
            )


            sigma = Association_transient.get_sigma_localy_catalog()

            # compute the final p value
            Association_transient.compute_pval()

            if catalog == "Pan-STARRS":
                Association_transient.parse_result(catalog)
            elif catalog == "HSC":
                Association_transient.parse_result_HSC()
            elif catalog == "AllWISE":
                Association_transient.parse_result(catalog)


            if do_plot_dist_pval:
                plotDir = os.path.join(outputdir)
                if not os.path.exists(plotDir):
                    mkdir_p(plotDir)
                # print("plotting the R(pavl) results")
                plotting.plot_dist_pval(Association_transient, plotDir)

            # add the metadatas
            Association_transient.table = Table(
                Association_transient.table,
                meta={
                    "RA": str(RA),
                    "DEC": str(DEC),
                    "pos_err": str(transient_error_radius),
                },
            )

            if len(Association_transient.table) == 0:
                if verbose:
                    print(
                        "No compatible object found for transient {}. Empty file saved.".format(
                            transients["ID"][i]
                        )
                    )

                # add an empty column here, otherwise astropy table fail to be read afterward
                Association_transient.table["empty"] = []

                Association_transient.table.write(
                    outputdir + "{}_{}.ecsv".format(transients["ID"][i], catalog),
                    format="ascii.ecsv",
                    overwrite=True,
                )

            else:
                Association_transient.table.write(
                    outputdir + "{}_{}.ecsv".format(transients["ID"][i], catalog),
                    format="ascii.ecsv",
                    overwrite=True,
                )

    return


#######  MAIN ######
if __name__ == "__main__":

    parser = argparse.ArgumentParser(
            description='Galclaim launching code.')

    
    parser.add_argument('--transient_file',
                        dest='transient_file',
                        required=True,
                        metavar='transient_file',
                        type=str,
                        help='Path to the transient file containing the transient localisation.')

    parser.add_argument('--catalogs',
                        dest='catalogs',
                        required=False,
                        metavar='catalogs',
                        type=str,
                        nargs='+',
                        default=[],
                        help='Catalogs to be used for the association. Available are Pan-STARRS, HSC and AllWISE.')
    
    parser.add_argument('--do_plots',
                        dest='do_plots',
                        required=False,
                        action="store_true",
                        help='Enable plots.')

    parser.add_argument('--verbose',
                        dest='verbose',
                        required=False,
                        action="store_true",
                        help='Enable verbose.')
    
    args = parser.parse_args()

    transient_file = Table.read(args.transient_file, format="ascii.ecsv")


    if not args.catalogs :
        args.catalogs = ["Pan-STARRS", "HSC", "AllWISE"]


    do_association(
        transient_file,
        outputdir="./results/",
        catalogs=args.catalogs,
        do_plot_dist_pval=args.do_plots,
        verbose=args.verbose,
    )
