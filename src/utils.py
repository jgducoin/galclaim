#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
from astropy.table import Table, vstack
from astropy.coordinates import SkyCoord
import astropy.units as u
from astroquery.mast import Catalogs
import numpy as np
from astroquery.irsa import Irsa
from astroquery.vizier import Vizier
import copy
import requests
from astropy.io import ascii
import os

"""
Created on Tue Apr 21 22:33:58 2020

@author: Ducoin Jean-Grégoire ducoin@iap.fr
         David Corre
"""


if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

hscapiurl = "https://catalogs.mast.stsci.edu/api/v0.1/hsc"

def add_object_type(table, catalog):
    """Adding information about object type"""
    type_list = []

    if catalog == 'Pan-STARRS':
        # Linear cut is advised for its simplicity in:
        # https://ui.adsabs.harvard.edu/abs/2016arXiv161205560C/abstract
        # But nn linear cuts are required for faint objects:
        # https://ui.adsabs.harvard.edu/abs/2014MNRAS.437..748F/abstract
        # So far only implemented linear cut
        maglim = 21
        for row in table:
            if not row['iKronMag'] or not row['iPSFMag']:
                type_list.append('unknown')
            elif row['iPSFMag'] - row['iKronMag'] > 0.05:
                type_list.append('galaxy')
            elif row['iPSFMag'] - row['iKronMag'] < 0.05:
                type_list.append('star')
            else:
                type_list.append('unknown')


    elif catalog == 'AllWISE':
        # https://ui.adsabs.harvard.edu/abs/2015MNRAS.448.1305K/abstract
        # W1mpro - J2MASS < 1.7 as a rough estimate
        # 78% of galaxies are selected, 1.8% cstellar contaminants
        # So does not mean that undefined are only stars
        # but when flagged as galaxy we are quite sure
        # BUT criteria tested on a bright sample, and we have usually faint objets
        # So this should be considered as a hint not a strict information

        # More complex criteria: https://arxiv.org/pdf/1607.01190.pdf
        for row in table:
            if row['w1sigmpro'] and row['j_msig_2mass']:
                if row['w1mpro'] - row['j_m_2mass'] < -1.7:
                    type_list.append('galaxy')
                else:
                    type_list.append('unknown')
            else:
                type_list.append('unknown')
    
    elif catalog == 'GLADE':
        #the GLADE catalog is supposed to only contain galaxy
        type_list.append('galaxy')
    
    else:
        for row in table:
            type_list.append('unknown')

    table['objType'] = type_list

    return table

def query_catalog(coords, radius, catalog, radius_min=20*u.arcsecond):
    """Query the provided position with the given radius and catalog"""

    if radius < radius_min:
        radius = radius_min

    if catalog == "Pan-STARRS":
        # Avoid duplicates using the flag below
        # Information taken from http://ps1images.stsci.edu/ps1_dr2_api.html
        constraints = {'primaryDetection':1, 'sort_by':'distance.ASC'}

        query = Catalogs.query_region(coords,
                                      radius=radius,
                                      catalog="Panstarrs",
                                      data_release="dr2",
                                      table="stack",
                                      primaryDetection = 1)
        if len(query) == 0:
            query = None
        else:
            query = add_object_type(query, catalog)
            query['distance_arcs'] = query['distance'] * 3600
            query.sort('distance_arcs')    

    elif catalog == 'AllWISE':
        # Need to sort these columns
        selcols = "designation,ra,dec,w1mag,w1mpro,w1gmag,ext_flg,var_flg,\
            ph_qual,xscprox,r_2mass,na,nb,w1snr,w1rchi2,rchi2,rchi2,w1sat,\
            n_2mass,tmass_key,j_m_2mass,h_m_2mass,k_m_2mass,sigra,sigdec,\
            sigradec,w1sigmpro,satnum,w1rchi2_pm,cc_flags,rel,w1nm,w1flux,\
            w1flg,w1mcor,w1magp,rho12,w1rsemi,w1ba,w1pa,w1gerr,w1gflg,\
            j_msig_2mass,h_msig_2mass,k_msig_2mass,w2mpro,w2sigmpro,w2snr,\
            w2rchi2,w3mpro,w3sigmpro,w3snr,w3rchi2,w4mpro,w4sigmpro,w4snr,\
            w4rchi2,w2gmag,w2gerr,w2gflg,w3gmag,w3gerr,w3gflg,w4gmag,w4gerr,\
            w4gflg"
        query = Irsa.query_region(coords,
                                  catalog='allwise_p3as_psd',
                                  selcols=selcols,
                                  radius=radius)
        if len(query) == 0:
            query = None
        else:
            query = add_object_type(query, catalog)
            
    elif catalog == 'GLADE':

        query = Vizier(columns=['all']).query_region(coords, 
                                    radius=radius,
                                    catalog='VII/291/gladep')
        
        if len(query) == 0:
            query = None
        else:
            query = query[0]
            query = add_object_type(query, catalog)

    return query


def count_galaxies_in_query(query_transient_box, query_random, RA, Dec, catalog):
    """
    This function count the number of galaxies in the outputed (and parsed) query from Pan-STARRS
    AllWISE and GLADE.
    
    Input parameters
    ----------------
    query_transient_box: astropy table
         Results of the query centered on the transient position
    query_random: astropy table
         Results of the query centered on the randomly drawned positions, 
         or results of the query with big_radius
    RA: float
        RA of the best transient localisation
    Dec: float
        Dec of the best transient localisation
    """
    
    # Number of objects in the transient error box query
    N_obj_transient_box = len(query_transient_box)

    nan_array = np.zeros(N_obj_transient_box)
    nan_array[:] = np.nan

    if catalog == 'Pan-STARRS':
        Ncounts_tot = {'gKronMag':copy.deepcopy(nan_array),
                   'rKronMag':copy.deepcopy(nan_array),
                   'iKronMag':copy.deepcopy(nan_array),
                   'zKronMag':copy.deepcopy(nan_array),
                   'yKronMag':copy.deepcopy(nan_array)}
        magerr_cols = [
                'gKronMagErr',
                'rKronMagErr',
                'iKronMagErr',
                'zKronMagErr',
                'yKronMagErr']
        RA_col = 'raMean'
        Dec_col = 'decMean'
        
    elif catalog == 'AllWISE':
        Ncounts_tot = {'w1mpro':copy.deepcopy(nan_array),
                       'w2mpro':copy.deepcopy(nan_array),
                       'w3mpro':copy.deepcopy(nan_array),
                       'w4mpro':copy.deepcopy(nan_array)}
        magerr_cols = [
                'w1sigmpro',
                'w2sigmpro',
                'w3sigmpro',
                'w4sigmpro']
        RA_col = 'ra'
        Dec_col = 'dec'
        
    elif catalog == 'GLADE':
        Ncounts_tot = {'Bmag':copy.deepcopy(nan_array),
                       'Jmag':copy.deepcopy(nan_array),
                       'Hmag':copy.deepcopy(nan_array),
                       'Kmag':copy.deepcopy(nan_array)}
        magerr_cols = [
                'e_Bmag',
                'e_Jmag',
                'e_Hmag',
                'e_Kmag']
        RA_col = 'RAJ2000'
        Dec_col = 'DEJ2000'
        
    if query_random is None:
        return Ncounts_tot

    else:
        for j in range(N_obj_transient_box):
            
            counter = 0
            
            for filt, N_list in Ncounts_tot.items():
                query = copy.deepcopy(query_random)

                mask_nan = np.isfinite(query[filt])

                query = query[mask_nan]
                #The following piece of code supress from the big_radius query the object which are
                #whithin the 30 arcsec from the best transient localisation. (creating a "donuts")
                #This ensure not to bias the backgroud population estimation with the true host

                pos_center = SkyCoord(RA*u.deg, Dec*u.deg, frame='icrs')
                pos_obj = SkyCoord(query[RA_col], query[Dec_col],
                                   unit=(u.degree, u.degree),
                                   frame='icrs')

                sep = pos_obj.separation(pos_center)

                mask_donuts = sep.to(u.arcsec).value >= 30
                
                mask_mag = query[filt] <= float(query_transient_box[filt][j])

                mask_star = query['objType'] != 'star'
                
                mask_bad_mag = ((query[magerr_cols[counter]] < 0.35) & \
                        (query[magerr_cols[counter]] > 0) & \
                        (query[magerr_cols[counter]].mask == False) & \
                        (query[filt].mask == False))
                        
                mask = mask_donuts & mask_mag & mask_star & mask_bad_mag

                Ncounts = len(query[mask])

                #if Ncounts == 0 that mean that this is a star and no galaxies have smaller mag
                #in this case the pval is not really meaningfull (this is a star!)
                #but as the classification is not 100% sure we want to compute the pval for this object also
                #so it as at least one object with <= mag, itself.
                #so put 1 as least
                Ncounts = max(1, Ncounts)

                #if there is no detection in a given band for the studied object
                if np.isnan(float(query_transient_box[filt][j])):
                    Ncounts = np.nan
                #we exept the GLADE catalog from this below error being finite criterium
                #indeed the GLADE catalog does not provide the error in B mag for most of the galaxies
                if np.isnan(float(query_transient_box[magerr_cols[counter]][j])) and catalog != 'GLADE':
                    Ncounts = np.nan

                Ncounts_tot[filt][j] = Ncounts
                counter += 1
                
        return Ncounts_tot


####HSC related functions

def hscmetadata(table="summary",release="v3",magtype="magaper2",baseurl=hscapiurl):
    """Return metadata for the specified catalog and table
    
    Parameters
    ----------
    table (string): summary, detailed, propermotions, or sourcepositions
    release (string): v3 or v2
    magtype (string): magaper2 or magauto (only applies to summary table)
    baseurl: base URL for the request
    
    Returns an astropy table with columns name, type, description
    """
    url = "{}/metadata".format(cat2url(table,release,magtype,baseurl=baseurl))
    r = requests.get(url)
    r.raise_for_status()
    v = r.json()
    # convert to astropy table
    tab = Table(rows=[(x['name'],x['type'],x['description']) for x in v],
               names=('name','type','description'))
    return tab


def cat2url(table="summary",release="v3",magtype="magaper2",baseurl=hscapiurl):
    """Return URL for the specified catalog and table
    
    Parameters
    ----------
    table (string): summary, detailed, propermotions, or sourcepositions
    release (string): v3 or v2
    magtype (string): magaper2 or magauto (only applies to summary table)
    baseurl: base URL for the request
    
    Returns a string with the base URL for this request
    """
    checklegal(table,release,magtype)
    if table == "summary":
        url = "{baseurl}/{release}/{table}/{magtype}".format(**locals())
    else:
        url = "{baseurl}/{release}/{table}".format(**locals())
    return url


def checklegal(table,release,magtype):
    """Checks if this combination of table, release and magtype is acceptable
    
    Raises a ValueError exception if there is problem
    """
    
    releaselist = ("v2", "v3")
    if release not in releaselist:
        raise ValueError("Bad value for release (must be one of {})".format(
            ', '.join(releaselist)))
    if release=="v2":
        tablelist = ("summary", "detailed")
    else:
        tablelist = ("summary", "detailed", "propermotions", "sourcepositions")
    if table not in tablelist:
        raise ValueError("Bad value for table (for {} must be one of {})".format(
            release, ", ".join(tablelist)))
    if table == "summary":
        magtypelist = ("magaper2", "magauto")
        if magtype not in magtypelist:
            raise ValueError("Bad value for magtype (must be one of {})".format(
                ", ".join(magtypelist)))


def mastQuery(request, url='https://mast.stsci.edu/api/v0/invoke'):
    """Perform a MAST query.

    Parameters
    ----------
    request (dictionary): The MAST request json object
    url (string): The service URL

    Returns the returned data content
    """
    
    # Encoding the request as a json string
    requestString = json.dumps(request)
    r = requests.post(url, data={'request': requestString})
    r.raise_for_status()
    return r.text


def resolve(name):
    """Get the RA and Dec for an object using the MAST name resolver
    
    Parameters
    ----------
    name (str): Name of object

    Returns RA, Dec tuple with position
    """

    resolverRequest = {'service':'Mast.Name.Lookup',
                       'params':{'input':name,
                                 'format':'json'
                                },
                      }
    resolvedObjectString = mastQuery(resolverRequest)
    resolvedObject = json.loads(resolvedObjectString)
    # The resolver returns a variety of information about the resolved object, 
    # however for our purposes all we need are the RA and Dec
    try:
        objRa = resolvedObject['resolvedCoordinate'][0]['ra']
        objDec = resolvedObject['resolvedCoordinate'][0]['decl']
    except IndexError as e:
        raise ValueError("Unknown object '{}'".format(name))
    return (objRa, objDec)


def create_table(query,verbose):
    """Transform the HSC query into an astropy table"""

    if query:
            
        # save typing a quoted list of columns
        columns = query[0:query.index('\r')].split(",")
        columns = [x.strip() for x in columns]
        columns = [x for x in columns if x and not x.startswith('#')]

        tab = ascii.read(query)
        
        if verbose:
            print("retrieved data and converted to {}-row astropy table".format(len(tab)))

    else:
        tab = None
    # clean up the output format    
    #for col in columns:
    #    if col[0:2] == 'A_':
    #        tab[col].format = "{:.3f}"
    
    #tab['MatchRA'].format = "{:.6f}"
    #tab['MatchDec'].format = "{:.6f}"
    #tab['StartMJD'].format = "{:.5f}"
    #tab['StopMJD'].format = "{:.5f}"
    
    return tab


def hscsearch(table="summary",release="v3",magtype="magaper2",format="csv",
              columns=None, baseurl=hscapiurl, verbose=False,
           **kw):
    """Do a general search of the HSC catalog (possibly without ra/dec/radius)
    
    Parameters
    ----------
    table (string): summary, detailed, propermotions, or sourcepositions
    release (string): v3 or v2
    magtype (string): magaper2 or magauto (only applies to summary table)
    format: csv, votable, json
    columns: list of column names to include (None means use defaults)
    baseurl: base URL for the request
    verbose: print info about request
    **kw: other parameters (e.g., 'numimages.gte':2).  Note this is required!
    """
    
    data = kw.copy()
    if not data:
        raise ValueError("You must specify some parameters for search")
    if format not in ("csv","votable","json"):
        raise ValueError("Bad value for format")
    url = "{}.{}".format(cat2url(table,release,magtype,baseurl=baseurl),format)
    if columns:
        # check that column values are legal
        # create a dictionary to speed this up
        dcols = {}
        for col in hscmetadata(table,release,magtype)['name']:
            dcols[col.lower()] = 1
        badcols = []
        for col in columns:
            if col.lower().strip() not in dcols:
                badcols.append(col)
        if badcols:
            raise ValueError('Some columns not found in table: {}'.format(', '.join(badcols)))
        # two different ways to specify a list of column values in the API
        # data['columns'] = columns
        data['columns'] = '[{}]'.format(','.join(columns))

    # either get or post works
    # r = requests.post(url, data=data)
    r = requests.get(url, params=data)

    if verbose:
        print(r.url)
    r.raise_for_status()
    if format == "json":
        return r.json()
    else:
        return r.text


def hsccone(ra,dec,radius,table="summary",release="v3",format="csv",magtype="magaper2",
            columns=None, baseurl=hscapiurl, verbose=False,
            **kw):
    """Do a cone search of the HSC catalog
    
    Parameters
    ----------
    ra (float): (degrees) J2000 Right Ascension
    dec (float): (degrees) J2000 Declination
    radius (float): (degrees) Search radius (<= 0.5 degrees)
    table (string): summary, detailed, propermotions, or sourcepositions
    release (string): v3 or v2
    magtype (string): magaper2 or magauto (only applies to summary table)
    format: csv, votable, json
    columns: list of column names to include (None means use defaults)
    baseurl: base URL for the request
    verbose: print info about request
    **kw: other parameters (e.g., 'numimages.gte':2)
    """
    
    data = kw.copy()
    data['ra'] = ra
    data['dec'] = dec
    data['radius'] = radius
    return hscsearch(table=table,release=release,format=format,magtype=magtype,
                     columns=columns,baseurl=baseurl,verbose=verbose,**data)


def hsccone_full(RA,DEC,
                 radius,
                 release="v3",
                 format="csv",
                 magtype="MagAuto",
                 columns=None,
                 verbose=False):
    """
    In the HSC summary catalog there is no extended flag, we have to fetch it from the detailed one.
    So in this function we do two query and add a columns in the summary output with the extended
    flag info for each objects
    
    Also we save the image ID for each object for later uses
    """
    #query the summary catalog
    query_summary = hsccone(ra=RA,dec=DEC,radius=radius,table="summary",release=release,format=format,magtype=magtype,
            columns=columns, verbose=verbose)
    
    #convert output to astropy table
    query_summary = create_table(query_summary,verbose)

    #query the detailed catalog
    query_detailed = hsccone(ra=RA,dec=DEC,radius=radius,table="detailed",release=release,format=format,magtype=magtype,
            columns=columns, verbose=verbose)
    #convert output to astropy table
    query_detailed = create_table(query_detailed,verbose)
    
    if query_summary:
        
        flag_col = []
        image_ID_col = []
        filter_col = []
        for ID in query_summary['MatchID']:
            #select the rows with the good ID
            mask_ID = query_detailed['MatchID'] == ID
    
            #add the extended flag to the dedicated list (take the most luminious occurence of the object)
            #also save the image ID
            query_detailed_temp = query_detailed[mask_ID]
            query_detailed_temp.sort('MagAuto')
            flag_col += [query_detailed_temp['Flags'][0]]
            image_ID_col += [np.array(query_detailed_temp['ImageID'])]
            filter_col += [np.array(query_detailed_temp['Filter'])]
            
        query_full = copy.deepcopy(query_summary)
        query_full['extFlag'] = flag_col
        query_full['ImageID'] = image_ID_col
        query_full['Filter'] = filter_col
    
    else:
        query_full = ''

    return query_full


def get_N_in_image(self,RA, Dec,release="v3",
                        format="csv",
                        magtype="magauto",
                        columns=None,
                        verbose=False, ):
    """
    This function count the number of galaxies in the image of each HSC objects compatible
    using the criterion : ext_flag = 1 to select extended object
    
    Then keep the number with only mag <= mag_object
    
    Input parameters
    ----------------
    query_transient_box: astropy table
         Results of the query centered on the transient position
    query_random: astropy table
         Results of the query centered on the randomly drawned positions, 
         or results of the query with big_radius
    """

    # Number of objects in the transient error box query
    N_obj_transient_box = len(self.query_transient_box)

    #Query a dummy big radius to be sure to get the whole image inside of it (image: 202*202 arcsec^2)
    dummy_radius = (202/3600)/np.cos(np.pi/4) + (1/3600) #~286 arcsec (radius > image diagonal)
    
    query_dummy = hsccone_full(RA=self.transient_RA,DEC=self.transient_Dec,radius=dummy_radius,release=release,format=format,magtype=magtype,
            columns=columns, verbose=verbose)
    
    Ncounts_tot = {}
    
    #fetch all filter apearance in the query_transient_box and suppress double
    all_filter_appearance = np.array([])
    
    for k in range(len(self.query_transient_box['Filter'])):
        all_filter_appearance = np.concatenate((all_filter_appearance,self.query_transient_box['Filter'][k]))
        
    uniq_filter = np.unique(all_filter_appearance)

    for filt in uniq_filter:
        Ncounts_tot[filt] = np.array([np.nan]*N_obj_transient_box)
        
    col_find_list = []
    for i in range(len(self.query_transient_box)):
        mask_extended = query_dummy['extFlag'] == 1
        
        #identifie the good image ID, take the first occurence of the filter (as this is already sorted in mag)
        list_filter, list_ind = np.unique(self.query_transient_box['Filter'][i],return_index=True)
        list_ImageID = self.query_transient_box['ImageID'][i][list_ind]

        for l in range(len(list_filter)):

            mask_image = []
            for m in range(len(query_dummy['ImageID'])):
                if list_ImageID[l] in query_dummy['ImageID'][m]:
                    mask_image += [True]
                else:
                    mask_image += [False]
          
            #the default column name for the mag is instrument_filtername so the following lines are a bit messy
            #in order to fetch instrument+filter name
            for col_name in self.query_transient_box.columns:
                if col_name.endswith(list_filter[l]):
                    col_find = col_name
                    col_find_list += [col_find]
                    #break
        
            mask_mag = query_dummy[col_find] <= self.query_transient_box[col_find ][i]

            mask_part = np.bitwise_and(mask_extended, mask_image)

            mask_tot = np.bitwise_and(mask_part, mask_mag)

            query_dummy_temp = copy.deepcopy(query_dummy[mask_tot])

            #The following piece of code supress from the big_radius query the object which are
            #whithin the 30 arcsec from the best transient localisation. (creating a "donuts")
            #This ensure not to bias the backgroud population estimation with the true host
            
            pos_center = SkyCoord(RA*u.deg, Dec*u.deg, frame='icrs')
            pos_obj = SkyCoord(query_dummy_temp['MatchRA']*u.deg, query_dummy_temp['MatchDec']*u.deg, frame='icrs')
            
            sep = pos_obj.separation(pos_center)
            
            mask_donuts = np.where( np.array(sep.arcsec) >= 30)
            
            query_dummy_temp = query_dummy_temp[mask_donuts]

            Ncounts_tot[list_filter[l]][i] = max(int(len(query_dummy_temp)),1)
            
    return Ncounts_tot, dummy_radius, col_find_list


def object_type_HSC(query_transient_box):
    """
    This function classifie as 'galaxy','star' or 'unknown' the compatible object
    using the flag present in the HSC catalog

    Input parameters
    ----------------
    """        
    
    type_list = []

    for row in query_transient_box:    
        if row['extFlag'] == 0:
            type_list.append('star')
        elif row['extFlag'] == 1:
            type_list.append('galaxy')
        else:
            type_list.append('unknown')

    return type_list


def download_GLADE(verbose=False):
    """
    This function download the GLADE catalog if not found locally

    """
    catalogFile = './catalogs/GLADE_2.4.txt'
    
    if not os.path.isfile(catalogFile):
        print("GLADE catalog not found locally, start the automatic download")
        url = 'http://glade.elte.hu/GLADE_2.4.txt'
        os.system("wget -O ./catalogs/GLADE_2.4.txt {}".format(url))  
    
    else:
        if verbose:
            print("GLADE catalog found locally")
        
    return
