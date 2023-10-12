# galclaim

GALaxy Chance of Local Alignment algorIthM
--------------------------------------------

The galclaim software is dedicated to identify association between astrophysical transient sources and host galaxy. This association is made by estimating the chance alignment between a given transient sky localisation and nearby galaxies.

The main dependencies are the following :

-numpy
-healpy
-astropy
-matplotlib


The current version of the code allows to use Pan-STARRS, HSC and AllWISE catalogs.

**User guide:**
First clone this git repertory:

`git clone https://gitlab.in2p3.fr/ducoin/galclaim.git`

We provide an example of input file format (Transient_sources.ecsv). You can simply check you installation by running on this example file :

`cd galclaim`

`python ./src/Association_pval.py --transient_file ./Transient_sources.ecsv --catalogs AllWISE Pan-STARRS HSC GLADE --do_plots`

Please check the results in the result directory, the output ecsv file is not empty when at least one compatible object is found.

If everything went well up to here you can now change the Transient_sources.ecsv file according to your sources.

While using the GLADE catalog, one can add the possibility to check if the redshift of the host candidates is compatible with a given reshift range. Example:

`python ./src/Association_pval.py --transient_file ./Transient_sources.ecsv --catalogs AllWISE Pan-STARRS HSC GLADE --do_plots --do_redshift 0.01 0.05`

When enabled the --do_redshift command add a column in the GLADE outputted files saying if the redshift is compatible with the provided range.
