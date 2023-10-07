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

`python ./src/Association_pval.py --transient_file ./Transient_sources.ecsv --catalogs AllWISE Pan-STARRS HSC --do_plots`

Please check the results in the result directory, the output ecsv file is not empty when at least one compatible object is found.

If everything went well up to here you can now change the Transient_sources.ecsv file according to your sources.
