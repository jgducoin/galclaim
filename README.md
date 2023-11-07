# galclaim

GALaxy Chance of Local Alignment algorIthM
--------------------------------------------

The galclaim software is dedicated to identify association between astrophysical transient sources and host galaxy. This association is made by estimating the chance alignment between a given transient sky localisation and nearby galaxies.

The main dependencies are the following :

-numpy

-healpy

-astropy

-matplotlib


The current version of the code allows to use Pan-STARRS, HSC, AllWISE and GLADE catalogs. The code Also pre-check for nearby bright galaxy using the RC3 catalog (https://heasarc.gsfc.nasa.gov/w3browse/all/rc3.html). When a nearby galaxy is found, a warning is raised to the user and the properties of the galaxy are saved in a dedicated output file.

**User guide:**

First clone this git repertory:

`git clone https://gitlab.in2p3.fr/ducoin/galclaim.git`

We provide an example of input file format (Transient_sources.ecsv). You can simply check you installation by running on this example file :

`cd galclaim`

`python ./src/Association_pval.py --transient_file ./Transient_sources.ecsv --catalogs AllWISE Pan-STARRS HSC GLADE`

Please check the results in the result directory, the output ecsv file is not empty when at least one compatible object is found.

If everything went well up to here you can now change the Transient_sources.ecsv file according to your sources.

While using the GLADE catalog, one can add the possibility to check if the redshift of the host candidates is compatible with a given reshift range. Example:

`python ./src/Association_pval.py --transient_file ./Transient_sources.ecsv --catalogs AllWISE Pan-STARRS HSC GLADE --do_plots --do_redshift 0.01 0.05`

When enabled the --do_redshift command add a column in the GLADE outputted files saying if the redshift is compatible with the provided range.

One can anable the verbose option to increase the number of information returned in the console by the software :

`python ./src/Association_pval.py --transient_file ./Transient_sources.ecsv --catalogs AllWISE --verbose`

One can enable plots displaying the computed pval in fuction of the radius for the found objects for each transient and each catalog. No plot is produced if no candidate is found (the code don't plot objects with pval>0.95). An example with the example input file can be done running :

`python ./src/Association_pval.py --transient_file ./Transient_sources.ecsv --catalogs Pan-STARRS HST --do_plots`

The obtained plots are stored in the result/plots directory.
