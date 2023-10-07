#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import copy
import numpy as np

"""
Created on Tue Jun 27 16:59:12 2023

@author: ducoin ducoin@iap.fr
"""

def define_plot_resolution():
    """
    This function allow to define the resolution of a matplotlib plot on a way
    wich is device independent. Put this before saving any of your plot to get
    homogeneous resolution.
    """

    fig = plt.gcf()  # get current figure

    DPI = fig.get_dpi()
    fig.set_size_inches(1920.0 / float(DPI), 1080.0 / float(DPI))

    return


def plot_dist_pval(self, outputdir):
    """
    This function plot the pval in fuction of the radius for the found objects
    (we don't plot objects with pval>0.95)
    """

    if self.table is None:
        if self.verbose:
            print("No plot done for this transient as the table is empty")
        return
    
    plt.ioff()

    Pan_STARRS_colors = {
        "gKronMag": "green",
        "rKronMag": "red",
        "iKronMag": "blue",
        "zKronMag": "maroon",
        "yKronMag": "yellow",
    }

    HSC_colors = [
        "teal",
        "lightskyblue",
        "grey",
        "black",
        "chartreuse",
        "darkblue",
        "indianred",
    ]

    AllWISE_colors = {
        "w1mpro" : "orange",
        "w2mpro" : "darksalmon",
        "w3mpro" : "red",
        "w4mpro" : "brown",
    }

    GLADE_colors = {
        "Bmag" : "BLUE",
        "Jmag" : "darksalmon",
        "Hmag" : "red",
        "Kmag" : "brown",
    }
    
    i = 0
    y_max_filter = []
    y_min_filter = []
    
    for filt, pval in self.pval_list.items():

        temp_table = copy.deepcopy(self.table)

        temp_table.sort("dist_to_center")

        if self.catalog == "Pan-STARRS":
            temp_col = []
            for i in range(len(temp_table[filt + "_pval"])):

                if temp_table[filt + "_pval"][i] == "None":
                    temp_col += [None]
                else:
                    temp_col += [float(temp_table[filt + "_pval"][i])]

            temp_table[filt + "_pval"] = temp_col

        mask_None = temp_table[filt + "_pval"] != None
        temp_table = temp_table[mask_None]
        mask_pval = temp_table[filt + "_pval"] <= 0.95
        temp_table = temp_table[mask_pval]

        if self.catalog == "Pan-STARRS":

            mask_not_star = temp_table["type"] != "star"
            temp_table_not_star = temp_table[mask_not_star]

            mask_gal = temp_table["type"] == "galaxy"
            temp_table_gal = temp_table[mask_gal]
            
            mask_unknown = temp_table["type"] == "unknown"
            temp_table_unknown = temp_table[mask_unknown]

        elif self.catalog == "HSC":

            mask_not_star = temp_table["type"] != "star"
            temp_table_not_star = temp_table[mask_not_star]
            
            mask_gal = temp_table["type"] == "galaxy"
            temp_table_gal = temp_table[mask_gal]
            
            mask_unknown = temp_table["type"] == "unknown"
            temp_table_unknown = temp_table[mask_unknown]   
            
        elif self.catalog == "AllWISE":

            mask_not_star = temp_table["type"] != "star"
            temp_table_not_star = temp_table[mask_not_star]
            
            mask_gal = temp_table["type"] == "galaxy"
            temp_table_gal = temp_table[mask_gal]
            
            mask_unknown = temp_table["type"] == "unknown"
            temp_table_unknown = temp_table[mask_unknown]       

        elif self.catalog == "GLADE":
            
            mask_not_star = temp_table["type"] != "star"
            temp_table_not_star = temp_table[mask_not_star]

            mask_gal = temp_table["type"] == "galaxy"
            temp_table_gal = temp_table[mask_gal]
            
            

        if self.catalog == "Pan-STARRS":

            #skip if there is no candidate in this filter
            if len(temp_table_not_star[filt + "_pval"]) :
                
                plt.plot(
                    temp_table_not_star["dist_to_center"],
                    temp_table_not_star[filt + "_pval"],
                    color=Pan_STARRS_colors[filt],
                    label=filt,
                )
    
                plt.scatter(
                    temp_table_gal["dist_to_center"],
                    temp_table_gal[filt + "_pval"],
                    color=Pan_STARRS_colors[filt],
                    marker=r"$\odot$",
                    s=150,
                )

                plt.scatter(
                    temp_table_unknown["dist_to_center"],
                    temp_table_unknown[filt + "_pval"],
                    color=Pan_STARRS_colors[filt],
                    marker='X',
                    s=120,
                )
                
                
                #register the max of each fiter for later ylim
                y_max_filter += [max(temp_table_not_star[filt + "_pval"])]
                y_min_filter += [min(temp_table_not_star[filt + "_pval"])]

        if self.catalog == "HSC":

            #skip if there is no candidate in this filter
            if len(temp_table_not_star[filt + "_pval"]) :
                
                plt.plot(
                    temp_table_not_star["dist_to_center"],
                    temp_table_not_star[filt + "_pval"],
                    color=HSC_colors[i],
                    label=filt,
                )
    
                plt.scatter(
                    temp_table_gal["dist_to_center"],
                    temp_table_gal[filt + "_pval"],
                    color=HSC_colors[i],
                    marker=r"$\odot$",
                    s=150,
                )
                
                plt.scatter(
                    temp_table_unknown["dist_to_center"],
                    temp_table_unknown[filt + "_pval"],
                    color=HSC_colors[i],
                    marker='X',
                    s=120,
                )          
                
                
                #register the max of each fiter for later ylim
                y_max_filter += [max(temp_table_not_star[filt + "_pval"])]
                y_min_filter += [min(temp_table_not_star[filt + "_pval"])]

        if self.catalog == "AllWISE":
            
            #skip if there is no candidate in this filter
            if len(temp_table_not_star[filt + "_pval"]) :
            
                plt.plot(
                    temp_table_not_star["dist_to_center"],
                    temp_table_not_star[filt + "_pval"],
                    color=AllWISE_colors[filt],
                    label=filt,
                )
    
                plt.scatter(
                    temp_table_gal["dist_to_center"],
                    temp_table_gal[filt + "_pval"],
                    color=AllWISE_colors[filt],
                    marker=r"$\odot$",
                    s=150,
                )
                
                plt.scatter(
                    temp_table_unknown["dist_to_center"],
                    temp_table_unknown[filt + "_pval"],
                    color=AllWISE_colors[filt],
                    marker="X",
                    s=120,
                )

                
                #register the max of each fiter for later ylim
                y_max_filter += [max(temp_table_not_star[filt + "_pval"])]
                y_min_filter += [min(temp_table_not_star[filt + "_pval"])]

        if self.catalog == "GLADE":
            
            #skip if there is no candidate in this filter
            if len(temp_table_not_star[filt + "_pval"]) :
    
                plt.scatter(
                    temp_table_gal["dist_to_center"],
                    temp_table_gal[filt + "_pval"],
                    color=GLADE_colors[filt],
                    marker=r"$\odot$",
                    s=150,
                    label=filt
                )
                
                #register the max of each fiter for later ylim
                y_max_filter += [max(temp_table_not_star[filt + "_pval"])]
                y_min_filter += [min(temp_table_not_star[filt + "_pval"])]


        i += 1



    #three last scatter for the legende
    plt.scatter(
        0,
        10,
        color='white',
        marker=r"$\odot$",
        s=150,
        label = '   ',
    )        
    
    plt.scatter(
        0,
        10,
        color='black',
        marker=r"$\odot$",
        s=150,
        label = 'galaxies',
    )

    plt.scatter(
        0,
        10,
        color='black',
        marker="X",
        s=120,
        label = 'unknown',
    )


    axes = plt.gca()
    axes.yaxis.grid(True, linestyle="--")

    axes.vlines(
        self.error_radius,
        0,
        10,
        linestyles="dashed",
        color="red",
        label="transient localisation error",
    )

    
    plt.legend()
    plt.xlabel("dist to the transient error box center [arcsec]", fontsize=30)
    plt.ylabel("pval", fontsize=30)
    plt.yscale("log")
    
    if len(y_min_filter):
        plt.ylim([0.9*min(y_min_filter), 1.1*max(y_max_filter)])
    define_plot_resolution()

    plt.savefig(outputdir +'plots/' + "R(pval)_" + self.transient + "_" + self.catalog + ".png")
    plt.close()
    
    return

