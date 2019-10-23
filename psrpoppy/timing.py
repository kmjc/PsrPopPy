#!/usr/bin/env python2

"""
Author: Tyler Cohen
Affiliation: New Mexico tech
Email: tyler.cohen@student.nmt.edu
Date: 2/13/18
Version: 1.0.0

Description: Calculate the timing precision distribution of an
MSP population model
"""

from __future__ import print_function
import populate
import cPickle
import numpy as np
import matplotlib as mlp
import matplotlib.ticker as mtick
import matplotlib.pyplot as plt

SAVEDIR = "/Users/tcohen/Documents/New_Mexico_Tech/Research/NANOGrav/plots/psrpoppy/"

def main(savefig=False, plotType='npulsar'):
    """
    Plot timing precision distribution from sigma_TOA population model

    Parameters
    __________
    savefig : bool, optional (Default Value = False)
           Save plot output to file in SAVEDIR
    plotType : string (Default Value = 'npulsar')
            Type of distribution to plot. Valid options are 'npulsar
            which plots histogram of number of pulsars with timing 
            precision X or 'cumulative' which plots cumulative 
            distribution of pulsars with timing precision better than X.
    """

    sigma_toa = []
    popfile = open("populate.model_lorimer12", 'r')
    pop = cPickle.load(popfile)
    
    # Compute sigma TOA for whole population
    for pulsar in pop.population:
        if not pulsar.dead:
            sigtoa = sigma_TOA(pulsar)
            sigma_toa.append(sigtoa)

    # Construct histogram
    histdata, bins = np.histogram(np.log10(np.array(sigma_toa)), bins=100)
    width = 0.9 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2   # bin centers are avg of adjacent bins
    print("Binwidth = " + str(width))
    
    # Plot histogram
    plt.ion()
    ax = plt.subplot(111)
    ax.bar(center, histdata, align='center', width=width,
           color='b', alpha = 0.5, label=r'$t_{int}=1800$s, $BW$=1GHz,' +
           ' $S_{sys}=$3Jy')
    ax.set_yscale('log')
    ax.set_xlabel(r"$\sigma_{TOA}$ (ms)")
    ax.set_ylabel("# of pulsars")

    ax.set_title("MSP Timing-Precision Distribution")
    legend = ax.legend(loc='upper right')

    if savefig=True:
        plt.savefig(SAVEDIR + "timing_precision_distro.pdf")
    plt.show()

    return sigma_toa, histdata, bins


def sigma_TOA(pulsar, n_pol=2, t_int=1800., bandwidth=1e9, s_sys=3.,
              verbose=False):
    """
    Calculate timing uncertainty of a model pulsar due to
    radiometer noise and pulse jitter noise given survey
    parameters

    Paramters
    _________
    pulsar : pulsar.Pulsar object
          Model-generate pulsar object whose timing precision
          to be calculated
    n_pol : int, optional (Default Value = 2)
          Number of polarization states
    t_int : float, optional (Default Value = 1800)
          Integration time in seconds
    bandwidth : fload, optional (Default Value = 1GHz)
          Receiver bandwidth in Hz
    s_sys : float, optional (Default Value = 3)
          System equivalent flux density in Jy
    verbose : bool, optional (Default Value = False)
          Print computed values

    Returns
    _______
    sigma_toa : float
            Uncertainty in timing precision in ms
    """
    
    width_ms = pulsar.period * pulsar.width_degree / 360.
    s_1400 = pulsar.s_1400() / 1000. # convert from mJy to Jy

    # Calculate uncertainty due to radiometer noise
    survey_params = np.sqrt(n_pol * t_int * bandwidth) / s_sys
    duty_cycle = width_ms / pulsar.period
    snr = survey_params * s_1400 * np.sqrt((1 - duty_cycle) / duty_cycle)
    sigma_rn = width_ms / snr

    # Calculate uncertainty due to jitter
    n_pulses = (t_int * 1000.) / pulsar.period
    sigma_j = width_ms / np.sqrt(n_pulses)

    # TOA uncertainty is quadrature sum
    sigma_toa = np.sqrt(sigma_rn**2 + sigma_j**2)

    if verbose == True:
        print("\nPeriod = ", pulsar.period, "ms")
        print("Width = ", width_ms, "ms")
        print("S_1400 = ", s_1400, "Jy")
        print("SNR = ", snr)
        print("Sigma_RN = ", sigma_rn, "ms")
        print("Sigma_j = ", sigma_j, "ms")
        print("Sigma_TOA = ", sigma_toa, "ms")

    return sigma_toa

if __name__ == '__main__':
   main()
