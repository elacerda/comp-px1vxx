#!/usr/bin/python
#
# Lacerda@Saco - 4/Aug/2014
#
from morph_type import get_morph, morph_number
import numpy as np
from pycasso import fitsQ3DataCube
from get_morfologia import get_morfologia
import sys
import matplotlib as mpl
from matplotlib import pyplot as plt
from get_morfologia import get_morfologia
import os

plot = True
#plot = False
debug = False
#debug = True

#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

mpl.rcParams['font.size']       = 16
mpl.rcParams['axes.labelsize']  = 22
mpl.rcParams['axes.titlesize']  = 26
mpl.rcParams['font.family']     = 'sans-serif'

CALIFAWorkDir = '/Users/lacerda/CALIFA/'
    
galaxiesListFile    = CALIFAWorkDir + 'listOf300GalPrefixes.txt'
baseCode            = 'Bgsd6e'
versionSuffix       = 'px1_q043.d14a'
#versionSuffix       = 'v20_q043.d14a'
superFitsDir        = '/Volumes/backupzeira/CALIFA/q043/px1/'
#superFitsDir        = CALIFAWorkDir + 'gal_fits/' + versionSuffix + '/'

#emLinesFitsDir      = CALIFAWorkDir + 'superfits/' + versionSuffix + '/'
emLinesFitsDir      = CALIFAWorkDir + 'rgb-gas/' + versionSuffix + '/'
imgDir              = CALIFAWorkDir + 'images/'

Zsun = 0.019
Lsun = 3.826e33 # erg/s
qCCM = {
    '4861' : 1.16427,
    '5007' : 1.12022,
    '6563' : 0.81775,
    '6583' : 0.81466,
}

f               = open(galaxiesListFile, 'r')
listOfPrefixes  = f.readlines()
f.close()

if debug:
    #listOfPrefixes = listOfPrefixes[0:20]        # Q&D tests ...
    listOfPrefixes = ['K0026\n']
    
N_gals = len(listOfPrefixes)
