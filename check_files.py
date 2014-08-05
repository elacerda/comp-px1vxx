#!/usr/bin/python
#
# Lacerda@Saco - 4/Aug/2014
#
import os
import numpy as np

def check_files(f1, f2, f3, f4):
    load = 4
    if not os.path.isfile(f1):
        load -= 1
        print '%s: file not found' % f1
    if not os.path.isfile(f2):
        load -= 1
        print '%s: file not found' % f2
    if not os.path.isfile(f3):
        load -= 1
        print '%s: file not found' % f3
    if not os.path.isfile(f4):
        load -= 1
        print '%s: file not found' % f4
    if load < 4:
        return False
    return True

if __name__ == '__main__':
    CALIFAWorkDir = '/Users/lacerda/CALIFA/'

    galaxiesListFile    = CALIFAWorkDir + 'listALL.csv'
    
    baseCode            = 'Bgsd6e'
    versionSuffix       = 'v20_q043.d14a'
    
    px1SuperFitsDir     = '/Volumes/backupzeira/CALIFA/q043/px1/'
    vxxSuperFitsDir     = '/Volumes/backupzeira/CALIFA/q043/v20/' + baseCode + '/'
    px1EmLinesFitsDir   = '/Volumes/backupzeira/CALIFA/q043/emlines/px1' + baseCode + '/'
    vxxEmLinesFitsDir   = '/Volumes/backupzeira/CALIFA/q043/emlines/v20' + baseCode + '/'
    imgDir              = CALIFAWorkDir + 'images/'
    
    f                   = open(galaxiesListFile, 'r')
    listOfPrefixes      = f.readlines()
    f.close()
            
    N_gals = len(listOfPrefixes)
    
    for iGal in np.arange(N_gals):
        galName         = listOfPrefixes[iGal][:-1]
        
        px1CALIFASuffix    = '_synthesis_eBR_px1_q043.d14a512.ps03.k1.mE.CCM.' + baseCode + '.fits'
        px1CALIFAFitsFile  = px1SuperFitsDir + galName + px1CALIFASuffix
        vxxCALIFASuffix    = '_synthesis_eBR_' + versionSuffix + '512.ps03.k1.mE.CCM.' + baseCode + '.fits'
        vxxCALIFAFitsFile  = vxxSuperFitsDir + galName + vxxCALIFASuffix

        px1EmLinesSuffix   = '_synthesis_eBR_px1_q043.d14a512.ps03.k1.mE.CCM.' + baseCode + '.EML.MC100.fits'
        px1EmLinesFitsFile = px1EmLinesFitsDir + galName + px1EmLinesSuffix
        vxxEmLinesSuffix   = '_synthesis_eBR_' + versionSuffix + '512.ps03.k1.mE.CCM.' + baseCode + '.EML.MC100.fits'
        vxxEmLinesFitsFile = vxxEmLinesFitsDir + galName + vxxEmLinesSuffix
                                
        if not check_files(px1CALIFAFitsFile, vxxCALIFAFitsFile, px1EmLinesFitsFile, vxxEmLinesFitsFile):
            continue