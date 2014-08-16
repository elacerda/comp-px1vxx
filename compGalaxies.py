#!/usr/bin/python
#
# Lacerda@Saco - 15/Aug/2014
#
from check_files import check_files
import numpy as np
from pycasso import fitsQ3DataCube
import sys
import matplotlib as mpl
from matplotlib import pyplot as plt
from get_morfologia import get_morfologia
from comp import CALIFACompare
import os

debug = False
#debug = True

CALIFAWorkDir = '/Users/lacerda/CALIFA/'
    
#galaxiesListFile    = CALIFAWorkDir + 'listOf300GalPrefixes.txt'
galaxiesListFile    = CALIFAWorkDir + 'listALL.csv'

baseCode            = 'Bgsd6e'
versionSuffix       = 'v20_q043.d14a'

px1SuperFitsDir     = '/Volumes/backupzeira/CALIFA/q043/px1/'
vxxSuperFitsDir     = '/Volumes/backupzeira/CALIFA/q043/v20/' + baseCode + '/'
px1EmLinesFitsDir   = '/Volumes/backupzeira/CALIFA/q043/emlines/px1' + baseCode + '/'
vxxEmLinesFitsDir   = '/Volumes/backupzeira/CALIFA/q043/emlines/v20' + baseCode + '/'
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
    listOfPrefixes = listOfPrefixes[0:10]        # Q&D tests ...
    #listOfPrefixes = ['K0026\n']
    
N_gals = len(listOfPrefixes)

RbinIni = 0.0
RbinFin = 2.0
RbinStep = 0.1
Rbin__r = np.arange(RbinIni, RbinFin + RbinStep, RbinStep)
RbinCenter__r = (Rbin__r[:-1] + Rbin__r[1:]) / 2.0
NRbins = len(RbinCenter__r)

#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

if __name__ == '__main__':
    at_flux_px1__rG = np.ma.zeros((NRbins, N_gals))
    at_flux_vxx__rG = np.ma.zeros((NRbins, N_gals))
    at_mass_px1__rG = np.ma.zeros((NRbins, N_gals))
    at_mass_vxx__rG = np.ma.zeros((NRbins, N_gals))
    alogZ_flux_px1__rG = np.ma.zeros((NRbins, N_gals))
    alogZ_flux_vxx__rG = np.ma.zeros((NRbins, N_gals))
    alogZ_mass_px1__rG = np.ma.zeros((NRbins, N_gals))
    alogZ_mass_vxx__rG = np.ma.zeros((NRbins, N_gals))
    A_V_px1__rG = np.ma.zeros((NRbins, N_gals))
    A_V_vxx__rG = np.ma.zeros((NRbins, N_gals))
    McorSD_px1__rG = np.ma.zeros((NRbins, N_gals))
    McorSD_vxx__rG = np.ma.zeros((NRbins, N_gals))
    
    for iGal in np.arange(N_gals):
        galName         = listOfPrefixes[iGal][:-1]
        #galName = 'K0277'
        
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

        # Just to remember, this K is not a PyCASSO DataCube.        
        K = CALIFACompare(px1File = px1CALIFAFitsFile, vxxFile = vxxCALIFAFitsFile, 
                          px1ELFile = px1EmLinesFitsFile, vxxELFile = vxxEmLinesFitsFile)

        # Setup elliptical-rings geometry
        pa_px1, ba_px1 = K.px1.getEllipseParams()
        K.px1.setGeometry(pa_px1, ba_px1)
        pa_vxx, ba_vxx = K.vxx.getEllipseParams()
        K.vxx.setGeometry(pa_vxx, ba_vxx)
        
        HLR_pix = K.prop('HLR_pix')
        at_flux__Vr = K.radialProfile(K.prop('at_flux__yx'), rad_scale = HLR_pix, weiProp__yx = K.prop('LobnSD__yx'))         
        at_flux_px1__rG[:, iGal] = at_flux__Vr['px1']
        at_flux_vxx__rG[:, iGal] = at_flux__Vr['vxx']

        at_mass__Vr = K.radialProfile(K.prop('at_mass__yx'), rad_scale = HLR_pix, weiProp__yx = K.prop('McorSD__yx')) 
        at_mass_px1__rG[:, iGal] = at_mass__Vr['px1']
        at_mass_vxx__rG[:, iGal] = at_mass__Vr['vxx']
        
        alogZ_flux__Vr = K.radialProfile(K.prop('alogZ_flux__yx'), rad_scale = HLR_pix, weiProp__yx = K.prop('LobnSD__yx'))
        alogZ_flux_px1__rG[:, iGal] = alogZ_flux__Vr['px1']
        alogZ_flux_vxx__rG[:, iGal] = alogZ_flux__Vr['vxx'] 
        
        alogZ_mass__Vr = K.radialProfile(K.prop('alogZ_mass__yx'), rad_scale = HLR_pix, weiProp__yx = K.prop('McorSD__yx'))
        alogZ_mass_px1__rG[:, iGal] = alogZ_mass__Vr['px1']
        alogZ_mass_px1__rG[:, iGal] = alogZ_mass__Vr['vxx']
        
        A_V__Vr = K.radialProfile(K.prop('A_V__yx'), rad_scale = HLR_pix, weiProp__yx = None)
        A_V_px1__rG[:, iGal] = A_V__Vr['px1']
        A_V_vxx__rG[:, iGal] = A_V__Vr['vxx'] 
        
        McorSD__Vr = K.radialProfile(K.prop('McorSD__yx'), rad_scale = HLR_pix, weiProp__yx = None)
        McorSD_px1__rG[:, iGal] = McorSD__Vr['px1']
        McorSD_vxx__rG[:, iGal] = McorSD__Vr['vxx']
        
        K.px1.close()
        K.px1.EL.close()
        K.vxx.close()
        K.vxx.EL.close()
        
    McorSD_diff = McorSD_px1__rG.flatten() - McorSD_vxx__rG.flatten()
    alogZ_flux_diff = alogZ_flux_px1__rG.flatten() - alogZ_flux_vxx__rG.flatten()
    at_flux_diff = at_flux_px1__rG.flatten() - at_flux_vxx__rG.flatten()
    A_V_diff = A_V_px1__rG.flatten() - A_V_vxx__rG.flatten()
    
    f = plt.figure()
    f.set_size_inches(8,6)
    ax = f.gca()
    ax.scatter(np.log10(McorSD_diff), alogZ_flux_diff)
    ax.set_xlabel(r'$\mathcal{M}_\star^{\mathtt{px1}}(R) - \mathcal{M}_\star^{\mathtt{v20}}(R)$')
    ax.set_ylabel(r'$\langle \log\ Z_\star \rangle_L^{\mathtt{px1}} - \langle \log\ Z_\star \rangle_L^{\mathtt{v20}}$')
    textbox = dict(boxstyle = 'round', facecolor = 'wheat', alpha = 0.)
    txt = '$Xmean: %.2f X\sigma: %.2f\nYmean: %.2f Y\sigma: %.2f$' % (McorSD_diff.mean(), McorSD_diff.std(),
                                                                      alogZ_flux_diff.mean(), alogZ_flux_diff.std()) 
    ax.text(0.95, 0.95, txt, fontsize = 10, transform = ax.transAxes,
            va = 'top', ha = 'right', bbox = textbox)
    ax.grid()
    f.savefig('McorSD-alogZ_flux.png')