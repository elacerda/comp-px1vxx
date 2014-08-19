#!/usr/bin/python
#
# Lacerda@Saco - 15/Aug/2014
#
# To understand how works CALIFACompare object, read comp.py.
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

mpl.rcParams['font.size']       = 16
mpl.rcParams['axes.labelsize']  = 22
mpl.rcParams['axes.titlesize']  = 26
mpl.rcParams['font.family']     = 'sans-serif'

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

def plotDiff(x, y, xlabel, ylabel, fileName = None): 
    f = plt.figure()
    f.set_size_inches(12,8)
    ax = f.gca()
    ax.scatter(x, y)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    textbox = dict(boxstyle = 'round', facecolor = 'wheat', alpha = 0.)
    txt = u'$X:\ \mathtt{mean}\ %.2f\ \ \sigma\ %.2f$\n$Y:\ \mathtt{mean}\ %.2f\ \ \sigma\ %.2f$' % (x.mean(), x.std(), y.mean(), y.std()) 
    ax.text(0.95, 0.95, txt, fontsize = 12, transform = ax.transAxes,
            va = 'top', ha = 'right', bbox = textbox)
    ax.grid()
    f.savefig(fileName)

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
    LobnSD_px1__rG = np.ma.zeros((NRbins, N_gals))
    LobnSD_vxx__rG = np.ma.zeros((NRbins, N_gals))
    
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

        LobnSD__Vr = K.radialProfile(K.prop('LobnSD__yx'), rad_scale = HLR_pix, weiProp__yx = None)
        LobnSD_px1__rG[:, iGal] = LobnSD__Vr['px1']
        LobnSD_vxx__rG[:, iGal] = LobnSD__Vr['vxx']
        
        K.px1.close()
        K.px1.EL.close()
        K.vxx.close()
        K.vxx.EL.close()
    
    delta_logMcorSD = np.log10(McorSD_px1__rG.flatten()) - np.log10(McorSD_vxx__rG.flatten())
    delta_logLobnSD = np.log10(LobnSD_px1__rG.flatten()) - np.log10(LobnSD_vxx__rG.flatten())
    delta_alogZ_flux = alogZ_flux_px1__rG.flatten() - alogZ_flux_vxx__rG.flatten()
    delta_alogZ_mass = alogZ_mass_px1__rG.flatten() - alogZ_mass_vxx__rG.flatten()
    delta_at_flux = at_flux_px1__rG.flatten() - at_flux_vxx__rG.flatten()
    delta_at_mass = at_mass_px1__rG.flatten() - at_mass_vxx__rG.flatten()
    delta_A_V = A_V_px1__rG.flatten() - A_V_vxx__rG.flatten()
    
    propPlotDelta = [
        dict(x = delta_logMcorSD, y = delta_alogZ_mass,
             xlabel = r'$\delta \log\mathcal{M}_\star(R)$',
             ylabel = r'$\delta \langle \log\ Z_\star \rangle_M(R)$',
             fileName = 'logMcorSD-alogZ_mass.png'),
        dict(x = delta_logMcorSD, y = delta_alogZ_flux,
             xlabel = r'$\delta \log\mathcal{M}_\star(R)$',
             ylabel = r'$\delta \langle \log\ Z_\star \rangle_L(R)$',
             fileName = 'logMcorSD-alogZ_flux.png'),
        dict(x = delta_logMcorSD, y = delta_at_mass,
             xlabel = r'$\delta \log\mathcal{M}_\star(R)$',
             ylabel = r'$\delta \langle \log\ t_\star \rangle_M(R)$',
             fileName = 'logMcorSD-at_mass.png'),
        dict(x = delta_logMcorSD, y = delta_at_flux,
             xlabel = r'$\delta \log\mathcal{M}_\star(R)$',
             ylabel = r'$\delta \langle \log\ t_\star \rangle_L(R)$',
             fileName = 'logMcorSD-at_flux.png'),
        dict(x = delta_logMcorSD, y = delta_A_V,
             xlabel = r'$\delta \log\mathcal{M}_\star(R)$',
             ylabel = r'$\delta A_V(R)$',
             fileName = 'logMcorSD-A_V.png'),
        dict(x = delta_logLobnSD, y = delta_alogZ_mass,
             xlabel = r'$\delta \log\mathcal{L}_\star(R)$',
             ylabel = r'$\delta \langle \log\ Z_\star \rangle_M(R)$',
             fileName = 'logLobnSD-alogZ_mass.png'),
        dict(x = delta_logLobnSD, y = delta_alogZ_flux,
             xlabel = r'$\delta \log\mathcal{L}_\star(R)$',
             ylabel = r'$\delta \langle \log\ Z_\star \rangle_L(R)$',
             fileName = 'logLobnSD-alogZ_flux.png'),
        dict(x = delta_logLobnSD, y = delta_at_mass,
             xlabel = r'$\delta \log\mathcal{L}_\star(R)$',
             ylabel = r'$\delta \langle \log\ t_\star \rangle_M(R)$',
             fileName = 'logLobnSD-at_mass.png'),
        dict(x = delta_logLobnSD, y = delta_at_flux,
             xlabel = r'$\delta \log\mathcal{L}_\star(R)$',
             ylabel = r'$\delta \langle \log\ t_\star \rangle_L(R)$',
             fileName = 'logLobnSD-at_flux.png'),
        dict(x = delta_logLobnSD, y = delta_A_V,
             xlabel = r'$\delta \log\mathcal{L}_\star(R)$',
             ylabel = r'$\delta A_V(R)$',
             fileName = 'logLobnSD-A_V.png'),
        dict(x = delta_logLobnSD, y = delta_logMcorSD,
             xlabel = r'$\delta \log\mathcal{L}_\star(R)$',
             ylabel = r'$\delta \log\mathcal{M}_\star(R)$',
             fileName = 'logLobnSD-logMcorSD.png'),
        dict(x = delta_at_flux, y = delta_alogZ_flux,
             xlabel = r'$\delta \langle \log\ t_\star \rangle_L(R)$',
             ylabel = r'$\delta \langle \log\ Z_\star \rangle_L(R)$',
             fileName = 'at_flux-alogZ_flux.png'),
        dict(x = delta_at_mass, y = delta_alogZ_mass,
             xlabel = r'$\delta \langle \log\ t_\star \rangle_M(R)$',
             ylabel = r'$\delta \langle \log\ Z_\star \rangle_M(R)$',
             fileName = 'at_mass-alogZ_mass.png'),
        dict(x = delta_at_flux, y = delta_A_V,
             xlabel = r'$\delta \langle \log\ t_\star \rangle_L(R)$',
             ylabel = r'$\delta A_V(R)$',
             fileName = 'at_flux-A_V.png'),
    ]
    
    for p in propPlotDelta:
        plotDiff(p['x'], p['y'], p['xlabel'], p['ylabel'], p['fileName'])
