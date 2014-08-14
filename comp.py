#!/usr/bin/python
#
# Lacerda@Saco - 4/Aug/2014
#
from check_files import check_files
import numpy as np
from pycasso import fitsQ3DataCube
import sys
import argparse as ap
import matplotlib as mpl
from matplotlib import pyplot as plt
from get_morfologia import get_morfologia
import os

mpl.rcParams['font.size']       = 16
mpl.rcParams['axes.labelsize']  = 22
mpl.rcParams['axes.titlesize']  = 26
mpl.rcParams['font.family']     = 'sans-serif'

Zsun = 0.019
Lsun = 3.826e33 # erg/s
qCCM = {
    '4861' : 1.16427,
    '5007' : 1.12022,
    '6563' : 0.81775,
    '6583' : 0.81466,
}

propStarlight = [
    dict(prop = 'at_flux__yx', title = r'$\langle \log\ t_\star \rangle_L (R)$', weiProp = 'LobnSD__yx'),
    dict(prop = 'at_mass__yx', title = r'$\langle \log\ t_\star \rangle_M (R)$', weiProp = 'McorSD__yx'),
    dict(prop = 'alogZ_flux__yx', title = r'$\langle \log\ Z_\star \rangle_L (R)$', weiProp = 'LobnSD__yx'),
    dict(prop = 'alogZ_mass__yx', title = r'$\langle \log\ Z_\star \rangle_M (R)$', weiProp = 'McorSD__yx'),
    dict(prop = 'v_d__yx', title = r'$\sigma_\star (R)$', weiProp = None),
    dict(prop = 'LobnSD__yx', title = r'$\mathcal{L}_\star$', weiProp = None),
    dict(prop = 'McorSD__yx', title = r'$\mu_\star$', weiProp = None),
]

propEL = [
    dict(prop = 'tau_V_neb__z', extensive = False, title = r'$\tau_V^{neb}$', weiProp = None),
    dict(prop = 'logZ_neb__z', extensive = False, title = r'$\log\ Z_{neb}$', weiProp = None),
    dict(prop = 'Ha_obs__z', extensive = True, title = r'$F_{H\alpha}^{obs}$', weiProp = None), 
    dict(prop = 'Hb_obs__z', extensive = True, title = r'$F_{H\beta}^{obs}$', weiProp = None),
    dict(prop = 'O3_obs__z', extensive = True, title = r'$F_{[OIII]}^{obs}$', weiProp = None),
    dict(prop = 'N2_obs__z', extensive = True, title = r'$F_{[NII]}^{obs}$', weiProp = None),
]

NStarlight = len(propStarlight)
NEL = len(propEL)
    
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

def parser_args(description):
    default = dict(outputImgSuffix = 'png',
                   vxx = 'v20',
                   outputDir = '.',
                   )
    
    parser = ap.ArgumentParser(description = description)
    parser.add_argument('--px1FitsFile', 
                        help = 'PyCASSO px1 FITS file',
                        metavar = 'FILENAME',
                        type = str,
                        default = None)
    parser.add_argument('--vxxFitsFile', 
                        help = 'PyCASSO vXX FITS file',
                        metavar = 'FILENAME',
                        type = str,
                        default = None)
    parser.add_argument('--px1EmLinesFitsFile', '-e',
                        help = 'Emission Lines px1 FITS File',
                        metavar = 'FILENAME',
                        type = str,
                        default = None)
    parser.add_argument('--vxxEmLinesFitsFile',
                        help = 'Emission Lines vXX FITS File',
                        metavar = 'FILENAME',
                        type = str,
                        default = None)
    parser.add_argument('--vxx',
                        help = 'Voronoi XX',
                        metavar = 'vXX',
                        type = str,
                        default = default['vxx'])
    parser.add_argument('--galaxyImgFile', '-g',
                        help = 'The image of the galaxy',
                        metavar = 'FILE',
                        type = str,
                        default = None)
    parser.add_argument('--outputImgSuffix', '-o',
                        help = 'Suffix of image file. Sometimes denote the image type. (Ex.: image.png)',
                        type = str,
                        default = default['outputImgSuffix'])
    parser.add_argument('--outputDir', '-d',
                        help = 'Image output directory',
                        metavar = 'DIR',
                        type = str,
                        default = default['outputDir'])

    return parser.parse_args()

def radialProfileWeighted(v__yx, w__yx, bins, rad_scale, func_radialProfile = None):
    v__r = None

    if func_radialProfile:
        w__r = func_radialProfile(w__yx, bin_r = bins, mode = 'sum', rad_scale = rad_scale)
        v_w__r = func_radialProfile(v__yx * w__yx, bin_r = bins, mode = 'sum', rad_scale = rad_scale)
        v__r = v_w__r / w__r

    return v__r

def Mpc_to_cm(dist):
    # google: 1 Mpc = 3.08567758e24 cm
    return dist * 3.08567758e24 

def tauV_to_AV(tauV):
    return tauV * np.log10(np.exp(1)) / 0.4

def flux_to_LSol(flux, distance):
    return 4. * np.pi * Mpc_to_cm(distance) ** 2.0 * flux / Lsun

def calc_Lint_Ha(Lobs__Lz, err_Lobs__Lz, tauVNeb__z,lines):
    i_Ha = lines.index('6563')
    i_Hb = lines.index('4861')
    
    q = qCCM['6563'] / (qCCM['4861'] - qCCM['6563'])
    
    eHa = np.ma.exp(qCCM['6563'] * tauVNeb__z)
    LobsHaHb = Lobs__Lz[i_Ha, :] / Lobs__Lz[i_Hb, :]

    Lint_Ha__z = Lobs__Lz[i_Ha, :] * eHa
    
    a = err_Lobs__Lz[i_Ha, :]
    b = q * LobsHaHb * err_Lobs__Lz[i_Hb, :]
    err_Lint_Ha__z = eHa * np.sqrt(a ** 2.0 + b ** 2.0)
    
    return Lint_Ha__z, err_Lint_Ha__z

def DrawHLRCircleInSDSSImage(ax, HLR_pix, pa, ba):
    from matplotlib.patches import Ellipse
    center , a , b_a , theta = np.array([ 256 , 256]) , HLR_pix * 512.0/75.0 , ba ,  pa*180/np.pi 
    e1 = Ellipse(center, height=2*a*b_a, width=2*a, angle=theta, fill=False, color='white',lw=2,ls='dotted')
    e2 = Ellipse(center, height=4*a*b_a, width=4*a, angle=theta, fill=False, color='white',lw=2,ls='dotted')
    ax.add_artist(e1)
    ax.add_artist(e2)

def plotRadialPropAxis(ax, x, prop__Vr, vxx):
    ax.plot(x, prop__Vr['px1'], 'b-', label = 'px1')
    ax.plot(x, prop__Vr['vxx'], 'r-', label = '%s' % vxx)
    leg = ax.legend()
    leg.get_frame().set_alpha(0.)

def plotNbyNprops(K, NRows, NCols, args):
    f, axArr = plt.subplots(NRows, NCols)
    f.set_size_inches(6 * NCols , 5 * NRows)
    
    for ax in f.axes:
        ax.set_axis_off()

    k = 0
    l = 0
    
    for i in range(0, NRows):
        for j in range(0, NCols):
            if i == 0 and j == 0:
                ax = axArr[0, 0]
                ax.set_axis_on()
                galimg = plt.imread(args.galaxyImgFile)[::-1,:,:]
                plt.setp(ax.get_xticklabels(), visible = False)
                plt.setp(ax.get_yticklabels(), visible = False)
                ax.imshow(galimg, origin = 'lower')
                pa_px1, ba_px1 = K.px1.getEllipseParams()
                DrawHLRCircleInSDSSImage(ax, K.px1.HLR_pix, pa_px1, ba_px1)
                ax.set_xticklabels([])
                ax.set_yticklabels([])
            elif k < NStarlight:
                p = propStarlight[k]
                ax = axArr[i, j] 
                ax.set_axis_on()
                ax.set_title(p['title'])
                prop__Vr = K.radialProp(p['prop'], weiProp = p['weiProp'])
                plotRadialPropAxis(ax, K.RbinCenter__r, prop__Vr, args.vxx)
                k += 1
            elif k >= NStarlight and l < NEL:
                p = propEL[l]
                ax = axArr[i, j] 
                ax.set_axis_on()
                ax.set_title(p['title'])
                prop__Vr = K.radialPropEL(p['prop'], weiProp = p['weiProp'], extensive = p['extensive'])
                plotRadialPropAxis(ax, K.RbinCenter__r, prop__Vr, args.vxx)
                l += 1
            else:
                continue
            
    f.savefig('%s.%s' % (K.px1.califaID, args.outputImgSuffix))

def plotImgRadProp(K, args):
    NCols = 4 

    for i in range(0, NStarlight + NEL):

        f, axArr = plt.subplots(1, NCols)
        f.set_size_inches(6 * NCols, 5)
        
        for ax in f.axes:
            ax.set_axis_off()

        ax = axArr[0]
        ax.set_axis_on()
        galimg = plt.imread(args.galaxyImgFile)[::-1,:,:]
        plt.setp(ax.get_xticklabels(), visible = False)
        plt.setp(ax.get_yticklabels(), visible = False)
        ax.imshow(galimg, origin = 'lower')
        pa_px1, ba_px1 = K.px1.getEllipseParams()
        DrawHLRCircleInSDSSImage(ax, K.px1.HLR_pix, pa_px1, ba_px1)
        ax.set_xticklabels([])
        ax.set_yticklabels([])

        if i < NStarlight:
            p = propStarlight[i]
            prop__yx = K.prop(p['prop'])
            prop__Vr = K.radialProp(p['prop'], weiProp = p['weiProp'])                
        else:
            k = NEL - i
            p = propEL[k]
            prop__yx = K.propELYX(p['prop'], extensive = p['extensive'])
            prop__Vr = K.radialPropEL(p['prop'], weiProp = p['weiProp'], extensive = p['extensive'])

        ax = axArr[1]
        ax.set_axis_on()
        im = ax.imshow(prop__yx['px1'], origin = 'lower')
        pa_px1, ba_px1 = K.px1.getEllipseParams()
        DrawHLRCircleInSDSSImage(ax, K.px1.HLR_pix, pa_px1, ba_px1)
        ax.set_title('%s px1' % p['title'])
        f.colorbar(ax = ax, mappable = im)

        ax = axArr[2]
        ax.set_axis_on()
        im = ax.imshow(prop__yx['vxx'], origin = 'lower')
        pa_vxx, ba_vxx = K.vxx.getEllipseParams()
        DrawHLRCircleInSDSSImage(ax, K.vxx.HLR_pix, pa_vxx, ba_vxx)
        ax.set_title('%s v20' % p['title'])
        f.colorbar(ax = ax, mappable = im)

        ax = axArr[3]
        ax.set_axis_on()
        ax.set_title('%s v20' % p['title'])
        plotRadialPropAxis(ax, K.RbinCenter__r, prop__Vr, args.vxx)

        f.tight_layout()
        f.savefig('%s-%s.%s' % (K.px1.califaID, p['prop'], args.outputImgSuffix))

class CALIFACompare(object):
    def __init__(self, px1File = None, vxxFile = None, px1ELFile = None, vxxELFile = None):
        self.px1 = None
        self.vxx = None
        
        self.px1 = fitsQ3DataCube(px1File)
        self.vxx = fitsQ3DataCube(vxxFile)

        if px1ELFile:
            self.px1.loadEmLinesDataCube(px1ELFile)
                
        if vxxELFile:
            self.vxx.loadEmLinesDataCube(vxxELFile)
            
        self.radial_init()

    def radial_init(self, RbinIni = 0.0, RbinFin = 2.0, RbinStep = 0.1):
        self.Rbin__r = np.arange(RbinIni, RbinFin + RbinStep, RbinStep)
        self.RbinCenter__r = (self.Rbin__r[:-1] + self.Rbin__r[1:]) / 2.0
        self.NRbins = len(self.RbinCenter__r)
            
    def prop(self, attrib):
        p = { 
             'px1' : self.px1.__getattribute__(attrib),
             'vxx' : self.vxx.__getattribute__(attrib)
        }
        return p  

    def propEL(self, attrib):
        p = { 
             'px1' : self.px1.EL.__getattribute__(attrib),
             'vxx' : self.vxx.EL.__getattribute__(attrib)
        }
        return p  

    def radial(self, prop, weiProp = None):
        HLR_pix = self.prop('HLR_pix')
        if weiProp:
            weiProp__yx = self.prop(weiProp)
            prop__r = {
                 'px1' : radialProfileWeighted(prop['px1'], weiProp__yx['px1'], 
                                               self.Rbin__r, HLR_pix['px1'], 
                                               self.px1.radialProfile),
                 'vxx' : radialProfileWeighted(prop['vxx'], weiProp__yx['vxx'], 
                                               self.Rbin__r, HLR_pix['vxx'],
                                               self.vxx.radialProfile),
            }
        else:
            prop__r = {
                 'px1' : self.px1.radialProfile(prop['px1'], self.Rbin__r, rad_scale = HLR_pix['px1']),
                 'vxx' : self.vxx.radialProfile(prop['vxx'], self.Rbin__r, rad_scale = HLR_pix['vxx']),
            }
        return prop__r

    def radialProp(self, attrib, weiProp = None):
        return self.radial(self.prop(attrib), weiProp)
    
    def propELYX(self, attrib, extensive):
        p__z = self.propEL(attrib)
        p__yx = {
            'px1' : self.px1.zoneToYX(p__z['px1'], extensive = extensive),
            'vxx' : self.vxx.zoneToYX(p__z['vxx'], extensive = extensive),
        }
        return p__yx

    def radialPropEL(self, attrib, weiProp = None, extensive = False):
        return self.radial(self.propELYX(attrib, extensive), weiProp)
    
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

if __name__ == '__main__':
    args = parser_args(sys.argv[0])
    if check_files(args.px1FitsFile, args.vxxFitsFile, args.px1EmLinesFitsFile, args.vxxEmLinesFitsFile):
        K = CALIFACompare(px1File = args.px1FitsFile, vxxFile = args.vxxFitsFile, 
                          px1ELFile = args.px1EmLinesFitsFile, vxxELFile = args.vxxEmLinesFitsFile)
                 
        # Setup elliptical-rings geometry
        pa_px1, ba_px1 = K.px1.getEllipseParams()
        K.px1.setGeometry(pa_px1, ba_px1)
        pa_vxx, ba_vxx = K.vxx.getEllipseParams()
        K.vxx.setGeometry(pa_vxx, ba_vxx)
          
        plotNbyNprops(K, NRows = 4, NCols = 4, args = args)
        plotImgRadProp(K, args)