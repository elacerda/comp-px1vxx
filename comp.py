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
    dict(prop = 'at_flux__yx', title = r'$\langle \log\ t_\star \rangle_L$', 
         weiProp = 'LobnSD__yx', log = False),
    dict(prop = 'at_mass__yx', title = r'$\langle \log\ t_\star \rangle_M$', 
         weiProp = 'McorSD__yx', log = False),
    dict(prop = 'alogZ_flux__yx', title = r'$\langle \log\ Z_\star \rangle_L$', 
         weiProp = 'LobnSD__yx', log = False),
    dict(prop = 'alogZ_mass__yx', title = r'$\langle \log\ Z_\star \rangle_M$', 
         weiProp = 'McorSD__yx', log = False),
    dict(prop = 'A_V__yx', title = r'$A_V$', weiProp = None, log = False),
    dict(prop = 'v_d__yx', title = r'$\sigma_\star$', weiProp = None, log = False),
    dict(prop = 'LobnSD__yx', title = r'$\log \mathcal{L}_\star$', weiProp = None, log = True),
    dict(prop = 'McorSD__yx', title = r'$\mathcal{M}_\star$', weiProp = None, log = False),
]

propEL = [
    dict(prop = 'tau_V_neb__z', extensive = False, title = r'$\tau_V^{neb}$', 
         weiProp = None, log = False),
    dict(prop = 'logZ_neb__z', extensive = False, title = r'$\log\ Z_{neb}$', 
         weiProp = None, log = False),
    dict(prop = 'Ha_obs__z', extensive = True, title = r'$\log F_{H\alpha}^{obs}$', 
         weiProp = None, log = True), 
    dict(prop = 'Hb_obs__z', extensive = True, title = r'$\log F_{H\beta}^{obs}$', 
         weiProp = None, log = True),
    dict(prop = 'O3_obs__z', extensive = True, title = r'$\log F_{[OIII]}^{obs}$', 
         weiProp = None, log = True),
    dict(prop = 'N2_obs__z', extensive = True, title = r'$\log F_{[NII]}^{obs}$', 
         weiProp = None, log = True),
]

NStarlight = len(propStarlight)
NEL = len(propEL)
Ntot = NStarlight + NEL
    
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
                        required = True,
                        default = None)
    parser.add_argument('--vxxFitsFile', 
                        help = 'PyCASSO vXX FITS file',
                        metavar = 'FILENAME',
                        type = str,
                        required = True,
                        default = None)
    parser.add_argument('--px1EmLinesFitsFile', '-e',
                        help = 'Emission Lines px1 FITS File',
                        metavar = 'FILENAME',
                        type = str,
                        required = True,
                        default = None)
    parser.add_argument('--vxxEmLinesFitsFile',
                        help = 'Emission Lines vXX FITS File',
                        metavar = 'FILENAME',
                        type = str,
                        required = True,
                        default = None)
    parser.add_argument('--galaxyImgFile', '-g',
                        help = 'The image of the galaxy',
                        metavar = 'FILE',
                        type = str,
                        required = True,
                        default = None)
    parser.add_argument('--vxx',
                        help = 'Voronoi XX',
                        metavar = 'vXX',
                        type = str,
                        default = default['vxx'])
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

def DrawHLRCircleInSDSSImage(ax, HLR_pix, pa, ba):
    from matplotlib.patches import Ellipse
    center , a , b_a , theta = np.array([ 256 , 256]) , HLR_pix * 512.0/75.0 , ba ,  pa*180/np.pi 
    e1 = Ellipse(center, height=2*a*b_a, width=2*a, angle=theta, fill=False, color='white',lw=2,ls='dotted')
    e2 = Ellipse(center, height=4*a*b_a, width=4*a, angle=theta, fill=False, color='white',lw=2,ls='dotted')
    ax.add_artist(e1)
    ax.add_artist(e2)
    
def DrawHLRCircle(ax, K, color='white', lw=1.5):
    from matplotlib.patches import Ellipse
    center , a , b_a , theta = np.array([ K.x0 , K.y0]) , K.HLR_pix , K.ba ,  K.pa*180/np.pi 
    e1 = Ellipse(center, height=2*a*b_a, width=2*a, angle=theta, fill=False, color=color,lw=lw,ls='dotted')
    e2 = Ellipse(center, height=4*a*b_a, width=4*a, angle=theta, fill=False, color=color,lw=lw,ls='dotted')
    ax.add_artist(e1)
    ax.add_artist(e2)

def plotNbyNprops(K, NRows, NCols, args):
    f, axArr = plt.subplots(NRows, NCols)
    f.set_size_inches(6 * NCols , 5 * NRows)
    
    for ax in f.axes:
        ax.set_axis_off()

    k = 0

    for i in range(0, NRows):
        for j in range(0, NCols):
            # Plot the SDSS Image
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
            else:
                if k < NStarlight:
                    p = propStarlight[k]
                    prop__yx = K.prop(p['prop'])
                elif k >= NStarlight and k < Ntot:
                    l = k - NStarlight
                    p = propEL[l]
                    prop__z = K.propEL(p['prop'])
                    prop__yx = K.zoneToYX(prop__z, extensive = p['extensive'])
                elif k >= Ntot:
                    continue

                if p['log']:
                    prop__yx = dict(px1 = np.log(prop__yx['px1']), vxx = np.log(prop__yx['vxx']))
                    
                ax = axArr[i, j]
                ax.set_axis_on()
                ax.set_title(p['title'])
                if p['weiProp']:
                    weiProp__yx = K.prop(p['weiProp'])
                HLR_pix = K.prop('HLR_pix')
                prop__Vr = K.radialProfile(prop__yx, rad_scale = HLR_pix, weiProp__yx = weiProp__yx)
                ax.plot(K.RbinCenter__r, prop__Vr['px1'], 'b-', label = 'px1')
                ax.plot(K.RbinCenter__r, prop__Vr['vxx'], 'r-', label = '%s' % args.vxx)
                diff__r = prop__Vr['px1'] - prop__Vr['vxx']
                textbox = dict(boxstyle = 'round', facecolor = 'wheat', alpha = 0.)
                txt = 'mean: %.2f\n$\sigma$: %.2f' % (diff__r.mean(), diff__r.std()) 
                ax.text(0.95, 0.95, txt, fontsize = 14, transform = ax.transAxes,
                        va = 'top', ha = 'right', bbox = textbox)
                if i == NRows - 1:
                    ax.set_xlabel('HLR') 
                #leg = ax.legend()
                #leg.get_frame().set_alpha(0.)
                k += 1
    f.suptitle('PX1 (blue) and V20 (red) Radial Profiles', fontsize = 24)            
    f.savefig('%s.%s' % (K.px1.califaID, args.outputImgSuffix))

def plotImgRadProp(K, args):
    NCols = 4 

    for i in range(0, Ntot):

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

        # Choose the property dictionary.
        if i < NStarlight:
            p = propStarlight[i]
            prop__yx = K.prop(p['prop'])
        else:
            k = i - NStarlight
            p = propEL[k]
            prop__z = K.propEL(p['prop'])
            prop__yx = K.zoneToYX(prop__z, extensive = p['extensive'])

        if p['log']:
            prop__yx = dict(px1 = np.log(prop__yx['px1']), vxx = np.log(prop__yx['vxx']))

        # Set vmax and vmin to the lower interval    
        vmax = prop__yx['px1'].max()
        vmin = prop__yx['px1'].min()

        if prop__yx['vxx'].max() < vmax:
            vmax = prop__yx['vxx'].max()
            
        if vmin < prop__yx['vxx'].min():
            vmin = prop__yx['vxx'].min()

        ax = axArr[1]
        ax.set_axis_on()
        im = ax.imshow(prop__yx['px1'], origin = 'lower', vmax = vmax, vmin = vmin)
        pa_px1, ba_px1 = K.px1.getEllipseParams()
        ax.set_title('%s px1' % p['title'])
        DrawHLRCircle(ax, K)
        f.colorbar(ax = ax, mappable = im)

        ax = axArr[2]
        ax.set_axis_on()
        im = ax.imshow(prop__yx['vxx'], origin = 'lower', vmax = vmax, vmin = vmin)
        pa_vxx, ba_vxx = K.vxx.getEllipseParams()
        ax.set_title('%s %s' % (p['title'], args.vxx))
        DrawHLRCircle(ax, K)
        f.colorbar(ax = ax, mappable = im)
        
        ax = axArr[3]
        ax.set_axis_on()
        ax.set_title('%s' % p['title'])
        HLR_pix = K.prop('HLR_pix')
        if p['weiProp']:
            weiProp__yx = K.prop(p['weiProp'])
        prop__Vr = K.radialProfile(prop__yx, rad_scale = HLR_pix, weiProp__yx = weiProp__yx)                
        ax.plot(K.RbinCenter__r, prop__Vr['px1'], 'b-', label = 'px1')
        ax.plot(K.RbinCenter__r, prop__Vr['vxx'], 'r-', label = '%s' % args.vxx)
        diff__r = prop__Vr['px1'] - prop__Vr['vxx']
        textbox = dict(boxstyle = 'round', facecolor = 'wheat', alpha = 0.)
        txt = 'mean: %.2f\n$\sigma$: %.2f' % (diff__r.mean(), diff__r.std()) 
        ax.text(0.95, 0.95, txt, fontsize = 14, transform = ax.transAxes,
                va = 'top', ha = 'right', bbox = textbox)
        #leg = ax.legend()
        #leg.get_frame().set_alpha(0.)
        
        f.tight_layout()
        f.savefig('%s-%s.%s' % (K.px1.califaID, p['prop'], args.outputImgSuffix))

class CALIFACompare(object):
    '''
    This object is build to manage two CALIFA Super FITS data cubes at same time. 
    One with all spectra sampled pixel by pixel, and another using Voronoi zones.
    The return of all properties is stored in a dictionary with px1 and vxx keys.
    '''
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
    
    # The Emission Lines properties are sampled by zones, so, this is why we have zoneToYX() too.
    def zoneToYX(self, prop__z, extensive):
        p__yx = {
            'px1' : self.px1.zoneToYX(prop__z['px1'], extensive = extensive),
            'vxx' : self.vxx.zoneToYX(prop__z['vxx'], extensive = extensive),
        }
        return p__yx

    def radialProfile(self, prop, rad_scale, weiProp__yx = None):
        if weiProp__yx:
            prop__r = {
                 'px1' : radialProfileWeighted(prop['px1'], weiProp__yx['px1'], 
                                               self.Rbin__r, rad_scale['px1'], 
                                               self.px1.radialProfile),
                 'vxx' : radialProfileWeighted(prop['vxx'], weiProp__yx['vxx'], 
                                               self.Rbin__r, rad_scale['vxx'],
                                               self.vxx.radialProfile),
            }
        else:
            prop__r = {
                 'px1' : self.px1.radialProfile(prop['px1'], self.Rbin__r, 
                                                rad_scale = rad_scale['px1']),
                 'vxx' : self.vxx.radialProfile(prop['vxx'], self.Rbin__r, 
                                                rad_scale = rad_scale['vxx']),
            }
        return prop__r
    
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