### MultiColorFits
### v2.1
### written by Phil Cigan
__author__ = "Phil Cigan <pcigan@gmu.edu>"
__version__ = "2.1.0"


#Some resources for now: 
#https://docs.enthought.com/traitsui/tutorials/traits_ui_scientific_app.html
#http://scipy-cookbook.readthedocs.io/items/EmbeddingInTraitsGUI.html
#http://www.scipy-lectures.org/advanced/traits/index.html#example
#https://gist.github.com/pierre-haessig/9838326
#
#(Examples)
#http://henrysmac.org/blog/2014/8/19/demo-of-enthoughts-traitsui-with-matplotlib-and-a-popup-menu.html
#Eric Tollerud's traitsUI fitting GUI  https://pythonhosted.org/PyModelFit/over.html#tutorial-examples
#
#Click interaction with matplotlib: http://scipy-cookbook.readthedocs.io/items/Matplotlib_Interactive_Plotting.html
#QT based toolbar  https://stackoverflow.com/questions/18468390/creating-matplotlib-toolbar-in-python-traits-editor



import numpy as np
try: import astropy.io.fits as pyfits; import astropy.wcs as pywcs; from astropy.wcs import WCS; 
except ImportError: import pyfits; import pywcs; from pywcs import WCS
from astropy.visualization import ZScaleInterval, LinearStretch, SqrtStretch, SquaredStretch, LogStretch, PowerStretch, PowerDistStretch, SinhStretch, AsinhStretch, ManualInterval
#from astropy import units as u
from scipy.stats import percentileofscore
import matplotlib
from pyface.api import FileDialog, OK
from pyface.qt import QtGui, QtCore


try:
    ### python3 and PyQt5
    matplotlib.use('Qt5Agg') 
    from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT
    from traitsui.qt4.editor import Editor 
    try: 
        #from traitsui.editor import Editor #--> No, non-qt4 version raises an error about 'set_size_policy'...
        from traitsui.basic_editor_factory import BasicEditorFactory
    except: 
        #from traitsui.qt4.editor import Editor 
        from traitsui.qt4.basic_editor_factory import BasicEditorFactory
    #These two for changing the cursor
    from PyQt5.QtCore import Qt #PyQt5 only for python3
    from PyQt5.QtWidgets import QApplication, QMenu #PyQt5 only for python3
except:
    ### PyQt4
    matplotlib.use('Qt4Agg') 
    from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT
    from traitsui.qt4.editor import Editor
    from traitsui.qt4.basic_editor_factory import BasicEditorFactory
    #These two for changing the cursor
    from PyQt4.QtCore import Qt
    from PyQt4.QtGui import QApplication, QMenu
    # Cursor shapes: http://ftp.ics.uci.edu/pub/centos0/ics-custom-build/BUILD/PyQt-x11-gpl-4.7.2/doc/html/qcursor.html

from matplotlib import pyplot as plt
from matplotlib.figure import Figure
from matplotlib.image import AxesImage
from matplotlib.axes import Axes
from matplotlib.widgets import AxesWidget
#from matplotlib.colors import to_hex #Not in mpl<2.0.  Use try statement below, otherwise define workaround

from traits.api import Any, Instance

from traits.api import HasTraits, Button, Instance, List, Str, Enum, Float, File, Any, Bool
from traits.api import Array, CInt, Int, CFloat, on_trait_change, Range, Dict, Button
from traitsui.api import View, Item, Group, VGroup, HSplit, CheckListEditor, HGroup, Handler, StatusItem, TextEditor, FileEditor, BooleanEditor, CodeEditor, Action, OKButton, CancelButton
from traitsui.ui_info import UIInfo
from traitsui.file_dialog import save_file

from traitsui.api import ColorTrait,ColorEditor

from skimage import color as ski_color
#from skimage.exposure import adjust_gamma  #Just use custom function, as this doesn't handle NaNs
from skimage.exposure import rescale_intensity
from skimage import filters

import colorsys

#For manual plots
import matplotlib.patheffects as PathEffects
from matplotlib import patches
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredEllipse, AnchoredSizeBar
try: from mpl_toolkits.axes_grid1.anchored_artists import AnchoredText  #Matplotlib <2.1
except: from matplotlib.offsetbox import AnchoredText                   #Matplotlib >=2.1
from matplotlib.offsetbox import AnchoredText
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredDrawingArea
from astropy.coordinates import SkyCoord

try: import reproject
except: print('multicolorfits: reproject package required for some options in task reproject2D()')
try: import kapteyn
except: print('multicolorfits: kapteyn package required for some options in task reproject2D()')


##########################################################################
##### ----- Put any custom color or colormap definitions here ----- #####
##########################################################################

#e.g.:  
#mycustomred='#C11B17'
#execfile('/path/to/custom/definitions/file.py-txt-whatever')


##########################################################################


scaling_fns={'linear':LinearStretch, 'sqrt':SqrtStretch, 'squared':SquaredStretch, 'log':LogStretch, 'power':PowerDistStretch, 'sinh':SinhStretch, 'asinh':AsinhStretch}

def adjust_gamma(array_in,gamma):
    """
    Replacement function for skimage.exposure.adjust_gamma, so that NaNs don't throw errors
    
    Parameters
    ----------
    array_in : array
        Input image array
    gamma : float  
        Gamma correction value
    
    Returns
    -------
    array
        Gamma-adjusted image values
    """
    return array_in**(float(gamma))

def drawProgressBar(percent, barlength = 20, prefix = '', suffix = ''):
    """
    This draws a percentage bar of the specified length to stdout\n
    
    Once you have calculated the compeltion percent in your task, call this function to update the display.\n
    
    (Alternatively, could use a module like tqdm or termcolor.colored)
    """
    #import sys
    sys.stdout.write("\r")
    ANSIred='\033[31m'; ANSIgreen='\033[32m'; ANSIyellow='\033[33m'; ANSIblue='\033[34m'; ANSIreset='\033[0m';
    progress = ''+ANSIgreen #Start with ANSI escape code for Green
    for i in range(barlength):
        if i < int(barlength * percent): progress += "=" #print (green) double dashes
        else: progress += ANSIblue+'-' #Change to ANSI escape code for blue, print single dashes
    progress+=ANSIreset #Reset style
    sys.stdout.write(prefix + "[ %s ] %s%.2f%%%s"%(progress,ANSIyellow,percent*100,ANSIreset) + suffix)
    sys.stdout.flush()
    
    #ANSI escape codes:
    #Foreground: black=30,red=31,green=32,yellow=33,blue=34,purple=35,cyan=36,white=37,reset=39
    #Background: black=40,red=41,green=42,yellow=43,blue=44,purple=45,cyan=46,white=47,reset=49
    #Style: NoEffect/ResetALL=0,Bold=1,Dim=2,Underline=4,Blink=5,Inverse=7,Hidden=8,NormalWidth=22
    #Mix & Match inside the \033[< >m  with commas ==> for example, \033[32,46,5m gives green,cyan bg, blink text

def nanpercofscore(array_in,score,**kwargs):
    """
    Usage is identical to scipy.stats.stats.percentileofscore(array_in,score), this just corrects for NaNs
    
    Parameters
    ----------
    array_in : array 
        input dataset
    score : float
        score for calculation
    **kwargs 
        keyword arguments to pass to scipy.stats.percentileofscore()
    
    Returns
    -------
    float, or array the length of score 
        percentile of the score
    """
    return percentileofscore(np.ma.masked_invalid(array_in).compressed(),score,**kwargs)

def force_hdr_to_2D(hdrin):
    """
    A simple function to take in a header and remove items related to 3D or 4D structure -- such as NAXIS3...
    
    Parameters
    ----------
    hdrin : astropy.io.fits.header
        Header object
    
    Returns
    -------
    astropy.io.fits.header
        Output 2D header
    """
    hdr2D=hdrin.copy()
    for item in ['CRVAL3','CRPIX3','CDELT3','CTYPE3','CUNIT3','NAXIS3', 'CRVAL4','CRPIX4','CDELT4','CUNIT4','CTYPE4', 'NAXIS4', 'PC03_01','PC04_01','PC03_02','PC04_02', 'PC01_03','PC02_03','PC03_03','PC04_03','PC01_04','PC02_04','PC03_04','PC04_04']: 
        #del hdr2D[item]
        try: hdr2D.remove(item) 
        except: pass
        hdr2D['NAXIS']=2
    return hdr2D

def force_hdr_floats(hdrin):
    """
    Helper function to force various header values to floats.
    Sometimes programs save header values as strings, which messes up the WCS...
    
    Parameters
    ----------
    hdrin : astropy.io.fits.header
        Header object
    
    Returns
    -------
    astropy.io.fits.header
        Output header
    """
    for fc in ['CRPIX1','CRPIX2','CRVAL1','CRVAL2','CDELT1','CDELT2','CD1_1','CD1_2','CD2_1','CD2_2', 'CROTA','CROTA2','EQUINOX']: 
        try: hdrin.set(fc,float(hdrin[fc])) 
        except: pass

def dec2sex(rain,decin,as_string=False,decimal_places=2):
    """
    Converts decimal coordinates to sexagesimal.
    
    Parameters
    ----------
    rain : float
        Input Right Ascension in decimal -- e.g.,  12.34567
    decin : float
        input Declination in decimal -- e.g. -34.56789
    as_string : bool
        Specifies whether to return output as a string (useful for making tables)
    decimal_places : int 
        Number of decimals places to use when as_string=True
    
    Returns
    -------
    list
        ['HH:MM:SS.ss', 'DD:MM:SS.ss']
    """
    rmins,rsec=divmod(24./360*rain*3600,60)
    rh,rmins=divmod(rmins,60)
    #dmins,dsec=divmod(decin*3600,60)
    #ddeg,dmins=divmod(dmins,60)
    ddeg=int(decin)
    if ddeg==0 and decin<0: ddeg=np.NZERO #Case of -0 degrees
    dmins=int(abs(decin-ddeg)*60)
    dsec=(abs((decin-ddeg)*60)-dmins)*60
    if as_string==True: 
        if ddeg==0 and decin<0: decdegstring='-00'
        else: decdegstring='{0:0>2d}'.format(ddeg)
        return ['{0:0>2d}:{1:0>2d}:{2:0>{4}.{3}f}'.format(int(rh),int(rmins),rsec, decimal_places,decimal_places+3), '{0}:{1:0>2d}:{2:0>{4}.{3}f}'.format(decdegstring,int(dmins),dsec, decimal_places,decimal_places+3)]
    else: return [int(rh),int(rmins),rsec],[ddeg,int(dmins),dsec]

def sex2dec(rain,decin):
    """
    Converts sexagesimal coordinates to decimal. HMS and DMS separated by colons (:)
    
    Parameters
    ----------
    rain : str 
        input Right Ascension as a sexagesimal string -- e.g.,  '03:45:6789'
    decin : str 
        input Declination as a sexagesimal string -- e.g.,  '-12:34:5678998765'
    
    Returns
    -------
    list
        [12.345678, -10.987654]
    """
    if ':' in rain: ra=[float(val)*360./24 for val in rain.split(':')]
    else: ra=[float(val)*360./24 for val in rain.split(' ')]
    raout=ra[0]+ra[1]/60.+ra[2]/3600.
    if ':' in decin: dec=[float(val) for val in decin.split(':')]
    else: dec=[float(val) for val in decin.split(' ')]
    if dec[0]<0 or '-' in decin.split(' ')[0]: decout=dec[0]-dec[1]/60.-dec[2]/3600.
    else: decout=dec[0]+dec[1]/60.+dec[2]/3600.
    return [raout,decout]

def getcdelts(hdrin,getrot=False):
    """
    Function to calculate CDELT1 and CDELT2 from the input header PCx_x cards.
    
    Parameters
    ----------
    hdrin : astropy.io.fits.header
        Header object
    getrot :bool
        Specifies whether to return the rotation value crota in the output
    
    Returns
    -------
    float, float (, float)
        CDELT1, CDELT2 (, CROTA)
    """
    try: 
        cdelt1=float(hdrin['CDELT1']); cdelt2=float(hdrin['CDELT2'])
        try: 
            #Checks for PC matrix, which will modify the CDELT values
            pc1_1=hdrin['PC1_1']; pc1_2=hdrin['PC1_2']; pc2_1=hdrin['PC2_1']; pc2_2=hdrin['PC2_2']; 
            cdelt1*=np.sqrt(pc1_1**2+pc1_2**2)*np.sign(pc1_1); 
            cdelt2*=np.sqrt(pc2_1**2+pc2_2**2)*np.sign(pc2_2);
            crota=np.arctan(pc1_2/pc1_1)*180./np.pi #CROTA(2) is in degrees, from North
        except: 
            try: crota=hdrin['CROTA2']
            except: 
                try: crota=hdrin['CROTA']  
                except: crota=0.
    except:
        try: 
            cd1_1=float(hdrin['CD1_1']); cd1_2=float(hdrin['CD1_2']); 
            cd2_1=float(hdrin['CD2_1']); cd2_2=float(hdrin['CD2_2']); 
        except: raise(Exception('Header does not contain CDELT2 or CD2_2 cards...'))
        cdelt1=np.sqrt(cd1_1**2+cd1_2**2)*np.sign(cd1_1); 
        cdelt2=np.sqrt(cd2_1**2+cd2_2**2)*np.sign(cd2_2);
        crota=np.arctan(cd1_2/cd1_1)*180./np.pi #CROTA(2) is in degrees, from North
    if getrot==False: return cdelt1,cdelt2
    else: return cdelt1,cdelt2,crota

def convsky2pix(headerin,rain,decin,precise=False,checksys=False,incoordsys='fk5', incoordequinox='J2000.0', forceimagesys=None, originindex=0):
    """
    Helper function to convert sky coordinates to pixel coordinates.  Now uses SkyCoord.
    
    Allowed SkyCoord frame systems: ['altaz', 'barycentrictrueecliptic', 'cirs', 'fk4', 'fk4noeterms', 'fk5', 'galactic', 'galactocentric', 'gcrs', 'geocentrictrueecliptic', 'hcrs', 'heliocentrictrueecliptic', 'icrs', 'itrs', 'precessedgeocentric', 'supergalactic']
    
    Parameters
    ----------
    headerin : astropy.io.fits.header
        Header object
    rain : float 
        Input Right Ascension, in decimal
    decin : float 
        Input Declination, in decimal
    precise : bool  
        False (default) to round to nearest integer pixel.  True to return fraction of pixel. 
    checksys : bool
        When True, checks that the input coordinate frame is the same as the header frame (e.g., fk5, icrs, etc.)
    incoordsys : str
        SkyCoord frame system, default = 'fk5'
    incoordequinox : str
        equinox, default = 'J2000.0'  
    forceimagesys : str
        SkyCoord frame system. User can specify a frame forceimagesys to use (in case no RADESYS in the header, or if it's known to be wrong...)
    originindex : int
        Default = 0.  Pixel origin to use for calculations (0 or 1)
    
    Returns
    -------
    list
        [X-pixel, Y-pixel]
    
    """
    if checksys==True:
        try: 
            try: headerframe=headerin['RADESYS']
            except: headerframe=headerin['RADECSYS']
        except: 
            if forceimagesys is not None: headerframe=forceimagesys
            else: raise(Exception('Input header has no valid RADESYS/RADECSYS frame, and forceimagesys not specified...'))
        if headerframe.lower()==incoordsys.lower(): pass
        else: 
            #Convert coordinates between RA/DEC systems
            incoordSC=SkyCoord(ra=rain, dec=decin, unit='deg', frame=incoordsys.lower(), equinox=incoordequinox).transform_to(headerframe.lower())
            rain=float(incoordSC.ra.value); decin=float(incoordSC.dec.value)
    
    try:
        wcstemp=pywcs.WCS(headerin)
        try:
            if headerin['NAXIS']==2: pixarr=wcstemp.wcs_world2pix([[rain,decin]],originindex)
            else: pixarr=wcstemp.wcs_world2pix([[rain,decin,0]],originindex) #wcstemp.wcs_sky2pix([[rain,decin,0]],0)
        except:
            if headerin['WCSAXES']==2: pixarr=wcstemp.wcs_world2pix([[rain,decin]],originindex)
            else: pixarr=wcstemp.wcs_world2pix([[rain,decin,0]],originindex) #wcstemp.wcs_sky2pix([[rain,decin,0]],0)
    except:
        wcstemp=pywcs.WCS(naxis=2)
        wcstemp.wcs.crpix=[headerin['CRPIX1'],headerin['CRPIX2']]
        try: wcstemp.wcs.cdelt=list(getcdelts(headerin))
        except: raise(Exception('Invalid WCS CDELTS'))
        wcstemp.wcs.crval=[headerin['CRVAL1'],headerin['CRVAL2']]
        wcstemp.wcs.ctype=[headerin['CTYPE1'],headerin['CTYPE2']]
        pixarr=wcstemp.wcs_world2pix([[rain,decin,0]],originindex) #wcstemp.wcs_sky2pix([[rain,decin,0]],0)
    if precise==True: return [pixarr[0][0],pixarr[0][1]]
    else: return [np.int(np.round(pixarr[0][0])),np.int(np.round(pixarr[0][1]))]

def convpix2sky(headerin,xin,yin,outcoordsys='same',outcoordequinox='J2000.0',forceimagesys=None, originindex=0):
    """
    Helper function to convert pixel coordinates to sky coordinates.  Now uses SkyCoord.
    
    Allowed SkyCoord frame systems: ['altaz', 'barycentrictrueecliptic', 'cirs', 'fk4', 'fk4noeterms', 'fk5', 'galactic', 'galactocentric', 'gcrs', 'geocentrictrueecliptic', 'hcrs', 'heliocentrictrueecliptic', 'icrs', 'itrs', 'precessedgeocentric', 'supergalactic']
    
    Parameters
    ----------
    headerin : astropy.io.fits.header
        Header object
    xin : float 
        Input x-axis pixel position
    yin : float 
        Input y-axis pixel position
    outcoordsys : str 
        'same' for frame in the header, or can alternatively specify a SkyCoord frame system.  
    outcoordequinox : str
        equinox, default = 'J2000.0'  
    forceimagesys : str
        SkyCoord frame system. User can specify a frame forceimagesys to use (in case no RADESYS in the header, or if it's known to be wrong...)
    originindex : int
        Default = 0.  Pixel origin to use for calculations (0 or 1)
    
    Returns
    -------
    list
        [RA_decimal, DEC_decimal]
    
    """
    try:
        wcstemp=pywcs.WCS(headerin)
        try:
            if headerin['NAXIS']==2: pixarr=wcstemp.wcs_pix2world([[xin,yin]],originindex)
            else: pixarr=wcstemp.wcs_pix2world([[xin,yin,0]],originindex) #wcstemp.wcs_pix2sky([[xin,yin,0]],0)
        except: 
            if headerin['WCSAXES']==2: pixarr=wcstemp.wcs_pix2world([[xin,yin]],originindex)
            else: pixarr=wcstemp.wcs_pix2world([[xin,yin,0]],originindex) #wcstemp.wcs_pix2sky([[xin,yin,0]],0)
    except: 
        wcstemp=pywcs.WCS(naxis=2)
        wcstemp.wcs.crpix=[headerin['CRPIX1'],headerin['CRPIX2']]
        wcstemp.wcs.crval=[headerin['CRVAL1'],headerin['CRVAL2']]
        wcstemp.wcs.ctype=[headerin['CTYPE1'],headerin['CTYPE2']]
        try: wcstemp.wcs.cdelt=list(getcdelts(headerin))
        except: raise(Exception('Invalid WCS CDELTS'))
        pixarr=wcstemp.wcs_pix2world([[xin,yin,0]],originindex) #wcstemp.wcs_pix2world([[xin,yin,0]],0)
    if outcoordsys=='same': return [pixarr[0][0],pixarr[0][1]]
    else: 
        try: 
            try: headerframe=headerin['RADESYS']
            except: headerframe=headerin['RADECSYS']
        except: 
            #User can specify a frame forceimagesys to use (in case no RADESYS in the header, or if it's known to be wrong...)
            if forceimagesys is not None: headerframe=forceimagesys
            else: raise(Exception('Input header has no valid RADESYS/RADECSYS frame, and forceimagesys not specified...'))
        outcoordSC=SkyCoord(ra=pixarr[0][0], dec=pixarr[0][1], unit='deg', frame=headerframe.lower(), equinox=outcoordequinox).transform_to(outcoordsys.lower())
        return [outcoordSC.ra.value,outcoordSC.dec.value]

def makesimpleheader(headerin,naxis=2,radesys=None,equinox=None,pywcsdirect=False):
    """
    Function to make a new 'simple header' from the WCS information in the input header.  
    
    Parameters
    ----------
    headerin : astropy.io.fits.header
        Header object
    naxis : int 
        Specifies how many axes the final header should have.  Default=2
    radesys :str 
        RA/DEC system to use (valid SkyCoord frame system, e.g. 'icrs')
    equinox : str
        Equinox to use for the output header
    pywcsdirect : bool
        True to create the header directly with astropy.wcs.WCS 
    
    Returns
    -------
    astropy.io.fits.header
        Output header
    """
    if type(headerin)==str:
        import astropy.io.fits as pyfits
        headerin=pyfits.getheader(headerin)
    if pywcsdirect==True: wcstemp=pywcs.WCS(header=headerin)
    else: 
        wcstemp=pywcs.WCS(naxis=naxis); 
        if naxis>2:
            wcstemp.wcs.crpix=[float(headerin['CRPIX1']), float(headerin['CRPIX2']), float(headerin['CRPIX3'])]
            wcstemp.wcs.crval=[float(headerin['CRVAL1']), float(headerin['CRVAL2']), float(headerin['CRVAL3'])]
            wcstemp.wcs.ctype=[headerin['CTYPE1'],headerin['CTYPE2'],headerin['CTYPE3']]
            try: wcstemp.wcs.cunit=[headerin['CUNIT1'],headerin['CUNIT2'],headerin['CUNIT3']]
            except: pass
            try: wcstemp.wcs.cdelt=list(getcdelts(headerin))+[headerin['CDELT3']]; 
            except: raise(Exception('Invalid WCS CDELTS'))
        else: 
            wcstemp.wcs.crpix=[float(headerin['CRPIX1']),float(headerin['CRPIX2'])]
            wcstemp.wcs.crval=[float(headerin['CRVAL1']),float(headerin['CRVAL2'])]
            wcstemp.wcs.ctype=[headerin['CTYPE1'],headerin['CTYPE2']]
            try: wcstemp.wcs.cunit=[headerin['CUNIT1'],headerin['CUNIT2']]
            except: pass
            try: wcstemp.wcs.cdelt=list(getcdelts(headerin)); 
            except: raise(Exception('Invalid WCS CDELTS'))
        try: crota=getcdelts(headerin,getrot=True)[-1] #degrees, from N
        except: raise(Exception('Invalid WCS params for CROTAx'))
        #if crota!=0.: wcstemp.wcs.crota=[crota]*2 #Header will include PC_x cards if crot not 0
        try: wcstemp.wcs.radesys=headerin['RADESYS']
        except: pass
        try: wcstemp.wcs.equinox=headerin['EQUINOX']
        except: pass
    if radesys is not None: wcstemp.wcs.radesys=radesys; #e.g. 'FK5', 'ICRS'. For manually forcing string, not true reprojection.
    if equinox is not None: wcstemp.wcs.equinox=equinox; #e.g. 2000.0
    simpleheader=wcstemp.to_header()
    if pywcsdirect==False: 
        if crota!=0.: simpleheader['CROTA2']=crota #Alternative method to just use (deprecated) CROTA2 card
    simpleheader['NAXIS']=naxis; 
    try: simpleheader['NAXIS1']=int(headerin['NAXIS1']); simpleheader['NAXIS2']=int(headerin['NAXIS2']);
    except: pass
    if naxis>2: 
        for card in ['NAXIS3','CRPIX3','CRVAL3','CDELT3','CTYPE3','CUNIT3', 'SPECSYS','ALTRVAL','ALTRPIX']:
            try: simpleheader[card]=headerin[card]
            except: pass 
    for card in ['CROTA','CROTA1','CROTA2','BSCALE','BZERO','ZSCALE','BMAJ','BMIN','BPA', 'JANSCALE','FLUXCONV', 
                 'WAVELEN','FREQ', 'RESTFRQ', 'LATPOLE','LONPOLE']:
        try: simpleheader[card]=float(headerin[card])
        except: pass
    for card in ['BUNIT','OBJECT','TELESCOP','ZUNITS','SPECSYS']:
        try: simpleheader[card]=headerin[card]
        except: pass
    return simpleheader

def getcdmatrix(hdrin,crot=None):
    """
    Calculate the CDn_m matrix from CDELTS and CROTA/CROTA2
    
    Parameters
    ----------
    hdrin : astropy.io.fits.header
        Header object
    crot : float  
        Rotation in degrees, if you know it but the header doesn't correctly have it
    
    Returns
    -------
    float, float, float, float
        CD1_1, CD1_2, CD2_1, CD2_2
    """
    try: cd1_1=hdrin['CD1_1']; cd1_2=hdrin['CD1_2']; cd2_1=hdrin['CD2_1']; cd2_2=hdrin['CD2_2']; 
    except: 
        if crot is None:
            try: crot=hdrin['CROTA2']
            except: 
                try:crot=hdrin['CROTA']
                except: crot=0.
        try: cdelt1=float(hdrin['CDELT1']); cdelt2=float(hdrin['CDELT2']) 
        except: raise(Exception('Header does not contain CDELT2 or CD2_2 cards...'))
        try: 
            pc1_1=hdrin['PC1_1']; pc1_2=hdrin['PC1_2']; pc2_1=hdrin['PC2_1']; pc2_2=hdrin['PC2_2']; 
            cd1_1=pc1_1*cdelt1; cd1_2=pc1_2*cdelt1; cd2_1=pc2_1*cdelt2; cd2_2=pc2_2*cdelt2
        except: 
            cd1_1=cdelt1*np.cos(crot*np.pi/180); cd1_2=cdelt1*np.sin(crot*np.pi/180); 
            cd2_1=cdelt2*-np.sin(crot*np.pi/180); cd2_2=cdelt2*np.cos(crot*np.pi/180)
        #Even if PC1_1 matrix used, CDELTs should still be included there...
    return cd1_1,cd1_2,cd2_1,cd2_2

def getdegperpix(hdrin):
    """
    Calculates degrees per pixel side.  Assumes input header CDELTs are in degrees.
    
    Parameters
    ----------
    hdrin : astropy.io.fits.header
        Header object
    Returns
    -------
    float
        Degrees per pixel side
    """
    cdelt1,cdelt2=getcdelts(hdrin)
    if abs(cdelt1)-abs(cdelt2)<1e-6: return abs(cdelt2)
    else: raise(Exception('CDELT1 and CDELT2 are significantly different!'))

def getasecperpix(hdrin):
    """
    Calculates arcseconds per pixel side.  Assumes input header CDELTs are in degrees.
    
    Parameters
    ----------
    hdrin : astropy.io.fits.header
        Header object
    Returns
    -------
    float
        Arcseconds per pixel side
    """
    return getdegperpix(hdrin)*3600

def getsteradperpix(hdrin):
    """
    Calculates steradians per pixel.  Assumes input header CDELTs are in degrees.
    
    Parameters
    ----------
    hdrin : astropy.io.fits.header
        Header object
    Returns
    -------
    float
        Steradians per pixel
    """
    cdelt1,cdelt2=getcdelts(hdrin)
    if abs(cdelt1)-abs(cdelt2)<1e-6: return (cdelt2*np.pi/180.)**2
    else: raise(Exception('CDELT1 and CDELT2 are significantly different!'))

def beampars_asec_fromhdr(hdrin): 
    """
    Returns beam parameters [BMAJ (arcsec), BMIN (arcsec), BPA (deg)] from input header
    
    Parameters
    ----------
    hdrin : astropy.io.fits.header
        Header object
    Returns
    -------
    list
        [BMAJ_asec, BMIN_asec, BPA_deg]
    """
    return [hdrin['BMAJ']*3600.,hdrin['BMIN']*3600.,hdrin['BPA']] #Return [BMAJ_asec,BMIN_asec,BPA_deg]

def pixperbeam_from_hdr(hdrin):
    """
    Calculates the number of pixels per beam from the beam parameters (BMAJ,BMIN) in the header.
    Beam area = 2*PI*sigma_maj*sigma_min   = 2*PI*FWHM_maj*FWHM_min/(sqrt(8*ln(2)))**2  = PI*FWHM1*FWHM2/(4*ln(2))
    That's in whatever units the FWHM are in, which is degrees in the case of hdrin['BMAJ'], so use CDELTS to get area in pixels
    Requires valid BMAJ,BMIN header cards, where BMAJ,BMIN are in degrees
    Note that the scaling factor is 1.13, not 2*pi, because BMAJ/BMIN are FWHM, not sigma
    
    Parameters
    ----------
    hdrin : astropy.io.fits.header
        Header object
    Returns
    -------
    float
        Pixels per beam
    """
    return np.pi/(4.*np.log(2)) * hdrin['BMAJ']*hdrin['BMIN']/getdegperpix(hdrin)**2

def reproject2D(mapin,hdrfrom,hdrto,scale=False,method='interp',interpdict={'order':1,'mode':'constant','cval':np.nan}, returnfootprint=False, parallel=True):
    """
    Function that reprojects a 2D map from the parameters in one header to the parameters in another header.
    
    Parameters
    ----------
    mapin : array 
        Input map / data array (np.ndarray). Usually would get this from astropy.io.fits.getdata(...)
    hdrfrom : astropy.io.fits.header
        Header for the original image, specifying the original pixel size, etc.
    hdrto : astropy.io.fits.header
        Header to reproject to (usually just a second image's header)
    scale : bool
        True if units are Flux [W/m^2, Jy, or similar].  False for brightness [W/m^2/sr, Jy/sr, or similar]
        --> Note that when convolving maps in beam units [e.g., Jy/beam], the reprojection will need scale=True because the beam sizes change.
    method : str 
        One of 'kapteyn' (kapteyn package, interpolation), 'interp' (reproject.py), or 'spi' (reproject package, spherical polygon intersection).  Drizzle not yet implemented
    interpdict : dict 
        For method='kapteyn'. Sets the interpol_dict (interpolation). interpdict={order:<>,mode:<>,cval:<>}
        --> order = spline order, 0 to 5; mode='constant','nearest','reflect','wrap'; cval = value outside bounds (NaN)
    returnfootprint : bool  
        True to also return the footprint of which reprojected pixels fell on original image grid
    parallel : bool 
        True to use parallel processing (method='spi' only)
    
    Returns
    -------
    array (, array)
        Reprojected map (, optional footprint)
    """
    if scale==True:
        # --> Jy_2 = Jy_1_reproj*([sr/pix_2]/[sr/pix_1])
        # Set to True for [Jy/pix, W/m2/pix, Jy/beam if beam has been convolved]
        mapin_brightness=mapin.copy()/getsteradperpix(hdrfrom)
    else: mapin_brightness=mapin.copy() 
    if method=='interp': 
        from reproject import reproject_interp
        map_reproj,map_footprint = reproject_interp((mapin_brightness,hdrfrom),hdrto)
    elif method=='spi': 
        from reproject import reproject_exact
        map_reproj,map_footprint = reproject_exact((mapin_brightness,hdrfrom),hdrto,parallel=parallel)
    else: 
        from kapteyn import maputils
        hdrfrom=hdrfrom.copy(); hdrto=hdrto.copy() #Otherwise it modifies original.  (Bad for reproject3D calls)
        if hdrto['NAXIS']>2: hdrto['NAXIS']=2
        if hdrfrom['NAXIS']>2: hdrfrom['NAXIS']=2
        map_in=maputils.FITSimage(externaldata=mapin_brightness,externalheader=hdrfrom)
        map_reproj=map_in.reproject_to(hdrto,interpol_dict=interpdict).dat; 
    if scale==True: map_reproj*=getsteradperpix(hdrto)
    if returnfootprint==True and method!='kapteyn': return map_reproj,map_footprint
    else: return map_reproj

def reproject3D(mapin,hdrfrom,hdrto,scale=False,method='kapteyn',parallel=True,returnfootprint=False, print_progress=False):
    """
    Function that reprojects a 3D cube from the parameters in one header to the parameters in another header.
    
    Calls reproject2D for each slice in the cube.  
    
    See reproject2D for argument desscriptions.
    """
    #Scale: True if units are Flux [W/m^2 or similar].  False for brightness [W/m^2/sr or similar]
    tmpcube=np.zeros(mapin.shape[-3]).astype(np.ndarray); tmpcube_foot=tmpcube.copy()
    if returnfootprint==True and method=='spi': 
        for zz in range(mapin.shape[-3]): 
            if print_progress==True: 
                drawProgressBar(float(zz)/mapin.shape[-3],prefix='  Reprojecting channels ',suffix='  %i of %i'%(zz+1,mapin.shape[-3]))
            tmpcube[zz],tmpcube_foot[zz]=reproject2D(mapin[zz,:,:],hdrfrom,hdrto, scale=scale, method=method, returnfootprint=True)
        if print_progress==True: drawProgressBar(1.,prefix='  Reprojecting channels ',suffix='  %i of %i    \n'%(zz+1,mapin.shape[-3]))
        return np.array([tmpcube[zz] for zz in range(mapin.shape[-3])]),np.array([tmpcube_foot[zz] for zz in range(mapin.shape[-3])])
    else: 
        for zz in range(mapin.shape[-3]): 
            if print_progress==True: 
                drawProgressBar(float(zz)/mapin.shape[-3],prefix='  Reprojecting channels ',suffix='  %i of %i'%(zz+1,mapin.shape[-3]))
            tmpcube[zz]=reproject2D(mapin[zz,:,:],hdrfrom,hdrto, scale=scale, method=method, returnfootprint=False)
        if print_progress==True: drawProgressBar(1.,prefix='  Reprojecting channels ',suffix='  %i of %i    \n'%(zz+1,mapin.shape[-3]))
        return np.array([tmpcube[zz] for zz in range(mapin.shape[-3])])

def cropfits2D(datain,hdrin,xbounds,ybounds,newref=None,savenew=False,overwrite=False):
    """
    Function to crop a 2D fits image to the specified pixel bounds.  
    
    Parameters
    ----------
    datain : array
        Input fits image array
    hdrin : astropy.io.fits.header 
        Input header
    xbounds :list 
        [min,max] x-axis pixel limits to use for new image slice
    ybounds :list 
        [min,max] y-axis pixel limits to use for new image slice
    newref : None or str 
        None, 'center', or 'origin'.  None (default) keeps the reference pixel sky coordinate the same.  
        'center' forces the reference pixel to be the new center, 'origin' forces it to the new origin.
    savenew : bool or str
        False (default) does not save, otherwise specify a save path (e.g. savenew='./mynewfits.fits') to save the crop to disk  
    overwrite : bool 
        Input option to astropy.io.fits.writeto(..., overwrite=False)
    
    Returns
    -------
    array, astropy.io.fits.header
        cropdata, crophdr
    """
    #xbounds, ybounds are pixel values formatted as list: [xlo,xhi] and [ylo,yhi]
    #savenew= (path+)filename if desired to save new cropped fits file
    if len(datain.shape)>3: raise(Exception('Data array has 4 or more dimensions - reduce them to 2!'))
    if len(datain.shape)==2: pass
    elif len(datain.shape)>2 and (datain.shape[i]==1 for i in range(len(datain.shape)-2)): datain=datain[0,:,:]
    else: raise(Exception('File is not flattened 2D!  Use cropfits3D.'))
    if True in [item<0 for item in np.concatenate([xbounds,ybounds])]: 
        raise(Exception('Specified X/Y crop bounds are negative - must be positive.'))
    try: cropdata=datain[ybounds[0]:ybounds[1]+1,xbounds[0]:xbounds[1]+1]
    except: cropdata=datain[int(np.round(ybounds[0])):int(np.round(ybounds[1]+1)), int(np.round(xbounds[0])):int(np.round(xbounds[1]+1))]
    crophdr=hdrin.copy()
    crophdr['NAXIS2'],crophdr['NAXIS1']=cropdata.shape
    if newref=='center':
        crophdr['CRPIX1']=(int(xbounds[1]-xbounds[0]))*.5; crophdr['CRPIX2']=(int(ybounds[1]-ybounds[0]))*.5;
        crophdr['CRVAL1'],crophdr['CRVAL2']=convpix2sky(hdrin,(int(xbounds[0])+int(xbounds[1]-xbounds[0])*.5)-1,(int(ybounds[0])+int(ybounds[1]-ybounds[0])*.5)-1) 
        ###Next, adjust CROTA2 based on offset - needed for big deviations near poles
        #radiff=crophdr['CRVAL1']-hdrin['CRVAL1'] #In this order because RA increases to left
        #if radiff<=-180.: radiff+=360.
        #elif radiff>=180.: radiff-=360.
        radiff=np.unwrap(np.deg2rad([0,hdrin['CRVAL1']-crophdr['CRVAL1']]))[1]*180./np.pi #Equivalent of above 3 lines
        crophdr['CROTA2']=radiff*np.sin(crophdr['CRVAL2']*np.pi/180) #sin because at equator, DEC=0
        #--> For large changes in CRPIX, there will be some slight offsets/position errors.
    elif newref=='origin':
        crophdr['CRPIX1']=0.; crophdr['CRPIX2']=0.;
        crophdr['CRVAL1'],crophdr['CRVAL2']=convpix2sky(hdrin,int(xbounds[0])-1,int(ybounds[0])-1); 
        ###Next, adjust CROTA2 based on offset - needed for big deviations near poles
        radiff=np.unwrap(np.deg2rad([0,hdrin['CRVAL1']-crophdr['CRVAL1']]))[1]*180./np.pi
        crophdr['CROTA2']=radiff*np.sin(crophdr['CRVAL2']*np.pi/180) #sin because at equator, DEC=0
        #print(xbounds[0],ybounds[0])
        ### There's still a fraction of an arcsec offset for 'center' and 'origin'....
    else: crophdr['CRPIX1']=hdrin['CRPIX1']-int(xbounds[0]); crophdr['CRPIX2']=hdrin['CRPIX2']-int(ybounds[0]); 
    crophdr['NAXIS']=2
    if type(savenew)==str: 
        try: 
            try: pyfits.writeto(savenew,cropdata,crophdr,overwrite=overwrite)
            except: pyfits.writeto(savenew,cropdata,crophdr,overwrite=overwrite)
        except: print('Could not save file to path supplied.')
    return cropdata,crophdr

def cropfits3D(datain,hdrin,xbounds,ybounds, zbounds=[0,None], newref=None, savenew=False, overwrite=False):
    """
    Function to crop a 3D fits cube to the specified pixel bounds.  
    
    zbounds (optional) are channel number boundaries [zmin,zmax], 0-indexed
    
    Other arguments the same as for cropfits2D
    """
    if len(datain.shape)>3:
        if (datain.shape[i]==1 for i in range(len(datain.shape)-1)): datain=datain[0,:,:,:]
        else: raise(Exception('Data array has 4 or more dimensions - reduce them to 3!'))
    elif len(datain.shape)<3: raise(Exception('File does not have 3 dimensions!  Use cropfits2D.'))
    else: pass
    if zbounds[1] is None: zbounds[1]=datain.shape[-3]+1
    cropdata=datain[int(zbounds[0]):int(zbounds[1]), int(np.round(ybounds[0])):int(np.round(ybounds[1]))+1, int(np.round(xbounds[0])):int(np.round(xbounds[1]))+1]
    crophdr=hdrin.copy()
    crophdr['NAXIS3'],crophdr['NAXIS2'],crophdr['NAXIS1']=cropdata.shape; crophdr['CRPIX3']-=int(zbounds[0])
    if newref=='center':
        crophdr['CRPIX1']=(int(xbounds[1]-xbounds[0]))*.5; crophdr['CRPIX2']=(int(ybounds[1]-ybounds[0]))*.5;
        crophdr['CRVAL1'],crophdr['CRVAL2']=convpix2sky(hdrin,(int(xbounds[0])+int(xbounds[1]-xbounds[0])*.5)-1,(int(ybounds[0])+int(ybounds[1]-ybounds[0])*.5)-1) 
        ###Next, adjust CROTA2 based on offset - needed for big deviations near poles
        radiff=np.unwrap(np.deg2rad([0,hdrin['CRVAL1']-crophdr['CRVAL1']]))[1]*180./np.pi
        crophdr['CROTA2']=radiff*np.sin(crophdr['CRVAL2']*np.pi/180) #sin because at equator, DEC=0
    elif newref=='origin':
        crophdr['CRPIX1']=0.; crophdr['CRPIX2']=0.;
        crophdr['CRVAL1'],crophdr['CRVAL2']=convpix2sky(hdrin,int(xbounds[0])-1,int(ybounds[0])-1); 
        ###Next, adjust CROTA2 based on offset - needed for big deviations near poles
        radiff=np.unwrap(np.deg2rad([0,hdrin['CRVAL1']-crophdr['CRVAL1']]))[1]*180./np.pi
        crophdr['CROTA2']=radiff*np.sin(crophdr['CRVAL2']*np.pi/180) #sin because at equator, DEC=0
        ### There's still a fraction of an arcsec offset for 'center' and 'origin'....
    else: crophdr['CRPIX1']=hdrin['CRPIX1']-int(xbounds[0]); crophdr['CRPIX2']=hdrin['CRPIX2']-int(ybounds[0]);
    crophdr['NAXIS']=3
    if type(savenew)==str: 
        try: 
            try: pyfits.writeto(savenew,cropdata,crophdr,overwrite=overwrite)
            except: pyfits.writeto(savenew,cropdata,crophdr,overwrite=overwrite)
        except: print('Could not save file to path supplied.')
    return cropdata,crophdr

def cropfits2D_coords(datain,hdrin,centerRADEC,radius_asec,newref=None,savenew=False,overwrite=False, return_cropcenterpix=False, coords_in='dec', checksys=False,incoordsys='fk5',incoordequinox='J2000.0', forceimagesys=None):
    """
    Function to crop a 2D fits image based on sky coordinates and specified width.  
    
    Parameters
    ----------
    datain : array 
        Input fits image array
    hdrin : astropy.io.fits.header 
        Input header
    centerRADEC : list 
        List of RA,DEC coords.  e.g. [12.345678,-10.98765]  
    radius_asec : float 
        Radius (box half-width) of new image, in arcseconds
    newref : str
        None, 'center', or 'origin'.  None keeps the reference pixel sky coordinate the same.  
        'center' forces the reference pixel to be the new center, 'origin' forces it to the new origin.
    savenew : bool or str
        False (default) does not save, otherwise specify a save path (e.g. savenew='./mynewfits.fits') to save the crop to disk 
    overwrite : bool
        Input option to astropy.io.fits.writeto(..., overwrite=False)
    return_cropcenterpix : bool  
        True will return the cropped image center pixel location.
    coords_in : str 
        'dec' for decimal (default) or 'sex' for sexagesimal.  
    checksys : bool
        Passes to convsky2pix().  When True, checks that the input coordinate frame is the same as the header frame (e.g., fk5, icrs, etc.)
    incoordsys : str
        Passes to convsky2pix().  SkyCoord frame system, default = 'fk5'
    incoordequinox : str
        Passes to convsky2pix().  equinox, default = 'J2000.0'  
    forceimagesys : str
        Passes to convsky2pix().  SkyCoord frame system. User can specify a frame forceimagesys to use (in case no RADESYS in the header, or if it's known to be wrong...)
    
    Returns
    -------
    array, astropy.io.fits.header
        cropdata, crophdr  
    """
    if 'sex' in coords_in.lower(): centerRADEC=sex2dec(*centerRADEC)
    centerpix=convsky2pix(hdrin,centerRADEC[0],centerRADEC[1],precise=True, checksys=checksys,incoordsys=incoordsys, incoordequinox=incoordequinox, forceimagesys=forceimagesys)
    pixextent=(radius_asec/3600)/getcdelts(hdrin)[1] #Total crop image width is center+/- radius_asec
    cropdat,crophdr=cropfits2D(datain,hdrin, [centerpix[0]-pixextent,centerpix[0]+pixextent],[centerpix[1]-pixextent,centerpix[1]+pixextent], newref=newref,savenew=savenew,overwrite=overwrite)
    if return_cropcenterpix==True: 
        centerpixcrop=convsky2pix(crophdr,centerRADEC[0],centerRADEC[1],precise=True)
        return cropdat,crophdr,centerpixcrop
    else: return cropdat,crophdr

def cropfits3D_coords(datain,hdrin,centerRADEC,radius_asec,zbounds=[0,None],newref=None,savenew=False, overwrite=False, coords_in='dec',  return_cropcenterpix=False):
    """
    Function to crop a 3D fits cube based on sky coordinates and specified width. 
    
    zbounds (optional) are channel number boundaries [zmin,zmax], 0-indexed
    
    Other arguments the same as for cropfits2D_coords()
    """
    if 'sex' in coords_in.lower(): centerRADEC=sex2dec(*centerRADEC)
    if zbounds[1] is None: zbounds[1]=datain.shape[0]+1
    centerpix=convsky2pix(hdrin,centerRADEC[0],centerRADEC[1],precise=True)
    pixextent=(radius_asec/3600)/getcdelts(hdrin)[1] #Total crop image width is center+/- radius_asec
    cropdat,crophdr=cropfits3D(datain,hdrin, [centerpix[0]-pixextent,centerpix[0]+pixextent],[centerpix[1]-pixextent,centerpix[1]+pixextent],zbounds=zbounds, newref=newref, savenew=savenew, overwrite=overwrite)
    if return_cropcenterpix==True: 
        centerpixcrop=convsky2pix(crophdr,centerRADEC[0],centerRADEC[1],precise=True)
        return cropdat,crophdr,centerpixcrop
    else: return cropdat,crophdr

def smooth_image(datain, sigma=3):
    """
    Simple wrapper to call skimage.filters.gaussian(), to perform simple image
    smoothing. This can be useful for e.g., matching another image resolution, 
    smoothing over noise to boost apparent dynamic range, etc.
    
    Parameters
    ----------
    datain : array 
        Input fits image array.  Can either be a 2D image or a 3D (multichannel)
        set of images.
    sigma : int or float
        The sigma value (in pixels) for Gaussian smoothing
    
    Returns
    -------
    array
        smoothed_image
    """
    if len(datain.shape)==2: return filters.gaussian(datain, sigma=sigma)
    else: return filters.gaussian(datain, sigma=sigma, multichannel=True)

def hex_to_rgb(hexstring):
    """
    Converts a hexadecimal string to RGB tuple 
    
    Parameters
    ----------
    hexstring : str
        Hexadecimal string such as '#FFFFFF'
    Returns
    -------
    tuple
        RGB tuple such as (256,256,256)
    """
    #From http://stackoverflow.com/a/214657
    hexstring = hexstring.lstrip('#')
    lv = len(hexstring)
    return tuple(int(hexstring[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))

def rgb_to_hex(rgb): 
    """
    Converts RGB tuple to a hexadecimal string
    
    Parameters
    ----------
    rgb : tuple
        RGB tuple such as (256,256,256)
    Returns
    -------
    str
        Hexadecimal string such as '#FFFFFF'
    """    
    return '#%02x%02x%02x'%rgb #input RGB as tuple.  e.g.: rgb_to_hex((255, 255, 255))

def hex_to_hsv(hexstring):
    """
    Convert a hexadecimal string to HSV tuple
    
    Parameters
    ----------
    hexstring : str
        Hexadecimal string such as '#3300FF'
    Returns
    -------
    tuple
        HSV tuple such as (0.7,1.,1.)
    """
    #See Wikipedia article for details -- https://en.wikipedia.org/wiki/HSL_and_HSV#From_HSV  .  Modified from colorsys.py,
    # HSV: Hue, Saturation, Value
    # H = position in the spectrum, S = color saturation ("purity"), V = color brightness
    r,g,b=np.array(hex_to_rgb(hexstring))/255. #Convert Hex to RGB fracs (i.e., 0..1 instead of 0..255)
    maxc = max(r,g,b); minc = min(r,g,b); 
    s = (maxc-minc) / maxc
    v = maxc
    if minc == maxc: return 0.0, 0.0, v
    rc = (maxc-r) / (maxc-minc);   gc = (maxc-g) / (maxc-minc);   bc = (maxc-b) / (maxc-minc)
    if r == maxc: h = bc-gc
    elif g == maxc: h = 2.0+rc-bc
    else: h = 4.0+gc-rc
    h = (h/6.0) % 1.0
    return h, s, v #All in fractions in range [0...1]

def hexinv(hexstring):
    """
    Convenience function to calculate the inverse color (opposite on the color wheel).
    e.g.: hexinv('#FF0000') = '#00FFFF'
    
    Parameters
    ----------
    hexstring : str
        Hexadecimal string such as '#FF0000'
    Returns
    -------
    str
        Hexadecimal string such as '#00FFFF'
    """
    if hexstring[0]=='#': hexstring=hexstring.lstrip('#')
    hexcolor=int(hexstring,16)
    color_comp=0xFFFFFF^hexcolor
    hexcolor_comp="#%06X"%color_comp
    return hexcolor_comp

try: from matplotlib.colors import to_hex #Converts matplotlib colors to hex/html.  Only in mpl>=2.0
except:
    def to_hex(c):
        #For mpl v<2.0
        import matplotlib.colors as mplc
        c_hex=mplc.rgb2hex(mplc.colorConverter.to_rgb(c)) #could use to_rgba to keep alpha
        return c_hex

def greyRGBize_image(datin,rescalefn='linear',scaletype='abs',min_max=[None,None], gamma=2.2, checkscale=False):
    """
    ### Takes an image and returns 3-frame [R,G,B] (vals from 0...1)
    
    Parameters
    ----------
    datin : array 
        Input 2D image data array
    rescalefn : func 
        Function to use for rescaling intensity.  imscale.linear/sqrt/squared/log/power/sinh/asinh
    scaletype : str 
        'abs' for absolute values, 'perc' for percentiles
    min_max : list
        [min,max] vals to use in rescale.  if scaletype='perc', list the percentiles to use, e.g. [1.,95.]
    gamma : float 
        Value for gamma correction.  For combining colorized frames, use default gamma=2.2.  For inverse, use gamma=(1./2.2)
    checkscale : bool  
        True to bring up plot to check the new image scale.
    
    Returns
    -------
    array
        Greyscale RGB image, shape=[ypixels,xpixels,3]
    """
    if 'per' in scaletype.lower(): 
        if min_max==[None,None]: min_max=[0.,100.]
        minval,maxval=np.percentile(np.ma.masked_invalid(datin).compressed(),min_max)
    else: 
        minval=[np.nanmin(datin) if min_max[0] is None else min_max[0]][0]
        maxval=[np.nanmax(datin) if min_max[1] is None else min_max[1]][0]
    #Used specified rescaling function
    datscaled=(scaling_fns[rescalefn]() + ManualInterval(vmin=minval,vmax=maxval))(datin)
    #datscaled=rescalefn(datin,vmin=minval,vmax=maxval)
    if gamma!=1: datscaled=adjust_gamma(datscaled,gamma)
    #Need to scale image between -1 and 1 if data type is float...
    datlinear=LinearStretch()(datscaled)
    #datlinear=imscale.linear(np.nan_to_num(datscaled))
    #Convert to RGB
    dat_greyRGB=ski_color.gray2rgb(datlinear)
    
    if checkscale is not False: 
        plt.clf(); plt.close('all')
        fig0=plt.figure(0); 
        ax1=fig0.add_subplot(121); 
        plt.imshow(datin,interpolation='nearest',origin='lower',cmap='gist_gray'); 
        plt.title('Input Image')
        ax2=fig0.add_subplot(122); 
        plt.imshow(dat_greyRGB**(1./gamma),interpolation='nearest',origin='lower'); 
        plt.title('Scaled Image')
        plt.show(); #plt.clf(); plt.close('all')

    return dat_greyRGB
  
def colorize_image(image, colorvals, colorintype='hsv',dtype=np.float64,gammacorr_color=1):
    """
    ### Add color of the given hue to an RGB greyscale image.
    
    Parameters
    ----------
    image : array
        Greyscale RGB image -- as would be output from greyRGBize_image()
    colorvals : str or list or tuple 
        color values to apply to image.  e.g., '#FF0000' if colorintype='hex'
    colorintype : str 
        'hsv' for [0..1,0..1,0..1],  'rgb' for [0..255,0..255,0..255], or 'hex' for '#XXXXXX'
    dtype : dtype 
        Defaults to standard numpy float, but may want to lower to e.g. float32 for large images (>~1000x1000)
    gammacorr_color : float
        To use color as-is, leave as 1 (default).  To gamma-correct color at this step (e.g., to match gamma for checking a scaled image), specify a factor
    
    Returns
    -------
    array
        Colorized RGB image, shape=[ypixels,xpixels,3]
    """
    if colorintype not in ['hsv', 'hsv_dict', 'rgb', 'hex']: raise Exception("  colorintype must be 'hsv', 'hsv_dict', 'rgb', or 'hex'")
    hsv = ski_color.rgb2hsv(image).astype(dtype)
    if colorintype.lower()=='rgb': colorvals=np.array(hex_to_hsv(rgb_to_hex(colorvals))).astype(dtype)
    elif colorintype.lower()=='hex': colorvals=np.array(hex_to_hsv(colorvals)).astype(dtype) #from custom_colormaps.py
    if colorintype.lower()=='hsv_dict': hue,saturation,v=colorvals['hue'],colorvals['sat'],colorvals['v'],
    else: hue,saturation,v=colorvals
    if gammacorr_color!=1: 
        hue,saturation,v = colorsys.rgb_to_hsv( *np.array( colorsys.hsv_to_rgb(hue,saturation,v) )**gammacorr_color )
    hsv[:, :, 2] *= v 
    hsv[:, :, 1] = saturation
    hsv[:, :, 0] = hue
    return ski_color.hsv2rgb(hsv).astype(dtype)

def combine_multicolor(im_list_colorized,gamma=2.2,inverse=False):
    """
    Combines input colorized RGB images [:,:,3] into one intensity-rescaled RGB image
    
    Parameters
    ----------
    im_list_colorized : list 
        List of colorized RGB images.  e.g., [ halpha_purple, co21_orange, sio54_teal ]
    gamma : float 
        Value used for gamma correction ^1/gamma.  Default=2.2.  
    inverse : bool  
        True will invert the scale so that white is the background
    
    Returns
    -------
    array
        Colorized RGB image (combined), shape=[ypixels,xpixels,3]
    """
    combined_RGB=LinearStretch()(np.nansum(im_list_colorized,axis=0))
    if inverse==True: RGB_maxints=tuple(1.-np.nanmax(combined_RGB[:,:,i]) for i in [0,1,2])
    else: RGB_maxints=tuple(np.nanmax(combined_RGB[:,:,i]) for i in [0,1,2])
    for i in [0,1,2]: 
        combined_RGB[:,:,i]=np.nan_to_num(rescale_intensity(combined_RGB[:,:,i], out_range=(0, combined_RGB[:,:,i].max()/np.max(RGB_maxints) )));
    combined_RGB=LinearStretch()(combined_RGB**(1./gamma)) #gamma correction
    if inverse==True: combined_RGB=1.-combined_RGB #gamma correction
    return combined_RGB

def plotsinglemulticolorRGB(multicolorin,hdrin,axtitle,savepath, xaxislabel='RA', yaxislabel='DEC', tickcolor='w', labelcolor='k', facecolor='w', minorticks=True, dpi=150):
    """
    Convenience function to plot and save a single multicolor RGB image.
    
    For tick/label colors: Any standard matplotlib color format, such as 'k', 'black', '#000000', or '0.0' 
    
    Parameters
    ----------
    multicolorin : array 
        Multicolor RGB image -- as would be output from combine_multicolor()
    hdrin : astropy.io.fits.header 
        Header to use for WCS information
    axtitle : str 
        String to use for plot title.  e.g. "My crazy 5-color image", or empty quotes "" for nothing.
    savepath : str 
        Path to save file to.  e.g., "./plots/mycrazyimage.pdf"
    xaxislabel: str
        Label to use for x-axis. Default='RA'
    yaxislabel: str
        Label to use for y-axis. Default='DEC'
    tickcolor : str  
        Color to use for ticks in plot, default = 'w'. 
    labelcolor : str 
        Color to use for ticklabels, default = 'k'
    facecolor : str 
        Color to use for figure facecolor (border around plot), default ='w'
    minorticks : bool 
        True to use minor ticks 
    dpi : int  
        Default=150. Dots per inch value to use for saved plot.
    """
    plt.rcParams.update({'font.family': 'serif','xtick.major.size':6,'ytick.major.size':6, \
                         'xtick.major.width':1.,'ytick.major.width':1., \
                         'xtick.direction':'in','ytick.direction':'in',
                         'xtick.top':True, 'ytick.right':True})
    fig1 = plt.figure(1)
    fig1.set_facecolor(facecolor)
    wcs = WCS(hdrin)
    ax1=fig1.add_subplot(111, projection=wcs)
    ax1.imshow(multicolorin,origin='lower',interpolation='nearest')
    ax1.set_title(axtitle,color=labelcolor,size=12)
    #ax1.tick_params(axis='both', colors=tickcolor); 
    ax1.coords.frame.set_color(tickcolor) #This way is particular to astropy wcs frames
    rapars = ax1.coords[0]
    decpars = ax1.coords[1]
    #rapars.set_ticks(spacing=10*u.arcmin, color='white', exclude_overlapping=True)
    #decpars.set_ticks(spacing=5*u.arcmin, color='white', exclude_overlapping=True)
    rapars.display_minor_ticks(minorticks); rapars.set_minor_frequency(5)
    decpars.display_minor_ticks(minorticks)
    rapars.set_major_formatter('hh:mm:ss')
    decpars.set_major_formatter('dd:mm:ss')
    rapars.set_separator(('$^\mathrm{H}$', "'", '"'))
    #decpars.ticklabels.set_rotation(45) #Rotate ticklabels
    #rapars.set_ticks(color=tickcolor); decpars.set_ticks(color=tickcolor)
    rapars.set_ticks(number=6,size=8,color=tickcolor); 
    decpars.set_ticks(number=6,size=8,color=tickcolor); 
    rapars.set_ticklabel(size=10,color=labelcolor); decpars.set_ticklabel(size=10,color=labelcolor); 
    ax1.set_xlabel(xaxislabel, color=labelcolor); 
    ax1.set_ylabel(yaxislabel, color=labelcolor)
    #ax1.coords.frame.set_linewidth(1.5)
    plt.savefig(savepath,bbox_inches='tight',dpi=dpi,facecolor=fig1.get_facecolor())
    plt.clf(); plt.close('all')

def comparemulticolorRGB_pureRGB(rgbin,multicolorin,hdrin,ax2title,suptitle,savepath,xaxislabel='RA', yaxislabel='DEC',tickcolor='w', labelcolor='k',facecolor='w', supy=.8, dpi=150, minorticks=True):
    """
    Convenience function to compare a pure RGB image with your multicolor RGB plot, and save to file.
    
    For tick/label colors: Any standard matplotlib color format, such as 'k', 'black', '#000000', or '0.0' 
    
    Parameters
    ----------
    rgbin : array 
        Input pure RGB image  -->  e.g.,   np.dstack( [redframe, gframe, bframe] ) 
    multicolorin : array
        Multicolor RGB image -- as would be output from combine_multicolor()
    hdrin : astropy.io.fits.header 
        Header to use for WCS information
    ax2title : str 
        String to use for 2nd axis plot title.  e.g. "My crazy 5-color image", or empty quotes "" for nothing.
    suptitle : str 
        String for super title (centered at top of plot, super title applies to whole figure)
    savepath : str 
        Path to save file to.  e.g., "./plots/mycrazyimage.pdf"
    xaxislabel: str
        Label to use for x-axes. Default='RA'
    yaxislabel: str
        Label to use for y-axes. Default='DEC'
    tickcolor : str  
        Color to use for ticks in plot, default = 'w'. 
    labelcolor : str 
        Color to use for ticklabels, default = 'k'
    facecolor : str 
        Color to use for figure facecolor (border around plot), default ='w'
    minorticks : bool 
        True to use minor ticks 
    dpi : int  
        Default=150. Dots per inch value to use for saved plot.
    supy : float 
        Position for suptitle (default = 0.8)
    """
    plt.rcParams.update({'font.family': 'serif','xtick.major.size':6,'ytick.major.size':6, \
                         'xtick.major.width':1.,'ytick.major.width':1., \
                         'xtick.direction':'in','ytick.direction':'in',
                         'xtick.top':True, 'ytick.right':True})
    fig1 = plt.figure(1,figsize=(10,8))
    fig1.set_facecolor(facecolor)
    wcs = WCS(hdrin)
    ax1=fig1.add_subplot(121, projection=wcs)
    #ax1.set_facecolor(facecolor)
    ax1.imshow(rgbin,origin='lower',interpolation='nearest')
    ax1.set_title('Simple RGB',color=labelcolor)
    ax2=fig1.add_subplot(122, projection=wcs)
    ax2.imshow(multicolorin,origin='lower',interpolation='nearest')
    for ax in [ax1,ax2]:
        rapars = ax.coords[0]
        decpars = ax.coords[1]
        #rapars.set_ticks(color=tickcolor); decpars.set_ticks(color=tickcolor)
        rapars.set_ticks(number=6,color=tickcolor); decpars.set_ticks(number=6,color=tickcolor); 
        rapars.set_ticklabel(size=8,color=labelcolor); 
        decpars.set_ticklabel(size=8,color=labelcolor); 
        #rapars.set_ticks(spacing=10*u.arcmin, color='white', exclude_overlapping=True)
        #decpars.set_ticks(spacing=5*u.arcmin, color='white', exclude_overlapping=True)
        rapars.display_minor_ticks(minorticks); #rapars.set_minor_frequency(10)
        decpars.display_minor_ticks(minorticks)
        rapars.set_major_formatter('hh:mm:ss')
        decpars.set_major_formatter('dd:mm:ss')
        rapars.set_separator(('$^\mathrm{H}$', "'", '"'))
        #decpars.ticklabels.set_rotation(45) #Rotate ticklabels
        #decpars.ticklabels.set_color(xkcdrust) #Ticklabel Color
        ax.set_xlabel(xaxislabel, color=labelcolor); 
        ax.set_ylabel(yaxislabel, color=labelcolor)
    plt.subplots_adjust(wspace=0.3)
    ax2.set_title(ax2title,color=labelcolor)
    for ax in [ax1, ax2]: 
        #ax.tick_params(axis='both', colors=tickcolor); 
        ax.coords.frame.set_color(tickcolor) #This way is particular to astropy wcs frames
        #ax.coords.frame.set_linewidth(1.5)
    plt.suptitle(suptitle,y=supy,color=labelcolor)
    plt.savefig(savepath,bbox_inches='tight',dpi=dpi,facecolor=fig1.get_facecolor())
    plt.clf(); plt.close('all')

def saveRGBfits(savepath,multicolorRGBdat,commonhdr,overwrite=True):
    """
    Convenience function to save out the multicolorRGB image to a fits file.
    
    Parameters
    ----------
    savepath : str 
        Path to save file to.  e.g., "./fits/mycrazy5colr.fits"
    multicolorRGBdat : array 
        Multicolor RGB image -- as would be output from combine_multicolor()
    commonhdr : astropy.io.fits.header 
        The common header with the WCS information for the output images
    overwrite : bool
        Passed to astropy.io.fits.writeto()
    """
    pyfits.writeto(savepath,np.swapaxes(np.swapaxes(multicolorRGBdat,0,2),2,1),commonhdr, overwrite=overwrite)

class Header_Editor(HasTraits):
    """
    FITS image header editor window
    """
    #info = Instance(UIInfo)
    contents=Str(multi_line=True, editor=CodeEditor() )
    
    #: The OK, Cancel and Load/Save Header buttons:
    ok = Button("OK")
    cancel = Button("Cancel")
    ##update_button = Button("Apply Update") #, action="update_header")
    #update = Action(name="Apply Update", action="update_header", tooltip='Apply changes to image header')
    loadheader = Action(name="Load from file", action="loadhdrfromfile", 
                 tooltip='Load header from text file')
    saveheader = Action(name="Save to file", action="savehdrtofile", 
                 tooltip='Save header to text file')
    
    traits_view=View(
        Item(name='contents',show_label=False),
        buttons=[loadheader, saveheader, OKButton, CancelButton],
        #buttons=[update, loadheader, saveheader, OKButton, CancelButton],
        title='FITS Header', resizable=True,
    )
    
    
    #def _contents_changed(self):
    #    #This will update the header after ANY changes (i.e., while typing)
    #    self.hdr=pyfits.Header.fromstring(self.contents) #No, this creates Header_Editor.hdr...
    
    #def _update_button_fired(self): 
    #    self.hdr=pyfits.Header.fromstring(self.contents) #No, this creates Header_Editor.hdr...
    
    #def update_header(self):
    #    self.hdr=pyfits.Header.fromstring(self.contents) #No, this creates Header_Editor.hdr...
    #    #--> Still cannot figure out how to make a pop-up window modify an
    #    #    object in the parent window...
    
    def loadhdrfromfile(self):
        hdrloadpath=open_file()
        with open(hdrloadpath, 'r') as headerinfile:
            self.contents=headerinfile.read()
    
    def savehdrtofile(self):
        #hdrstring=self.hdr.tostring(sep='\n')
        hdrsavepath=save_file()
        f=open(hdrsavepath,'w')
        f.write( self.contents )
        f.close()

class _MPLFigureEditor(Editor):
    """
    GUI Figure Editor 
    """
    scrollable  = True
    def init(self, parent):
        self.control = self._create_canvas(parent)
        self.set_tooltip()
    def update_editor(self): pass
    def _create_canvas(self, parent):
        """ Create the MPL canvas. """
        # matplotlib commands to create a canvas
        frame = QtGui.QWidget()
        mpl_canvas = FigureCanvas(self.value)
        mpl_canvas.setParent(frame)
        #mpl_toolbar = NavigationToolbar2QTAgg(mpl_canvas,frame)
        mpl_toolbar = NavigationToolbar2QT(mpl_canvas,frame)

        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(mpl_canvas)
        vbox.addWidget(mpl_toolbar)
        frame.setLayout(vbox)

        return frame

class MPLFigureEditor(BasicEditorFactory): klass = _MPLFigureEditor

class MPLInitHandler(Handler):
    ui_info = Instance(UIInfo)
    def init(self, info):
        """
        This method gets called after the controls have all been
        created but before they are displayed.
        """
        self.ui_info = info
        self.ui_info.object.setup_mpl_events()
        return True


##### The left side panels #####

class ControlPanel(HasTraits):
    """
    This is the control panel where the various parameters for the images are specified in the GUI.
    (The left-side panel)
    """
    
    gamma=Float(2.2)
    
    AutoRefresh=Bool(True)
    
    fitsfile    = File(filter=[u"*.fit*"],editor = FileEditor(dialog_style='open'))
    #fitsfile    = File(filter=[u"*.fits",u"*.*"])
    image_figure = Instance(Figure, ())
    image = Array()
    image_axes = Instance(Axes)
    image_axesimage = Instance(AxesImage)
    image_xsize = Int(256); image_ysize = Int(256)
    #hdrstring=Str(multiline=True)
    
    datamin= Float(0.0,auto_set=False,enter_set=True) #Say, in mJy
    datamax= Float(1.0,auto_set=False,enter_set=True) #auto_set=input set on each keystroke, enter_set=set after Enter
    percent_min=Range(value=0.0, low=0.0, high=100.); 
    percent_max=Range(value=100.0, low=0.0, high=100.) #Percentile of data values for rescaling
    minmaxbutton=Button('Min/Max')
    zscalebutton=Button('Zscale')
    updatebutton=Button('Update Preview Plot')
    
    image_scale=Str('linear')
    scale_dropdown=Enum(['linear','sqrt','squared','log','power','sinh','asinh'])
    smooth_image=Bool(False)
    smooth_sigma=Float(3.0, auto_set=False,enter_set=True)
    
    imagecolor=Str('#FFFFFF') 
    imagecolor_picker=ColorTrait((255,255,255))
    
    #plotbeam_button=Button('Add Beam (FWHM)')
    
    plotbutton_individual = Button(u"Plot Single Image")
    plotbutton_inverted_individual = Button(u"Plot Inverted Single Image")
    clearbutton_individual = Button(u"Clear Single Image")
    
    headerwindow=Instance(HasTraits)
    open_headerwindow=Button('Edit Header')
    applyheader=Button('Apply Header Changes')
    
    status_string_left=Str('')
    status_string_right=Str('')
    
    def __init__(self):
        self._init_params() #Set placeholder things like the WCS, tick color, map units...
        self.image = self._fresh_image() #Sets a blank image
        self.image_axes = self.image_figure.add_subplot(111,aspect=1)
        self.image_axesimage = self.image_axes.imshow(self.image, cmap='gist_gray', origin='lower', interpolation='nearest')
        self.image_axes.axis('off')
        #self.hdrstring='COMMENT  No image loaded'
        self.hdr = pyfits.Header.fromstring('COMMENT  No image loaded')
    

    view = View(
      HSplit(
        VGroup( 
            Item("fitsfile", label=u"Select 2D FITS file", style='simple', show_label=True ), #,height=100),
            
            HGroup(
              VGroup(
                     Item('plotbutton_individual', tooltip=u"Plot the single image",show_label=False),
                     Item('plotbutton_inverted_individual', \
                        tooltip=u"Plot the single inverted image", show_label=False),
                     Item('clearbutton_individual', tooltip=u"Clear the single image",show_label=False),
                     Item('_'), Item('50'), #Spacer of 50 pixels
                     
                     Item(name='open_headerwindow',show_label=False, \
                        tooltip='Open header for editing in new window.\n' + 
                                'To apply changes, click "Apply Header Changes".\n' + 
                                'Clicking "Edit Header" again before applying\n' +
                                '  will discard/cancel changes.'),
                     Item(name='applyheader', show_label=False, \
                        tooltip='Apply the modified header to the loaded .fits image'), 
                     Item('_'),  Item('50'), #Spacer of 50 pixels
                      
                     Item('imagecolor',label='Image Color',show_label=True, \
                      tooltip='Color of ticks: standard name float[0..1], or #hex', \
                      editor=TextEditor(auto_set=False, enter_set=True,)),
                     Item('imagecolor_picker',label='Pick Color',show_label=True,editor=ColorEditor()),
                     
                     ),
              Item('image_figure', editor=MPLFigureEditor(), show_label=False, width=300, height=300, resizable=True),
            ),
            HGroup(Item('datamin', tooltip=u"Minimum data val for scaling", show_label=True),
                   Item('datamax', tooltip=u"Maximum data val for scaling", show_label=True)),
            Item('percent_min', tooltip=u"Min. percentile for scaling", show_label=True),
            Item('percent_max', tooltip=u"Max. percentile for scaling", show_label=True),
            HGroup(Item('minmaxbutton', tooltip=u"Reset to data min/max", show_label=False),
                   Item('zscalebutton', tooltip=u"Compute scale min/max from zscale algorithm",  
                        show_label=False),
                   Item('50'), #Spacer of 50 pixels
                   Item('updatebutton', tooltip=u"Update with the new min/max values", show_label=False),
                   ), #Hgroup of Min/Max, Zscale, update
            Item('scale_dropdown',label='Scale',show_label=True),
            HGroup(Item('smooth_image', style='custom', label='Smooth image', tooltip='Click to perform Gaussian smoothing'),
                   Item('smooth_sigma', label='Smoothing \u03C3 (pixels)', 
                        tooltip=u"Number of pixels for sigma in Gaussian smoothing" ),
                  ),
                
        ), #End of Left column
        
      ), #End of HSplit
      resizable=True, #height=0.75, width=0.75, #title=u"Multi-Color Image Combiner", 
      handler=MPLInitHandler,
      statusbar = [StatusItem(name = 'status_string_left', width = 0.5), \
                   StatusItem(name = 'status_string_right', width = 0.5)]
    ) #End of View
    
    def _init_params(self):
        self.in_use=False
        plt.rcParams.update({'font.family': 'serif','xtick.major.size':6,'ytick.major.size':6, \
                             'xtick.major.width':1.,'ytick.major.width':1., \
                             'xtick.direction':'in','ytick.direction':'in'})
        try: plt.rcParams.update({'xtick.top':True, 'ytick.right':True}) #apparently not in mpl v<2.0...
        except: pass #Make a workaround for mpl<2.0 later...
        self.datamin_initial=0.; self.datamax_initial=1.; 
        self.datamin=0.; self.datamax=1. #This will be the displayed value of the scaling min/max
    
    def _fresh_image(self): 
        blankdata=np.zeros([100,100]); blankdata[-1,-1]=1
        return blankdata
    
    def _fitsfile_changed(self):
        self.data,self.hdr=pyfits.getdata(self.fitsfile,header=True)
        force_hdr_floats(self.hdr) #Ensure that WCS cards such as CDELT are floats instead of strings
        
        naxis = int(self.hdr['NAXIS'])
        if naxis > 2:
            #print('Dropping Extra axes')
            self.hdr=force_hdr_to_2D(self.hdr)
            try: self.data=self.data[0,0,:,:]
            except: self.data=self.data[0,:,:]
            self.status_string_right = 'Dropped extra axes'
        
        self.datamax_initial = np.nanmax(self.data).item()
        self.datamin_initial = np.nanmin(self.data).item()
        self.datamax = np.nanmax(self.data).item()
        self.datamin = np.nanmin(self.data).item()
        
        self.in_use=True
        
    def _headerwindow_default(self):
        return Header_Editor()

    def _open_headerwindow_fired(self):
        self.headerwindow.contents=self.hdr.tostring(sep='\n')
        self.headerwindow.edit_traits()
    #def _open_headerwindow_changed(self):
    #    self.hdr=pyfits.Header.fromstring(self.headerwindow.contents,sep='\n')

    def _applyheader_fired(self):
        self.hdr=pyfits.Header.fromstring(self.headerwindow.contents,sep='\n')
        
    #@on_trait_change('imagecolor')
    #def update_imagecolor(self): 
    def _imagecolor_changed(self): 
        try: 
            #Catch case when you've predefined a color variable in hex string format, e.g., mynewred='#C11B17'
            #--> Need to do this first, otherwise traits throws a fit up the stack even despite the try/except check
            globals()[self.imagecolor] #This check should catch undefined inputs
            self.imagecolor_picker=hex_to_rgb(globals()[self.imagecolor])
            self.status_string_right = 'Image color changed to '+self.imagecolor
        except:
            try: 
                self.imagecolor=to_hex(self.imagecolor)
                self.status_string_right = 'Image color changed to '+to_hex(self.imagecolor)
            except: self.status_string_right = "Color name %s not recognized.  Must be standard mpl.colors string, float[0..1] or #hex string"%(self.imagecolor) 
        try: self.imagecolor_picker=hex_to_rgb(to_hex(self.imagecolor)) #update the picker color...
        except: pass
        ### self.image_greyRGB and self.image_colorRGB may not yet be instantiated if the color is changed before clicking 'plot'
        try: self.image_colorRGB=colorize_image(self.image_greyRGB, self.imagecolor, colorintype='hex', gammacorr_color=self.gamma)
        except: pass
        try: self.image_axesimage.set_data(self.image_colorRGB**(1./self.gamma))
        except: pass
        self.in_use=True
        self.image_figure.canvas.draw()
    
    #@on_trait_change('imagecolor_picker')
    #def update_imagecolorpicker(self):
    def _imagecolor_picker_changed(self):
        #print(self.tickcolor_picker.name())
        self.imagecolor=self.imagecolor_picker.name()
    
    
    #@on_trait_change('percent_min')
    #def update_scalepercmin(self):
    def _percent_min_changed(self):
        self.datamin = np.nanpercentile(self.data,self.percent_min)
        if self.AutoRefresh==True:
            #Disable automatic plot updating for large file sizes, to avoid hangs
            self.data_scaled=(scaling_fns[self.image_scale]() + ManualInterval(vmin=self.datamin,vmax=self.datamax))(self.data)
            self.image_greyRGB=ski_color.gray2rgb(adjust_gamma(self.data_scaled,self.gamma))
            self.image_colorRGB=colorize_image(self.image_greyRGB, self.imagecolor, colorintype='hex', gammacorr_color=self.gamma)
            self.image_axesimage.set_data(self.image_colorRGB**(1./self.gamma))
            self.image_figure.canvas.draw()
            self.status_string_right = "Updated scale using percentiles"
    #@on_trait_change('percent_max')
    #def update_scalepercmax(self):
    def _percent_max_changed(self):
        self.datamax = np.nanpercentile(self.data,self.percent_max)
        if self.AutoRefresh==True:
            #Disable automatic plot updating for large file sizes, to avoid hangs
            self.data_scaled=(scaling_fns[self.image_scale]() + ManualInterval(vmin=self.datamin,vmax=self.datamax))(self.data)
            self.image_greyRGB=ski_color.gray2rgb(adjust_gamma(self.data_scaled,self.gamma))
            self.image_colorRGB=colorize_image(self.image_greyRGB, self.imagecolor, colorintype='hex', gammacorr_color=self.gamma)
            self.image_axesimage.set_data(self.image_colorRGB**(1./self.gamma))
            self.image_figure.canvas.draw()
            self.status_string_right = "Updated scale using percentiles"
    
    def _datamin_changed(self): 
        self.percent_min=np.round(nanpercofscore(self.data.ravel(),self.datamin,kind='strict'),2)
    def _datamax_changed(self): 
        self.percent_max=np.round(nanpercofscore(self.data.ravel(),self.datamax,kind='strict'),2)
    #@on_trait_change('datamin')
    #def update_datamin(self): 
    #    self.percent_min=np.round(nanpercofscore(self.data.ravel(),self.datamin,kind='strict'),2)
    #@on_trait_change('datamax')
    #def update_datamax(self): 
    #    self.percent_max=np.round(nanpercofscore(self.data.ravel(),self.datamax,kind='strict'),2)
    
    #@on_trait_change('scale_dropdown')
    #def update_image_scale(self):
    def _scale_dropdown_changed(self):
        self.image_scale=self.scale_dropdown
        #self.norm=ImageNormalize(self.sregion,stretch=scaling_fns[self.image_scale]() )
        self.data_scaled=(scaling_fns[self.image_scale]() + ManualInterval(vmin=self.datamin,vmax=self.datamax))(self.data)
        #*** Instead, should I just integrate my imscale class here instead of astropy? ...
        self.image_greyRGB=ski_color.gray2rgb(adjust_gamma(self.data_scaled,self.gamma))
        self.image_colorRGB=colorize_image(self.image_greyRGB, self.imagecolor, colorintype='hex', gammacorr_color=self.gamma)
        
        self.image_axesimage.set_data(self.image_colorRGB**(1./self.gamma))
        
        self.in_use=True
        
        self.image_figure.canvas.draw()
        self.status_string_right = 'Image scale function changed to '+self.image_scale
    
    def _minmaxbutton_fired(self): 
        self.datamin=self.datamin_initial; self.datamax=self.datamax_initial
        #self.image_axesimage.norm.vmin=self.datamin
        #self.image_axesimage.norm.vmax=self.datamax
        self.percent_min=np.round(nanpercofscore(self.data.ravel(),self.datamin,kind='strict'),2)
        self.percent_max=np.round(nanpercofscore(self.data.ravel(),self.datamax,kind='strict'),2)
        #self.image_figure.canvas.draw()
        self.status_string_right = "Scale reset to min/max"
    
    def _zscalebutton_fired(self):
        tmpZscale=ZScaleInterval().get_limits(self.data)
        self.datamin=float(tmpZscale[0])
        self.datamax=float(tmpZscale[1])
        self.percent_min=np.round(nanpercofscore(self.data.ravel(),self.datamin,kind='strict'),2)
        self.percent_max=np.round(nanpercofscore(self.data.ravel(),self.datamax,kind='strict'),2)
        #self.image_figure.canvas.draw()
        self.status_string_right = "Min/max determined by zscale"
    
    def _smooth_image_changed(self):
        try: self.data
        except: self.status_string_right = "No fits file loaded yet!"; return
        if self.smooth_image == True:
            self.image_colorRGB=filters.gaussian(self.image_colorRGB, sigma=self.smooth_sigma, )#multichannel=True)
        else: 
            self.image_colorRGB=colorize_image(self.image_greyRGB, self.imagecolor, colorintype='hex', gammacorr_color=self.gamma)
            self.image_colorRGB=colorize_image(self.image_greyRGB, self.imagecolor, colorintype='hex', gammacorr_color=self.gamma)
        self.image_axesimage.set_data(self.image_colorRGB**(1./self.gamma))
        self.image_figure.canvas.draw()
        self.status_string_right = 'Image smoothing turned '+['On' if self.smooth_image == True else 'Off'][0]
    def _smooth_sigma_changed(self):
        try: self.data
        except: self.status_string_right = "No fits file loaded yet!"; return
        if self.smooth_image == True:
            self.image_colorRGB=colorize_image(self.image_greyRGB, self.imagecolor, colorintype='hex', gammacorr_color=self.gamma)
            self.image_colorRGB=filters.gaussian(self.image_colorRGB, sigma=self.smooth_sigma, )#multichannel=True)
            self.image_axesimage.set_data(self.image_colorRGB**(1./self.gamma))
            self.image_figure.canvas.draw()
            self.status_string_right = 'Image smoothed with \u03C3 = %.2f pixels'%(self.smooth_sigma)
    
    def _plotbutton_individual_fired(self):
        try: self.data
        except: self.status_string_right = "No fits file loaded yet!"; return
        #self.image=self.data
        ###Using this command is preferable, as long as the projection doesn't need to be updated...
        #  The home zoom button will work, but no WCS labels because projection wasn't set during init.  
        #Scale the data to [0,1] range
        self.data_scaled=(scaling_fns[self.image_scale]() + ManualInterval(vmin=self.datamin,vmax=self.datamax))(self.data)
        #Convert scale[0,1] image to greyscale RGB image
        self.image_greyRGB=ski_color.gray2rgb(adjust_gamma(self.data_scaled,self.gamma))
        self.image_colorRGB=colorize_image(self.image_greyRGB, self.imagecolor, colorintype='hex', gammacorr_color=self.gamma)
        self.image_axesimage.set_data(self.image_colorRGB**(1./self.gamma))
        ###Using this set instead properly updates the axes labels to WCS, but the home zoom button won't work
        #self.image_figure.clf()
        #self.image_axes = self.image_figure.add_subplot(111,aspect=1)#,projection=self.wcs)
        #self.image_axesimage = self.image_axes.imshow(self.image, cmap=self.image_cmap,origin='lower',interpolation='nearest', norm=self.norm)
        
        self.percent_min=np.round(nanpercofscore(self.data.ravel(),self.datamin,kind='strict'),2)
        self.percent_max=np.round(nanpercofscore(self.data.ravel(),self.datamax,kind='strict'),2)
        
        self.in_use=True
        
        #self.update_radecpars()
        self.image_figure.canvas.draw()
        self.status_string_right = "Plot updated"

    def _plotbutton_inverted_individual_fired(self):
        try: self.data
        except: self.status_string_right = "No fits file loaded yet!"; return
        self.data_scaled=(scaling_fns[self.image_scale]() + ManualInterval(vmin=self.datamin,vmax=self.datamax))(self.data)
        self.image_greyRGB=ski_color.gray2rgb(adjust_gamma(self.data_scaled,self.gamma))
        self.image_colorRGB=colorize_image(self.image_greyRGB, hexinv(self.imagecolor), colorintype='hex', gammacorr_color=self.gamma)
        self.image_axesimage.set_data(1.-self.image_colorRGB**(1./self.gamma)) 
        #self.image_axesimage.set_data(combine_multicolor([self.image_colorRGB,],gamma=self.gamma,inverse=True))  
        self.percent_min=np.round(nanpercofscore(self.data.ravel(),self.datamin,kind='strict'),2)
        self.percent_max=np.round(nanpercofscore(self.data.ravel(),self.datamax,kind='strict'),2)
        self.in_use=True
        self.image_figure.canvas.draw()
        self.status_string_right = "Plot updated"
    
    def _updatebutton_fired(self):
        try: self.data
        except: self.status_string_right = "No fits file loaded yet!"; return
        self.data_scaled=(scaling_fns[self.image_scale]() + ManualInterval(vmin=self.datamin,vmax=self.datamax))(self.data)
        self.image_greyRGB=ski_color.gray2rgb(adjust_gamma(self.data_scaled,self.gamma))
        self.image_colorRGB=colorize_image(self.image_greyRGB, self.imagecolor, colorintype='hex', gammacorr_color=self.gamma)
        self.image_axesimage.set_data(self.image_colorRGB**(1./self.gamma))
        self.image_figure.canvas.draw()
        self.status_string_right = "Plot updated"
        
    def _clearbutton_individual_fired(self): 
        try: del self.data, self.data_scaled, self.image_greyRGB; self.image_colorRGB #In case clear already pressed once
        except: pass
        self.in_use=False
        self.image_figure.clf()
        self.image = self._fresh_image()
        self.hdr = pyfits.Header.fromstring('COMMENT  No image loaded')
        self.image_axes = self.image_figure.add_subplot(111,aspect=1)
        self.image_axesimage = self.image_axes.imshow(self.image, cmap='gist_gray',origin='lower',interpolation='nearest')
        self.image_axes.axis('off')
        self.image_figure.canvas.draw()
        self.status_string_right = "Plot cleared"

    def setup_mpl_events(self):
        self.image_axeswidget = AxesWidget(self.image_axes)
        self.image_axeswidget.connect_event('motion_notify_event', self.image_on_motion)
        self.image_axeswidget.connect_event('figure_leave_event', self.on_cursor_leave)
        self.image_axeswidget.connect_event('figure_enter_event', self.on_cursor_enter)
        self.image_axeswidget.connect_event('button_press_event', self.image_on_click)

    def image_on_motion(self, event):
        if event.xdata is None or event.ydata is None: return
        x = int(np.round(event.xdata))
        y = int(np.round(event.ydata))
        if ((x >= 0) and (x < self.image.shape[1]) and
            (y >= 0) and (y < self.image.shape[0])):
            imval = self.image[y, x]
            self.status_string_left = "x,y={},{}  {:.5g}".format(x, y, imval)
        else:  self.status_string_left = ""
    
    def image_on_click(self, event):
        if event.xdata is None or event.ydata is None or event.button is not 1: return #Covers when click outside of main plot
        #print(event)
        x = int(np.round(event.xdata)) #xdata is the actual pixel position.  xy is in 'display space', i.e. pixels in the canvas
        y = int(np.round(event.ydata))
        #xwcs,ywcs=self.wcs.wcs_pix2world([[x,y]],0)[0]; #print(xwcs,ywcs)
        if ((x >= 0) and (x < self.image.shape[1]) and
            (y >= 0) and (y < self.image.shape[0])):
            imval = self.image[y, x]
            #self.status_string_right = "x,y=[{},{}], RA,DEC=[{}, {}], value = {:.5g}".format(x, y,xwcs,ywcs, imval)
            self.status_string_right = "x,y[{},{}] = {:.3f},{:.3f}  {:.5g}".format(x, y,event.xdata,event.ydata, imval)
        else: self.status_string_right = ""
        ## left-click: event.button = 1, middle-click: event.button=2, right-click: event.button=3.  
        ## For double-click, event.dblclick = False for first click, True on second
        #print(event.button, event.dblclick)

    def on_cursor_leave(self, event):
        QApplication.restoreOverrideCursor()
        self.status_string_left = ''

    def on_cursor_enter(self, event): QApplication.setOverrideCursor(Qt.CrossCursor)


##### The main viewer #####

class multicolorfits_viewer(HasTraits):
    """The main window. (Right-side panel.)
    Has instructions for creating and destroying the app.
    """
    
    panel1 = Instance(ControlPanel)
    panel2 = Instance(ControlPanel)
    panel3 = Instance(ControlPanel)
    panel4 = Instance(ControlPanel)
    
    figure_combined = Instance(Figure,())
    image = Array()
    image_axes = Instance(Axes)
    image_axesimage = Instance(AxesImage)
    image_xsize = Int(256); image_ysize = Int(256)
    
    gamma = Float(2.2)
    
    AutoRefresh=Bool(True)
    #count_changes = Int(0)
    
    tickcolor=Str('0.9')#,auto_set=False,enter_set=True) #Apparently need to set to TextEditor explicitly below...
    tickcolor_picker=ColorTrait((230,230,230))
    showminorticks=Bool(True)
    sexdec=Enum('Sexagesimal','Decimal')
    sexdec_x_formatter=Str('hh:mm:ss.ss')
    sexdec_y_formatter=Str('dd:mm:ss.ss')
    sex_x_formatter=Str('hh:mm:ss.ss')
    sex_y_formatter=Str('dd:mm:ss.ss')
    dec_x_formatter=Str('d.dddddd')
    dec_y_formatter=Str('d.dddddd')
    
    plotbutton_combined = Button(u"Plot Combined Image")
    plotbutton_inverted_combined = Button(u"Plot Inverted Combined Image")
    clearbutton_combined = Button(u"Clear Combined Image")
    save_the_image = Button(u"Save Image")
    save_the_fits = Button(u"Save RGB Fits")
    print_params = Button(u"Print Params")
    
    status_string_left = Str('')
    status_string_right = Str('')
    
    def _panel1_default(self):
        return ControlPanel()#figure=self.figure)
    def _panel2_default(self):
        return ControlPanel()#figure=self.figure)
    def _panel3_default(self):
        return ControlPanel()#figure=self.figure)
    def _panel4_default(self):
        return ControlPanel()#figure=self.figure)
    
    def __init__(self):
        super(multicolorfits_viewer, self).__init__()
        
        self._init_params() #Set placeholder things like the WCS, tick color, map units...
        self.image = self._fresh_image() #Sets a blank image
        self.image_axes = self.figure_combined.add_subplot(111,aspect=1)
        self.image_axesimage = self.image_axes.imshow(self.image, cmap='gist_gray', origin='lower', interpolation='nearest')
        self.image_axes.text(50,50,'Currently, all images must share a common pixel grid before \nloading, as this GUI does not yet reproject on the fly. \nThe provided mcf.reproject2D() function can be used to achieve this.',color='w',ha='center')
        self.image_axes.set_xlabel(self.xlabel); self.image_axes.set_ylabel(self.ylabel)
        self.image_axes.tick_params(axis='both',color=self.tickcolor) #colors=... also sets label color
        try: self.image_axes.coords.frame.set_color(self.tickcolor) #Updates the frame color.  .coords won't exist until WCS set
        except: 
            [self.image_axes.spines[s].set_color(self.tickcolor) for s in ['top','bottom','left','right']]
    
    view = View(HGroup(
                    Item('AutoRefresh', style='custom', label='Auto-Refresh Plots (un-check box to increase speed for large files)'),
                    Item('50'), #Spacer of 50 pixels
                    Item("gamma",label=u"Gamma",show_label=True), 
                ), #Top Panel HGroup
                Item('_'),
                
                HSplit(
                  Group(
                    Group(Item('panel1', style="custom",show_label=False),label='Image 1'),
                    Group(Item('panel2', style="custom",show_label=False),label='Image 2'),
                    Group(Item('panel3', style="custom",show_label=False),label='Image 3'),
                    Group(Item('panel4', style="custom",show_label=False),label='Image 4'),
                  orientation='horizontal',layout='tabbed',springy=True),
                
                VGroup(
                    HGroup(
                      VGroup(
                        Item('tickcolor',label='Tick Color',show_label=True, \
                            tooltip='Color of ticks: standard name float[0..1], or #hex', \
                            editor=TextEditor(auto_set=False, enter_set=True,)),
                        Item('tickcolor_picker',label='Pick',show_label=True,editor=ColorEditor()), 
                        Item('showminorticks', style='custom', label='Display minor ticks'),
                        ),
                      VGroup(
                        Item('sexdec',label='Coordinate Style',tooltip=u'Display coordinates in sexagesimal or decimal', \
                            show_label=True),
                        #auto_set=input set on each keystroke, enter_set=set after Enter
                        Item('sexdec_x_formatter',label='X-axis tick format',tooltip=u'Format for the x tick labels', \
                            show_label=True, editor=TextEditor(auto_set=False, enter_set=True,)),
                        Item('sexdec_y_formatter',label='Y-axis tick format',tooltip=u'Format for the y tick labels', \
                            show_label=True, editor=TextEditor(auto_set=False, enter_set=True,)),
                        ),
                      ),
                    Item('figure_combined', editor=MPLFigureEditor(),show_label=False, width=900, height=800, resizable=True),
                  HGroup(
                    
                    Item('plotbutton_combined', tooltip=u"Plot the image",show_label=False),
                    Item('plotbutton_inverted_combined', tooltip=u"Plot the inverted image",show_label=False),
                    Item('clearbutton_combined',tooltip=u'Clear the combined figure',show_label=False),
                    Item("save_the_image", tooltip=u"Save current image. Mileage may vary...",show_label=False),
                    Item("save_the_fits", tooltip=u"Save RGB frames as single fits file with header.",show_label=False),
                    Item("print_params", tooltip=u"Print out current settings for use in manual image scripting.",show_label=False),
                  ), #HGroup
                ), #VGroup
                show_labels=False,),
           resizable=True,
           height=0.75, width=0.75,
           statusbar = [StatusItem(name = 'status_string_left', width = 0.5),
                        StatusItem(name = 'status_string_right', width = 0.5)],
           title=u"Fits Multi-Color Combiner",handler=MPLInitHandler ) #View
    
    def _init_params(self):
        plt.rcParams.update({'font.family': 'serif','xtick.major.size':6,'ytick.major.size':6, \
                             'xtick.major.width':1.,'ytick.major.width':1., \
                             'xtick.direction':'in','ytick.direction':'in'})
        try: plt.rcParams.update({'xtick.top':True, 'ytick.right':True}) #apparently not in mpl v<2.0...
        except: pass #Make a workaround for mpl<2.0 later...
        self.datamin_initial=0.; self.datamax_initial=1.; 
        self.datamin=0.; self.datamax=1. #This will be the displayed value of the scaling min/max
        #self.mapunits='Pixel Value'
        #self.tickcolor='0.5'#'white', 'black', '0.5'
        self.wcs=WCS(); 
        self.xlabel='x'; self.ylabel='y'; 
    
    def _fresh_image(self): 
        #self.norm=ImageNormalize(self.image,stretch=scaling_fns['linear']() )
        blankdata=np.zeros([100,100]); blankdata[-1,-1]=1
        return blankdata
    
    def update_radecpars(self):
        try: 
            if '-' in self.hdr['CTYPE1']: self.xlabel=self.hdr['CTYPE1'].split('-')[0] 
            else: self.xlabel=['CTYPE1'][:4]
        except: self.xlabel='x'
        try: 
            if '-' in self.hdr['CTYPE2']: self.ylabel=self.hdr['CTYPE2'].split('-')[0]
            else: self.ylabel=['CTYPE2'][:4]
        except: self.ylabel='y'
        self.image_axes.set_xlabel(self.xlabel); 
        self.image_axes.set_ylabel(self.ylabel)
        self.rapars = self.image_axes.coords[0]
        self.decpars = self.image_axes.coords[1]
        self.rapars.set_ticks(color=self.tickcolor); self.decpars.set_ticks(color=self.tickcolor)
        self.rapars.set_ticks(number=6); self.decpars.set_ticks(number=6); 
        #self.rapars.set_ticklabel(size=8); self.decpars.set_ticklabel(size=8); #size here means the tick length
        ##self.rapars.set_ticks(spacing=10*u.arcmin, color='white', exclude_overlapping=True)
        ##self.decpars.set_ticks(spacing=5*u.arcmin, color='white', exclude_overlapping=True)
        self.rapars.display_minor_ticks(self.showminorticks)#(True); #self.rapars.set_minor_frequency(10)
        self.decpars.display_minor_ticks(self.showminorticks)#(True)
        if self.sexdec=='Sexagesimal':
            self.rapars.set_major_formatter(self.sex_x_formatter)#('hh:mm:ss.ss')
            self.decpars.set_major_formatter(self.sex_y_formatter)#('dd:mm:ss.ss')
            #self.rapars.set_separator(('$^\mathrm{H}$', "'", '"'))
            self.rapars.set_separator(('$^\mathrm{H}$ ', "' ", '" '))
            #format_xcoord=lambda x,y: '{}i$^\mathrm{H}${}{}{}'.format(x[0],x[1],"'",x[2],'"')
            #self.image_axes.format_coord=format_xcoord
        else: 
            self.rapars.set_major_formatter(self.dec_x_formatter)#('d.dddddd')
            self.decpars.set_major_formatter(self.dec_y_formatter)#('d.dddddd')
            #self.rapars.set_separator(('$\degree$ ','')) #<-- Nope, currently this only affects sexagesimal
        ##self.decpars.ticklabels.set_rotation(45) #Rotate ticklabels
        ##self.decpars.ticklabels.set_color(xkcdrust) #Ticklabel Color
    
    @on_trait_change('gamma')
    def update_gamma(self):
        self.panel1.gamma=self.gamma; self.panel2.gamma=self.gamma; 
        self.panel3.gamma=self.gamma; self.panel4.gamma=self.gamma; 
        self.status_string_right = 'Gamma changed to '+str(self.gamma)
    
    @on_trait_change('AutoRefresh')
    #def _my_boolean_trait_changed(self): self.count_changes += 1 # Check for changes in AutoRefresh
    def update_AutoRefresh(self):
        self.panel1.AutoRefresh=self.AutoRefresh; self.panel2.AutoRefresh=self.AutoRefresh; 
        self.panel3.AutoRefresh=self.AutoRefresh; self.panel4.AutoRefresh=self.AutoRefresh; 
        self.status_string_right = 'AutoRefresh plots set to '+str(self.AutoRefresh)
    
    @on_trait_change('tickcolor')
    def update_tickcolor(self): 
        try: 
            #Catch case when you've predefined a color variable in hex string format, e.g., mynewred='#C11B17'
            #--> Need to do this first, otherwise traits throws a fit up the stack even despite the try/except check
            globals()[self.tickcolor] #This check should catch undefined inputs
            self.image_axes.tick_params(axis='both',color=globals()[self.tickcolor])
            self.image_axes.coords.frame.set_color(self.tickcolor)
            self.tickcolor_picker=hex_to_rgb(globals()[self.tickcolor])
            self.status_string_right = 'Tick color changed to '+to_hex(self.tickcolor)
        except:
            try: 
                self.tickcolor=to_hex(self.tickcolor)
                try: self.update_radecpars()
                except: 
                    self.image_axes.tick_params(axis='both',color=to_hex(self.tickcolor)); 
                    self.image_axes.coords.frame.set_color(to_hex(self.tickcolor));
                self.status_string_right = 'Tick color changed to '+to_hex(self.tickcolor)
            except: self.status_string_right = "Color name %s not recognized.  Must be standard mpl.colors string, float[0..1] or #hex string"%(self.tickcolor) 
        try: self.tickcolor_picker=hex_to_rgb(to_hex(self.tickcolor)) #update the picker color...
        except: pass
        self.figure_combined.canvas.draw()
    
    @on_trait_change('tickcolor_picker')
    def update_tickcolorpicker(self):
        #print(self.tickcolor_picker.name())
        self.tickcolor=self.tickcolor_picker.name()
    
    @on_trait_change('showminorticks')
    def update_showminorticks(self):
        self.update_radecpars()
        self.figure_combined.canvas.draw()
        self.status_string_right = 'Minor ticks turned '+['On' if self.showminorticks is True else 'Off'][0]
    
    @on_trait_change('sexdec')
    def update_sexdec(self): 
        self.update_radecpars()
        if self.sexdec=='Sexagesimal': 
            self.sexdec_x_formatter=self.sex_x_formatter
            self.sexdec_y_formatter=self.sex_y_formatter
        else: 
            self.sexdec_x_formatter=self.dec_x_formatter
            self.sexdec_y_formatter=self.dec_y_formatter
        self.figure_combined.canvas.draw()
        self.status_string_right = 'Coordinate style changed to '+self.sexdec
    
    @on_trait_change('sexdec_x_formatter')
    def update_sexdec_x_formatter(self): 
        if self.sexdec=='Sexagesimal': self.sex_x_formatter=self.sexdec_x_formatter
        else: self.dec_x_formatter=self.sexdec_x_formatter
        self.figure_combined.canvas.draw()
        self.status_string_right = 'X-axis tick format changed to '+self.sex_x_formatter
    @on_trait_change('sexdec_y_formatter')
    def update_sexdec_y_formatter(self): 
        if self.sexdec=='Sexagesimal': self.sex_y_formatter=self.sexdec_y_formatter
        else: self.dec_y_formatter=self.sexdec_y_formatter
        self.figure_combined.canvas.draw()
        self.status_string_right = 'Y-axis tick format changed to '+self.sex_y_formatter
    
    def _plotbutton_combined_fired(self):
        try: self.panel1.data
        except: self.status_string_right = "No fits file loaded yet!"; return
        #self.image=self.panel1.data
        self.wcs=WCS(self.panel1.hdr)
        self.hdr=self.panel1.hdr
        
        self.combined_RGB=combine_multicolor( [pan.image_colorRGB for pan in [self.panel1,self.panel2,self.panel3,self.panel4] if pan.in_use==True], gamma=self.gamma)
        
        ###Using this command is preferable, as long as the projection doesn't need to be updated...
        #  The home zoom button will work, but no WCS labels because projection wasn't set during init.  
        #self.image_axesimage.set_data(self.data)
        ###Using this set instead properly updates the axes labels to WCS, but the home zoom button won't work
        self.figure_combined.clf()
        self.image_axes = self.figure_combined.add_subplot(111,aspect=1,projection=self.wcs)
        self.image_axesimage = self.image_axes.imshow(self.combined_RGB, origin='lower', interpolation='nearest')
        
        self.update_radecpars()
        self.figure_combined.canvas.draw()
        self.status_string_right = "Plot updated"
    
    def _plotbutton_inverted_combined_fired(self):
        try: self.panel1.data
        except: self.status_string_right = "No fits file loaded yet!"; return
        self.wcs=WCS(self.panel1.hdr)
        self.hdr=self.panel1.hdr
        self.combined_RGB=combine_multicolor( [pan.image_colorRGB for pan in [self.panel1,self.panel2,self.panel3,self.panel4] if pan.in_use==True], inverse=True, gamma=self.gamma, )
        self.figure_combined.clf()
        self.image_axes = self.figure_combined.add_subplot(111,aspect=1,projection=self.wcs)
        self.image_axesimage = self.image_axes.imshow(self.combined_RGB, origin='lower', interpolation='nearest')
        self.update_radecpars()
        self.figure_combined.canvas.draw()
        self.status_string_right = "Plot updated"
        
    def _clearbutton_combined_fired(self): 
        try: del self.combined_RGB #If clear already pressed once, data will already have been deleted...
        except: pass
        self.in_use=False
        self.figure_combined.clf()
        self.image = self._fresh_image()
        self.image_axes = self.figure_combined.add_subplot(111,aspect=1)
        self.image_axesimage = self.image_axes.imshow(self.image, cmap='gist_gray', origin='lower', interpolation='nearest')
        self.xlabel='x'; self.ylabel='y'
        self.image_axes.set_xlabel(self.xlabel); self.image_axes.set_ylabel(self.ylabel)
        self.image_axes.tick_params(axis='both',color=self.tickcolor)
        try: self.image_axes.coords.frame.set_color(self.tickcolor)
        except: self.tickcolor_picker=hex_to_rgb(to_hex(self.tickcolor))
        self.figure_combined.canvas.draw()
        self.status_string_right = "Plot cleared"
        
    def setup_mpl_events(self):
        self.image_axeswidget = AxesWidget(self.image_axes)
        self.image_axeswidget.connect_event('motion_notify_event', self.image_on_motion)
        self.image_axeswidget.connect_event('figure_leave_event', self.on_cursor_leave)
        self.image_axeswidget.connect_event('figure_enter_event', self.on_cursor_enter)
        self.image_axeswidget.connect_event('button_press_event', self.image_on_click)

    def image_on_motion(self, event):
        if event.xdata is None or event.ydata is None: return
        x = int(np.round(event.xdata))
        y = int(np.round(event.ydata))
        if ((x >= 0) and (x < self.image.shape[1]) and
            (y >= 0) and (y < self.image.shape[0])):
            imval = self.image[y, x]
            self.status_string_left = "x,y={},{}  {:.5g}".format(x, y, imval)
        else:  self.status_string_left = ""
    
    def image_on_click(self, event):
        if event.xdata is None or event.ydata is None or event.button is not 1: return #Covers when click outside of main plot
        #print(event)
        x = int(np.round(event.xdata)) #xdata is the actual pixel position.  xy is in 'display space', i.e. pixels in the canvas
        y = int(np.round(event.ydata))
        xwcs,ywcs=self.wcs.wcs_pix2world([[x,y]],0)[0]; #print(xwcs,ywcs)
        if ((x >= 0) and (x < self.image.shape[1]) and
            (y >= 0) and (y < self.image.shape[0])):
            imval = self.image[y, x]
            self.status_string_right = "x,y=[{},{}], RA,DEC=[{}, {}], value = {:.5g}".format(x, y,xwcs,ywcs, imval)
            #self.status_string_right = "x,y[{},{}] = {:.3f},{:.3f}  {:.5g}".format(x, y,event.xdata,event.ydata, imval)
        else: self.status_string_right = ""
        ## left-click: event.button = 1, middle-click: event.button=2, right-click: event.button=3.  
        ## For double-click, event.dblclick = False for first click, True on second
        #print(event.button, event.dblclick)

    def on_cursor_leave(self, event):
        QApplication.restoreOverrideCursor()
        self.status_string_left = ''

    def on_cursor_enter(self, event): QApplication.setOverrideCursor(Qt.CrossCursor)
    
    def _save_the_image_fired(self):
        dlg = FileDialog(action='save as')
        if dlg.open() == OK: 
            self.figure_combined.canvas.draw()
            self.figure_combined.savefig(dlg.path,size=(800,800),dpi=300,bbox_inches='tight')
    
    def _save_the_fits_fired(self): 
        #Generate a generic header with correct WCS and comments about the colors that made it
        #... come back and finish this later...
        dlg = FileDialog(action='save as')
        if dlg.open() == OK: pyfits.writeto(dlg.path, np.swapaxes(np.swapaxes(self.combined_RGB,0,2),2,1), self.hdr)
    
    def _print_params_fired(self): 
        print('\n\nRGB Image plot params:')
        pan_i=0
        for pan in [self.panel1,self.panel2,self.panel3,self.panel4]: 
            pan_i+=1
            if pan.in_use==True: 
                print('image%i: '%(pan_i))
                print('    vmin = %.3e , vmax = %.3e, scale = %s'%(pan.datamin, pan.datamax, pan.image_scale))
                print("    image color = '%s'"%(pan.imagecolor))
        print("gamma = %.1f , tick color = '%s'\n"%(self.gamma,self.tickcolor))
        #print('vmin=%.3e, vmax=%.3e'%(self.panel1.datamin,self.panel1.datamax))
        
def mcf_gui():
    """
    Call this function to run the GUI. (Keeping as legacy)
    """
    multicolorfits_viewer().configure_traits()

def gui():
    """
    Call this function to run the GUI. 
    """
    multicolorfits_viewer().configure_traits()

if __name__ == '__main__':
    #ControlPanel().configure_traits()
    multicolorfits_viewer().configure_traits()



### TODO wishlist ###

##
## GUI
#* Put fields for including titles - overall title on main image, colored interior titles for individual images.
#* streamline the circular dependencies (datamin/datamax and percent_min/percent_max)
#* incorporate ability to regrid the images -- OR new version that loads in data arrays that haven't yet been saved to .fits
#* Options for plot text font (size, dropdown menu for font type which queries those available)
#* Include options for tick length/width
#* Include options for displaying coordinate grid (thin lines over the image) - color and maybe transparency to control on/off
#* Option for GLON/GLAT or others
#* Button to plot final image in greyscale (for printing), and maybe for colorblindness varieties?
#* Simple on-the-fly smoothing option (to decrease appearance of noise)
#* Save/restore state, similar to DS9 backup option
#* Possibly enable significantly downsampled image thumbnails for a 'fast preview mode', to speed up interactive limit adj.
#* Add option to select between for ICRS (RA/DEC) and Galactic (GLON/GLAT)
#* rcParams editor

