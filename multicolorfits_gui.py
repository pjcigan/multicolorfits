"""
Some resources for now: traits(ui) and chaco --> Though probably don't want to use Chaco asmy normal code uses matplotlib
https://docs.enthought.com/traitsui/tutorials/traits_ui_scientific_app.html
http://scipy-cookbook.readthedocs.io/items/EmbeddingInTraitsGUI.html
http://www.scipy-lectures.org/advanced/traits/index.html#example
https://gist.github.com/pierre-haessig/9838326

(Examples)
http://henrysmac.org/blog/2014/8/19/demo-of-enthoughts-traitsui-with-matplotlib-and-a-popup-menu.html
Eric Tollerud's traitsUI fitting GUI  https://pythonhosted.org/PyModelFit/over.html#tutorial-examples

Click interaction with matplotlib: http://scipy-cookbook.readthedocs.io/items/Matplotlib_Interactive_Plotting.html

QT based toolbar  https://stackoverflow.com/questions/18468390/creating-matplotlib-toolbar-in-python-traits-editor
"""


import numpy as np
try: import astropy.io.fits as pyfits; import astropy.wcs as pywcs; from astropy.wcs import WCS
except ImportError: import pyfits; import pywcs; from pywcs import WCS
from astropy.visualization import ZScaleInterval, LinearStretch, SqrtStretch, SquaredStretch, LogStretch, PowerStretch, PowerDistStretch, SinhStretch, AsinhStretch, ManualInterval
#from astropy import units as u
from scipy.stats import percentileofscore
import matplotlib
from pyface.api import FileDialog, OK
from pyface.qt import QtGui, QtCore

matplotlib.use('Qt4Agg') 
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT

from matplotlib import pyplot as plt
from matplotlib.figure import Figure
from matplotlib.image import AxesImage
from matplotlib.axes import Axes
from matplotlib.widgets import AxesWidget
#from matplotlib.colors import to_hex #Not in mpl<2.0.  Use try statement below, otherwise define workaround

from traits.api import Any, Instance
from traitsui.qt4.editor import Editor
from traitsui.qt4.basic_editor_factory import BasicEditorFactory

from traits.api import HasTraits, Button, Instance, List, Str, Enum, Float, File, Any
from traits.api import Array, CInt, Int, CFloat, on_trait_change, Range, Dict, Button
from traitsui.api import View, Item, Group, VGroup, HSplit, CheckListEditor, HGroup, Handler, StatusItem, TextEditor
from traitsui.ui_info import UIInfo

from traitsui.api import ColorTrait,ColorEditor

#These two for changing the cursor
from PyQt4.QtCore import Qt
from PyQt4.QtGui import QApplication, QMenu
# Cursor shapes: http://ftp.ics.uci.edu/pub/centos0/ics-custom-build/BUILD/PyQt-x11-gpl-4.7.2/doc/html/qcursor.html

from skimage import color as ski_color
from skimage.exposure import adjust_gamma, rescale_intensity, adjust_sigmoid

import colorsys


##########################################################################
##### ----- Put any custom color or colormap definitions here ----- #####
##########################################################################

#e.g.:  
#mycustomred='#C11B17'
#execfile('/path/to/custom/definitions/file.py-txt-whatever')


##########################################################################


scaling_fns={'linear':LinearStretch, 'sqrt':SqrtStretch, 'squared':SquaredStretch, 'log':LogStretch, 'power':PowerDistStretch, 'sinh':SinhStretch, 'asinh':AsinhStretch}

def getcdelts(hdrin,getrot=False):
    try: 
        cdelt1=float(hdrin['CDELT1']); cdelt2=float(hdrin['CDELT2'])
        try: 
            #Checks for PC matrix, which will modify the CDELT values
            pc1_1=hdrin['PC1_1']; pc1_2=hdrin['PC1_2']; pc2_1=hdrin['PC2_1']; pc2_2=hdrin['PC2_2']; 
            cdelt1*=np.sqrt(pc1_1**2+pc1_2**2)*np.sign(pc1_1); cdelt2*=np.sqrt(pc2_1**2+pc2_2**2)*np.sign(pc2_2)
            crota=np.arctan(pc1_2/pc1_1)*180./np.pi #CROTA(2) is in degrees, from North
        except: 
            try: crota=hdrin['CROTA2']
            except: 
                try: crota=hdrin['CROTA']  
                except: crota=0.
    except:
        try: cd1_1=float(hdrin['CD1_1']); cd1_2=float(hdrin['CD1_2']); cd2_1=float(hdrin['CD2_1']); cd2_2=float(hdrin['CD2_2']); 
        except: raise(Exception('Header does not contain CDELT2 or CD2_2 cards...'))
        cdelt1=np.sqrt(cd1_1**2+cd1_2**2)*np.sign(cd1_1); cdelt2=np.sqrt(cd2_1**2+cd2_2**2)*np.sign(cd2_2)
        crota=np.arctan(cd1_2/cd1_1)*180./np.pi #CROTA(2) is in degrees, from North
    if getrot==False: return cdelt1,cdelt2
    else: return cdelt1,cdelt2,crota

def force_hdr_to_2D(hdrin):
    hdr2D=hdrin.copy()
    for item in ['CRVAL3','CRPIX3','CDELT3','CTYPE3','CUNIT3','NAXIS3','CRVAL4','CRPIX4','CDELT4','CUNIT4','CTYPE4', 'NAXIS4', 'PC03_01','PC04_01','PC03_02','PC04_02','PC01_03','PC02_03','PC03_03','PC04_03','PC01_04','PC02_04','PC03_04','PC04_04']: 
        #del hdr2D[item]
        try: hdr2D.remove(item) 
        except: pass
        hdr2D['NAXIS']=2
    return hdr2D

def force_hdr_floats(hdrin):
    for fc in ['CRPIX1','CRPIX2','CRVAL1','CRVAL2','CDELT1','CDELT2','CD1_1','CD1_2','CD2_1','CD2_2','CROTA','CROTA2','EQUINOX']: 
        try: hdrin.set(fc,float(hdrin[fc])) #Sometimes people save header values as strings, which messes up the WCS...
        except: pass

def hex_to_rgb(hexstring):
    #From http://stackoverflow.com/a/214657
    hexstring = hexstring.lstrip('#')
    lv = len(hexstring)
    return tuple(int(hexstring[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))

def hex_to_hsv(hexstring):
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
    #convenience function to calculate the inverse color (opposite on the color wheel)
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
    
def colorize_image(image, colorvals, colorintype='hsv',dtype=np.float64,gammacorr_color=1):
    ### Add color of the given hue to an RGB greyscale image.
    #colorintype='hsv' [0..1,0..1,0..1],  'rgb' [0..255,0..255,0..255], or 'hex' '#XXXXXX'
    #dtype defaults tostandard numpy float, but may want to lower to e.g. float32 for large images (>~1000x1000)
    #gammacorr_color=(float).  To use color as-is, leave as 1.  To gamma-correct color (e.g., to match gamma for scaled image), specify factor
    if colorintype not in ['hsv', 'hsv_dict', 'rgb', 'hex']: raise Exception("  colorintype must be 'hsv', 'hsv_dict', 'rgb', or 'hex'")
    hsv = ski_color.rgb2hsv(image).astype(dtype)
    if colorintype.lower()=='rgb': colorvals=ski_color.rgb2hsv(colorvals).astype(dtype)
    elif colorintype.lower()=='hex': colorvals=np.array(hex_to_hsv(colorvals)).astype(dtype) #from custom_colormaps.py
    if colorintype.lower()=='hsv_dict': hue,saturation,v=colorvals['hue'],colorvals['sat'],colorvals['v'],
    else: hue,saturation,v=colorvals
    if gammacorr_color!=1: 
        hue,saturation,v=colorsys.rgb_to_hsv(*np.array(colorsys.hsv_to_rgb(hue,saturation,v))**gammacorr_color)
    hsv[:, :, 2] *= v 
    hsv[:, :, 1] = saturation
    hsv[:, :, 0] = hue
    return ski_color.hsv2rgb(hsv).astype(dtype)

def combine_multicolor(im_list_colorized,gamma=2.2,inverse=False):
    #Combines input colorized RGB images [:,:,3] into one intensity-rescaled 
    #im_list_colorized: list of colorized images.  e.g., [halpha_purple,co21_orange,sio54_teal]
    #gamma = value used for gamma correction ^1/gamma [default=2.2].  If inverse=True, will automatically use ^gamma instead.
    #inverse=True  will invert the scale so that white is the background
    combined_RGB=LinearStretch()(np.nansum(im_list_colorized,axis=0))
    RGB_maxints=tuple(np.nanmax(combined_RGB[:,:,i]) for i in [0,1,2])
    for i in [0,1,2]: 
        combined_RGB[:,:,i]=np.nan_to_num(rescale_intensity(combined_RGB[:,:,i], out_range=(0, combined_RGB[:,:,i].max()/np.max(RGB_maxints) )));
    combined_RGB=LinearStretch()(combined_RGB**(1./gamma)) #gamma correction
    if inverse==True: combined_RGB=1.-combined_RGB #gamma correction
    return combined_RGB

class _MPLFigureEditor(Editor):
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
    """This is the control panel where the various parameters for the images are specified
    """
    
    gamma=2.2
    
    fitsfile    = File(filter=[u"*.fits"])
    image_figure = Instance(Figure, ())
    image = Array()
    image_axes = Instance(Axes)
    image_axesimage = Instance(AxesImage)
    image_xsize = Int(256); image_ysize = Int(256)
    
    datamin= Float(0.0,auto_set=False,enter_set=True) #Say, in mJy
    datamax= Float(1.0,auto_set=False,enter_set=True) #auto_set=input set on each keystroke, enter_set=set after Enter
    perc_min=Range(value=0.0, low=0.0, high=100.); 
    perc_max=Range(value=100.0, low=0.0, high=100.) #Percentile of data values for rescaling
    minmaxbutton=Button('Min/Max')
    zscalebutton=Button('Zscale')
    
    image_scale=Str('linear')
    scale_dropdown=Enum(['linear','sqrt','squared','log','power','sinh','asinh'])
    
    imagecolor=Str('#FFFFFF') 
    imagecolor_picker=ColorTrait((255,255,255))
    
    #plotbeam_button=Button('Add Beam (FWHM)')
    
    plotbutton_individual = Button(u"Plot Single")
    plotbutton_inverted_individual = Button(u"Plot Inverted Single")
    clearbutton_individual = Button(u"Clear Single")
    
    status_string_left=Str('')
    status_string_right=Str('')
    
    def __init__(self):
        self._init_params() #Set placeholder things like the WCS, tick color, map units...
        self.image = self._fresh_image() #Sets a blank image
        self.image_axes = self.image_figure.add_subplot(111,aspect=1)
        self.image_axesimage = self.image_axes.imshow(self.image, cmap='gist_gray',origin='lower',interpolation='nearest')
        self.image_axes.axis('off')
    
    view = View(
      HSplit(
        VGroup( 
            Item("fitsfile", label=u"Select 2D FITS file", show_label=True), #,height=100),
            
            HGroup(
              VGroup(
                     Item('plotbutton_individual', tooltip=u"Plot the single image",show_label=False),
                     Item('plotbutton_inverted_individual', tooltip=u"Plot the single inverted image",show_label=False),
                     Item('clearbutton_individual', tooltip=u"Clear the single image",show_label=False),
                     Item('_'),
                     
                     Item('imagecolor',label='Image Color',show_label=True, \
                      tooltip='Color of ticks: standard name float[0..1], or #hex', \
                      editor=TextEditor(auto_set=False, enter_set=True,)),
                     Item('imagecolor_picker',label='Pick',show_label=True,editor=ColorEditor()),
                     
                     ),
              Item('image_figure', editor=MPLFigureEditor(), show_label=False, width=300, height=300,resizable=True),
            ),
            HGroup(Item('datamin', tooltip=u"Minimum data val for scaling", show_label=True),
                   Item('datamax', tooltip=u"Maximum data val for scaling", show_label=True)),
            Item('perc_min', tooltip=u"Min. percentile for scaling", show_label=True),
            Item('perc_max', tooltip=u"Max. percentile for scaling", show_label=True),
            HGroup(Item('minmaxbutton', tooltip=u"Reset to data min/max", show_label=False),
                   Item('zscalebutton', tooltip=u"Compute scale min/max from zscale algorithm", show_label=False)),
            Item('scale_dropdown',label='Scale',show_label=True),
                
        ), #End of Left column
        
      ), #End of HSplit
      resizable=True, #height=0.75, width=0.75, #title=u"Multi-Color Image Combiner", 
      handler=MPLInitHandler,
      statusbar = [StatusItem(name = 'status_string_left', width = 0.5), StatusItem(name = 'status_string_right', width = 0.5)]
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
        
        self.datamax_initial = np.asscalar(np.nanmax(self.data))
        self.datamin_initial = np.asscalar(np.nanmin(self.data))
        self.datamax = np.asscalar(np.nanmax(self.data))
        self.datamin = np.asscalar(np.nanmin(self.data))
        
        self.in_use=True
        
        
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
                self.status_string_right = 'Image color changed to '+self.imagecolor
            except: self.status_string_right = "Color name %s not recognized.  Must be standard mpl.colors string, float[0..1] or #hex string"%(self.imagecolor) 
        try: self.imagecolor_picker=hex_to_rgb(to_hex(self.imagecolor)) #update the picker color...
        except: pass
        ### self.image_greyRGB and self.image_colorRGB may not yet be instantiated if the color is changed before clicking 'plot'
        try: self.image_colorRGB=colorize_image(self.image_greyRGB,self.imagecolor,colorintype='hex',gammacorr_color=self.gamma)
        except: pass
        try: self.image_axesimage.set_data(self.image_colorRGB**(1./self.gamma))
        except: pass
        self.in_use=True
        self.image_figure.canvas.draw()
    
    #@on_trait_change('imagecolor_picker')
    #def update_imagecolorpicker(self):
    def _imagecolor_picker_changed(self):
        #print self.tickcolor_picker.name()
        self.imagecolor=self.imagecolor_picker.name()
    
    
    #@on_trait_change('perc_min')
    #def update_scalepercmin(self):
    def _perc_min_changed(self):
        self.datamin = np.nanpercentile(self.data,self.perc_min)
        self.data_scaled=(scaling_fns[self.image_scale]() + ManualInterval(vmin=self.datamin,vmax=self.datamax))(self.data)
        self.image_greyRGB=ski_color.gray2rgb(adjust_gamma(self.data_scaled,self.gamma))
        self.image_colorRGB=colorize_image(self.image_greyRGB,self.imagecolor,colorintype='hex',gammacorr_color=self.gamma)
        self.image_axesimage.set_data(self.image_colorRGB**(1./self.gamma))
        self.image_figure.canvas.draw()
        self.status_string_right = "Updated scale using percentiles"
    #@on_trait_change('perc_max')
    #def update_scalepercmax(self):
    def _perc_max_changed(self):
        self.datamax = np.nanpercentile(self.data,self.perc_max)
        self.data_scaled=(scaling_fns[self.image_scale]() + ManualInterval(vmin=self.datamin,vmax=self.datamax))(self.data)
        self.image_greyRGB=ski_color.gray2rgb(adjust_gamma(self.data_scaled,self.gamma))
        self.image_colorRGB=colorize_image(self.image_greyRGB,self.imagecolor,colorintype='hex',gammacorr_color=self.gamma)
        self.image_axesimage.set_data(self.image_colorRGB**(1./self.gamma))
        self.image_figure.canvas.draw()
        self.status_string_right = "Updated scale using percentiles"
    
    ### Very slow to update datamin and datamax as well as percs... Can comment these if desired and just hit plot after datamin
    
    #@on_trait_change('datamin')
    #def update_datamin(self): self.perc_min=np.round(percentileofscore(self.data.ravel(),self.datamin,kind='strict'),2)
    def _datamin_changed(self): self.perc_min=np.round(percentileofscore(self.data.ravel(),self.datamin,kind='strict'),2)
    #@on_trait_change('datamax')
    #def update_datamax(self): self.perc_max=np.round(percentileofscore(self.data.ravel(),self.datamax,kind='strict'),2)
    def _datamax_changed(self): self.perc_max=np.round(percentileofscore(self.data.ravel(),self.datamax,kind='strict'),2)
    
    #@on_trait_change('scale_dropdown')
    #def update_image_scale(self):
    def _scale_dropdown_changed(self):
        self.image_scale=self.scale_dropdown
        #self.norm=ImageNormalize(self.sregion,stretch=scaling_fns[self.image_scale]() )
        self.data_scaled=(scaling_fns[self.image_scale]() + ManualInterval(vmin=self.datamin,vmax=self.datamax))(self.data)
        #*** Instead, should I just integrate my imscale class here instead of astropy? ...
        self.image_greyRGB=ski_color.gray2rgb(adjust_gamma(self.data_scaled,self.gamma))
        self.image_colorRGB=colorize_image(self.image_greyRGB,self.imagecolor,colorintype='hex',gammacorr_color=self.gamma)
        
        self.image_axesimage.set_data(self.image_colorRGB**(1./self.gamma))
        
        self.in_use=True
        
        self.image_figure.canvas.draw()
        self.status_string_right = 'Image scale function changed to '+self.image_scale
    
    def _minmaxbutton_fired(self): 
        self.datamin=self.datamin_initial; self.datamax=self.datamax_initial
        #self.image_axesimage.norm.vmin=self.datamin
        #self.image_axesimage.norm.vmax=self.datamax
        self.perc_min=np.round(percentileofscore(self.data.ravel(),self.datamin,kind='strict'),2)
        self.perc_max=np.round(percentileofscore(self.data.ravel(),self.datamax,kind='strict'),2)
        #self.image_figure.canvas.draw()
        self.status_string_right = "Scale reset to min/max"
    
    def _zscalebutton_fired(self):
        tmpZscale=ZScaleInterval().get_limits(self.data)
        self.datamin=float(tmpZscale[0])
        self.datamax=float(tmpZscale[1])
        self.perc_min=np.round(percentileofscore(self.data.ravel(),self.datamin,kind='strict'),2)
        self.perc_max=np.round(percentileofscore(self.data.ravel(),self.datamax,kind='strict'),2)
        #self.image_figure.canvas.draw()
        self.status_string_right = "Min/max determined by zscale"
    
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
        self.image_colorRGB=colorize_image(self.image_greyRGB,self.imagecolor,colorintype='hex',gammacorr_color=self.gamma)
        self.image_axesimage.set_data(self.image_colorRGB**(1./self.gamma))
        ###Using this set instead properly updates the axes labels to WCS, but the home zoom button won't work
        #self.image_figure.clf()
        #self.image_axes = self.image_figure.add_subplot(111,aspect=1)#,projection=self.wcs)
        #self.image_axesimage = self.image_axes.imshow(self.image, cmap=self.image_cmap,origin='lower',interpolation='nearest', norm=self.norm)
        
        self.perc_min=np.round(percentileofscore(self.data.ravel(),self.datamin,kind='strict'),2)
        self.perc_max=np.round(percentileofscore(self.data.ravel(),self.datamax,kind='strict'),2)
        
        self.in_use=True
        
        #self.update_radecpars()
        self.image_figure.canvas.draw()
        self.status_string_right = "Plot updated"

    def _plotbutton_inverted_individual_fired(self):
        try: self.data
        except: self.status_string_right = "No fits file loaded yet!"; return
        self.data_scaled=(scaling_fns[self.image_scale]() + ManualInterval(vmin=self.datamin,vmax=self.datamax))(self.data)
        self.image_greyRGB=ski_color.gray2rgb(adjust_gamma(self.data_scaled,self.gamma))
        self.image_colorRGB=colorize_image(self.image_greyRGB,hexinv(self.imagecolor),colorintype='hex',gammacorr_color=self.gamma)
        self.image_axesimage.set_data(1.-self.image_colorRGB**(1./2.2)) 
        self.image_axesimage.set_data(combine_multicolor([self.image_colorRGB,],gamma=self.gamma,inverse=True))  
        self.perc_min=np.round(percentileofscore(self.data.ravel(),self.datamin,kind='strict'),2)
        self.perc_max=np.round(percentileofscore(self.data.ravel(),self.datamax,kind='strict'),2)
        self.in_use=True
        self.image_figure.canvas.draw()
        self.status_string_right = "Plot updated"
    
    def _clearbutton_individual_fired(self): 
        try: del self.data, self.data_scaled, self.image_greyRGB; self.image_colorRGB #In case clear already pressed once
        except: pass
        self.in_use=False
        self.image_figure.clf()
        self.image = self._fresh_image()
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
        #print event
        x = int(np.round(event.xdata)) #xdata is the actual pixel position.  xy is in 'display space', i.e. pixels in the canvas
        y = int(np.round(event.ydata))
        #xwcs,ywcs=self.wcs.wcs_pix2world([[x,y]],0)[0]; #print xwcs,ywcs
        if ((x >= 0) and (x < self.image.shape[1]) and
            (y >= 0) and (y < self.image.shape[0])):
            imval = self.image[y, x]
            #self.status_string_right = "x,y=[{},{}], RA,DEC=[{}, {}], value = {:.5g}".format(x, y,xwcs,ywcs, imval)
            self.status_string_right = "x,y[{},{}] = {:.3f},{:.3f}  {:.5g}".format(x, y,event.xdata,event.ydata, imval)
        else: self.status_string_right = ""
        ## left-click: event.button = 1, middle-click: event.button=2, right-click: event.button=3.  
        ## For double-click, event.dblclick = False for first click, True on second
        #print event.button, event.dblclick

    def on_cursor_leave(self, event):
        QApplication.restoreOverrideCursor()
        self.status_string_left = ''

    def on_cursor_enter(self, event): QApplication.setOverrideCursor(Qt.CrossCursor)


##### The main viewer #####

class multicolorfits_viewer(HasTraits):
    """The main window. Has instructions for creating and destroying the app.
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
    
    tickcolor=Str('0.9')#,auto_set=False,enter_set=True) #Apparently need to set to TextEditor explicitly below...
    tickcolor_picker=ColorTrait((230,230,230))
    sexdec=Enum('Sexagesimal','Decimal')
    
    plotbutton_combined = Button(u"Plot Combined")
    plotbutton_inverted_combined = Button(u"Plot Inverted Combined")
    clearbutton_combined = Button(u"Clear Combined")
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
        self.image_axesimage = self.image_axes.imshow(self.image, cmap='gist_gray',origin='lower',interpolation='nearest')
        self.image_axes.set_xlabel(self.xlabel); self.image_axes.set_ylabel(self.ylabel)
        self.image_axes.tick_params(axis='both',color=self.tickcolor) #colors=... also sets label color
        try: self.image_axes.coords.frame.set_color(self.tickcolor) #Updates the frame color.  .coords won't exist until WCS set
        except: 
            [self.image_axes.spines[s].set_color(self.tickcolor) for s in ['top','bottom','left','right']]
    
    view = View(Item("gamma",label=u"Gamma",show_label=True), 
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
                      Item('tickcolor',label='Tick Color',show_label=True, \
                          tooltip='Color of ticks: standard name float[0..1], or #hex', \
                          editor=TextEditor(auto_set=False, enter_set=True,)),
                      Item('tickcolor_picker',label='Pick',show_label=True,editor=ColorEditor()), 
                      Item('sexdec',label='Coordinate Style',tooltip=u'Display coordinates in sexagesimal or decimal', \
                           show_label=True),
                      ),
                    Item('figure_combined', editor=MPLFigureEditor(),show_label=False, width=900, height=800,resizable=True),
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
        self.rapars = self.image_axes.coords[0]
        self.decpars = self.image_axes.coords[1]
        self.rapars.set_ticks(color=self.tickcolor); self.decpars.set_ticks(color=self.tickcolor)
        self.rapars.set_ticks(number=6); self.decpars.set_ticks(number=6); 
        #self.rapars.set_ticklabel(size=8); self.decpars.set_ticklabel(size=8); #size here means the tick length
        ##self.rapars.set_ticks(spacing=10*u.arcmin, color='white', exclude_overlapping=True)
        ##self.decpars.set_ticks(spacing=5*u.arcmin, color='white', exclude_overlapping=True)
        self.rapars.display_minor_ticks(True); #self.rapars.set_minor_frequency(10)
        self.decpars.display_minor_ticks(True)
        if self.sexdec=='Sexagesimal':
            self.rapars.set_major_formatter('hh:mm:ss.ss')
            self.decpars.set_major_formatter('dd:mm:ss.ss')
            #self.rapars.set_separator(('$^\mathrm{H}$', "'", '"'))
            self.rapars.set_separator(('H ', "' ", '" '))
            #format_xcoord=lambda x,y: '{}i$^\mathrm{H}${}{}{}'.format(x[0],x[1],"'",x[2],'"')
            #self.image_axes.format_coord=format_xcoord
        else: 
            self.rapars.set_major_formatter('d.dddddd')
            self.decpars.set_major_formatter('d.dddddd')
        ##self.decpars.ticklabels.set_rotation(45) #Rotate ticklabels
        ##self.decpars.ticklabels.set_color(xkcdrust) #Ticklabel Color
    
    @on_trait_change('tickcolor')
    def update_tickcolor(self): 
        try: 
            #Catch case when you've predefined a color variable in hex string format, e.g., mynewred='#C11B17'
            #--> Need to do this first, otherwise traits throws a fit up the stack even despite the try/except check
            globals()[self.tickcolor] #This check should catch undefined inputs
            self.image_axes.tick_params(axis='both',color=globals()[self.tickcolor])
            self.image_axes.coords.frame.set_color(self.tickcolor)
            self.tickcolor_picker=hex_to_rgb(globals()[self.tickcolor])
            self.status_string_right = 'Tick color changed to '+self.tickcolor
        except:
            try: 
                self.tickcolor=to_hex(self.tickcolor)
                try: self.update_radecpars()
                except: 
                    self.image_axes.tick_params(axis='both',color=to_hex(self.tickcolor)); 
                    self.image_axes.coords.frame.set_color(to_hex(self.tickcolor));
                self.status_string_right = 'Tick color changed to '+self.tickcolor
            except: self.status_string_right = "Color name %s not recognized.  Must be standard mpl.colors string, float[0..1] or #hex string"%(self.tickcolor) 
        try: self.tickcolor_picker=hex_to_rgb(to_hex(self.tickcolor)) #update the picker color...
        except: pass
        self.figure_combined.canvas.draw()
    
    @on_trait_change('tickcolor_picker')
    def update_tickcolorpicker(self):
        #print self.tickcolor_picker.name()
        self.tickcolor=self.tickcolor_picker.name()
    
    @on_trait_change('sexdec')
    def update_sexdec(self): 
        self.update_radecpars()
        self.figure_combined.canvas.draw()
        self.status_string_right = 'Coordinate style changed to '+self.sexdec
    
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
        self.image_axesimage = self.image_axes.imshow(self.combined_RGB, origin='lower',interpolation='nearest')
        
        self.update_radecpars()
        self.figure_combined.canvas.draw()
        self.status_string_right = "Plot updated"
    
    def _plotbutton_inverted_combined_fired(self):
        try: self.panel1.data
        except: self.status_string_right = "No fits file loaded yet!"; return
        self.wcs=WCS(self.panel1.hdr)
        self.hdr=self.panel1.hdr
        self.combined_RGB=combine_multicolor( [pan.image_colorRGB for pan in [self.panel1,self.panel2,self.panel3,self.panel4] if pan.in_use==True], gamma=self.gamma, inverse=True)
        self.figure_combined.clf()
        self.image_axes = self.figure_combined.add_subplot(111,aspect=1,projection=self.wcs)
        self.image_axesimage = self.image_axes.imshow(self.combined_RGB, origin='lower',interpolation='nearest')
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
        self.image_axesimage = self.image_axes.imshow(self.image, cmap='gist_gray', origin='lower',interpolation='nearest')
        self.xlabel='x'; self.ylabel='y'
        self.image_axes.set_xlabel(self.xlabel); self.image_axes.set_ylabel(self.ylabel)
        self.image_axes.tick_params(axis='both',color=self.tickcolor)
        self.image_axes.coords.frame.set_color(self.tickcolor)
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
        #print event
        x = int(np.round(event.xdata)) #xdata is the actual pixel position.  xy is in 'display space', i.e. pixels in the canvas
        y = int(np.round(event.ydata))
        xwcs,ywcs=self.wcs.wcs_pix2world([[x,y]],0)[0]; #print xwcs,ywcs
        if ((x >= 0) and (x < self.image.shape[1]) and
            (y >= 0) and (y < self.image.shape[0])):
            imval = self.image[y, x]
            self.status_string_right = "x,y=[{},{}], RA,DEC=[{}, {}], value = {:.5g}".format(x, y,xwcs,ywcs, imval)
            #self.status_string_right = "x,y[{},{}] = {:.3f},{:.3f}  {:.5g}".format(x, y,event.xdata,event.ydata, imval)
        else: self.status_string_right = ""
        ## left-click: event.button = 1, middle-click: event.button=2, right-click: event.button=3.  
        ## For double-click, event.dblclick = False for first click, True on second
        #print event.button, event.dblclick

    def on_cursor_leave(self, event):
        QApplication.restoreOverrideCursor()
        self.status_string_left = ''

    def on_cursor_enter(self, event): QApplication.setOverrideCursor(Qt.CrossCursor)
    
    def _save_the_image_fired(self):
        dlg = FileDialog(action='save as')
        if dlg.open() == OK: plt.savefig(dlg.path,size=(800,800),dpi=300)
    
    def _save_the_fits_fired(self): 
        #Generate a generic header with correct WCS and comments about the colors that made it
        #... come back and finish this later...
        dlg = FileDialog(action='save as')
        if dlg.open() == OK: pyfits.writeto(dlg.path,np.swapaxes(np.swapaxes(self.combined_RGB,0,2),2,1),self.hdr)
    
    def _print_params_fired(self): 
        print('\n\nRGB Image plot params:')
        pan_i=0
        for pan in [self.panel1,self.panel2,self.panel3,self.panel4]: 
            pan_i+=1
            if pan.in_use==True: 
                print('image%i: '%(pan_i))
                print('    vmin = %.3e , vmax = %.3e, scale = %s'%(pan.datamin,pan.datamax,pan.image_scale))
                print("    image color = '%s'"%(pan.imagecolor))
        print("gamma = %.1f , tick color = '%s'\n"%(self.gamma,self.tickcolor))
        #print('vmin=%.3e, vmax=%.3e'%(self.panel1.datamin,self.panel1.datamax))
        

if __name__ == '__main__':
    #ControlPanel().configure_traits()
    multicolorfits_viewer().configure_traits()



### TODO wishlist ###

#* Put fields for including titles - overall title on main image, colored interior titles for individual images.
#* streamline the circular dependencies (datamin/datamax and perc_min/perc_max)
#* incorporate ability to regrid the images -- OR new version that loads in data arrays that haven't yet been saved to .fits
#* Options for plot text font (size, dropdown menu for font type which queries those available)
#* Include options for tick length/width
#* Include options for displaying coordinate grid (thin lines over the image) - color and maybe transparency to control on/off
#* update manual plotting  param printout button 



