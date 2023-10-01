# 2.1

2023-09-27

* Updated PyQt imports to use PyQt6 by default, to work with py 3.10 and M processor Macs.  PyQt5 is still there as a backup.  PyQt4 should still work for Python 3.4 and lower at the moment, but this may be phased out in future releases.

2022-09-04

* Added buttons in the left side panels to edit image headers
    - This may be useful for tweaking things (e.g. bump reference pixel values to shift an HST image alignment), to declutter HISTORY and COMMENT cards, etc.  You can also try deleting offending header cards if something is preventing proper image combination.
    - Currently, after making edits in the popup window the user must also click the Apply Header Changes button to apply the changes.  Clicking Edit Header again before applying will discard the changes. 
    - Headers can also be saved to / loaded from text files
* Fixed the Save Image button functionality 
* Added options for setting xaxislabel and yaxislabel in the convenience functions plotsinglemulticolorRGB() and comparemulticolorRGB_pureRGB()
* Added fields for x/y axis ticklabel formats.  (See [the astropy documentation](https://docs.astropy.org/en/stable/visualization/wcsaxes/ticks_labels_grid.html#tick-label-format for options) for available formats and syntax)
* Combined image x/y axis labels now automatically generate from header CTYPE1/2 values (when present)
* Added a button for turning minor ticks on/off
* Added simple gaussian smoothing function.  This can be useful for reducing noise, matching resolutions, etc.


# 2.0

2019-09-12

* Updated to use in python 3 or python 2. (Tested in 3.6.8 and 2.7.6)
* Updated to PyQt5 
* Preferred usage is now either 
    1. run multicolorfits.py from terminal, or 
    2. import and run GUI with multicolorfits.mcf_gui()  
* Added several functions to aid scripting (outside of the GUI)
    - reprojection to specified image header (requires reproject package or kapteyn package)
    - cropping by specifying coordinates and desired size
    - color conversions
    - some header wcs handling
    - Quick-look plotting functions to save out final combined images, and to compare to the simple RGB case
* Added an option to disable auto-refresh/update of the preview plots in the left panel after adjusting min/max values.  This will increase speed for large files, but with this option users will need to click the 'Plot' button (or the new 'Update plot' button) after each change to see the results.
* improved handling of image NaNs


# 1.2 

2018-12-12

* Roughly fixed incompatibility with matplotlib 1.5. 
    - defined a replacement to_hex() function
    - separated out rcParams keywords that don't exist in mpl 1.5


# 1.1 

2018-12-11

* Implemented the 'inverse plot' buttons


# 1.0

2018-12-10

* Initial version on github


