# 2.0

2019-09-12

* Updated to use in python 3 or python 2. (Tested in 3.6.8 and 2.7.6)
* Updated to PyQt5 
* Preferred usage is now either 1. run multicolorfits.py from terminal, or 2. import and run GUI with multicolorfits.mcf_gui()  
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


