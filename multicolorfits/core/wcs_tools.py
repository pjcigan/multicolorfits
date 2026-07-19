"""
FITS header / WCS / coordinate helpers, and pixel- and coordinate-based
image cropping utilities.
"""

import numpy as np
import warnings
import astropy.io.fits as pyfits
import astropy.wcs as pywcs
from astropy.coordinates import SkyCoord

__all__ = [
    'force_header_2d',
    'force_header_3d',
    'force_header_floats',
    'normalize_wcs_meta',
    'deg2dms',
    'dms2deg',
    'deg2hour',
    'hour2deg',
    'dec2sex',
    'sex2dec',
    'angular_distance',
    'get_cdelts',
    'sky_to_pixel',
    'pixel_to_sky',
    'make_simple_header',
    'get_cd_matrix',
    'deg_per_pixel',
    'arcsec_per_pixel',
    'sterad_per_pixel',
    'squeeze_image',
    'header_coord_grids',
    'beam_params_arcsec',
    'pixels_per_beam',
    'crop_image',
    'crop_cube',
    'crop_image_sky',
    'crop_cube_sky',
]


def _axis_cards(axes):
    """WCS card names belonging to the given FITS axis numbers.

    Ported/adapted from skyplothelper for robust header stripping.
    """
    cards = []
    for ax in axes:
        for base in ('NAXIS', 'CTYPE', 'CRVAL', 'CRPIX', 'CDELT', 'CUNIT', 'CROTA'):
            cards.append('%s%d' % (base, ax))
        for i in (1, 2, 3, 4):
            cards += ['PC%d_%d' % (ax, i), 'PC%d_%d' % (i, ax),
                      'CD%d_%d' % (ax, i), 'CD%d_%d' % (i, ax),
                      'PC%02d_%02d' % (ax, i), 'PC%02d_%02d' % (i, ax)]
    return cards


def _has_higher_axis_cards(header, axes=(3, 4)):
    """True if a nominally 2D header still carries higher-axis WCS cards."""
    return any(card in header for card in _axis_cards(axes))


_FLOAT_WCS_CARDS = (
    'CRPIX1', 'CRPIX2', 'CRVAL1', 'CRVAL2', 'CDELT1', 'CDELT2',
    'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2', 'PC1_1', 'PC1_2', 'PC2_1', 'PC2_2',
    'CROTA', 'CROTA2', 'EQUINOX', 'LONPOLE', 'LATPOLE',
)


def normalize_wcs_meta(hdrin):
    """Coerce integer WCS metadata cards so astropy.wcs does not FITS-fix them.

    Some pipelines (e.g. Chandra/HST) write ``WCSAXES`` or ``NAXIS*`` as floats;
    astropy emits ``FITSFixedWarning`` when it has to correct those on WCS build.
    """
    for key in ('NAXIS', 'NAXIS1', 'NAXIS2', 'NAXIS3', 'NAXIS4', 'WCSAXES'):
        if key not in hdrin:
            continue
        try:
            hdrin[key] = int(float(hdrin[key]))
        except (TypeError, ValueError):
            pass
    try:
        if int(hdrin.get('NAXIS', 2)) == 2:
            hdrin['WCSAXES'] = 2
    except (TypeError, ValueError):
        pass
    return hdrin


def force_header_2d(hdrin):
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
    hdr2D = hdrin.copy()
    for item in _axis_cards((3, 4)):
        if item in hdr2D:
            del hdr2D[item]
    hdr2D['NAXIS'] = 2
    if 'WCSAXES' in hdr2D:
        hdr2D['WCSAXES'] = 2
    return normalize_wcs_meta(hdr2D)


def force_header_3d(hdrin):
    """Return copy of header with all 4th-axis cards removed (NAXIS=3).

    Ported/adapted from skyplothelper for parity with modern header tools.
    """
    hdr3D = hdrin.copy()
    for item in _axis_cards((4,)):
        if item in hdr3D:
            del hdr3D[item]
    hdr3D['NAXIS'] = 3
    if 'WCSAXES' in hdr3D:
        hdr3D['WCSAXES'] = 3
    return normalize_wcs_meta(hdr3D)


def force_header_floats(hdrin):
    """
    Helper function to force various header values to floats.
    Sometimes programs save header values as strings, which messes up the WCS...

    Parameters
    ----------
    hdrin : astropy.io.fits.header
        Header object
    """
    for fc in _FLOAT_WCS_CARDS:
        try:
            hdrin.set(fc, float(hdrin[fc]))
        except Exception:
            pass
    return hdrin


def deg2dms(valin, return_type=list, str_delimiter=':', str_decimal_places=2,
            str_zero_pad=True):
    """Convert decimal degrees to DMS components.

    Ported/adapted from skyplothelper to fix sign handling and 60-second
    rollover edge-cases.
    """
    sign = -1 if valin < 0 else 1
    absval = abs(valin)
    ddeg = int(absval)
    dmins_f = (absval - ddeg) * 60
    dmins = int(dmins_f)
    dsec = (dmins_f - dmins) * 60

    if return_type is str or return_type == 'str':
        dsec_r = round(dsec, str_decimal_places)
        if dsec_r >= 60:
            dsec_r -= 60
            dmins += 1
        if dmins >= 60:
            dmins -= 60
            ddeg += 1
        if str_zero_pad:
            deg_pad = '0>2'
            sec_pad = '0>%d' % (3 + str_decimal_places if str_decimal_places > 0 else 2)
        else:
            deg_pad = ''
            sec_pad = ''
        prefix = '-' if sign < 0 else ''
        return '{6}{0:{7}d}{3}{1:{7}d}{3}{2:{5}.{4}f}'.format(
            ddeg, dmins, dsec_r, str_delimiter, str_decimal_places, sec_pad,
            prefix, deg_pad)

    if dsec >= 60 - 1e-9:
        dsec = 0.0
        dmins += 1
    if dmins >= 60:
        dmins -= 60
        ddeg += 1
    if sign < 0:
        if ddeg != 0:
            ddeg = -ddeg
        elif dmins != 0:
            dmins = -dmins
        else:
            dsec = -dsec
    return return_type([int(ddeg), int(dmins), dsec])


def dms2deg(valin):
    """Convert DMS string ('dd:mm:ss.s') or [d, m, s] list to decimal degrees.

    Ported/adapted from skyplothelper.
    """
    if isinstance(valin, str):
        stripped = valin.strip()
        sign = -1 if stripped.startswith('-') else 1
        body = stripped.lstrip('+-')
        comps = [float(v) for v in body.replace('d', ':').replace('m', ':')
                 .replace('s', '').split(':')]
        return sign * (abs(comps[0]) + comps[1] / 60. + comps[2] / 3600.)
    comps = list(valin)
    sign = 1
    for c in comps:
        if c < 0:
            sign = -1
            break
        if c > 0:
            break
    return sign * (abs(comps[0]) + abs(comps[1]) / 60. + abs(comps[2]) / 3600.)


def deg2hour(valin):
    """Convert decimal degrees to [hours, minutes, seconds] list."""
    rmins, rsec = divmod(24. / 360 * valin * 3600, 60)
    rh, rmins = divmod(rmins, 60)
    return [int(rh), int(rmins), rsec]


def hour2deg(valin):
    """Convert HMS string ('hh:mm:ss.s') or [h, m, s] list to decimal degrees."""
    if isinstance(valin, str):
        cleaned = valin.lower().replace('h', ':').replace('m', ':').replace('s', '')
        comps = [float(v) * 360. / 24 for v in cleaned.split(':')]
    else:
        comps = [v * 360. / 24 for v in valin]
    return comps[0] + comps[1] / 60. + comps[2] / 3600.


def dec2sex(rain, decin, as_string=False, decimal_places=2):
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
    hms = deg2hour(float(rain))
    dms = deg2dms(float(decin))
    if as_string:
        ras = '{0:0>2d}:{1:0>2d}:{2:0>{4}.{3}f}'.format(
            int(hms[0]), int(hms[1]), hms[2], decimal_places, decimal_places + 3)
        decs = deg2dms(float(decin), return_type=str, str_decimal_places=decimal_places)
        return [ras, decs]
    return hms, dms


def sex2dec(rain, decin):
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
    raout = hour2deg(rain)
    dec_clean = str(decin).lower().replace('d', ':').replace('m', ':').replace('s', '')
    if ':' in dec_clean:
        dec = [float(val) for val in dec_clean.split(':')]
    else:
        dec = [float(val) for val in str(decin).split(' ')]
    if dec[0] < 0 or str(decin).strip().startswith('-'):
        decout = dec[0] - dec[1] / 60. - dec[2] / 3600.
    else:
        decout = dec[0] + dec[1] / 60. + dec[2] / 3600.
    return [raout, decout]


def angular_distance(coords1_deg, coords2_deg, pythag_approx=False, returncomponents=False,
                    input_precision=np.float64):
    """Angular distance between two sky positions using the Vincenty formula.

    Ported/adapted from skyplothelper.
    """
    coords1_rad = np.array(coords1_deg, dtype=input_precision) * (np.pi / 180)
    coords2_rad = np.array(coords2_deg, dtype=input_precision) * (np.pi / 180)
    if coords1_rad.ndim == 1:
        lon1, lat1 = coords1_rad[0], coords1_rad[1]
        lon2, lat2 = coords2_rad[0], coords2_rad[1]
    else:
        lon1, lat1 = coords1_rad[:, 0], coords1_rad[:, 1]
        lon2, lat2 = coords2_rad[:, 0], coords2_rad[:, 1]
    lat_mean = 0.5 * (lat1 + lat2)
    if pythag_approx:
        dRA_rad = (lon2 - lon1) * np.cos(lat_mean)
        dDEC_rad = lat2 - lat1
        sep_rad = np.sqrt(dRA_rad**2 + dDEC_rad**2)
        if returncomponents:
            return np.degrees(sep_rad), np.degrees(dRA_rad), np.degrees(dDEC_rad)
        return np.degrees(sep_rad)
    sin_dlon = np.sin(lon1 - lon2)
    cos_dlon = np.cos(lon1 - lon2)
    sin_lat1, cos_lat1 = np.sin(lat1), np.cos(lat1)
    sin_lat2, cos_lat2 = np.sin(lat2), np.cos(lat2)
    A = cos_lat2 * sin_dlon
    B = cos_lat1 * sin_lat2 - sin_lat1 * cos_lat2 * cos_dlon
    sep_rad = np.arctan2(np.hypot(A, B), sin_lat1 * sin_lat2 + cos_lat1 * cos_lat2 * cos_dlon)
    if returncomponents:
        dRA_rad = (lon2 - lon1) * np.cos(lat_mean)
        dDEC_rad = lat2 - lat1
        return np.degrees(sep_rad), np.degrees(dRA_rad), np.degrees(dDEC_rad)
    return np.degrees(sep_rad)


def get_cdelts(hdrin, getrot=False):
    """
    Function to calculate CDELT1 and CDELT2 from the input header PCx_x cards.

    Parameters
    ----------
    hdrin : astropy.io.fits.header
        Header object
    getrot : bool
        Specifies whether to return the rotation value crota in the output

    Returns
    -------
    float, float (, float)
        CDELT1, CDELT2 (, CROTA)
    """
    try:
        cdelt1 = float(hdrin['CDELT1']); cdelt2 = float(hdrin['CDELT2'])
        try:
            # Checks for PC matrix, which will modify the CDELT values
            pc1_1 = hdrin['PC1_1']; pc1_2 = hdrin['PC1_2']; pc2_1 = hdrin['PC2_1']; pc2_2 = hdrin['PC2_2']
            cdelt1 *= np.sqrt(pc1_1**2 + pc1_2**2) * np.sign(pc1_1)
            cdelt2 *= np.sqrt(pc2_1**2 + pc2_2**2) * np.sign(pc2_2)
            crota = np.degrees(np.arctan2(pc1_2, pc1_1))  # robust quadrant handling
        except Exception:
            try: crota = hdrin['CROTA2']
            except Exception:
                try: crota = hdrin['CROTA']
                except Exception: crota = 0.
    except Exception:
        try:
            cd1_1 = float(hdrin['CD1_1']); cd1_2 = float(hdrin['CD1_2'])
            cd2_1 = float(hdrin['CD2_1']); cd2_2 = float(hdrin['CD2_2'])
        except Exception:
            raise Exception('Header does not contain CDELT2 or CD2_2 cards...')
        cdelt1 = np.sqrt(cd1_1**2 + cd1_2**2) * np.sign(cd1_1)
        cdelt2 = np.sqrt(cd2_1**2 + cd2_2**2) * np.sign(cd2_2)
        crota = np.degrees(np.arctan2(cd1_2, cd1_1))  # robust quadrant handling
    if getrot is False: return cdelt1, cdelt2
    else: return cdelt1, cdelt2, crota


def sky_to_pixel(headerin, rain, decin, precise=False, checksys=False, incoordsys='fk5',
                incoordequinox='J2000.0', forceimagesys=None, originindex=0):
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
    if checksys is True:
        try:
            try: headerframe = headerin['RADESYS']
            except Exception: headerframe = headerin['RADECSYS']
        except Exception:
            if forceimagesys is not None: headerframe = forceimagesys
            else: raise Exception('Input header has no valid RADESYS/RADECSYS frame, and forceimagesys not specified...')
        if headerframe.lower() == incoordsys.lower(): pass
        else:
            # Convert coordinates between RA/DEC systems
            incoordSC = SkyCoord(ra=rain, dec=decin, unit='deg', frame=incoordsys.lower(),
                                 equinox=incoordequinox).transform_to(headerframe.lower())
            rain = float(incoordSC.ra.value); decin = float(incoordSC.dec.value)

    try:
        wcstemp = pywcs.WCS(headerin)
        try:
            if headerin['NAXIS'] == 2: pixarr = wcstemp.wcs_world2pix([[rain, decin]], originindex)
            else: pixarr = wcstemp.wcs_world2pix([[rain, decin, 0]], originindex)
        except Exception:
            if headerin['WCSAXES'] == 2: pixarr = wcstemp.wcs_world2pix([[rain, decin]], originindex)
            else: pixarr = wcstemp.wcs_world2pix([[rain, decin, 0]], originindex)
    except Exception:
        wcstemp = pywcs.WCS(naxis=2)
        wcstemp.wcs.crpix = [headerin['CRPIX1'], headerin['CRPIX2']]
        try: wcstemp.wcs.cdelt = list(get_cdelts(headerin))
        except Exception: raise Exception('Invalid WCS CDELTS')
        wcstemp.wcs.crval = [headerin['CRVAL1'], headerin['CRVAL2']]
        wcstemp.wcs.ctype = [headerin['CTYPE1'], headerin['CTYPE2']]
        pixarr = wcstemp.wcs_world2pix([[rain, decin, 0]], originindex)
    if precise is True: return [pixarr[0][0], pixarr[0][1]]
    else: return [int(np.round(pixarr[0][0])), int(np.round(pixarr[0][1]))]


def pixel_to_sky(headerin, xin, yin, outcoordsys='same', outcoordequinox='J2000.0',
                forceimagesys=None, originindex=0):
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
        wcstemp = pywcs.WCS(headerin)
        try:
            if headerin['NAXIS'] == 2: pixarr = wcstemp.wcs_pix2world([[xin, yin]], originindex)
            else: pixarr = wcstemp.wcs_pix2world([[xin, yin, 0]], originindex)
        except Exception:
            if headerin['WCSAXES'] == 2: pixarr = wcstemp.wcs_pix2world([[xin, yin]], originindex)
            else: pixarr = wcstemp.wcs_pix2world([[xin, yin, 0]], originindex)
    except Exception:
        wcstemp = pywcs.WCS(naxis=2)
        wcstemp.wcs.crpix = [headerin['CRPIX1'], headerin['CRPIX2']]
        wcstemp.wcs.crval = [headerin['CRVAL1'], headerin['CRVAL2']]
        wcstemp.wcs.ctype = [headerin['CTYPE1'], headerin['CTYPE2']]
        try: wcstemp.wcs.cdelt = list(get_cdelts(headerin))
        except Exception: raise Exception('Invalid WCS CDELTS')
        pixarr = wcstemp.wcs_pix2world([[xin, yin, 0]], originindex)
    if outcoordsys == 'same':
        return [pixarr[0][0], pixarr[0][1]]
    else:
        try:
            try: headerframe = headerin['RADESYS']
            except Exception: headerframe = headerin['RADECSYS']
        except Exception:
            # User can specify a frame forceimagesys to use (in case no RADESYS in the header, or if it's known to be wrong...)
            if forceimagesys is not None: headerframe = forceimagesys
            else: raise Exception('Input header has no valid RADESYS/RADECSYS frame, and forceimagesys not specified...')
        outcoordSC = SkyCoord(ra=pixarr[0][0], dec=pixarr[0][1], unit='deg', frame=headerframe.lower(),
                              equinox=outcoordequinox).transform_to(outcoordsys.lower())
        return [outcoordSC.ra.value, outcoordSC.dec.value]


def make_simple_header(headerin, naxis=2, radesys=None, equinox=None, pywcsdirect=False):
    """
    Function to make a new 'simple header' from the WCS information in the input header.

    Parameters
    ----------
    headerin : astropy.io.fits.header
        Header object
    naxis : int
        Specifies how many axes the final header should have.  Default=2
    radesys : str
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
    if type(headerin) == str:
        headerin = pyfits.getheader(headerin)
    if pywcsdirect is True:
        wcstemp = pywcs.WCS(header=headerin)
    else:
        wcstemp = pywcs.WCS(naxis=naxis)
        if naxis > 2:
            wcstemp.wcs.crpix = [float(headerin['CRPIX1']), float(headerin['CRPIX2']), float(headerin['CRPIX3'])]
            wcstemp.wcs.crval = [float(headerin['CRVAL1']), float(headerin['CRVAL2']), float(headerin['CRVAL3'])]
            wcstemp.wcs.ctype = [headerin['CTYPE1'], headerin['CTYPE2'], headerin['CTYPE3']]
            try: wcstemp.wcs.cunit = [headerin['CUNIT1'], headerin['CUNIT2'], headerin['CUNIT3']]
            except Exception: pass
            try: wcstemp.wcs.cdelt = list(get_cdelts(headerin)) + [headerin['CDELT3']]
            except Exception: raise Exception('Invalid WCS CDELTS')
        else:
            wcstemp.wcs.crpix = [float(headerin['CRPIX1']), float(headerin['CRPIX2'])]
            wcstemp.wcs.crval = [float(headerin['CRVAL1']), float(headerin['CRVAL2'])]
            wcstemp.wcs.ctype = [headerin['CTYPE1'], headerin['CTYPE2']]
            try: wcstemp.wcs.cunit = [headerin['CUNIT1'], headerin['CUNIT2']]
            except Exception: pass
            try: wcstemp.wcs.cdelt = list(get_cdelts(headerin))
            except Exception: raise Exception('Invalid WCS CDELTS')
        try: crota = get_cdelts(headerin, getrot=True)[-1]  # degrees, from N
        except Exception: raise Exception('Invalid WCS params for CROTAx')
        try: wcstemp.wcs.radesys = headerin['RADESYS']
        except Exception: pass
        try: wcstemp.wcs.equinox = headerin['EQUINOX']
        except Exception: pass
    if radesys is not None: wcstemp.wcs.radesys = radesys  # e.g. 'FK5', 'ICRS'. For manually forcing string, not true reprojection.
    if equinox is not None: wcstemp.wcs.equinox = equinox  # e.g. 2000.0
    simpleheader = wcstemp.to_header()
    if pywcsdirect is False:
        if crota != 0.: simpleheader['CROTA2'] = crota  # Alternative method to just use (deprecated) CROTA2 card
    simpleheader['NAXIS'] = naxis
    try:
        simpleheader['NAXIS1'] = int(headerin['NAXIS1']); simpleheader['NAXIS2'] = int(headerin['NAXIS2'])
    except Exception: pass
    if naxis > 2:
        for card in ['NAXIS3', 'CRPIX3', 'CRVAL3', 'CDELT3', 'CTYPE3', 'CUNIT3', 'SPECSYS', 'ALTRVAL', 'ALTRPIX']:
            try: simpleheader[card] = headerin[card]
            except Exception: pass
    for card in ['CROTA', 'CROTA1', 'CROTA2', 'BSCALE', 'BZERO', 'ZSCALE', 'BMAJ', 'BMIN', 'BPA',
                 'JANSCALE', 'FLUXCONV', 'WAVELEN', 'FREQ', 'RESTFRQ', 'LATPOLE', 'LONPOLE']:
        try: simpleheader[card] = float(headerin[card])
        except Exception: pass
    for card in ['BUNIT', 'OBJECT', 'TELESCOP', 'ZUNITS', 'SPECSYS']:
        try: simpleheader[card] = headerin[card]
        except Exception: pass
    return simpleheader


def get_cd_matrix(hdrin, crot=None):
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
    try:
        cd1_1 = hdrin['CD1_1']; cd1_2 = hdrin['CD1_2']; cd2_1 = hdrin['CD2_1']; cd2_2 = hdrin['CD2_2']
    except Exception:
        if crot is None:
            try: crot = hdrin['CROTA2']
            except Exception:
                try: crot = hdrin['CROTA']
                except Exception: crot = 0.
        try: cdelt1 = float(hdrin['CDELT1']); cdelt2 = float(hdrin['CDELT2'])
        except Exception: raise Exception('Header does not contain CDELT2 or CD2_2 cards...')
        try:
            pc1_1 = hdrin['PC1_1']; pc1_2 = hdrin['PC1_2']; pc2_1 = hdrin['PC2_1']; pc2_2 = hdrin['PC2_2']
            cd1_1 = pc1_1 * cdelt1; cd1_2 = pc1_2 * cdelt1; cd2_1 = pc2_1 * cdelt2; cd2_2 = pc2_2 * cdelt2
        except Exception:
            cd1_1 = cdelt1 * np.cos(crot * np.pi / 180); cd1_2 = cdelt1 * np.sin(crot * np.pi / 180)
            cd2_1 = cdelt2 * -np.sin(crot * np.pi / 180); cd2_2 = cdelt2 * np.cos(crot * np.pi / 180)
        # Even if PC1_1 matrix used, CDELTs should still be included there...
    return cd1_1, cd1_2, cd2_1, cd2_2


def deg_per_pixel(hdrin):
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
    cdelt1, cdelt2 = get_cdelts(hdrin)
    if abs(abs(cdelt1) - abs(cdelt2)) < 1e-6:
        return abs(cdelt2)
    warnings.warn('CDELT1 and CDELT2 are significantly different; returning mean(abs(CDELT1), abs(CDELT2))')
    return np.mean([abs(cdelt1), abs(cdelt2)])


def arcsec_per_pixel(hdrin):
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
    return deg_per_pixel(hdrin) * 3600


def sterad_per_pixel(hdrin):
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
    return np.radians(deg_per_pixel(hdrin)) ** 2


def squeeze_image(data, header=None, verbose=True):
    """Squeeze a FITS image to 2D, cleaning up the header to match.

    Ported/adapted from skyplothelper for robust degenerate-axis handling.
    """
    data = np.asarray(data)
    if data.ndim == 2:
        if header is None:
            return data, None
        hdr_2d = header.copy()
        if _has_higher_axis_cards(hdr_2d) or int(hdr_2d.get('NAXIS', 2)) > 2:
            hdr_2d = force_header_2d(hdr_2d)
        else:
            hdr_2d = normalize_wcs_meta(hdr_2d)
        return data, hdr_2d
    if data.ndim < 2:
        raise ValueError('Data must be at least 2D, got %dD' % data.ndim)
    non_spatial = data.shape[:-2]
    non_degenerate = [i for i, s in enumerate(non_spatial) if s > 1]
    if non_degenerate:
        raise ValueError('Cannot squeeze to 2D: non-degenerate extra axes present: %r' % (data.shape,))
    data_2d = data.squeeze()
    while data_2d.ndim > 2:
        data_2d = data_2d[0]
    if verbose and data.ndim > 2:
        squeezed_axes = data.ndim - 2
        print('squeeze_image: %s -> %s (%d degenerate axis%s removed)' %
              (data.shape, data_2d.shape, squeezed_axes, 'es' if squeezed_axes > 1 else ''))
    if header is None:
        return data_2d, None
    hdr_2d = force_header_2d(header)
    return data_2d, hdr_2d


def header_coord_grids(hdr_or_wcs, shape=None, x=None, y=None, return_1d=False):
    """Per-pixel world (lon/lat) coordinates for an image grid.

    Ported/adapted from skyplothelper for internal coordinate-aware helpers.
    """
    if isinstance(hdr_or_wcs, pywcs.WCS):
        wcs = hdr_or_wcs
        if shape is not None:
            ny, nx = int(shape[0]), int(shape[1])
        elif wcs.pixel_shape is not None:
            nx, ny = int(wcs.pixel_shape[0]), int(wcs.pixel_shape[1])
        elif x is not None and y is not None:
            nx, ny = len(np.atleast_1d(x)), len(np.atleast_1d(y))
        else:
            raise ValueError('WCS input needs shape=(ny, nx) or x/y vectors')
    else:
        wcs = pywcs.WCS(hdr_or_wcs)
        nx, ny = int(hdr_or_wcs['NAXIS1']), int(hdr_or_wcs['NAXIS2'])
    wcs2d = wcs.celestial if getattr(wcs, 'naxis', 2) > 2 else wcs
    if x is None:
        x = np.arange(nx)
    if y is None:
        y = np.arange(ny)
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if return_1d:
        cy = y[len(y) // 2] if len(y) else 0.0
        cx = x[len(x) // 2] if len(x) else 0.0
        lon1d, _ = wcs2d.pixel_to_world_values(x, np.full_like(x, cy))
        _, lat1d = wcs2d.pixel_to_world_values(np.full_like(y, cx), y)
        return np.asarray(lon1d), np.asarray(lat1d)
    xx, yy = np.meshgrid(x, y)
    lon, lat = wcs2d.pixel_to_world_values(xx, yy)
    return np.asarray(lon), np.asarray(lat)


def beam_params_arcsec(hdrin):
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
    return [hdrin['BMAJ'] * 3600., hdrin['BMIN'] * 3600., hdrin['BPA']]


def pixels_per_beam(hdrin):
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
    return np.pi / (4. * np.log(2)) * hdrin['BMAJ'] * hdrin['BMIN'] / deg_per_pixel(hdrin)**2


def crop_image(datain, hdrin, xbounds, ybounds, newref=None, savenew=False, overwrite=False):
    """
    Function to crop a 2D fits image to the specified pixel bounds.

    Parameters
    ----------
    datain : array
        Input fits image array
    hdrin : astropy.io.fits.header
        Input header
    xbounds : list
        [min,max] x-axis pixel limits to use for new image slice
    ybounds : list
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
    if len(datain.shape) > 3: raise Exception('Data array has 4 or more dimensions - reduce them to 2!')
    if len(datain.shape) == 2: pass
    elif len(datain.shape) > 2 and (datain.shape[i] == 1 for i in range(len(datain.shape) - 2)): datain = datain[0, :, :]
    else: raise Exception('File is not flattened 2D!  Use crop_cube.')
    if True in [item < 0 for item in np.concatenate([xbounds, ybounds])]:
        raise Exception('Specified X/Y crop bounds are negative - must be positive.')
    try: cropdata = datain[ybounds[0]:ybounds[1] + 1, xbounds[0]:xbounds[1] + 1]
    except Exception: cropdata = datain[int(np.round(ybounds[0])):int(np.round(ybounds[1] + 1)), int(np.round(xbounds[0])):int(np.round(xbounds[1] + 1))]
    crophdr = hdrin.copy()
    crophdr['NAXIS2'], crophdr['NAXIS1'] = cropdata.shape
    if newref == 'center':
        crophdr['CRPIX1'] = (int(xbounds[1] - xbounds[0])) * .5; crophdr['CRPIX2'] = (int(ybounds[1] - ybounds[0])) * .5
        crophdr['CRVAL1'], crophdr['CRVAL2'] = pixel_to_sky(hdrin, (int(xbounds[0]) + int(xbounds[1] - xbounds[0]) * .5) - 1, (int(ybounds[0]) + int(ybounds[1] - ybounds[0]) * .5) - 1)
        # Next, adjust CROTA2 based on offset - needed for big deviations near poles
        radiff = np.unwrap(np.deg2rad([0, hdrin['CRVAL1'] - crophdr['CRVAL1']]))[1] * 180. / np.pi
        crophdr['CROTA2'] = radiff * np.sin(crophdr['CRVAL2'] * np.pi / 180)  # sin because at equator, DEC=0
        # --> For large changes in CRPIX, there will be some slight offsets/position errors.
    elif newref == 'origin':
        crophdr['CRPIX1'] = 0.; crophdr['CRPIX2'] = 0.
        crophdr['CRVAL1'], crophdr['CRVAL2'] = pixel_to_sky(hdrin, int(xbounds[0]) - 1, int(ybounds[0]) - 1)
        # Next, adjust CROTA2 based on offset - needed for big deviations near poles
        radiff = np.unwrap(np.deg2rad([0, hdrin['CRVAL1'] - crophdr['CRVAL1']]))[1] * 180. / np.pi
        crophdr['CROTA2'] = radiff * np.sin(crophdr['CRVAL2'] * np.pi / 180)  # sin because at equator, DEC=0
        # There's still a fraction of an arcsec offset for 'center' and 'origin'....
    else:
        crophdr['CRPIX1'] = hdrin['CRPIX1'] - int(xbounds[0]); crophdr['CRPIX2'] = hdrin['CRPIX2'] - int(ybounds[0])
    crophdr['NAXIS'] = 2
    if type(savenew) == str:
        try: pyfits.writeto(savenew, cropdata, crophdr, overwrite=overwrite)
        except Exception: print('Could not save file to path supplied.')
    return cropdata, crophdr


def crop_cube(datain, hdrin, xbounds, ybounds, zbounds=[0, None], newref=None, savenew=False, overwrite=False):
    """
    Function to crop a 3D fits cube to the specified pixel bounds.

    zbounds (optional) are channel number boundaries [zmin,zmax], 0-indexed

    Other arguments the same as for crop_image
    """
    if len(datain.shape) > 3:
        if (datain.shape[i] == 1 for i in range(len(datain.shape) - 1)): datain = datain[0, :, :, :]
        else: raise Exception('Data array has 4 or more dimensions - reduce them to 3!')
    elif len(datain.shape) < 3:
        raise Exception('File does not have 3 dimensions!  Use crop_image.')
    else:
        pass
    if zbounds[1] is None: zbounds[1] = datain.shape[-3] + 1
    cropdata = datain[int(zbounds[0]):int(zbounds[1]), int(np.round(ybounds[0])):int(np.round(ybounds[1])) + 1, int(np.round(xbounds[0])):int(np.round(xbounds[1])) + 1]
    crophdr = hdrin.copy()
    crophdr['NAXIS3'], crophdr['NAXIS2'], crophdr['NAXIS1'] = cropdata.shape; crophdr['CRPIX3'] -= int(zbounds[0])
    if newref == 'center':
        crophdr['CRPIX1'] = (int(xbounds[1] - xbounds[0])) * .5; crophdr['CRPIX2'] = (int(ybounds[1] - ybounds[0])) * .5
        crophdr['CRVAL1'], crophdr['CRVAL2'] = pixel_to_sky(hdrin, (int(xbounds[0]) + int(xbounds[1] - xbounds[0]) * .5) - 1, (int(ybounds[0]) + int(ybounds[1] - ybounds[0]) * .5) - 1)
        radiff = np.unwrap(np.deg2rad([0, hdrin['CRVAL1'] - crophdr['CRVAL1']]))[1] * 180. / np.pi
        crophdr['CROTA2'] = radiff * np.sin(crophdr['CRVAL2'] * np.pi / 180)  # sin because at equator, DEC=0
    elif newref == 'origin':
        crophdr['CRPIX1'] = 0.; crophdr['CRPIX2'] = 0.
        crophdr['CRVAL1'], crophdr['CRVAL2'] = pixel_to_sky(hdrin, int(xbounds[0]) - 1, int(ybounds[0]) - 1)
        radiff = np.unwrap(np.deg2rad([0, hdrin['CRVAL1'] - crophdr['CRVAL1']]))[1] * 180. / np.pi
        crophdr['CROTA2'] = radiff * np.sin(crophdr['CRVAL2'] * np.pi / 180)  # sin because at equator, DEC=0
    else:
        crophdr['CRPIX1'] = hdrin['CRPIX1'] - int(xbounds[0]); crophdr['CRPIX2'] = hdrin['CRPIX2'] - int(ybounds[0])
    crophdr['NAXIS'] = 3
    if type(savenew) == str:
        try: pyfits.writeto(savenew, cropdata, crophdr, overwrite=overwrite)
        except Exception: print('Could not save file to path supplied.')
    return cropdata, crophdr


def crop_image_sky(datain, hdrin, centerRADEC, radius_asec, newref=None, savenew=False, overwrite=False,
                      return_cropcenterpix=False, coords_in='dec', checksys=False, incoordsys='fk5',
                      incoordequinox='J2000.0', forceimagesys=None):
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
        Passes to sky_to_pixel().  When True, checks that the input coordinate frame is the same as the header frame (e.g., fk5, icrs, etc.)
    incoordsys : str
        Passes to sky_to_pixel().  SkyCoord frame system, default = 'fk5'
    incoordequinox : str
        Passes to sky_to_pixel().  equinox, default = 'J2000.0'
    forceimagesys : str
        Passes to sky_to_pixel().  SkyCoord frame system. User can specify a frame forceimagesys to use (in case no RADESYS in the header, or if it's known to be wrong...)

    Returns
    -------
    array, astropy.io.fits.header
        cropdata, crophdr
    """
    if 'sex' in coords_in.lower(): centerRADEC = sex2dec(*centerRADEC)
    centerpix = sky_to_pixel(hdrin, centerRADEC[0], centerRADEC[1], precise=True, checksys=checksys,
                            incoordsys=incoordsys, incoordequinox=incoordequinox, forceimagesys=forceimagesys)
    pixextent = (radius_asec / 3600) / get_cdelts(hdrin)[1]  # Total crop image width is center+/- radius_asec
    cropdat, crophdr = crop_image(datain, hdrin, [centerpix[0] - pixextent, centerpix[0] + pixextent],
                                  [centerpix[1] - pixextent, centerpix[1] + pixextent],
                                  newref=newref, savenew=savenew, overwrite=overwrite)
    if return_cropcenterpix is True:
        centerpixcrop = sky_to_pixel(crophdr, centerRADEC[0], centerRADEC[1], precise=True)
        return cropdat, crophdr, centerpixcrop
    else:
        return cropdat, crophdr


def crop_cube_sky(datain, hdrin, centerRADEC, radius_asec, zbounds=[0, None], newref=None, savenew=False,
                      overwrite=False, coords_in='dec', return_cropcenterpix=False):
    """
    Function to crop a 3D fits cube based on sky coordinates and specified width.

    zbounds (optional) are channel number boundaries [zmin,zmax], 0-indexed

    Other arguments the same as for crop_image_sky()
    """
    if 'sex' in coords_in.lower(): centerRADEC = sex2dec(*centerRADEC)
    if zbounds[1] is None: zbounds[1] = datain.shape[0] + 1
    centerpix = sky_to_pixel(hdrin, centerRADEC[0], centerRADEC[1], precise=True)
    pixextent = (radius_asec / 3600) / get_cdelts(hdrin)[1]  # Total crop image width is center+/- radius_asec
    cropdat, crophdr = crop_cube(datain, hdrin, [centerpix[0] - pixextent, centerpix[0] + pixextent],
                                  [centerpix[1] - pixextent, centerpix[1] + pixextent], zbounds=zbounds,
                                  newref=newref, savenew=savenew, overwrite=overwrite)
    if return_cropcenterpix is True:
        centerpixcrop = sky_to_pixel(crophdr, centerRADEC[0], centerRADEC[1], precise=True)
        return cropdat, crophdr, centerpixcrop
    else:
        return cropdat, crophdr
