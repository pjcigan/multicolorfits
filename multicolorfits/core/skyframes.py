"""
Celestial frame conversion helpers (for example equatorial → Galactic).

Build a target header in another celestial frame from an existing FITS header
and reproject data onto it so a survey image can be shown with Galactic north
up (or another frame).

Requires the optional ``reproject`` extra for the reprojecting entry points.
"""

import numpy as np
import astropy.io.fits as pyfits
import astropy.wcs as pywcs
from astropy.coordinates import SkyCoord

from .wcs_tools import get_cdelts, tidy_header
from .reproject_tools import reproject_image, _require

__all__ = [
    'header_frame_name',
    'convert_header_frame',
    'reproject_to_frame',
    'reproject_to_galactic',
    'optimal_common_header',
    'make_rotated_header',
    'make_north_up_header',
    'reproject_to_rotation',
    'reproject_north_up',
]

# frame name -> (CTYPE1 4-char prefix, CTYPE2 4-char prefix)
_FRAME_CTYPES = {
    'galactic': ('GLON', 'GLAT'),
    'icrs': ('RA--', 'DEC-'),
    'fk5': ('RA--', 'DEC-'),
    'fk4': ('RA--', 'DEC-'),
    'ecliptic': ('ELON', 'ELAT'),
}
# frame name -> astropy SkyCoord frame string
_ASTROPY_FRAMES = {
    'galactic': 'galactic',
    'icrs': 'icrs',
    'fk5': 'fk5',
    'fk4': 'fk4',
    'ecliptic': 'barycentricmeanecliptic',
}


def header_frame_name(hdrin, default='fk5'):
    """
    Determine the celestial frame of a FITS header, as an astropy frame name.

    Uses CTYPE prefixes (GLON/ELON) first, then RADESYS/RADECSYS for
    equatorial frames.

    Parameters
    ----------
    hdrin : astropy.io.fits.header
        Header object
    default : str
        Frame to assume for equatorial headers without a RADESYS card.

    Returns
    -------
    str
        e.g. 'fk5', 'icrs', 'galactic', 'barycentricmeanecliptic'
    """
    ctype1 = str(hdrin.get('CTYPE1', ''))
    if ctype1.startswith('GLON'):
        return 'galactic'
    if ctype1.startswith('ELON'):
        return _ASTROPY_FRAMES['ecliptic']
    try:
        try:
            radesys = hdrin['RADESYS']
        except KeyError:
            radesys = hdrin['RADECSYS']
        return str(radesys).strip().lower()
    except KeyError:
        return default


def _projection_from_header(hdrin):
    """Extract the 3-letter projection code (e.g. 'TAN') from CTYPE1."""
    ctype1 = str(hdrin.get('CTYPE1', ''))
    if '-' in ctype1:
        proj = ctype1.split('-')[-1].strip()
        if len(proj) == 3:
            return proj
    return 'TAN'


def convert_header_frame(hdrin, frame='galactic', projection=None, fit_footprint=True):
    """
    Build a new 2D header in a different celestial frame, covering the same
    sky area as the input header with the same pixel scale, with the new
    frame's north axis up (no rotation).

    Use with reproject_to_frame() (or mcf.reproject_image) to actually resample
    image data onto the new grid.

    Parameters
    ----------
    hdrin : astropy.io.fits.header
        Input header (must contain a valid celestial WCS and NAXIS1/2).
    frame : str
        Target frame: 'galactic', 'icrs', 'fk5', 'fk4', or 'ecliptic'.
    projection : str or None
        3-letter FITS projection code for the output (e.g. 'TAN', 'SIN').
        None (default) keeps the projection of the input header.
    fit_footprint : bool
        True (default) sizes the output so the entire (rotated) input
        footprint fits.  False keeps the input NAXIS1/NAXIS2 dimensions.

    Returns
    -------
    astropy.io.fits.header
        New header in the target frame.
    """
    frame = frame.lower()
    if frame not in _FRAME_CTYPES:
        raise ValueError('frame must be one of %s' % (sorted(_FRAME_CTYPES),))
    target_frame = _ASTROPY_FRAMES[frame]

    nx = int(hdrin['NAXIS1'])
    ny = int(hdrin['NAXIS2'])
    if projection is None:
        projection = _projection_from_header(hdrin)

    wcs_in = pywcs.WCS(hdrin).celestial
    src_frame = header_frame_name(hdrin)

    # Center of the image in the source frame, converted to the target frame
    center_world = wcs_in.wcs_pix2world([[(nx - 1) / 2., (ny - 1) / 2.]], 0)[0]
    c_center = SkyCoord(center_world[0], center_world[1], unit='deg', frame=src_frame)
    c_center_t = c_center.transform_to(target_frame)
    lon0 = float(c_center_t.spherical.lon.deg)
    lat0 = float(c_center_t.spherical.lat.deg)

    cdelt1, cdelt2 = get_cdelts(hdrin)

    if fit_footprint:
        # Transform corner + edge-midpoint pixels and find the on-sky extent
        # around the new center, so the whole rotated footprint fits.
        xs = np.array([0, nx - 1, 0, nx - 1, (nx - 1) / 2., (nx - 1) / 2., 0, nx - 1])
        ys = np.array([0, 0, ny - 1, ny - 1, 0, ny - 1, (ny - 1) / 2., (ny - 1) / 2.])
        world = wcs_in.wcs_pix2world(np.column_stack([xs, ys]), 0)
        c_pts = SkyCoord(world[:, 0], world[:, 1], unit='deg', frame=src_frame).transform_to(target_frame)
        dlon, dlat = c_center_t.spherical_offsets_to(c_pts)
        half_lon = float(np.max(np.abs(dlon.deg)))
        half_lat = float(np.max(np.abs(dlat.deg)))
        nx_new = 2 * int(np.ceil(half_lon / abs(cdelt1))) + 1
        ny_new = 2 * int(np.ceil(half_lat / abs(cdelt2))) + 1
    else:
        nx_new, ny_new = nx, ny

    ctype1_prefix, ctype2_prefix = _FRAME_CTYPES[frame]
    wcs_out = pywcs.WCS(naxis=2)
    wcs_out.wcs.crpix = [(nx_new + 1) / 2., (ny_new + 1) / 2.]  # FITS 1-indexed center
    wcs_out.wcs.crval = [lon0, lat0]
    wcs_out.wcs.cdelt = [-abs(cdelt1), abs(cdelt2)]  # lon increases to the left
    wcs_out.wcs.ctype = ['%s-%s' % (ctype1_prefix, projection),
                         '%s-%s' % (ctype2_prefix, projection)]
    if frame in ('icrs', 'fk5', 'fk4'):
        wcs_out.wcs.radesys = frame.upper()
        if frame == 'fk5':
            wcs_out.wcs.equinox = 2000.0
        elif frame == 'fk4':
            wcs_out.wcs.equinox = 1950.0

    hdrout = wcs_out.to_header()
    hdrout['NAXIS'] = 2
    hdrout['NAXIS1'] = nx_new
    hdrout['NAXIS2'] = ny_new
    # Carry over non-WCS metadata that remains valid after reprojection
    for card in ['BUNIT', 'OBJECT', 'TELESCOP', 'BMAJ', 'BMIN', 'BPA']:
        try:
            hdrout[card] = hdrin[card]
        except KeyError:
            pass
    return hdrout


def reproject_to_frame(datain, hdrin, frame='galactic', method='interp', scale=False,
                       projection=None, fit_footprint=True, returnfootprint=False):
    """
    Reproject a 2D image into a different celestial frame (e.g. Galactic),
    with the new frame's north up.

    Requires the optional 'reproject' package (pip install multicolorfits[reproject]).

    Parameters
    ----------
    datain : array
        Input 2D image array.
    hdrin : astropy.io.fits.header
        Input header with valid WCS.
    frame : str
        Target frame: 'galactic', 'icrs', 'fk5', 'fk4', or 'ecliptic'.
    method : str
        Passed to reproject_image: 'interp' (default), 'spi', or 'kapteyn'.
    scale : bool
        Passed to reproject_image: True for flux units (Jy/pix etc.), False for
        brightness units (Jy/sr etc.).
    projection : str or None
        Output projection code; None keeps the input projection.
    fit_footprint : bool
        Size the output canvas to fit the full rotated footprint (default).
    returnfootprint : bool
        Also return the reprojection footprint array.

    Returns
    -------
    array, astropy.io.fits.header (, array)
        Reprojected data, the new header (, optional footprint)
    """
    hdrto = convert_header_frame(hdrin, frame=frame, projection=projection,
                                 fit_footprint=fit_footprint)
    result = reproject_image(datain, hdrin, hdrto, scale=scale, method=method,
                         returnfootprint=returnfootprint)
    if returnfootprint:
        return result[0], hdrto, result[1]
    return result, hdrto


def reproject_to_galactic(datain, hdrin, **kwargs):
    """
    Convenience alias for reproject_to_frame(..., frame='galactic').
    Returns (reprojected_data, new_header).
    """
    return reproject_to_frame(datain, hdrin, frame='galactic', **kwargs)


def optimal_common_header(images, frame=None, projection='TAN', **kwargs):
    """
    Compute an optimal common header covering a set of images, using
    reproject.mosaicking.find_optimal_celestial_wcs.  Useful as the 'master
    header' that all layers get reprojected onto before combining.

    Requires the optional 'reproject' package.

    Parameters
    ----------
    images : list
        List of (data, header) tuples (or of astropy HDU objects).
    frame : str or None
        Target frame name ('galactic', 'icrs', 'fk5', ...).  None keeps the
        frame chosen automatically by reproject.
    projection : str
        Output projection code (default 'TAN').
    **kwargs
        Passed through to find_optimal_celestial_wcs (e.g. resolution=...,
        auto_rotate=True).

    Returns
    -------
    astropy.io.fits.header
        Header (with NAXIS1/NAXIS2 set) covering all the input images.
    """
    _require('reproject', "pip install reproject  (or pip install multicolorfits[reproject])")
    from reproject.mosaicking import find_optimal_celestial_wcs

    if frame is not None:
        from astropy.coordinates import frame_transform_graph
        frame_cls = frame_transform_graph.lookup_name(_ASTROPY_FRAMES.get(frame.lower(), frame.lower()))
        if frame_cls is None:
            raise ValueError('Unknown frame: %s' % frame)
        kwargs['frame'] = frame_cls()

    wcs_out, shape_out = find_optimal_celestial_wcs(images, projection=projection, **kwargs)
    hdrout = wcs_out.to_header()
    hdrout['NAXIS'] = 2
    hdrout['NAXIS1'] = int(shape_out[1])
    hdrout['NAXIS2'] = int(shape_out[0])
    return hdrout


def _frame_key(hdrin):
    """Short frame name used by convert_header_frame ('galactic', 'fk5', ...)."""
    name = header_frame_name(hdrin)
    if name in ('galactic', 'icrs', 'fk5', 'fk4'):
        return name
    if 'ecliptic' in name:
        return 'ecliptic'
    return 'fk5'


def _copy_nonwcs(hdrout, hdrin):
    for card in ['BUNIT', 'OBJECT', 'TELESCOP', 'BMAJ', 'BMIN', 'BPA']:
        if card in hdrin and card not in hdrout:
            hdrout[card] = hdrin[card]
    return hdrout


def make_rotated_header(hdrin, rotation_deg=0.0, oversample=1.0, fit_footprint=True,
                        frame='auto', allow_flip=False, warn_flip=True):
    """
    Build a 2D header with a chosen on-sky rotation, in the image's own frame.

    ``rotation_deg=0`` (the default) is north-up in that frame: Galactic
    headers stay Galactic (+b up), equatorial stay equatorial. This does not
    silently convert to ICRS. Pass ``frame='galactic'`` (etc.) only when you
    also want a frame change — that path composes with convert_header_frame.

    The output is east-left (negative longitude scale) unless ``allow_flip``
    is set and the input WCS is mirrored. Pixel scale is
    ``|CDELT| / oversample``.

    Requires no optional packages (the header only). Resample with
    reproject_to_rotation / reproject_image.
    """
    hdrin = tidy_header(hdrin, warn_flip=warn_flip)
    src_key = _frame_key(hdrin)
    want = None if frame in (None, 'auto') else str(frame).lower()
    if want is not None and want != src_key:
        hdrin = convert_header_frame(hdrin, frame=want, fit_footprint=fit_footprint)
        src_key = want

    nx = int(hdrin['NAXIS1'])
    ny = int(hdrin['NAXIS2'])
    oversample = float(oversample or 1.0)
    if oversample <= 0:
        raise ValueError('oversample must be > 0')
    cdelt1, cdelt2 = get_cdelts(hdrin)
    s1 = abs(float(cdelt1)) / oversample
    s2 = abs(float(cdelt2)) / oversample
    flipped = False
    try:
        from .wcs_tools import wcs_is_flipped
        flipped = bool(wcs_is_flipped(hdrin))
    except Exception:
        flipped = False
    lon_sign = 1.0 if (allow_flip and flipped) else -1.0

    wcs_in = pywcs.WCS(hdrin).celestial
    src_frame = header_frame_name(hdrin)
    cx = (nx - 1) / 2.0
    cy = (ny - 1) / 2.0
    center_world = wcs_in.wcs_pix2world([[cx, cy]], 0)[0]
    c_center = SkyCoord(center_world[0], center_world[1], unit='deg', frame=src_frame)
    lon0 = float(c_center.spherical.lon.deg)
    lat0 = float(c_center.spherical.lat.deg)
    th = np.radians(float(rotation_deg or 0.0))

    if fit_footprint:
        xs = np.array([0, nx - 1, 0, nx - 1, cx, cx, 0, nx - 1])
        ys = np.array([0, 0, ny - 1, ny - 1, 0, ny - 1, cy, cy])
        world = wcs_in.wcs_pix2world(np.column_stack([xs, ys]), 0)
        c_pts = SkyCoord(world[:, 0], world[:, 1], unit='deg', frame=src_frame)
        dlon, dlat = c_center.spherical_offsets_to(c_pts)
        east = np.asarray(dlon.deg, dtype=float)
        north = np.asarray(dlat.deg, dtype=float)
        # Rotate offsets so +y is north rotated by rotation_deg (CCW, CROTA2 sense).
        east_r = east * np.cos(th) + north * np.sin(th)
        north_r = -east * np.sin(th) + north * np.cos(th)
        half_x = float(np.max(np.abs(east_r)))
        half_y = float(np.max(np.abs(north_r)))
        nx_new = max(nx, 2 * int(np.ceil(half_x / s1)) + 1)
        ny_new = max(ny, 2 * int(np.ceil(half_y / s2)) + 1)
    else:
        nx_new, ny_new = nx, ny

    projection = _projection_from_header(hdrin)
    ctype1_prefix, ctype2_prefix = _FRAME_CTYPES[src_key]
    wcs_out = pywcs.WCS(naxis=2)
    wcs_out.wcs.crpix = [(nx_new + 1) / 2.0, (ny_new + 1) / 2.0]
    wcs_out.wcs.crval = [lon0, lat0]
    wcs_out.wcs.cdelt = [lon_sign * s1, s2]
    wcs_out.wcs.ctype = ['%s-%s' % (ctype1_prefix, projection),
                         '%s-%s' % (ctype2_prefix, projection)]
    if abs(th) > 1e-12:
        # CROTA2-style PC: latitude axis rotated from +y by rotation_deg.
        c, s = np.cos(th), np.sin(th)
        wcs_out.wcs.pc = [[c, s], [-s, c]]
    if src_key in ('icrs', 'fk5', 'fk4'):
        wcs_out.wcs.radesys = src_key.upper()
        if src_key == 'fk5':
            wcs_out.wcs.equinox = 2000.0
        elif src_key == 'fk4':
            wcs_out.wcs.equinox = 1950.0

    hdrout = wcs_out.to_header()
    hdrout['NAXIS'] = 2
    hdrout['NAXIS1'] = int(nx_new)
    hdrout['NAXIS2'] = int(ny_new)
    return _copy_nonwcs(hdrout, hdrin)


def make_north_up_header(hdrin, **kwargs):
    """Alias for make_rotated_header(..., rotation_deg=0)."""
    kwargs = dict(kwargs)
    kwargs['rotation_deg'] = 0.0
    return make_rotated_header(hdrin, **kwargs)


def reproject_to_rotation(datain, hdrin, rotation_deg=0.0, oversample=1.0,
                          fit_footprint=True, frame='auto', allow_flip=False,
                          method='interp', order=1, scale=False,
                          returnfootprint=False):
    """
    Reproject an image onto a rotated (default: north-up) grid in its own frame.

    See make_rotated_header. ``order`` is the reproject_interp spline order;
    ``method='exact'`` (or ``order=0``) uses reproject_exact. Requires the
    optional reproject extra.

    Returns (data, header) or (data, header, footprint).
    """
    hdrto = make_rotated_header(
        hdrin, rotation_deg=rotation_deg, oversample=oversample,
        fit_footprint=fit_footprint, frame=frame, allow_flip=allow_flip)
    result = reproject_image(datain, hdrin, hdrto, scale=scale, method=method,
                             order=order, returnfootprint=returnfootprint)
    if returnfootprint:
        return result[0], hdrto, result[1]
    return result, hdrto


def reproject_north_up(datain, hdrin, **kwargs):
    """Alias for reproject_to_rotation(..., rotation_deg=0)."""
    kwargs = dict(kwargs)
    kwargs['rotation_deg'] = 0.0
    return reproject_to_rotation(datain, hdrin, **kwargs)
