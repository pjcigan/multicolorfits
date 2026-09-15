"""
Reprojection utilities.

The heavy optional dependencies (reproject, kapteyn) are imported lazily,
only when actually needed, so that the core package imports cleanly without
them.  A clear ImportError with install instructions is raised if the
requested method's package is missing.
"""

import numpy as np

from .scaling import draw_progress_bar
from .wcs_tools import sterad_per_pixel

__all__ = ['reproject_image', 'reproject_cube']


def _require(module_name, extra_hint):
    """Import an optional dependency or raise a helpful error."""
    import importlib
    try:
        return importlib.import_module(module_name)
    except ImportError as exc:
        raise ImportError(
            "multicolorfits: the '%s' package is required for this reprojection "
            "method. Install it with: %s" % (module_name, extra_hint)
        ) from exc


def reproject_image(mapin, hdrfrom, hdrto, scale=False, method='interp', order=1,
                interpdict={'order': 1, 'mode': 'constant', 'cval': np.nan},
                returnfootprint=False, parallel=True):
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
        One of 'interp' (reproject, spline order ``order``), 'exact' / 'spi'
        (reproject spherical-polygon intersection), or 'kapteyn'.
        ``order=0`` also selects the exact method.
    order : int
        Spline order for ``method='interp'`` (1 bilinear, 3 cubic). Ignored
        for exact / kapteyn.
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
    if scale is True:
        # --> Jy_2 = Jy_1_reproj*([sr/pix_2]/[sr/pix_1])
        # Set to True for [Jy/pix, W/m2/pix, Jy/beam if beam has been convolved]
        mapin_brightness = mapin.copy() / sterad_per_pixel(hdrfrom)
    else:
        mapin_brightness = mapin.copy()
    if method in ('exact', 'spi') or (method == 'interp' and int(order or 0) == 0):
        method = 'spi'
    if method == 'interp':
        reproject = _require('reproject', "pip install reproject  (or pip install multicolorfits[reproject])")
        map_reproj, map_footprint = reproject.reproject_interp(
            (mapin_brightness, hdrfrom), hdrto, order=int(order))
    elif method == 'spi':
        reproject = _require('reproject', "pip install reproject  (or pip install multicolorfits[reproject])")
        map_reproj, map_footprint = reproject.reproject_exact((mapin_brightness, hdrfrom), hdrto, parallel=parallel)
    else:
        kapteyn = _require('kapteyn', "see https://www.astro.rug.nl/software/kapteyn/ for install instructions")
        maputils = kapteyn.maputils
        hdrfrom = hdrfrom.copy(); hdrto = hdrto.copy()  # Otherwise it modifies original.  (Bad for reproject_cube calls)
        if hdrto['NAXIS'] > 2: hdrto['NAXIS'] = 2
        if hdrfrom['NAXIS'] > 2: hdrfrom['NAXIS'] = 2
        map_in = maputils.FITSimage(externaldata=mapin_brightness, externalheader=hdrfrom)
        map_reproj = map_in.reproject_to(hdrto, interpol_dict=interpdict).dat
    if scale is True:
        map_reproj *= sterad_per_pixel(hdrto)
    if returnfootprint is True and method != 'kapteyn':
        return map_reproj, map_footprint
    else:
        return map_reproj


def reproject_cube(mapin, hdrfrom, hdrto, scale=False, method='interp', parallel=True,
                returnfootprint=False, print_progress=False):
    """
    Reproject a 3D cube from one WCS header onto another, plane by plane.

    Calls :func:`reproject_image` for each spectral (or other non-spatial)
    slice.  For long cubes, set ``print_progress=True`` to show a stdout
    progress bar via :func:`~multicolorfits.draw_progress_bar` (zero
    dependencies).  Scripts that prefer `tqdm <https://tqdm.github.io/>`_
    can loop over planes themselves and call ``reproject_image`` — see
    that helper's docstring for an example.

    Parameters
    ----------
    mapin, hdrfrom, hdrto, scale, method, parallel, returnfootprint
        See :func:`reproject_image`.
    print_progress : bool, optional
        If True, update an ANSI progress bar on stdout after each plane
        (default False).

    Returns
    -------
    array or (array, array)
        Reprojected cube (, optional footprint cube when
        ``returnfootprint=True`` and ``method='spi'``).
    """
    # Scale: True if units are Flux [W/m^2 or similar].  False for brightness [W/m^2/sr or similar]
    tmpcube = np.zeros(mapin.shape[-3]).astype(np.ndarray); tmpcube_foot = tmpcube.copy()
    if returnfootprint is True and method == 'spi':
        for zz in range(mapin.shape[-3]):
            if print_progress is True:
                draw_progress_bar(float(zz) / mapin.shape[-3], prefix='  Reprojecting channels ', suffix='  %i of %i' % (zz + 1, mapin.shape[-3]))
            tmpcube[zz], tmpcube_foot[zz] = reproject_image(mapin[zz, :, :], hdrfrom, hdrto, scale=scale, method=method, returnfootprint=True)
        if print_progress is True:
            draw_progress_bar(1., prefix='  Reprojecting channels ', suffix='  %i of %i    \n' % (zz + 1, mapin.shape[-3]))
        return (np.array([tmpcube[zz] for zz in range(mapin.shape[-3])]),
                np.array([tmpcube_foot[zz] for zz in range(mapin.shape[-3])]))
    else:
        for zz in range(mapin.shape[-3]):
            if print_progress is True:
                draw_progress_bar(float(zz) / mapin.shape[-3], prefix='  Reprojecting channels ', suffix='  %i of %i' % (zz + 1, mapin.shape[-3]))
            tmpcube[zz] = reproject_image(mapin[zz, :, :], hdrfrom, hdrto, scale=scale, method=method, returnfootprint=False)
        if print_progress is True:
            draw_progress_bar(1., prefix='  Reprojecting channels ', suffix='  %i of %i    \n' % (zz + 1, mapin.shape[-3]))
        return np.array([tmpcube[zz] for zz in range(mapin.shape[-3])])
