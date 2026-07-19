"""
High-level end-to-end pipeline helpers for scripts.

These wrap core stretch → colorize → combine steps without requiring a GUI::

    rgb = mcf.combine_layers([d1, d2], colors=['#C11B17', '#4CC417'],
                             colorspace='lab')
    mcf.save_combined(rgb, header, 'out.fits')
"""

import numpy as np
import astropy.io.fits as pyfits

from .core.colorize import to_grey_rgb, colorize_image, combine_multicolor, smooth_image
from .core.colormix import combine_multicolor_colorspace, COLORSPACES, BLENDS
from .core.io_output import save_rgb_fits, annotate_provenance_header
from .core.stack_tools import align_stack

__all__ = [
    'combine_layers',
    'combine_from_files',
    'save_combined',
]


def combine_layers(data_list, colors, stretches=None, min_maxes=None,
                   colorspace='rgb', blend='screen', weights=None,
                   gamma=2.2, inverse=False, smooth_sigmas=None,
                   preview=False, dtype=None):
    """
    Stretch, colorize, and combine a list of 2D arrays into one RGB image.

    Parameters
    ----------
    data_list : list of array
        Greyscale 2D images (same shape if not pre-aligned).
    colors : list of str
        Hex colors, one per layer (e.g. ``['#C11B17', '#4CC417']``).
    stretches : list of str or None
        Per-layer stretch names; ``None`` uses ``'linear'`` for all.
    min_maxes : list of [vmin, vmax] or None
        Per-layer scaling limits; ``None`` uses each array's full range.
    colorspace : str
        ``'rgb'`` (classic channel sum) or ``'lab'`` / ``'hsv'`` / ``'hsl'``.
    blend : str
        Lightness blend for Lab/HSV/HSL: ``'screen'``, ``'sum'``, ``'max'``,
        ``'mean'``.
    weights : list or None
        Optional per-layer chroma weights (Lab/HSV/HSL).
    gamma : float
        Gamma for to_grey_rgb / colorize (default 2.2).
    inverse : bool
        Legacy white-background mode via ``1 - image`` (prefer a white
        ``combine_background`` in session / alpha APIs for new work).
    smooth_sigmas : list of float or None
        Optional Gaussian sigma per colorized layer.
    preview : bool
        If True, run the compose in float32 (lower memory; slightly reduced
        precision).  Re-render with ``preview=False`` for final float64 export.
        Ignored when ``dtype`` is set.
    dtype : numpy dtype or None
        Explicit working dtype; overrides ``preview``.

    Returns
    -------
    array
        Combined display-ready RGB image in [0, 1].

    Examples
    --------
    >>> rgb = combine_layers(  # doctest: +SKIP
    ...     [d1, d2], colors=['#C11B17', '#4CC417'], colorspace='lab')
    """
    n = len(data_list)
    if len(colors) != n:
        raise ValueError('Need one color per data layer')
    if stretches is None:
        stretches = ['linear'] * n
    if min_maxes is None:
        min_maxes = [[None, None]] * n
    if smooth_sigmas is None:
        smooth_sigmas = [None] * n

    if dtype is None:
        dtype = np.float32 if preview else np.float64

    colorized = []
    for i, (data, color, stretch, mm, sigma) in enumerate(
            zip(data_list, colors, stretches, min_maxes, smooth_sigmas)):
        grey = to_grey_rgb(data, rescalefn=stretch, min_max=[mm[0], mm[1]],
                                gamma=gamma, dtype=dtype)
        col = colorize_image(grey, color, colorintype='hex', gammacorr_color=gamma,
                             dtype=dtype)
        if sigma is not None and sigma > 0:
            col = smooth_image(col, sigma=sigma)
        colorized.append(col)

    if colorspace == 'rgb':
        return combine_multicolor(colorized, gamma=gamma, inverse=inverse, dtype=dtype)
    if colorspace not in COLORSPACES:
        raise ValueError("colorspace must be 'rgb' or one of %s" % (COLORSPACES,))
    if blend not in BLENDS:
        raise ValueError("blend must be one of %s" % (BLENDS,))
    return combine_multicolor_colorspace(
        colorized, colorspace=colorspace, blend=blend, weights=weights,
        gamma=gamma, inverse=inverse, dtype=dtype)


def combine_from_files(fits_paths, colors, align_reference=0, align_frame=None,
                       align_optimal=False, **combine_kwargs):
    """
    Load FITS files, optionally align them to a common grid, then combine.

    Parameters
    ----------
    fits_paths : list of str
        Paths to 2D FITS files.
    colors : list of str
        Hex colors, one per file.
    align_reference : int
        When aligning, index of the reference layer (default 0).
    align_frame : str or None
        Celestial frame for the master grid (``'galactic'``, …).  ``None``
        keeps the reference header's native frame.
    align_optimal : bool
        If True, build an optimal common header instead of using one reference.
    **combine_kwargs
        Forwarded to :func:`combine_layers` (stretches, colorspace, gamma, …).

    Returns
    -------
    array, header
        Combined RGB image and the common WCS header.
    """
    images = []
    for path in fits_paths:
        data, hdr = pyfits.getdata(path, header=True)
        images.append((np.asarray(data, dtype=float), hdr))
    if len(images) > 1:
        aligned = align_stack(
            images,
            reference=None if align_optimal else align_reference,
            optimal=align_optimal,
            frame=align_frame,
        )
        data_list = [a[0] for a in aligned]
        hdr = aligned[0][1]
    else:
        data_list = [images[0][0]]
        hdr = images[0][1]
    combined = combine_layers(data_list, colors, **combine_kwargs)
    return combined, hdr


def save_combined(savepath, combined, header, session_info=None, overwrite=True,
                  annotate=True):
    """
    Save a combined RGB image to FITS, optionally annotating the header with
    multicolorfits provenance (version, timestamp, layer colors/scales).

    Parameters
    ----------
    savepath : str
    combined : array
        RGB image from combine_layers / render_combined.
    header : astropy.io.fits.header
    session_info : dict or McfSession or None
        If a McfSession, panel/compose settings are recorded in HISTORY cards.
        If a dict, used as-is for annotate_provenance_header.
    overwrite : bool
    annotate : bool
        Add HISTORY/COMMENT cards (default True).
    """
    hdr = header.copy() if header is not None else pyfits.Header()
    if annotate:
        if session_info is not None and hasattr(session_info, 'to_dict'):
            session_info = session_info.to_dict()
        annotate_provenance_header(hdr, session_info or {})
    save_rgb_fits(savepath, combined, hdr, overwrite=overwrite)
