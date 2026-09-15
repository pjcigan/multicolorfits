"""
Helpers for aligning a stack of FITS images onto a common pixel grid before
colorizing and combining.

Typical workflow when layers do not yet share a grid::

    aligned = mcf.align_stack([(data1, hdr1), (data2, hdr2)], reference=0)
    # aligned[i] is (data_reproj, common_header) ready for to_grey_rgb

For Galactic display, pass ``frame='galactic'`` so the master header is built
in that frame first (see :func:`~multicolorfits.core.skyframes.convert_header_frame`).
"""

import numpy as np

from .reproject_tools import reproject_image
from .skyframes import convert_header_frame, optimal_common_header

__all__ = [
    'reproject_stack_to_header',
    'reproject_stack_to_reference',
    'align_stack',
    'downsample_for_preview',
    'prep_layers',
]


def reproject_stack_to_header(images, hdrto, method='interp', scale=False,
                              return_footprints=False, order=1):
    """
    Reproject every (data, header) pair onto a single target header.

    Parameters
    ----------
    images : list of (array, header) tuples
    hdrto : astropy.io.fits.header
        Target WCS grid (the 'master' header all layers will share).
    method, scale : passed to reproject_image
    return_footprints : bool
        If True, return list of footprint arrays alongside data.

    Returns
    -------
    list of array (, list of footprints)
        One reprojected data array per input image, all on hdrto's grid.
    """
    out = []
    footprints = []
    for data, hdr in images:
        result = reproject_image(data, hdr, hdrto, scale=scale, method=method,
                                 order=order, returnfootprint=return_footprints)
        if return_footprints:
            out.append(result[0])
            footprints.append(result[1])
        else:
            out.append(result)
    if return_footprints:
        return out, footprints
    return out


def reproject_stack_to_reference(images, reference=0, method='interp', scale=False,
                                 frame=None, projection=None, fit_footprint=True,
                                 return_footprints=False, order=1):
    """
    Reproject all images onto the grid of a reference layer.

    Parameters
    ----------
    images : list of (array, header) tuples
    reference : int
        Index of the layer whose header defines the target grid (default 0).
    frame : str or None
        If set (e.g. 'galactic'), the reference header is first converted to
        that celestial frame before reprojection (see convert_header_frame).
    projection, fit_footprint : passed to convert_header_frame when frame set
    method, scale : passed to reproject_image
    return_footprints : bool

    Returns
    -------
    list of array, header (, footprints)
        Reprojected data arrays and the common target header.
    """
    if not images:
        raise ValueError('images list is empty')
    if reference < 0 or reference >= len(images):
        raise ValueError('reference index out of range')
    _, hdr_ref = images[reference]
    if frame is not None:
        hdr_ref = convert_header_frame(hdr_ref, frame=frame, projection=projection,
                                       fit_footprint=fit_footprint)
    result = reproject_stack_to_header(images, hdr_ref, method=method, scale=scale,
                                       return_footprints=return_footprints, order=order)
    if return_footprints:
        data_list, footprints = result
        return data_list, hdr_ref, footprints
    return result, hdr_ref


def align_stack(images, reference=None, optimal=False, frame=None, **kwargs):
    """
    High-level alignment: pick a master header, reproject every layer onto it,
    return ``[(data, hdr), ...]`` ready for the color pipeline.

    Exactly one of ``reference`` or ``optimal`` must be used to choose the
    master grid:

    * ``reference=N`` -- use image N's header (optionally converted to
      ``frame`` first).
    * ``optimal=True`` -- compute the smallest common header with
      ``optimal_common_header()`` (requires reproject).

    Parameters
    ----------
    images : list of (array, header)
    reference : int or None
        Index of reference layer.  Mutually exclusive with ``optimal=True``.
    optimal : bool
        Use reproject.mosaicking.find_optimal_celestial_wcs.
    frame : str or None
        Target celestial frame for the master header ('galactic', 'icrs', ...).
    **kwargs
        Passed through to reproject_image (method, scale, ...).

    Returns
    -------
    list of (array, header)

    Examples
    --------
    ::

        import multicolorfits as mcf
        # Needs pip install "multicolorfits[reproject]"
        aligned = mcf.align_stack(
            [(data_a, hdr_a), (data_b, hdr_b)], reference=0)
        # aligned = mcf.align_stack(..., optimal=True, frame='galactic')
    """
    if optimal and reference is not None:
        raise ValueError('Specify either reference= or optimal=True, not both')
    if optimal:
        wcs_kwargs = {k: kwargs.pop(k) for k in list(kwargs) if k in ('projection', 'resolution', 'auto_rotate')}
        reproj_kwargs = {k: kwargs.pop(k) for k in list(kwargs) if k in ('method', 'scale', 'return_footprints', 'order')}
        hdr_master = optimal_common_header(images, frame=frame, **wcs_kwargs)
        data_list = reproject_stack_to_header(images, hdr_master, **reproj_kwargs)
    else:
        ref = 0 if reference is None else reference
        data_list, hdr_master = reproject_stack_to_reference(
            images, reference=ref, frame=frame, **kwargs)
    return list(zip(data_list, [hdr_master] * len(data_list)))


def downsample_for_preview(arr, max_size=512, order=1):
    """
    Downsample a 2D or RGB image so its longest side is at most ``max_size``.

    Intended for interactive GUI / notebook previews only.  Final exports and
    publications should use the full-resolution arrays.

    Uses ``skimage.transform.resize`` with anti-aliasing when shrinking;
    returns the input unchanged if it already fits.

    Parameters
    ----------
    arr : array
        2D ``(ny, nx)`` or RGB ``(ny, nx, 3)`` array.
    max_size : int
        Maximum pixels along the longest axis.
    order : int
        Spline order for resize (``1`` = bilinear, default).

    Returns
    -------
    array
        Downsampled array, same dimensionality as the input.
    """
    arr = np.asarray(arr)
    if arr.ndim not in (2, 3):
        raise ValueError('arr must be 2D or RGB (ny, nx, 3)')
    ny, nx = arr.shape[:2]
    longest = max(ny, nx)
    if longest <= max_size:
        return arr
    scale = max_size / float(longest)
    new_shape = (max(1, int(round(ny * scale))), max(1, int(round(nx * scale))))
    from skimage.transform import resize
    if arr.ndim == 2:
        return resize(arr, new_shape, order=order, mode='reflect',
                      anti_aliasing=True, preserve_range=True)
    return resize(arr, new_shape + (3,), order=order, mode='reflect',
                  anti_aliasing=True, preserve_range=True)


def prep_layers(images, north_up=True, rotation_deg=0.0, oversample=1.0,
                crop='overlap', frame='auto', method='interp', order=1,
                scale=False, blank_zeros=False, fit_footprint=True, pad=0):
    """
    Tidy headers, optionally rotate to north-up (or an arbitrary angle),
    reproject onto one grid, and crop to the overlapping finite pixels.

    Parameters
    ----------
    images : list of (array, header)
    north_up : bool
        Build the target grid with make_rotated_header. False keeps the
        reference (first) grid after tidying, unless ``rotation_deg`` or
        ``oversample`` is set.
    rotation_deg : float
        Extra rotation from north-up in the image frame (0 = north up).
    oversample : float
        Shrink the target pixel scale by this factor (2 = half the pixel size).
    crop : {'overlap', 'none'}
        Crop to pixels finite in every layer after reprojection.
    frame : str
        'auto' keeps each layer's frame on a target built from the first
        image. A named frame ('galactic', ...) changes the target frame.
    blank_zeros : bool
        Treat exact zeros as missing before reprojection.
    pad : int
        Passed to crop_to_overlap.

    Returns
    -------
    list of (array, header)
        Every header is the common target (after the optional crop).

    Examples
    --------
    ::

        import multicolorfits as mcf
        prepared = mcf.prep_layers(
            [(data_a, hdr_a), (data_b, hdr_b)],
            north_up=True, oversample=2, crop='overlap', order=1)
    """
    from .skyframes import make_rotated_header
    from .wcs_tools import blank_missing, crop_to_overlap, tidy_header

    if not images:
        raise ValueError('prep_layers: images list is empty')
    cleaned = []
    for data, hdr in images:
        cleaned.append((blank_missing(data, zeros=blank_zeros), tidy_header(hdr)))

    rotate = bool(north_up) or abs(float(rotation_deg or 0.0)) > 0 or float(oversample or 1) != 1.0
    if rotate:
        hdr_tgt = make_rotated_header(
            cleaned[0][1], rotation_deg=rotation_deg, oversample=oversample,
            fit_footprint=fit_footprint, frame=frame)
        arrays = reproject_stack_to_header(
            cleaned, hdr_tgt, method=method, scale=scale, order=order)
        common = hdr_tgt
    else:
        frame_arg = None if frame in (None, 'auto') else frame
        arrays, common = reproject_stack_to_reference(
            cleaned, reference=0, method=method, scale=scale, order=order,
            frame=frame_arg, fit_footprint=fit_footprint)
    if str(crop).lower() in ('overlap', 'auto', 'true', '1'):
        arrays, common = crop_to_overlap(arrays, common, pad=pad)
    return [(arr, common) for arr in arrays]
