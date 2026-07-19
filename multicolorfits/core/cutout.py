"""
Transparent cutouts and background deblending for multicolor composites.

Two sibling operations sharing one matte toolkit (percentile threshold ramp,
fade gamma, smoothing, geometric mattes, crop/resize):

* :func:`make_transparent_cutout` — *estimate* a matte for an RGB composite
  whose sky is part of the image (opaque object fading to transparent sky).
  Useful for slides, posters, README marks, and map stamps (e.g. the M74
  example in the skyplothelper markers tutorial).
* :func:`deblend_background` — *solve* for the alpha of a composite that was
  flattened onto a known solid background color (un-premultiply; exact
  round-trip with :func:`flatten_rgba`).

Typical use::

    rgba = mcf.make_transparent_cutout(combined_rgb, size=400, alpha_gamma=0.5)
    mcf.save_transparent_cutout(rgba, 'galaxy.png')

    # From a session (uses current compose settings, transparent background):
    rgba = session.export_transparent_cutout(size=256, savepath='stamp.png')

    # Recover RGBA from a JPEG/PNG flattened onto white:
    rgba = mcf.deblend_background(flat_rgb, background='white')
"""

from __future__ import annotations

import os

import numpy as np

__all__ = [
    'luminance',
    'alpha_from_intensity',
    'make_transparent_cutout',
    'save_transparent_cutout',
    'make_checkerboard',
    'preview_cutout_on_backgrounds',
    'batch_transparent_cutouts',
    'flatten_rgba',
    'deblend_background',
]


# --------------------------------------------------------------------------- helpers

def luminance(rgb):
    """
    Rec. 709 luma of an RGB (or greyscale) image in [0, 1].

    Parameters
    ----------
    rgb : array
        Shape ``(ny, nx)`` or ``(ny, nx, 3[+])``.

    Returns
    -------
    array
        2D luminance map in [0, 1].
    """
    arr = np.asarray(rgb, dtype=float)
    if arr.ndim == 2:
        return np.clip(arr, 0.0, 1.0)
    if arr.ndim != 3 or arr.shape[-1] < 3:
        raise ValueError('rgb must be (ny, nx) or (ny, nx, 3[+])')
    r, g, b = arr[..., 0], arr[..., 1], arr[..., 2]
    return np.clip(0.2126 * r + 0.7152 * g + 0.0722 * b, 0.0, 1.0)


def _finite_percentile(values, p):
    flat = np.asarray(values, dtype=float).ravel()
    flat = flat[np.isfinite(flat)]
    if flat.size == 0:
        return 0.0
    return float(np.percentile(flat, p))


def alpha_from_intensity(intensity, lo=55.0, hi=99.3, gamma=0.5, invert=False,
                         smooth=0.0):
    """
    Build a soft alpha matte from a 2D intensity map.

    Parameters
    ----------
    intensity : array (ny, nx)
        Brightness map in arbitrary units (typically luma of an RGB composite).
    lo, hi : float
        Percentiles of *intensity* mapped to alpha 0 and 1.  Defaults (~55 / 99.3)
        leave most sky transparent while keeping the bright object body opaque.
        Think of *hi* as the opacity **threshold**: everything at or above it
        is fully opaque; below it, alpha ramps down to 0 at the *lo* level.
    gamma : float
        Shape of the fade below the *hi* threshold (any value > 0).
        ``gamma < 1`` (e.g. 0.5) keeps more of the object body opaque —
        better over light slide backgrounds.  ``gamma = 1`` is a linear
        ramp.  ``gamma > 1`` (e.g. 1.4) tightens the matte, so only the
        brightest structure survives at high opacity.
    invert : bool
        If True, treat *dark* as signal (white-paper / emission-on-white look).
    smooth : float
        Gaussian sigma (pixels) applied to the intensity map *before* the
        percentile ramp.  A few pixels of smoothing makes the matte follow
        coherent structure (whole spiral arms) instead of per-pixel noise, so
        extended features stay opaque while the sky still cuts to zero.
        Requires scipy when > 0.

    Returns
    -------
    array (ny, nx) float in [0, 1]
    """
    inten = np.nan_to_num(np.asarray(intensity, dtype=float), nan=0.0)
    if smooth and float(smooth) > 0:
        try:
            from scipy.ndimage import gaussian_filter
        except ImportError as exc:  # pragma: no cover
            raise ImportError('smooth requires scipy') from exc
        inten = gaussian_filter(inten, sigma=float(smooth))
    if invert:
        inten = inten.max() - inten
    a_lo = _finite_percentile(inten, lo)
    a_hi = _finite_percentile(inten, hi)
    if not (a_hi > a_lo):
        a_hi = a_lo + 1e-9
    alpha = np.clip((inten - a_lo) / (a_hi - a_lo), 0.0, 1.0)
    g = float(gamma) if gamma is not None else 1.0
    if g <= 0:
        raise ValueError('gamma must be > 0')
    if g != 1.0:
        alpha = np.power(alpha, g)
    return alpha.astype(float)


def _saturation_map(rgb3):
    """HSV-style saturation ``(max - min) / max`` in [0, 1] (0 where max == 0)."""
    mx = np.nanmax(rgb3, axis=-1)
    mn = np.nanmin(rgb3, axis=-1)
    return np.divide(mx - mn, mx, out=np.zeros_like(mx, dtype=float), where=mx > 0)


def _sky_distance_map(rgb3, sky_color='black'):
    """
    Per-pixel color distance from *sky_color*, normalized to [0, 1].

    Reduces to a luminance-like map for a black sky, and to an inverted one
    for white — dark features on a light background survive the matte.
    """
    from matplotlib.colors import to_rgb
    sky = np.asarray(to_rgb(sky_color), dtype=float).reshape(1, 1, 3)
    dist = np.sqrt(np.nansum((rgb3 - sky) ** 2, axis=-1))
    return np.clip(dist / np.sqrt(3.0), 0.0, 1.0)


def _intensity_map(rgb, alpha_source='luma', layer=None, sky_color='black'):
    """
    Resolve the intensity map used for the alpha matte.

    alpha_source : 'luma' | 'max' | 'mean' | 'sat' | 'dist' | 'layer'
        ``'sat'`` — HSV saturation (colorful pixels survive; note that
        near-white *and* near-black pixels both fade).
        ``'dist'`` — color distance from *sky_color* (general form of luma:
        works for white or colored backgrounds, e.g. deblended composites).
        ``'layer'`` requires *layer* — a 2D array, or an integer channel index
        into *rgb* (0=R, 1=G, 2=B), or an RGBA alpha channel (3) if present.
    """
    arr = np.asarray(rgb, dtype=float)
    src = str(alpha_source or 'luma').strip().lower()
    if src in ('luma', 'lum', 'luminance', 'brightness'):
        return luminance(arr)
    if src == 'max':
        if arr.ndim == 2:
            return arr
        return np.nanmax(arr[..., :3], axis=-1)
    if src == 'mean':
        if arr.ndim == 2:
            return arr
        return np.nanmean(arr[..., :3], axis=-1)
    if src in ('sat', 'saturation'):
        if arr.ndim == 2:
            raise ValueError("alpha_source='sat' requires an RGB image")
        return _saturation_map(arr[..., :3])
    if src in ('dist', 'skydist', 'background'):
        if arr.ndim == 2:
            raise ValueError("alpha_source='dist' requires an RGB image")
        return _sky_distance_map(arr[..., :3], sky_color=sky_color)
    if src in ('layer', 'channel', 'band'):
        if layer is None:
            raise ValueError("alpha_source='layer' requires layer= (array or channel index)")
        if isinstance(layer, (int, np.integer)):
            if arr.ndim < 3 or int(layer) >= arr.shape[-1]:
                raise ValueError('layer index %s out of range for shape %s' % (layer, arr.shape))
            return np.asarray(arr[..., int(layer)], dtype=float)
        lay = np.asarray(layer, dtype=float)
        if lay.ndim != 2:
            raise ValueError('layer array must be 2D')
        if lay.shape != arr.shape[:2]:
            raise ValueError('layer shape %s does not match image %s' % (lay.shape, arr.shape[:2]))
        return lay
    raise ValueError("alpha_source must be 'luma', 'max', 'mean', 'sat', 'dist', "
                     "or 'layer' (got %r)" % alpha_source)


def _apply_matte(alpha, matte='none', soft_edge=0.0):
    """Optional circular/elliptical window and Gaussian feather on *alpha*."""
    alpha = np.asarray(alpha, dtype=float)
    ny, nx = alpha.shape
    kind = str(matte or 'none').strip().lower()
    if kind in ('circle', 'circular', 'ellipse', 'elliptical'):
        yy, xx = np.mgrid[0:ny, 0:nx].astype(float)
        cy, cx = (ny - 1) / 2.0, (nx - 1) / 2.0
        # Inscribe circle (or fill ellipse to the frame).
        if kind.startswith('circ'):
            r = min(nx, ny) / 2.0
            dist = np.hypot((xx - cx) / max(r, 1e-9), (yy - cy) / max(r, 1e-9))
        else:
            dist = np.hypot((xx - cx) / max(nx / 2.0, 1e-9), (yy - cy) / max(ny / 2.0, 1e-9))
        # Soft rim: 1 inside, 0 outside, linear falloff in last ~3% of radius.
        rim = 0.03
        window = np.clip((1.0 - dist) / rim, 0.0, 1.0)
        alpha = alpha * window
    elif kind not in ('none', '', 'full', 'free'):
        raise ValueError("matte must be 'none', 'circle', or 'ellipse' (got %r)" % matte)

    soft = float(soft_edge or 0.0)
    if soft > 0:
        try:
            from scipy.ndimage import gaussian_filter
        except ImportError as exc:  # pragma: no cover
            raise ImportError('soft_edge requires scipy') from exc
        # Soften only the matte edge without brightening empty sky much:
        # blur alpha then re-clamp; small sigma in pixels.
        alpha = gaussian_filter(alpha, sigma=soft)
        alpha = np.clip(alpha, 0.0, 1.0)
    return alpha


def _auto_crop(rgba, pad=0.05, alpha_threshold=1e-3):
    """Tight bbox around pixels with alpha > threshold, plus fractional padding."""
    rgba = np.asarray(rgba)
    alpha = rgba[..., 3]
    mask = alpha > float(alpha_threshold)
    if not np.any(mask):
        return rgba
    rows = np.any(mask, axis=1)
    cols = np.any(mask, axis=0)
    r0, r1 = int(np.argmax(rows)), int(len(rows) - np.argmax(rows[::-1]))
    c0, c1 = int(np.argmax(cols)), int(len(cols) - np.argmax(cols[::-1]))
    h, w = r1 - r0, c1 - c0
    pad_f = max(0.0, float(pad))
    pr = int(np.ceil(h * pad_f))
    pc = int(np.ceil(w * pad_f))
    ny, nx = alpha.shape
    r0 = max(0, r0 - pr)
    c0 = max(0, c0 - pc)
    r1 = min(ny, r1 + pr)
    c1 = min(nx, c1 + pc)
    return rgba[r0:r1, c0:c1]


def _resize_rgba(rgba, size, order=1):
    """Downsample (or leave) so the longest side is at most *size* pixels."""
    if size is None:
        return rgba
    size = int(size)
    if size < 1:
        raise ValueError('size must be a positive int or None')
    rgba = np.asarray(rgba, dtype=float)
    ny, nx = rgba.shape[:2]
    longest = max(ny, nx)
    if longest <= size:
        return rgba
    scale = size / float(longest)
    new_shape = (max(1, int(round(ny * scale))), max(1, int(round(nx * scale))), rgba.shape[-1])
    from skimage.transform import resize
    out = resize(rgba, new_shape, order=order, mode='reflect',
                 anti_aliasing=True, preserve_range=True)
    return np.clip(out, 0.0, 1.0)


def _crop_by_sky(arr, wcs, sky_box):
    """
    Crop an image array to a sky rectangle using *wcs*.

    Parameters
    ----------
    arr : array (ny, nx, ...)
    wcs : astropy.wcs.WCS
    sky_box : sequence of 4 floats
        ``(ra0, dec0, ra1, dec1)`` in degrees (order-independent corners).

    Returns
    -------
    array
        Cropped view (copy) of *arr*.
    """
    if sky_box is None:
        return arr
    box = [float(x) for x in sky_box]
    if len(box) != 4:
        raise ValueError('sky_box must be (ra0, dec0, ra1, dec1) in degrees')
    ra0, dec0, ra1, dec1 = box
    corners = np.array([
        [ra0, dec0], [ra0, dec1], [ra1, dec0], [ra1, dec1],
    ], dtype=float)
    try:
        pix = np.asarray(wcs.all_world2pix(corners, 0), dtype=float)
    except Exception:
        pix = np.asarray(wcs.wcs_world2pix(corners, 0), dtype=float)
    xs, ys = pix[:, 0], pix[:, 1]
    ny, nx = arr.shape[:2]
    c0 = max(0, int(np.floor(min(xs))))
    c1 = min(nx, int(np.ceil(max(xs))) + 1)
    r0 = max(0, int(np.floor(min(ys))))
    r1 = min(ny, int(np.ceil(max(ys))) + 1)
    if c1 <= c0 or r1 <= r0:
        raise ValueError('sky_box does not intersect the image')
    return np.asarray(arr[r0:r1, c0:c1]).copy()


# --------------------------------------------------------------------------- main API

def make_transparent_cutout(rgb, size=None, alpha_lo=55.0, alpha_hi=99.3,
                            alpha_gamma=0.5, alpha_source='luma', layer=None,
                            alpha_smooth=0.0, sky_color='black',
                            crop='auto', pad=0.05, matte='none', soft_edge=0.0,
                            invert=False, existing_alpha='replace',
                            wcs=None, sky_box=None):
    """
    Build a transparent-sky RGBA cutout from an RGB (or RGBA) composite.

    Parameters
    ----------
    rgb : array
        ``(ny, nx, 3)`` or ``(ny, nx, 4)`` in [0, 1] (values outside are clipped).
    size : int or None
        If set, downsample so the longest side is at most this many pixels
        (after cropping).  ``None`` keeps native resolution.
    alpha_lo, alpha_hi : float
        Percentiles of the intensity map mapped to fully transparent / opaque.
        ``alpha_hi`` acts as the opacity threshold: pixels at or above it are
        fully opaque, with alpha ramping down to zero at ``alpha_lo``.
    alpha_gamma : float
        Fade shape, any value > 0 (``< 1`` = bolder object body, ``> 1`` =
        tighter matte around the brightest structure; default 0.5 matches
        the skyplothelper M74 recipe).
    alpha_source : {'luma', 'max', 'mean', 'sat', 'dist', 'layer'}
        How to derive the intensity map for the matte.  ``'dist'`` uses the
        color distance from *sky_color* (use for non-black skies, e.g.
        white-background or deblended composites); ``'sat'`` uses HSV
        saturation.
    layer : array or int or None
        Extra band (2D) or channel index when ``alpha_source='layer'``.
    alpha_smooth : float
        Gaussian sigma (pixels) applied to the intensity map before the
        percentile ramp.  A few px makes the matte trace coherent structure
        (whole spiral arms opaque) instead of pixel noise; the RGB data
        itself is untouched.
    sky_color : color
        Background color for ``alpha_source='dist'`` (any matplotlib color).
    crop : {'auto', 'none'} or False
        ``'auto'`` trims to the non-transparent bbox (+ *pad*).
    pad : float
        Fractional padding around the auto-crop bbox (0.05 = 5%).
    matte : {'none', 'circle', 'ellipse'}
        Optional geometric window on top of the luminance matte.
    soft_edge : float
        Gaussian sigma (pixels) applied to alpha after the matte (0 = off).
    invert : bool
        Treat dark pixels as signal (white-background composites).
    existing_alpha : {'replace', 'multiply', 'keep'}
        If *rgb* already has an alpha channel: replace it, multiply into the
        new matte, or keep the original and only apply size/crop/matte.
    wcs, sky_box
        Optional WCS-aware pre-crop.  ``sky_box=(ra0, dec0, ra1, dec1)`` in
        degrees; requires an :class:`~astropy.wcs.WCS` matching *rgb*.

    Returns
    -------
    array
        Float RGBA ``(ny, nx, 4)`` in [0, 1].
    """
    arr = np.clip(np.nan_to_num(np.asarray(rgb, dtype=float)), 0.0, 1.0)
    if arr.ndim != 3 or arr.shape[-1] not in (3, 4):
        raise ValueError('rgb must have shape (ny, nx, 3) or (ny, nx, 4)')
    if sky_box is not None:
        if wcs is None:
            raise ValueError('sky_box requires wcs=')
        arr = _crop_by_sky(arr, wcs, sky_box)
        if layer is not None and not isinstance(layer, (int, np.integer)):
            layer = _crop_by_sky(np.asarray(layer), wcs, sky_box)

    rgb3 = arr[..., :3]
    mode = str(existing_alpha or 'replace').strip().lower()
    if arr.shape[-1] == 4 and mode == 'keep':
        alpha = np.clip(arr[..., 3], 0.0, 1.0)
    else:
        inten = _intensity_map(arr, alpha_source=alpha_source, layer=layer,
                               sky_color=sky_color)
        alpha = alpha_from_intensity(inten, lo=alpha_lo, hi=alpha_hi,
                                     gamma=alpha_gamma, invert=invert,
                                     smooth=alpha_smooth)
        if arr.shape[-1] == 4 and mode == 'multiply':
            alpha = alpha * np.clip(arr[..., 3], 0.0, 1.0)

    alpha = _apply_matte(alpha, matte=matte, soft_edge=soft_edge)
    rgba = np.dstack([rgb3, alpha])

    do_crop = crop not in (False, None, 'none', 'full', '')
    if do_crop:
        rgba = _auto_crop(rgba, pad=pad)

    rgba = _resize_rgba(rgba, size)
    return np.clip(rgba, 0.0, 1.0).astype(float)


def flatten_rgba(rgba, background='white', gamma=1.0):
    """
    Composite an RGBA image onto a solid background color (the inverse of
    :func:`deblend_background`).

    Per channel: ``C = alpha * F + (1 - alpha) * B``, optionally computed in
    gamma-decoded linear space when ``gamma != 1``.

    Parameters
    ----------
    rgba : array (ny, nx, 4)
        Foreground image, float [0, 1].
    background : color
        Solid background (any matplotlib color spec).
    gamma : float
        Encoding gamma of the inputs.  ``1.0`` blends the stored values
        directly; e.g. ``2.2`` decodes to linear light first and re-encodes
        the result.

    Returns
    -------
    array (ny, nx, 3) float in [0, 1]
    """
    from matplotlib.colors import to_rgb

    arr = np.clip(np.nan_to_num(np.asarray(rgba, dtype=float)), 0.0, 1.0)
    if arr.ndim != 3 or arr.shape[-1] != 4:
        raise ValueError('rgba must have shape (ny, nx, 4)')
    g = float(gamma)
    if g <= 0:
        raise ValueError('gamma must be > 0')
    bg = np.asarray(to_rgb(background), dtype=float).reshape(1, 1, 3)
    fg = arr[..., :3]
    alpha = arr[..., 3:4]
    if g != 1.0:
        fg, bg = fg ** g, bg ** g
    out = alpha * fg + (1.0 - alpha) * bg
    if g != 1.0:
        out = out ** (1.0 / g)
    return np.clip(out, 0.0, 1.0)


def _deblend_alpha(C, B):
    """
    Closed-form minimal alpha for ``C = alpha*F + (1-alpha)*B`` with F in
    the RGB gamut (both inputs in linear blend space).

    Per channel, ``F = B + (C - B)/alpha`` stays in [0, 1] exactly when
    ``alpha >= (C-B)/(1-B)`` (pixel brighter than background) or
    ``alpha >= (B-C)/B`` (darker).  The per-pixel minimum is the max over
    channels — fully vectorized, no iterative search.  Degenerate
    denominators (B at 0 or 1) are guarded: that direction's requirement
    reduces to ``alpha >= |C - B|``.
    """
    diff = C - B
    denom_up = np.where(B < 1.0, 1.0 - B, 1.0)
    denom_dn = np.where(B > 0.0, B, 1.0)
    need = np.where(diff >= 0.0, diff / denom_up, -diff / denom_dn)
    return np.clip(np.max(need, axis=-1), 0.0, 1.0)


def deblend_background(rgb, background='white', gamma=1.0, tol=1.5 / 255,
                       alpha_smooth=0.0, alpha_lo=None, alpha_hi=None,
                       alpha_gamma=1.0, fill='background',
                       matte='none', soft_edge=0.0,
                       crop='none', pad=0.05, size=None):
    """
    Recover an RGBA foreground from a composite flattened onto a known solid
    background color (un-premultiply) — the inverse of :func:`flatten_rgba`
    and the principled sibling of :func:`make_transparent_cutout`.

    Solves ``C = alpha * F + (1 - alpha) * B`` per pixel for the **minimum
    alpha** whose foreground ``F`` stays inside the RGB gamut (closed form;
    see :func:`_deblend_alpha`).  Works for any background color: black
    reduces to ``max(C)``, white to ``max(1 - C)``, and dark features on a
    light background survive.

    With the default settings the reconstruction is exact: flattening the
    returned RGBA back onto *background* (same *gamma*) reproduces the
    input.  The matte-shaping knobs below share the vocabulary — and the
    implementation — of :func:`make_transparent_cutout`; using them trades
    exactness for presentation control.

    Parameters
    ----------
    rgb : array (ny, nx, 3)
        Flattened composite, float [0, 1] (e.g. a JPEG/PNG loaded and scaled).
    background : color
        The known solid background it was flattened onto (any matplotlib
        color spec).
    gamma : float
        Encoding gamma assumed for the blend.  ``1.0`` (default) treats
        stored values as blend-space; use e.g. ``2.2`` if the composite was
        blended in linear light.
    tol : float
        Snap alpha to exactly 0 where the pixel is within *tol* of the
        background in every channel (default 1.5/255).  This absorbs 8-bit
        quantization / JPEG noise that would otherwise leave a faint
        near-transparent haze over the whole frame.
    alpha_smooth : float
        Gaussian sigma (pixels) applied to the recovered alpha map before
        the foreground is re-solved (gamut-clipped) — same trick as in
        :func:`make_transparent_cutout`: mild smoothing trades exact
        reconstruction near edges for a less noisy matte.  Requires scipy
        when > 0.
    alpha_lo, alpha_hi : float or None
        Optional percentile ramp re-applied to the *recovered* alpha, with
        the same meaning as in :func:`make_transparent_cutout`:
        ``alpha_hi`` is the full-opacity threshold, fading to transparent at
        ``alpha_lo``.  ``None`` (default) keeps the physical alpha.
    alpha_gamma : float
        Fade shape (> 0).  With a ramp, shapes it exactly as in
        :func:`make_transparent_cutout`; without one, applied as a plain
        power on the physical alpha (``< 1`` bolder body, ``> 1`` tighter
        matte).
    fill : {'background', 'black', 'white'} or color
        Foreground color to store where alpha == 0 (invisible; affects only
        later edits).  Default keeps the background color there.
    matte : {'none', 'circle', 'ellipse'}
        Optional geometric window (same as :func:`make_transparent_cutout`).
    soft_edge : float
        Gaussian feather (pixels) on the finished alpha (same as
        :func:`make_transparent_cutout`).
    crop : {'none', 'auto'} or False
        ``'auto'`` trims to the non-transparent bbox (+ *pad*).  Default
        ``'none'`` — deblending is a recovery operation, so the frame is
        preserved unless cropping is requested.
    pad : float
        Fractional padding around the auto-crop bbox.
    size : int or None
        Downsample so the longest side is at most this many pixels.

    Returns
    -------
    array (ny, nx, 4) float in [0, 1]

    Notes
    -----
    The problem is underdetermined (4 unknowns, 3 equations); minimum alpha
    is the standard, maximally-transparent convention (as used by
    color-keying / un-premultiply tools).  The recovered *F* touches the
    gamut boundary in at least one channel wherever ``0 < alpha < 1`` — that
    is inherent, not a bug.
    """
    from matplotlib.colors import to_rgb

    arr = np.clip(np.nan_to_num(np.asarray(rgb, dtype=float)), 0.0, 1.0)
    if arr.ndim != 3 or arr.shape[-1] < 3:
        raise ValueError('rgb must have shape (ny, nx, 3[+])')
    arr = arr[..., :3]
    g = float(gamma)
    if g <= 0:
        raise ValueError('gamma must be > 0')

    bg = np.asarray(to_rgb(background), dtype=float).reshape(1, 1, 3)
    C = arr ** g if g != 1.0 else arr
    B = bg ** g if g != 1.0 else bg

    alpha = _deblend_alpha(C, B)

    if tol and float(tol) > 0:
        alpha[np.max(np.abs(arr - bg), axis=-1) <= float(tol)] = 0.0

    if alpha_smooth and float(alpha_smooth) > 0:
        try:
            from scipy.ndimage import gaussian_filter
        except ImportError as exc:  # pragma: no cover
            raise ImportError('alpha_smooth requires scipy') from exc
        alpha = np.clip(gaussian_filter(alpha, sigma=float(alpha_smooth)),
                        0.0, 1.0)

    # Un-premultiply with the physical (optionally smoothed) alpha
    # (gamut-clipped; exact when alpha is the un-smoothed closed form).
    a3 = alpha[..., None]
    F = np.divide(C - (1.0 - a3) * B, a3,
                  out=np.broadcast_to(B, C.shape).copy(), where=a3 > 0)
    F = np.clip(F, 0.0, 1.0)
    if g != 1.0:
        F = F ** (1.0 / g)

    # Aesthetic matte shaping — same machinery as make_transparent_cutout,
    # applied on top of the physical alpha (breaks exact reconstruction).
    if alpha_lo is not None or alpha_hi is not None:
        alpha = alpha_from_intensity(
            alpha,
            lo=alpha_lo if alpha_lo is not None else 0.0,
            hi=alpha_hi if alpha_hi is not None else 100.0,
            gamma=alpha_gamma)
    elif alpha_gamma is not None and float(alpha_gamma) != 1.0:
        if float(alpha_gamma) <= 0:
            raise ValueError('alpha_gamma must be > 0')
        alpha = np.power(alpha, float(alpha_gamma))

    fill_kind = str(fill or 'background').strip().lower()
    if fill_kind != 'background':
        fill_rgb = np.asarray(to_rgb(fill), dtype=float)
        F = np.where(alpha[..., None] > 0, F, fill_rgb.reshape(1, 1, 3))

    alpha = _apply_matte(alpha, matte=matte, soft_edge=soft_edge)
    rgba = np.dstack([F, alpha])

    if crop not in (False, None, 'none', 'full', ''):
        rgba = _auto_crop(rgba, pad=pad)
    rgba = _resize_rgba(rgba, size)
    return np.clip(rgba, 0.0, 1.0).astype(float)


def save_transparent_cutout(rgba, savepath, overwrite=True):
    """
    Write an RGBA cutout to PNG or TIFF (format from the file extension).

    Arrays in this package use FITS / matplotlib ``origin='lower'`` (row 0
    at the bottom).  Files are written **vertically flipped** so typical
    image viewers (and HTML ``<img>``) show the same sky orientation as a
    matplotlib ``imshow(..., origin='lower')`` plot.  In-memory arrays
    returned by :func:`make_transparent_cutout` are *not* flipped.

    Parameters
    ----------
    rgba : array (ny, nx, 4)
        Float [0, 1] or uint8 (origin lower-left).
    savepath : str
        ``.png`` or ``.tif`` / ``.tiff``.
    overwrite : bool

    Returns
    -------
    str
        Absolute path written.
    """
    path = os.path.abspath(os.path.expanduser(savepath))
    if (not overwrite) and os.path.exists(path):
        raise FileExistsError(path)
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)

    arr = np.asarray(rgba)
    if arr.ndim != 3 or arr.shape[-1] != 4:
        raise ValueError('rgba must have shape (ny, nx, 4)')
    if np.issubdtype(arr.dtype, np.floating):
        arr8 = (np.clip(arr, 0, 1) * 255).astype(np.uint8)
    else:
        arr8 = np.clip(arr, 0, 255).astype(np.uint8)
    # FITS origin is lower-left; PNG/TIFF viewers treat row 0 as the top.
    arr8 = arr8[::-1]

    ext = os.path.splitext(path)[1].lower()
    try:
        from PIL import Image
        img = Image.fromarray(arr8, mode='RGBA')
        if ext in ('.tif', '.tiff'):
            img.save(path, format='TIFF')
        else:
            # PNG is the default for .png and any other extension.
            img.save(path, format='PNG')
    except ImportError:
        # matplotlib fallback (always available in this project).
        from matplotlib import image as mpl_image
        mpl_image.imsave(path, arr8)

    return path


def make_checkerboard(ny, nx, cell=16, light=0.85, dark=0.55):
    """
    Soft grey checkerboard background for previewing transparency.

    Returns
    -------
    array (ny, nx, 3) float
    """
    yy, xx = np.mgrid[0:ny, 0:nx]
    board = ((yy // cell) + (xx // cell)) % 2
    grey = np.where(board, float(light), float(dark))
    return np.dstack([grey, grey, grey])


def preview_cutout_on_backgrounds(rgba, backgrounds=None, figsize=None):
    """
    Show the cutout over several backgrounds (default: checkerboard, dark, light).

    Returns a matplotlib Figure (caller may ``plt.show()`` / ``savefig``).
    """
    import matplotlib.pyplot as plt

    rgba = np.asarray(rgba, dtype=float)
    if rgba.ndim != 3 or rgba.shape[-1] != 4:
        raise ValueError('rgba must have shape (ny, nx, 4)')
    ny, nx = rgba.shape[:2]
    if backgrounds is None:
        backgrounds = [
            ('checkerboard', make_checkerboard(ny, nx)),
            ('dark', np.full((ny, nx, 3), 0.08)),
            ('light', np.full((ny, nx, 3), 0.94)),
        ]

    n = len(backgrounds)
    fig, axes = plt.subplots(1, n, figsize=figsize or (3.2 * n, 3.2))
    if n == 1:
        axes = [axes]
    rgb = rgba[..., :3]
    alpha = rgba[..., 3:4]
    for ax, (label, bg) in zip(axes, backgrounds):
        bg = np.asarray(bg, dtype=float)
        if bg.ndim == 2:
            bg = np.dstack([bg, bg, bg])
        if bg.shape[:2] != (ny, nx):
            from skimage.transform import resize
            bg = resize(bg, (ny, nx, 3), order=1, mode='edge',
                        anti_aliasing=False, preserve_range=True)
        comp = rgb * alpha + bg * (1.0 - alpha)
        ax.imshow(np.clip(comp, 0, 1), origin='lower')
        ax.set_title(str(label), fontsize=10)
        ax.set_xticks([])
        ax.set_yticks([])
    fig.suptitle('Transparent cutout preview', fontsize=12)
    fig.tight_layout()
    return fig


def batch_transparent_cutouts(images, **kwargs):
    """
    Apply :func:`make_transparent_cutout` to a list of RGB images with shared settings.

    Parameters
    ----------
    images : sequence of arrays
    **kwargs
        Passed to :func:`make_transparent_cutout`.

    Returns
    -------
    list of RGBA arrays
    """
    return [make_transparent_cutout(im, **kwargs) for im in images]
