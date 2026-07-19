"""
Color conversion and image colorization pipeline: the heart of multicolorfits.

Workflow:  to_grey_rgb() -> colorize_image() -> combine_multicolor()
"""

import colorsys

import numpy as np
from skimage import filters
from astropy.visualization import ManualInterval

from .scaling import stretch_functions, adjust_gamma, SIGNED_STRETCHES, rescale_image

__all__ = [
    'hex_to_rgb',
    'rgb_to_hex',
    'hex_to_hsv',
    'rgb_to_hsv',
    'hsv_to_rgb',
    'hex_complement',
    'to_hex',
    'to_grey_rgb',
    'colorize_image_direct_rgb',
    'colorize_image',
    'combine_multicolor',
    'smooth_image',
]


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
        RGB tuple such as (255,255,255)
    """
    # From http://stackoverflow.com/a/214657
    hexstring = hexstring.lstrip('#')
    lv = len(hexstring)
    return tuple(int(hexstring[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))


def rgb_to_hex(rgb):
    """
    Converts RGB tuple to a hexadecimal string

    Parameters
    ----------
    rgb : tuple
        RGB tuple such as (255,255,255)

    Returns
    -------
    str
        Hexadecimal string such as '#ffffff'
    """
    return '#%02x%02x%02x' % tuple(rgb)  # input RGB as tuple.  e.g.: rgb_to_hex((255, 255, 255))


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
    # See Wikipedia article for details -- https://en.wikipedia.org/wiki/HSL_and_HSV#From_HSV  .  Modified from colorsys.py,
    # HSV: Hue, Saturation, Value
    # H = position in the spectrum, S = color saturation ("purity"), V = color brightness
    r, g, b = np.array(hex_to_rgb(hexstring)) / 255.  # Convert Hex to RGB fracs (i.e., 0..1 instead of 0..255)
    maxc = max(r, g, b); minc = min(r, g, b)
    s = (maxc - minc) / maxc
    v = maxc
    if minc == maxc: return 0.0, 0.0, v
    rc = (maxc - r) / (maxc - minc); gc = (maxc - g) / (maxc - minc); bc = (maxc - b) / (maxc - minc)
    if r == maxc: h = bc - gc
    elif g == maxc: h = 2.0 + rc - bc
    else: h = 4.0 + gc - rc
    h = (h / 6.0) % 1.0
    return h, s, v  # All in fractions in range [0...1]


def rgb_to_hsv(rgb):
    """
    Vectorized RGB to HSV conversion,
    ~3-10x faster than skimage.color.rgb2hsv for large arrays

    Parameters
    ----------
    rgb : numpy.ndarray
        RGB image with shape (H, W, 3) and values in [0, 1]

    Returns
    -------
    numpy.ndarray
        HSV image with same shape, values in [0, 1]
    """
    rgb = np.asarray(rgb, dtype=np.float32)

    r, g, b = rgb[..., 0], rgb[..., 1], rgb[..., 2]

    max_rgb = np.maximum(np.maximum(r, g), b)
    min_rgb = np.minimum(np.minimum(r, g), b)
    delta = max_rgb - min_rgb

    h = np.zeros_like(max_rgb)
    s = np.zeros_like(max_rgb)
    v = max_rgb  # Value is just the max RGB

    # Saturation (avoid division by zero)
    mask_nonzero = max_rgb != 0
    s[mask_nonzero] = delta[mask_nonzero] / max_rgb[mask_nonzero]

    # Hue (avoid division by zero)
    mask_delta = delta != 0

    # Case 1: Red is maximum
    mask_r = (max_rgb == r) & mask_delta
    h[mask_r] = (60 * ((g[mask_r] - b[mask_r]) / delta[mask_r]) + 360) % 360

    # Case 2: Green is maximum
    mask_g = (max_rgb == g) & mask_delta
    h[mask_g] = (60 * ((b[mask_g] - r[mask_g]) / delta[mask_g]) + 120) % 360

    # Case 3: Blue is maximum
    mask_b = (max_rgb == b) & mask_delta
    h[mask_b] = (60 * ((r[mask_b] - g[mask_b]) / delta[mask_b]) + 240) % 360

    h = h / 360.0

    return np.stack([h, s, v], axis=-1)


def hsv_to_rgb(hsv):
    """
    Vectorized HSV → RGB conversion for image-shaped arrays.

    Faster than ``skimage.color.hsv2rgb`` on large arrays while matching its
    [0, 1] channel convention.

    Parameters
    ----------
    hsv : numpy.ndarray
        HSV image with shape ``(..., 3)`` and values in [0, 1]
        (H as a fraction of a turn, not degrees).

    Returns
    -------
    numpy.ndarray
        RGB image of the same shape, values in [0, 1].
    """
    hsv = np.asarray(hsv, dtype=np.float32)

    h, s, v = hsv[..., 0], hsv[..., 1], hsv[..., 2]

    h = h * 360.0

    c = v * s  # Chroma
    h_prime = h / 60.0
    x = c * (1 - np.abs((h_prime % 2) - 1))
    m = v - c

    r = np.zeros_like(h)
    g = np.zeros_like(h)
    b = np.zeros_like(h)

    mask0 = (0 <= h_prime) & (h_prime < 1)
    mask1 = (1 <= h_prime) & (h_prime < 2)
    mask2 = (2 <= h_prime) & (h_prime < 3)
    mask3 = (3 <= h_prime) & (h_prime < 4)
    mask4 = (4 <= h_prime) & (h_prime < 5)
    mask5 = (5 <= h_prime) & (h_prime < 6)

    # Sector 0: (c, x, 0)
    r[mask0] = c[mask0]; g[mask0] = x[mask0]; b[mask0] = 0
    # Sector 1: (x, c, 0)
    r[mask1] = x[mask1]; g[mask1] = c[mask1]; b[mask1] = 0
    # Sector 2: (0, c, x)
    r[mask2] = 0; g[mask2] = c[mask2]; b[mask2] = x[mask2]
    # Sector 3: (0, x, c)
    r[mask3] = 0; g[mask3] = x[mask3]; b[mask3] = c[mask3]
    # Sector 4: (x, 0, c)
    r[mask4] = x[mask4]; g[mask4] = 0; b[mask4] = c[mask4]
    # Sector 5: (c, 0, x)
    r[mask5] = c[mask5]; g[mask5] = 0; b[mask5] = x[mask5]

    r += m
    g += m
    b += m

    rgb = np.stack([r, g, b], axis=-1)

    return np.clip(rgb, 0, 1)


def hex_complement(hexstring):
    """
    Return the complementary hex color (bitwise inverse on the RGB cube).

    Examples
    --------
    >>> hex_complement('#FF0000')
    '#00FFFF'
    >>> hex_complement('#00FF00')
    '#FF00FF'

    Parameters
    ----------
    hexstring : str
        Hexadecimal color such as ``'#FF0000'`` (leading ``#`` optional).

    Returns
    -------
    str
        Complementary hex color such as ``'#00FFFF'``.

    Notes
    -----
    Used by the legacy ``inverse`` display path (replace each layer color with
    its complement before compositing).  Prefer a non-black
    ``combine_background`` when you want a white paper look without flipping hues.
    """
    if hexstring[0] == '#': hexstring = hexstring.lstrip('#')
    hexcolor = int(hexstring, 16)
    color_comp = 0xFFFFFF ^ hexcolor
    hexcolor_comp = "#%06X" % color_comp
    return hexcolor_comp


try:
    from matplotlib.colors import to_hex  # Converts matplotlib colors to hex/html.  Only in mpl>=2.0
except ImportError:
    def to_hex(c):
        # For mpl v<2.0
        import matplotlib.colors as mplc
        return mplc.rgb2hex(mplc.colorConverter.to_rgb(c))  # could use to_rgba to keep alpha


def _intensity_to_grey_rgb(plane, dtype=None):
    """
    Expand a 2D intensity map to a greyscale RGB cube ``(ny, nx, 3)``.

    Same contract as ``skimage.color.gray2rgb`` (writable float array in [0, 1]
    with identical R=G=B), but uses ``broadcast_to`` + copy instead of skimage.
    """
    plane = np.clip(np.asarray(plane), 0.0, 1.0)
    if dtype is not None:
        plane = plane.astype(dtype, copy=False)
    # Materialize: downstream may write, and broadcast views share memory across
    # channels (unsafe if a channel is mutated independently).
    out = np.broadcast_to(plane[..., np.newaxis], plane.shape + (3,)).copy()
    return out


def to_grey_rgb(datin, rescalefn='linear', scaletype='abs', min_max=[None, None], gamma=2.2, checkscale=False, dtype=None):
    """
    Stretch a 2D intensity array and expand it to a greyscale RGB cube in [0, 1].

    This is the first step of the classic scripting pipeline before
    :func:`colorize_image` and :func:`combine_multicolor`.

    Examples
    --------
    >>> grey = to_grey_rgb(data, rescalefn='asinh', min_max=[0., 1.2], gamma=2.2)
    >>> grey.shape  # doctest: +SKIP
    (ny, nx, 3)

    Parameters
    ----------
    datin : array
        Input 2D image data.
    rescalefn : str
        Intensity stretch: ``'linear'``, ``'sqrt'``, ``'squared'``, ``'log'``,
        ``'power'``, ``'sinh'``, ``'asinh'``, or signed-data ``'symlog'`` /
        ``'symmetric_log'`` (the latter requires ``pip install pysymlog``).
    scaletype : str
        ``'abs'`` interprets ``min_max`` as data values; ``'perc'`` as percentiles.
    min_max : list
        ``[min, max]`` for the stretch.  With ``scaletype='perc'``, use e.g.
        ``[1., 95.]``.
    gamma : float
        Gamma applied after the stretch.  Use ``2.2`` when combining colorized
        frames; use ``1/2.2`` only for legacy inverse workflows.
    checkscale : bool
        If True, show a matplotlib comparison of input vs scaled preview.
    dtype : numpy dtype or None
        Output dtype.  Default preserves float64 behaviour.  Pass
        ``numpy.float32`` for a lower-memory interactive preview.

    Returns
    -------
    array
        Greyscale RGB image with shape ``(ny, nx, 3)`` and values in [0, 1].
    """
    if 'per' in scaletype.lower():
        if min_max == [None, None]: min_max = [0., 100.]
        minval, maxval = np.percentile(np.ma.masked_invalid(datin).compressed(), min_max)
    else:
        minval = [np.nanmin(datin) if min_max[0] is None else min_max[0]][0]
        maxval = [np.nanmax(datin) if min_max[1] is None else min_max[1]][0]
    # Use specified rescaling function (astropy stretches — keep as-is).
    if rescalefn.lower() in SIGNED_STRETCHES:
        datscaled = rescale_image(datin, rescalefn=rescalefn,
                                  vmin=minval, vmax=maxval)
    else:
        datscaled = (stretch_functions[rescalefn]() + ManualInterval(vmin=minval, vmax=maxval))(datin)
    if gamma != 1:
        datscaled = adjust_gamma(datscaled, gamma)
    # ManualInterval already maps into ~[0, 1]; clip replaces LinearStretch().
    dat_greyRGB = _intensity_to_grey_rgb(datscaled, dtype=dtype)

    if checkscale is not False:
        from matplotlib import pyplot as plt
        plt.clf(); plt.close('all')
        fig0 = plt.figure(0)
        ax1 = fig0.add_subplot(121)
        plt.imshow(datin, interpolation='nearest', origin='lower', cmap='gist_gray')
        plt.title('Input Image')
        ax2 = fig0.add_subplot(122, sharex=ax1, sharey=ax1)
        plt.imshow(dat_greyRGB**(1. / gamma), interpolation='nearest', origin='lower')
        plt.title('Scaled Image')
        plt.show()

    return dat_greyRGB


def colorize_image_direct_rgb(grey_rgb_image, hex_color, brightness_factor=1.0, gammacorr_color=1):
    """
    Tint a greyscale RGB cube by multiplying channels with a solid hex color.

    Equivalent to the HSV path in :func:`colorize_image` for greyscale inputs
    (``hsv_to_rgb(h, s, V)`` is linear in ``V``, so tinting reduces to
    ``intensity * (color_rgb ** gammacorr_color)``), but avoids a full
    RGB↔HSV round-trip.  Used automatically by :func:`colorize_image` for
    hex colors.

    Parameters
    ----------
    grey_rgb_image : numpy.ndarray
        Greyscale RGB image (R = G = B), typically from :func:`to_grey_rgb`.
    hex_color : str
        Hex color such as ``'#FF0000'``.
    brightness_factor : float
        Multiplier on the tinted result (analogous to HSV value).
    gammacorr_color : float
        Gamma applied to the color itself (matches ``colorize_image``).
        Use ``1`` for no color gamma.

    Returns
    -------
    numpy.ndarray
        Colorized RGB image in [0, 1], same shape as the input.
    """
    # Convert hex to RGB
    if isinstance(hex_color, str):
        hex_color = hex_color.lstrip('#')
        color_rgb = np.array([int(hex_color[i:i + 2], 16) for i in (0, 2, 4)]) / 255.0
    else:
        color_rgb = np.array(hex_color)

    if gammacorr_color != 1:
        color_rgb = color_rgb ** gammacorr_color

    # Take the intensity from the first channel (assuming grayscale)
    intensity = grey_rgb_image[..., 0:1]  # Keep dimensions for broadcasting

    # Apply color and brightness in one vectorized operation
    colorized = intensity * color_rgb[np.newaxis, np.newaxis, :] * brightness_factor

    return np.clip(colorized, 0, 1).astype(np.float32)


def colorize_image(image, colorvals, colorintype='hex', dtype=np.float64, gammacorr_color=1):
    """
    Add color of the given hue to an RGB greyscale image.

    Parameters
    ----------
    image : array
        Greyscale RGB image -- as would be output from to_grey_rgb()
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
    if colorintype not in ['hsv', 'hsv_dict', 'rgb', 'hex']:
        raise Exception("  colorintype must be 'hsv', 'hsv_dict', 'rgb', or 'hex'")

    # Guard the common misuse colorize_image(raw_2d_data, color): the 1st
    # argument must be a GREY RGB CUBE (ny, nx, 3) from to_grey_rgb(), not a
    # raw 2D intensity array.
    _img = np.asarray(image)
    if _img.ndim == 2:
        raise TypeError(
            "colorize_image: the 1st argument must be a greyscale RGB cube "
            "with shape (ny, nx, 3), got a 2D array with shape %r. Stretch "
            "your raw data first: grey = to_grey_rgb(data, rescalefn='asinh', "
            "min_max=[lo, hi]); layer = colorize_image(grey, '#FF0000')."
            % (_img.shape,))
    if _img.ndim != 3 or _img.shape[-1] != 3:
        raise TypeError(
            "colorize_image: the 1st argument must be a greyscale RGB cube "
            "with shape (ny, nx, 3), got shape %r. Build it with "
            "to_grey_rgb(data, ...)." % (_img.shape,))

    # Fast color conversion to HSV
    if colorintype.lower() == 'rgb':
        r, g, b = np.array(colorvals) / 255.0
        hue, saturation, v = colorsys.rgb_to_hsv(r, g, b)
    elif colorintype.lower() == 'hex':
        if isinstance(colorvals, str):
            # For a string hex color, MUCH faster to bypass HSV entirely.
            # The color's gamma folds directly into the RGB fast path
            # (hsv_to_rgb is linear in V), so this stays equivalent to the
            # HSV method for greyscale input while being ~3-4x faster.
            return colorize_image_direct_rgb(image, colorvals, gammacorr_color=gammacorr_color).astype(dtype)
        else:
            hue, saturation, v = hex_to_hsv(colorvals)
    elif colorintype.lower() == 'hsv_dict':
        hue, saturation, v = colorvals['hue'], colorvals['sat'], colorvals['v']
    else:  # hsv
        hue, saturation, v = colorvals

    # Gamma correction for color
    if gammacorr_color != 1:
        rgb_corrected = np.array(colorsys.hsv_to_rgb(hue, saturation, v)) ** gammacorr_color
        hue, saturation, v = colorsys.rgb_to_hsv(*rgb_corrected)

    hsv = rgb_to_hsv(image).astype(dtype)

    # Vectorized HSV channel updates
    hsv[..., 0] = hue
    hsv[..., 1] = saturation
    hsv[..., 2] *= v

    return hsv_to_rgb(hsv).astype(dtype)


def combine_multicolor(im_list_colorized, gamma=2.2, inverse=False, dtype=np.float64):
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
    dtype : numpy dtype
        Working/output dtype (default numpy.float64).  Pass numpy.float32 for an
        opt-in low-memory / faster "preview" mode (~half the memory traffic) at
        reduced numeric precision.

    Returns
    -------
    array
        Colorized RGB image (combined), shape=[ypixels,xpixels,3]
    """
    # Guard the common misuse combine_multicolor(single_image) or a list of raw
    # 2D arrays: the argument is a LIST of colorized (ny, nx, 3) layers from
    # colorize_image().
    if isinstance(im_list_colorized, np.ndarray):
        raise TypeError(
            "combine_multicolor: pass a LIST of colorized layers, not a single "
            "array (got an ndarray with shape %r). Each layer is an (ny, nx, 3) "
            "result of colorize_image(); e.g. combine_multicolor([layer_a, "
            "layer_b])." % (im_list_colorized.shape,))
    if len(im_list_colorized) == 0:
        raise ValueError("combine_multicolor: got an empty list of layers.")
    _first = np.asarray(im_list_colorized[0])
    if _first.ndim != 3 or _first.shape[-1] != 3:
        raise TypeError(
            "combine_multicolor: each list element must be a colorized RGB "
            "layer with shape (ny, nx, 3), got an element with shape %r. "
            "Colorize each band first: colorize_image(to_grey_rgb(data), "
            "'#FF0000')." % (_first.shape,))

    # Normalize dtype to a scalar type so typed zeros work across numpy versions
    dtype = np.dtype(dtype).type
    # Sum the colorized layers (NaN treated as 0, matching np.nansum) without
    # materializing an (N, ny, nx, 3) stack.  This is numerically identical to
    # the previous np.nansum + skimage implementation but ~1.5-2x faster on
    # large images, since it avoids the list->array copy, per-channel
    # skimage.rescale_intensity calls, and repeated nan_to_num passes.
    combined_RGB = np.asarray(im_list_colorized[0], dtype=dtype)
    combined_RGB = np.where(np.isnan(combined_RGB), dtype(0), combined_RGB)  # copy + NaN->0
    for im in im_list_colorized[1:]:
        im = np.asarray(im, dtype=dtype)
        combined_RGB += np.where(np.isnan(im), dtype(0), im)
    np.clip(combined_RGB, 0.0, 1.0, out=combined_RGB)  # == LinearStretch()

    channel_maxes = combined_RGB.max(axis=(0, 1))
    channel_mins = combined_RGB.min(axis=(0, 1))
    if inverse is True:
        max_of_maxints = float(np.max(1. - channel_maxes))
    else:
        max_of_maxints = float(np.max(channel_maxes))
    # Guard against divide-by-zero when all channels saturate (e.g. fully
    # overlapping inverse-colorized images) -- no rescale needed in that case.
    if max_of_maxints > 0:
        # Per channel, the old code did rescale_intensity(ch, out_range=(0, cmax/M)),
        # i.e. (ch - cmin)/(cmax - cmin) * (cmax / M).  Vectorize across channels.
        denom = channel_maxes - channel_mins
        valid = (channel_maxes > 0) & (denom > 0)
        with np.errstate(divide='ignore', invalid='ignore'):
            scale = np.where(valid, (channel_maxes / max_of_maxints) / denom, 0.0).astype(dtype)
        offset = np.where(valid, channel_mins, 0.0).astype(dtype)
        combined_RGB = (combined_RGB - offset) * scale
        combined_RGB = np.nan_to_num(combined_RGB)
    else:
        combined_RGB = np.nan_to_num(combined_RGB)
    np.power(combined_RGB, dtype(1. / gamma), out=combined_RGB)
    np.clip(combined_RGB, 0.0, 1.0, out=combined_RGB)  # == LinearStretch()
    if inverse is True: combined_RGB = (1. - combined_RGB).astype(dtype)
    return combined_RGB.astype(dtype, copy=False)


def smooth_image(datain, sigma=3):
    """
    Apply isotropic Gaussian smoothing to a 2D or multichannel image.

    Useful for matching a coarser neighbour's resolution or reducing noise
    before colorizing.

    Parameters
    ----------
    datain : array
        2D image or 3D multichannel cube (channels last).
    sigma : int or float
        Gaussian sigma in pixels.

    Returns
    -------
    array
        Smoothed image, same shape as ``datain``.
    """
    if len(datain.shape) == 2:
        return filters.gaussian(datain, sigma=sigma)
    else:
        try:
            return filters.gaussian(datain, sigma=sigma, channel_axis=-1)
        except TypeError:
            # older scikit-image (<0.19) used the multichannel keyword
            return filters.gaussian(datain, sigma=sigma, multichannel=True)
