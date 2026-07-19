"""
Color-space compositing beyond simple RGB addition.

The classic multicolorfits pipeline combines colorized layers by summing RGB
channels (``combine_multicolor``).  That is fast and familiar, but bright
overlaps can saturate toward white and dilute individual layer hues.

This module provides alternatives in perceptually motivated color spaces:

* ``'lab'`` — **supported.**  Layers are converted to CIE L\*a\*b\*; chromatic
  (a\*, b\*) components are a lightness-weighted average and L\* uses a
  selectable blend.  Lab is approximately perceptually uniform, so transitions
  are smoother and white clipping is reduced.
  ``colorspace='lab', blend='screen'`` is the recommended general alternative
  to classic RGB; ``blend='max'`` best preserves each layer's hue.
* ``'hsv'`` / ``'hsl'`` — **experimental.**  Hues are combined as 2D chroma
  vectors (``x = sat·cos(hue)``, ``y = sat·sin(hue)``), weighted by each
  layer's brightness, summed, then converted back
  (``h = atan2``, ``s = |vector|``).  Lightness/value is blended separately.
  Simple hue interpolation fails for distant hues (for example red→blue via
  green); vector summation does not (see https://stackoverflow.com/a/53328189).
  These modes can shift the global cast for palettes unbalanced around the hue
  wheel, so they remain experimental.

Complementary hues of equal brightness neutralize toward grey in any additive
light model — their chroma vectors cancel.  For paint-like subtractive mixes
(yellow + blue → green), use RYB/CMYK via ``combine_multicolor_alpha`` or
choose non-complementary layer colors (see ``multicolorfits.palettes``).

Also included:

* ``mix_colors_hex`` — mix plain colors (not images) for palette design.
* ``greyscale_image`` / ``simulate_colorblindness`` — accessibility checks on
  a finished composite.

All functions are pure NumPy / scikit-image.
"""

import numpy as np
from skimage import color as ski_color

from .colorize import (
    hex_to_rgb,
    rgb_to_hex,
    rgb_to_hsv,
    hsv_to_rgb,
)

__all__ = [
    'rgb_to_hsl',
    'hsl_to_rgb',
    'rgb_to_lab',
    'lab_to_rgb',
    'combine_multicolor_colorspace',
    'mix_colors_hex',
    'greyscale_image',
    'simulate_colorblindness',
]

COLORSPACES = ('lab', 'hsv', 'hsl')
BLENDS = ('screen', 'sum', 'max', 'mean')


# ---------------------------------------------------------------------------
# HSL conversions (vectorized; parallel to the HSV pair in colorize.py)
# ---------------------------------------------------------------------------

def rgb_to_hsl(rgb):
    """
    Vectorized RGB to HSL conversion.

    Parameters
    ----------
    rgb : numpy.ndarray
        RGB image with shape (..., 3) and values in [0, 1]

    Returns
    -------
    numpy.ndarray
        HSL image with same shape: channels (hue, saturation, lightness),
        all in [0, 1]
    """
    rgb = np.asarray(rgb, dtype=np.float64)
    r, g, b = rgb[..., 0], rgb[..., 1], rgb[..., 2]

    maxc = np.maximum(np.maximum(r, g), b)
    minc = np.minimum(np.minimum(r, g), b)
    delta = maxc - minc

    l = (maxc + minc) / 2.

    s = np.zeros_like(l)
    denom = 1. - np.abs(2. * l - 1.)
    mask_s = denom > 1e-12
    s[mask_s] = delta[mask_s] / denom[mask_s]

    h = np.zeros_like(l)
    mask_delta = delta != 0
    mask_r = (maxc == r) & mask_delta
    h[mask_r] = (60 * ((g[mask_r] - b[mask_r]) / delta[mask_r]) + 360) % 360
    mask_g = (maxc == g) & mask_delta
    h[mask_g] = (60 * ((b[mask_g] - r[mask_g]) / delta[mask_g]) + 120) % 360
    mask_b = (maxc == b) & mask_delta
    h[mask_b] = (60 * ((r[mask_b] - g[mask_b]) / delta[mask_b]) + 240) % 360

    return np.stack([h / 360., s, l], axis=-1)


def hsl_to_rgb(hsl):
    """
    Vectorized HSL to RGB conversion.  Channels (hue, saturation, lightness)
    in [0, 1] -> RGB in [0, 1].
    """
    hsl = np.asarray(hsl, dtype=np.float64)
    h, s, l = hsl[..., 0], hsl[..., 1], hsl[..., 2]

    c = (1. - np.abs(2. * l - 1.)) * s  # chroma
    h_prime = (h % 1.0) * 6.0
    x = c * (1. - np.abs((h_prime % 2) - 1.))
    m = l - c / 2.

    r = np.zeros_like(h)
    g = np.zeros_like(h)
    b = np.zeros_like(h)

    masks = [((0 <= h_prime) & (h_prime < 1)), ((1 <= h_prime) & (h_prime < 2)),
             ((2 <= h_prime) & (h_prime < 3)), ((3 <= h_prime) & (h_prime < 4)),
             ((4 <= h_prime) & (h_prime < 5)), ((5 <= h_prime) & (h_prime <= 6))]
    # (r, g, b) per 60-degree sector
    sector_vals = [(c, x, 0), (x, c, 0), (0, c, x), (0, x, c), (x, 0, c), (c, 0, x)]
    for mask, (rv, gv, bv) in zip(masks, sector_vals):
        r[mask] = rv[mask] if isinstance(rv, np.ndarray) else rv
        g[mask] = gv[mask] if isinstance(gv, np.ndarray) else gv
        b[mask] = bv[mask] if isinstance(bv, np.ndarray) else bv

    rgb = np.stack([r + m, g + m, b + m], axis=-1)
    return np.clip(rgb, 0, 1)


# ---------------------------------------------------------------------------
# CIELAB wrappers
# ---------------------------------------------------------------------------

def rgb_to_lab(rgb):
    """
    Convert an sRGB image (values in [0,1], shape (...,3)) to CIE L*a*b*.
    L* in [0,100]; a*,b* roughly in [-128,127].  Thin wrapper around
    skimage.color.rgb2lab with NaN protection.
    """
    return ski_color.rgb2lab(np.clip(np.nan_to_num(np.asarray(rgb, dtype=float)), 0, 1))


def lab_to_rgb(lab):
    """
    Convert a CIE L*a*b* image back to sRGB in [0,1].  Out-of-gamut colors
    are clipped.
    """
    return np.clip(ski_color.lab2rgb(np.asarray(lab, dtype=float)), 0, 1)


# ---------------------------------------------------------------------------
# Compositing
# ---------------------------------------------------------------------------

def _blend_stack(stack, blend):
    """Combine a (N, ...) stack of [0,1] lightness arrays into one."""
    if blend == 'sum':
        return np.clip(np.sum(stack, axis=0), 0, 1)
    elif blend == 'screen':
        return 1. - np.prod(1. - stack, axis=0)
    elif blend == 'max':
        return np.max(stack, axis=0)
    elif blend == 'mean':
        return np.mean(stack, axis=0)
    raise ValueError("blend must be one of %s" % (BLENDS,))


def combine_multicolor_colorspace(im_list_colorized, colorspace='lab', blend='screen',
                                  weights=None, gamma=2.2, inverse=False, dtype=np.float64):
    """
    Combine colorized RGB layers in a perceptual color space instead of by
    direct RGB channel addition.

    Drop-in alternative to combine_multicolor(): takes the same input (a list
    of colorized images from colorize_image(), i.e. the standard pipeline
    where to_grey_rgb applied ``adjust_gamma(scaled, gamma)``), and
    returns a display-ready [0,1] RGB image.

    The ``'lab'`` color space is a stable, supported option (the recommended
    general alternative to classic RGB).  ``'hsv'`` / ``'hsl'`` remain
    EXPERIMENTAL -- they can shift the overall color cast for palettes that are
    unbalanced around the hue wheel.

    Parameters
    ----------
    im_list_colorized : list
        List of colorized RGB images, e.g. [ halpha_purple, co21_orange ]
    colorspace : str
        'lab' (default, STABLE) -- combine chroma (a*,b*) as a lightness-weighted
        average in CIELAB; perceptually smooth, resists white-saturation.
        'hsv' or 'hsl' (EXPERIMENTAL) -- combine hues as brightness-weighted
        chroma vectors (x=s*cos(h), y=s*sin(h)); saturation falls out of the
        vector length.
    blend : str
        How the lightness channels combine across layers:
        'screen' (default) -- 1 - prod(1 - L_i).  Accumulates brightness like
        addition but saturates gracefully instead of clipping.
        'sum' -- clipped sum (closest to classic combine_multicolor behavior).
        'max' -- maximum (each pixel shows its brightest layer).
        'mean' -- average (best for mixing plain colors, dims composites).
    weights : list or None
        Optional per-layer weight factors for the chroma combination.
    gamma : float
        The gamma value used in the to_grey_rgb/colorize steps (default 2.2).
        Layers are converted to display-referred values (raised to
        ``1/gamma``) before the color-space math, so the output needs no
        further correction.
        Use gamma=1 if your inputs are already display-ready RGB.
    inverse : bool
        True inverts the lightness channel (white background) while keeping
        the chroma -- without flipping layer hues via ``hex_complement``.
    dtype : numpy dtype
        Output dtype.

    Returns
    -------
    array
        Combined display-ready RGB image, shape=[ypixels,xpixels,3] in [0,1]

    Notes
    -----
    Complementary hues at equal brightness cancel to grey by construction
    (additive light).  See the module docstring for discussion.
    """
    if colorspace not in COLORSPACES:
        raise ValueError("colorspace must be one of %s" % (COLORSPACES,))
    if blend not in BLENDS:
        raise ValueError("blend must be one of %s" % (BLENDS,))
    n = len(im_list_colorized)
    if n == 0:
        raise ValueError('im_list_colorized is empty')
    if weights is None:
        weights = np.ones(n)
    weights = np.asarray(weights, dtype=float)
    if weights.shape != (n,):
        raise ValueError('weights must have one value per image')

    # Display-referred layers: undo the gamma applied in to_grey_rgb/colorize so
    # perceptual mixing happens in linearish display space.
    disp = []
    for im in im_list_colorized:
        arr = np.clip(np.nan_to_num(np.asarray(im, dtype=float)), 0, 1)
        disp.append(arr ** (1. / gamma) if gamma != 1 else arr)

    eps = 1e-9

    if colorspace == 'lab':
        labs = [ski_color.rgb2lab(d) for d in disp]
        # skimage uses L* in [0, 100]; normalize to [0, 1] for blending.
        L = np.stack([lab[..., 0] / 100. for lab in labs])   # (N, ny, nx)
        a = np.stack([lab[..., 1] for lab in labs])
        b = np.stack([lab[..., 2] for lab in labs])
        # Weight chroma by each layer's lightness so faint overlays contribute less.
        w = L * weights[:, None, None]
        wsum = w.sum(axis=0)
        safe = wsum > eps
        a_out = np.zeros_like(wsum)
        b_out = np.zeros_like(wsum)
        a_out[safe] = (w * a).sum(axis=0)[safe] / wsum[safe]
        b_out[safe] = (w * b).sum(axis=0)[safe] / wsum[safe]
        L_out = _blend_stack(L, blend)
        if inverse:
            # Invert lightness only (white background) without flipping hues.
            L_out = 1. - L_out
        out = ski_color.lab2rgb(np.stack([L_out * 100., a_out, b_out], axis=-1))
    else:
        # HSV/HSL: sum chroma as 2D vectors so distant hues don't interpolate
        # through unrelated intermediates (e.g. red→blue via green).
        conv = rgb_to_hsv if colorspace == 'hsv' else rgb_to_hsl
        back = hsv_to_rgb if colorspace == 'hsv' else hsl_to_rgb
        hsx = [np.asarray(conv(d), dtype=np.float64) for d in disp]
        hue = np.stack([c[..., 0] for c in hsx])
        sat = np.stack([c[..., 1] for c in hsx])
        light = np.stack([c[..., 2] for c in hsx])
        w = light * weights[:, None, None]
        wsum = w.sum(axis=0)
        ang = 2. * np.pi * hue
        X = (w * sat * np.cos(ang)).sum(axis=0)
        Y = (w * sat * np.sin(ang)).sum(axis=0)
        safe = wsum > eps
        s_out = np.zeros_like(wsum)
        s_out[safe] = np.clip(np.hypot(X, Y)[safe] / wsum[safe], 0, 1)
        h_out = (np.arctan2(Y, X) / (2. * np.pi)) % 1.0
        l_out = _blend_stack(light, blend)
        if inverse:
            l_out = 1. - l_out
        out = back(np.stack([h_out, s_out, l_out], axis=-1))

    # Mixing is display-referred (gamma canceled per layer).  Restore a
    # classic-style gamma tone so session/GUI gamma changes the combined look.
    # (Same helper idea as paintmix._apply_relative_compose_gamma; inlined to
    # avoid a colormix↔paintmix import cycle.)
    g = float(gamma) if gamma is not None else 2.2
    ref = 2.2
    if np.isfinite(g) and g > 0 and abs(g - ref) >= 1e-9:
        out = np.power(np.clip(out, 0, 1), g / ref)
    return np.clip(out, 0, 1).astype(dtype)


def mix_colors_hex(hex_list, colorspace='lab', blend='mean', weights=None):
    """
    Mix solid colors (not images) in a chosen color space.

    Useful when designing palettes or predicting how overlapping layer colors
    will look before loading FITS data.

    Parameters
    ----------
    hex_list : list of str
        Hex colors such as ``['#FF0000', '#FFFF00']``.
    colorspace : str
        ``'lab'``, ``'hsv'``, or ``'hsl'`` (same semantics as
        :func:`combine_multicolor_colorspace`).
    blend : str
        Lightness blend: ``'mean'`` (default), ``'screen'``, ``'sum'``, or
        ``'max'``.
    weights : list or None
        Optional per-color weights.

    Returns
    -------
    str
        Mixed color as an uppercase hex string.

    Examples
    --------
    >>> mix_colors_hex(['#FF0000', '#FFFF00'], colorspace='hsv')  # doctest: +SKIP
    '#FF8800'
    """
    ims = [np.array(hex_to_rgb(h)).reshape(1, 1, 3) / 255. for h in hex_list]
    out = combine_multicolor_colorspace(ims, colorspace=colorspace, blend=blend,
                                        weights=weights, gamma=1.0)
    r, g, b = (np.clip(out[0, 0], 0, 1) * 255).round().astype(int)
    return rgb_to_hex((int(r), int(g), int(b))).upper()


# ---------------------------------------------------------------------------
# Accessibility / print checks
# ---------------------------------------------------------------------------

def greyscale_image(rgb, method='lab'):
    """
    Convert a finished RGB composite to greyscale, e.g. to check how it will
    look in a B&W print.

    Parameters
    ----------
    rgb : array
        RGB image (...,3) with values in [0,1]
    method : str
        'lab' (default) uses the CIELAB L* channel (perceptual lightness);
        'luma' uses the Rec.709 luma weights (0.2126R + 0.7152G + 0.0722B).

    Returns
    -------
    array
        Greyscale image with the same (...,3) shape (all channels equal).
    """
    arr = np.clip(np.nan_to_num(np.asarray(rgb, dtype=float)), 0, 1)
    if method == 'lab':
        L = ski_color.rgb2lab(arr)[..., 0] / 100.
    elif method == 'luma':
        L = 0.2126 * arr[..., 0] + 0.7152 * arr[..., 1] + 0.0722 * arr[..., 2]
    else:
        raise ValueError("method must be 'lab' or 'luma'")
    return np.stack([L] * 3, axis=-1)


def _srgb_to_linear(c):
    return np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)


def _linear_to_srgb(c):
    c = np.clip(c, 0, 1)
    return np.where(c <= 0.0031308, 12.92 * c, 1.055 * c ** (1. / 2.4) - 0.055)


# Machado, Oliveira & Fernandes (2009), severity 1.0 simulation matrices,
# applied in linear RGB space.
_CVD_MATRICES = {
    'protanopia': np.array([[0.152286, 1.052583, -0.204868],
                            [0.114503, 0.786281, 0.099216],
                            [-0.003882, -0.048116, 1.051998]]),
    'deuteranopia': np.array([[0.367322, 0.860646, -0.227968],
                              [0.280085, 0.672501, 0.047413],
                              [-0.011820, 0.042940, 0.968881]]),
    'tritanopia': np.array([[1.255528, -0.076749, -0.178779],
                            [-0.078411, 0.930809, 0.147602],
                            [0.004733, 0.691367, 0.303900]]),
}
_CVD_ALIASES = {'protan': 'protanopia', 'deutan': 'deuteranopia', 'tritan': 'tritanopia'}


def simulate_colorblindness(rgb, kind='deuteranopia'):
    """
    Simulate how an RGB composite appears with a color vision deficiency,
    to help pick layer colors that stay distinguishable for all viewers.

    Uses the Machado et al. (2009) severity-1.0 transformation matrices in
    linear RGB space.

    Parameters
    ----------
    rgb : array
        RGB image (...,3) with values in [0,1]
    kind : str
        'protanopia', 'deuteranopia' (default, the most common), or
        'tritanopia'.  Short forms 'protan'/'deutan'/'tritan' also accepted.

    Returns
    -------
    array
        Simulated RGB image, same shape, values in [0,1]
    """
    key = _CVD_ALIASES.get(kind.lower(), kind.lower())
    if key not in _CVD_MATRICES:
        raise ValueError("kind must be one of %s" % (sorted(_CVD_MATRICES),))
    arr = np.clip(np.nan_to_num(np.asarray(rgb, dtype=float)), 0, 1)
    linear = _srgb_to_linear(arr)
    simulated = np.einsum('ij,...j->...i', _CVD_MATRICES[key], linear)
    return np.clip(_linear_to_srgb(np.clip(simulated, 0, 1)), 0, 1)
