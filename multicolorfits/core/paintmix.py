"""
Subtractive color spaces (RYB, CMYK) and background / transparency compositing.

This module complements the RGB / HSV / HSL / Lab mixers in ``colormix.py``.
RYB and CMYK modes, plus the background control, are available in both GUIs
and via :func:`combine_multicolor_alpha` / :func:`combine_colorized_layers`.

1. **RYB (red–yellow–blue)** — a subtractive "artist's wheel" model following
   Sugita & Takahashi (2015/2017).  Black RGB maps to RYB ``(1, 1, 1)`` (full
   pigment) and white RGB maps to ``(0, 0, 0)`` (blank canvas).  Mixing in RYB
   reproduces familiar paint results (for example yellow + blue → green) rather
   than additive light (yellow + blue → white).

2. **CMYK (cyan–magenta–yellow–black)** — a print-ink model using standard
   device-independent RGB ↔ CMYK transforms (no UCR/GCR).  Overlapping layers
   accumulate ink coverage, then convert back to RGB.  Hue mixes differ from
   RYB (for example yellow + blue → black, not green).

3. **Background compositing.**  Classic ``combine_multicolor`` places the
   result on a black background.  The legacy ``inverse=True`` option applies
   ``1 - image``, which also replaces each hue with its complement.  A clearer
   approach treats each colorized layer as premultiplied (color × intensity)
   and composites it over a chosen background — white for publications, a
   custom color for slides, or ``None`` for a transparent PNG that fades to
   empty where there is no signal.  Bright regions keep the same additive
   colors as black-background RGB; only faint regions fade to the background.

All functions are pure NumPy.  Entry points are re-exported at the package root
(``import multicolorfits as mcf``).
"""

import numpy as np

from .colorize import hex_to_rgb, rgb_to_hex

__all__ = [
    'rgb_to_ryb',
    'ryb_to_rgb',
    'hex_to_ryb',
    'ryb_to_hex',
    'rgb_to_cmyk',
    'cmyk_to_rgb',
    'hex_to_cmyk',
    'cmyk_to_hex',
    'layer_coverage',
    'composite_over_background',
    'combine_multicolor_alpha',
    'ALPHA_SPACES',
]

ALPHA_SPACES = ('rgb', 'ryb', 'cmyk')  # valid ``mode`` values for combine_multicolor_alpha


# ---------------------------------------------------------------------------
# RYB <-> RGB (Sugita & Takahashi 2015/2017), fully vectorized
# ---------------------------------------------------------------------------

def rgb_to_ryb(rgb):
    """
    Vectorized sRGB -> RYB conversion (Sugita & Takahashi 2017).

    Parameters
    ----------
    rgb : numpy.ndarray
        Array with last axis of length 3, values in [0, 1].

    Returns
    -------
    numpy.ndarray
        RYB array of the same shape, values in [0, 1].  Note the subtractive
        convention: RGB black -> RYB (1,1,1); RGB white -> RYB (0,0,0).
    """
    rgb = np.clip(np.asarray(rgb, dtype=np.float64), 0.0, 1.0)
    r, g, b = rgb[..., 0], rgb[..., 1], rgb[..., 2]

    # Sugita & Takahashi: peel white/black so the remaining RGB is "pure" pigment.
    white = np.minimum(np.minimum(r, g), b)   # shared white component
    black = 1.0 - np.maximum(np.maximum(r, g), b)   # shared black component

    r_ = r - white
    g_ = g - white
    b_ = b - white

    # Map the residual RGB pigment into RYB axes (red, yellow, blue).
    min_rg = np.minimum(r_, g_)
    ry = r_ - min_rg
    yy = (g_ + min_rg) / 2.0
    by = (g_ + b_ - min_rg) / 2.0

    ryb = np.stack([ry, yy, by], axis=-1)
    max_rgb = np.maximum(np.maximum(r_, g_), b_)
    max_ryb = np.max(ryb, axis=-1)

    # Rescale so peak RYB intensity matches the residual RGB peak, then add black.
    with np.errstate(divide='ignore', invalid='ignore'):
        n = np.where(max_ryb > 0, max_ryb / max_rgb, 0.0)
        scale = np.where(n > 0, 1.0 / np.where(n > 0, n, 1.0), 1.0)
    ryb = ryb * scale[..., None] + black[..., None]
    return np.clip(ryb, 0.0, 1.0)


def ryb_to_rgb(ryb):
    """
    Vectorized RYB -> sRGB conversion (Sugita & Takahashi 2017).

    Parameters
    ----------
    ryb : numpy.ndarray
        Array with last axis of length 3, values in [0, 1].

    Returns
    -------
    numpy.ndarray
        sRGB array of the same shape, values in [0, 1].
    """
    ryb = np.clip(np.asarray(ryb, dtype=np.float64), 0.0, 1.0)
    R, Y, B = ryb[..., 0], ryb[..., 1], ryb[..., 2]

    black = np.minimum(np.minimum(R, Y), B)   # I_b
    white = 1.0 - np.maximum(np.maximum(R, Y), B)   # I_w

    r_ = R - black
    y_ = Y - black
    b_ = B - black

    min_yb = np.minimum(y_, b_)
    rr = r_ + y_ - min_yb
    gg = y_ + min_yb
    bb = 2.0 * (b_ - min_yb)

    rgb = np.stack([rr, gg, bb], axis=-1)
    max_ryb = np.maximum(np.maximum(r_, y_), b_)
    max_rgb = np.max(rgb, axis=-1)

    with np.errstate(divide='ignore', invalid='ignore'):
        n = np.where(max_rgb > 0, max_rgb / max_ryb, 0.0)
        scale = np.where(n > 0, 1.0 / np.where(n > 0, n, 1.0), 1.0)
    rgb = rgb * scale[..., None] + white[..., None]
    return np.clip(rgb, 0.0, 1.0)


def hex_to_ryb(hexstring):
    """
    Convert a hex color string to an RYB triplet in [0, 1].

    Parameters
    ----------
    hexstring : str
        Color such as ``'#FF0000'``.

    Returns
    -------
    numpy.ndarray
        Shape ``(3,)`` RYB values (subtractive convention).
    """
    rgb = np.array(hex_to_rgb(hexstring), dtype=float) / 255.
    return rgb_to_ryb(rgb)


def ryb_to_hex(ryb):
    """
    Convert an RYB triplet in [0, 1] to a hex color string.

    Parameters
    ----------
    ryb : array-like
        Length-3 RYB values.

    Returns
    -------
    str
        Uppercase hex color such as ``'#FF8800'``.
    """
    rgb = ryb_to_rgb(np.asarray(ryb, dtype=float))
    vals = tuple(int(round(c * 255)) for c in np.clip(rgb, 0, 1))
    return rgb_to_hex(vals).upper()


# ---------------------------------------------------------------------------
# CMYK <-> RGB (standard device-independent transforms), fully vectorized
# ---------------------------------------------------------------------------

def rgb_to_cmyk(rgb):
    """
    Vectorized sRGB -> CMYK conversion (standard, no UCR/GCR).

    Parameters
    ----------
    rgb : numpy.ndarray
        Array with last axis of length 3, values in [0, 1].

    Returns
    -------
    numpy.ndarray
        CMYK array of the same shape with last axis length 4, values in [0, 1].
        White RGB -> (0,0,0,0); black RGB -> (0,0,0,1).
    """
    rgb = np.clip(np.asarray(rgb, dtype=np.float64), 0.0, 1.0)
    r, g, b = rgb[..., 0], rgb[..., 1], rgb[..., 2]
    k = 1.0 - np.maximum(np.maximum(r, g), b)
    denom = np.maximum(1.0 - k, 1e-12)
    c = np.where(k < 1.0 - 1e-12, (1.0 - r - k) / denom, 0.0)
    m = np.where(k < 1.0 - 1e-12, (1.0 - g - k) / denom, 0.0)
    y = np.where(k < 1.0 - 1e-12, (1.0 - b - k) / denom, 0.0)
    return np.stack([c, m, y, k], axis=-1)


def cmyk_to_rgb(cmyk):
    """
    Vectorized CMYK -> sRGB conversion.

    Parameters
    ----------
    cmyk : numpy.ndarray
        Array with last axis of length 4, values in [0, 1].

    Returns
    -------
    numpy.ndarray
        sRGB array of the same shape with last axis length 3, values in [0, 1].
    """
    cmyk = np.clip(np.asarray(cmyk, dtype=np.float64), 0.0, 1.0)
    c, m, y, k = cmyk[..., 0], cmyk[..., 1], cmyk[..., 2], cmyk[..., 3]
    one_minus_k = 1.0 - k
    rgb = np.stack([(1.0 - c) * one_minus_k,
                    (1.0 - m) * one_minus_k,
                    (1.0 - y) * one_minus_k], axis=-1)
    return np.clip(rgb, 0.0, 1.0)


def hex_to_cmyk(hexstring):
    """
    Convert a hex color string to a CMYK quadruplet in [0, 1].

    Parameters
    ----------
    hexstring : str
        Color such as ``'#00FFFF'``.

    Returns
    -------
    numpy.ndarray
        Shape ``(4,)`` — cyan, magenta, yellow, key (black).
    """
    rgb = np.array(hex_to_rgb(hexstring), dtype=float) / 255.
    return rgb_to_cmyk(rgb)


def cmyk_to_hex(cmyk):
    """
    Convert a CMYK quadruplet in [0, 1] to a hex color string.

    Parameters
    ----------
    cmyk : array-like
        Length-4 CMYK values (key/black is folded into RGB).

    Returns
    -------
    str
        Uppercase hex color.
    """
    rgb = cmyk_to_rgb(np.asarray(cmyk, dtype=float))
    vals = tuple(int(round(c * 255)) for c in np.clip(rgb, 0, 1))
    return rgb_to_hex(vals).upper()


# ---------------------------------------------------------------------------
# Alpha / background compositing
# ---------------------------------------------------------------------------

def _to_rgb_triplet(color):
    """Parse 'white' / '#rrggbb' / (r,g,b) -> float RGB triplet in [0,1]."""
    from matplotlib.colors import to_rgb
    return np.array(to_rgb(color), dtype=np.float64)


def layer_coverage(rgb):
    """
    Per-pixel coverage (HSV-style "value") of a colorized layer.

    Returns the maximum of the R, G, and B channels.  For a single-hue layer
    produced by :func:`~multicolorfits.core.colorize.colorize_image` this is
    zero on the black background and one at full colored signal — the natural
    alpha for premultiplied compositing.
    """
    arr = np.clip(np.nan_to_num(np.asarray(rgb, dtype=np.float64)), 0.0, 1.0)
    return arr.max(axis=-1)


def composite_over_background(rgb, background=(1., 1., 1.), alpha=None):
    """
    Composite a black-background (premultiplied-alpha) colorized/combined image
    over a solid background color, or produce a transparent RGBA image.

    This is the clean "inverse" (white-background) operation: a normal additive
    composite is already premultiplied alpha (color x value with black where
    there is no signal), so bright regions keep exactly their RGB colors while
    faint / empty regions fade to ``background`` instead of to black.

    Parameters
    ----------
    rgb : array
        Combined/colorized RGB image, shape (...,3), values in [0,1], on a black
        background (as returned by ``combine_multicolor(..., inverse=False)`` or
        ``colorize_image``).
    background : color or None
        Any matplotlib color ('white', '#204060', (0.1,0.1,0.1), ...) to
        composite onto, or ``None`` to return a transparent RGBA image
        (shape (...,4)) that fades to nothing at zero signal.
    alpha : array or None
        Optional explicit per-pixel coverage in [0,1].  Defaults to the image
        'value' (max channel), which is the correct premultiplied coverage for
        colorized layers.

    Returns
    -------
    array
        RGB image (...,3) when ``background`` is a color, or RGBA (...,4) when
        ``background is None``.
    """
    premult = np.clip(np.nan_to_num(np.asarray(rgb, dtype=np.float64)), 0.0, 1.0)
    a = layer_coverage(premult) if alpha is None else \
        np.clip(np.nan_to_num(np.asarray(alpha, dtype=np.float64)), 0.0, 1.0)

    if background is None:
        # Un-premultiply to straight color for a proper transparent PNG.
        with np.errstate(divide='ignore', invalid='ignore'):
            straight = np.where(a[..., None] > 0, premult / a[..., None], 0.0)
        return np.concatenate([np.clip(straight, 0, 1), a[..., None]], axis=-1)

    bg = _to_rgb_triplet(background)
    out = premult + (1.0 - a)[..., None] * bg
    return np.clip(out, 0.0, 1.0)


def _apply_relative_compose_gamma(rgb, gamma, ref_gamma=2.2):
    """
    Restore a creative gamma control when mixing was done in display-referred
    space (where per-layer ``**(1/gamma)`` cancels the build gamma).

    Classic ``combine_multicolor`` mixes in gamma-space then applies
    ``**(1/gamma)`` once, so changing gamma changes the look.  RYB / CMYK /
    Lab undo gamma *before* mixing, which made the GUI gamma spinner a no-op.
    Applying ``**(gamma/ref)`` after the mix restores the same direction
    (higher gamma → darker midtones) while leaving ``gamma=ref`` unchanged.
    """
    if gamma is None:
        return rgb
    g = float(gamma)
    ref = float(ref_gamma)
    if not np.isfinite(g) or g <= 0 or abs(g - ref) < 1e-9:
        return rgb
    arr = np.asarray(rgb, dtype=np.float64)
    out = np.clip(arr, 0.0, 1.0)
    if out.ndim >= 1 and out.shape[-1] == 4:
        rgb3 = np.power(out[..., :3], g / ref)
        return np.concatenate([rgb3, out[..., 3:4]], axis=-1)
    return np.power(out, g / ref)


def combine_multicolor_alpha(im_list_colorized, background='white', mode='rgb',
                             gamma=2.2, weights=None, dtype=np.float64):
    """
    Combine colorized layers by compositing over a chosen background.

    Drop-in alternative to :func:`~multicolorfits.core.colorize.combine_multicolor`
    and :func:`~multicolorfits.core.colormix.combine_multicolor_colorspace` when
    you want:

    * a white (or custom) background without the legacy ``1 - image`` inverse
      (which also flips hues);
    * paint-like mixing with ``mode='ryb'`` (yellow + blue → green);
    * print-ink mixing with ``mode='cmyk'``;
    * transparent PNG output with ``background=None`` (RGBA that fades out at
      zero signal).

    Parameters
    ----------
    im_list_colorized : list
        Colorized RGB layers from :func:`~multicolorfits.core.colorize.colorize_image`
        (black background).
    background : color or None
        Solid background (``'white'``, hex, name, or RGB tuple), or ``None``
        for a transparent RGBA result.
    mode : str
        ``'rgb'`` (additive), ``'ryb'`` (subtractive / paint), or ``'cmyk'``
        (print ink).  Same keyword vocabulary as ``combine_colorized_layers``
        and ``ComposeState.combine_mode``.
    gamma : float
        Gamma used when layers were built; values are converted to
        display-referred space (``** (1/gamma)``) before compositing.
        Use ``1.0`` if inputs are already display-ready.
    weights : list or None
        Optional per-layer weights on each layer's coverage / pigment.
    dtype : numpy dtype
        Output dtype.

    Returns
    -------
    array
        RGB ``(..., 3)`` if ``background`` is a color, else RGBA ``(..., 4)``.

    Examples
    --------
    >>> out = combine_multicolor_alpha(  # doctest: +SKIP
    ...     [col_ir, col_opt], background='white', mode='ryb', gamma=2.2)
    """
    if mode not in ALPHA_SPACES:
        raise ValueError("mode must be one of %s" % (ALPHA_SPACES,))
    n = len(im_list_colorized)
    if n == 0:
        raise ValueError('im_list_colorized is empty')
    if weights is None:
        weights = np.ones(n)
    weights = np.asarray(weights, dtype=float)
    if weights.shape != (n,):
        raise ValueError('weights must have one value per image')

    # Display-referred layers, black background.
    disp = []
    for im in im_list_colorized:
        arr = np.clip(np.nan_to_num(np.asarray(im, dtype=np.float64)), 0.0, 1.0)
        disp.append(arr ** (1. / gamma) if gamma != 1 else arr)

    eps = 1e-9
    # Per-layer coverage (value) and straight (full-brightness) color.
    a_list, c_list = [], []
    for d, w in zip(disp, weights):
        val = layer_coverage(d)
        a = np.clip(val * w, 0.0, 1.0)
        with np.errstate(divide='ignore', invalid='ignore'):
            c = np.where(val[..., None] > eps,
                         d / np.maximum(val[..., None], eps), 0.0)
        a_list.append(a)
        c_list.append(np.clip(c, 0.0, 1.0))

    a_stack = np.stack(a_list, axis=0)
    c_stack = np.stack(c_list, axis=0)
    # Combined coverage via "screen" (order-independent, saturates gracefully).
    alpha_tot = 1.0 - np.prod(1.0 - a_stack, axis=0)

    if mode == 'rgb':
        # Additive premultiplied colors (identical hues to combine_multicolor).
        premult = np.clip(np.sum(a_stack[..., None] * c_stack, axis=0), 0.0, 1.0)
    else:
        # Subtractive: sum pigments in RYB or CMYK, then convert back to RGB.
        to_ink = rgb_to_ryb if mode == 'ryb' else rgb_to_cmyk
        from_ink = ryb_to_rgb if mode == 'ryb' else cmyk_to_rgb
        pigment = np.sum(a_stack[..., None] * to_ink(c_stack), axis=0)
        paint_white = from_ink(np.clip(pigment, 0.0, 1.0))
        # paint_white is the color as seen on white; premultiply by coverage so
        # it composites correctly over an arbitrary background / transparency.
        premult = paint_white * alpha_tot[..., None]

    out = composite_over_background(premult, background=background, alpha=alpha_tot)
    # RGB additive already encodes gamma via the layer build + this path's
    # display convert; subtractive modes cancel gamma during the mix — restore
    # a classic-style gamma tone so the GUI/session gamma control has effect.
    if mode in ('ryb', 'cmyk'):
        out = _apply_relative_compose_gamma(out, gamma)
    return out.astype(dtype)
