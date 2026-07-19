"""
Composite-mode-aware color-combination swatch (overlapping tinted circles) plus
the shared layer-combining dispatch used by both the real image and the swatch.

The swatch renders one circle per channel in its color, arranged so they overlap,
with each circle a synthetic full-signal "band" image pushed through the *exact
same* colorize + combine pipeline as the displayed image.  Overlap regions
therefore show the honest mixed colors for the current mode (additive RGB,
subtractive RYB paint, Lab, HSV), background, gamma, and inverse setting -- like
Figure 1 of Sugita & Takahashi (2017).
"""

import numpy as np

from .scaling import adjust_gamma
from .colorize import colorize_image, hex_complement, combine_multicolor
from .colormix import combine_multicolor_colorspace
from .paintmix import combine_multicolor_alpha, composite_over_background

__all__ = ['combine_colorized_layers', 'swatch_layout', 'combo_swatch']


def combine_colorized_layers(colorized, mode='rgb', blend='screen', gamma=2.2,
                             background='black', inverse=False, dtype=np.float64):
    """
    Combine a list of already-colorized RGB layers into a final image, applying a
    compositing mode / blend / gamma / inverse / background.

    Single source of truth for the compositing math, shared by the session's
    ``render_combined`` (real data) and ``combo_swatch`` (synthetic circles), so
    the legend's overlap colors always match the displayed image.

    Parameters
    ----------
    colorized : list of arrays
        Colorized RGB layers (as from ``colorize_image``).
    mode : str
        'rgb' (additive), 'ryb' (subtractive paint), 'cmyk' (print ink), or
        'lab'/'hsv'/'hsl'.
    blend : str
        Lightness blend for the lab/hsv/hsl modes ('screen'/'sum'/'max'/'mean').
    gamma : float
    background : color or None or 'black'
        'black' => classic additive (no flatten); any other color composites the
        result over it; None (or 'transparent'/'none') => transparent RGBA output.
    inverse : bool
        Legacy complementary-color inverse (ignored when a non-black background is used).
    dtype : numpy dtype

    Returns
    -------
    array
        Combined RGB [...,3] (or RGBA [...,4] for a transparent background).
    """
    bg_key = str(background).strip().lower() if background is not None else 'transparent'
    use_bg = bg_key not in ('black', '')
    bg_arg = None if (background is None or bg_key in ('transparent', 'none')) else background
    eff_inverse = False if use_bg else inverse

    if mode == 'ryb':
        return combine_multicolor_alpha(colorized, background=bg_arg, mode='ryb',
                                        gamma=gamma, dtype=dtype)
    if mode == 'cmyk':
        return combine_multicolor_alpha(colorized, background=bg_arg, mode='cmyk',
                                        gamma=gamma, dtype=dtype)
    if mode == 'rgb':
        combined = combine_multicolor(colorized, gamma=gamma, inverse=eff_inverse,
                                      dtype=dtype)
    else:
        combined = combine_multicolor_colorspace(colorized, colorspace=mode, blend=blend,
                                                 gamma=gamma, inverse=eff_inverse, dtype=dtype)
    if use_bg:
        combined = composite_over_background(combined, background=bg_arg).astype(dtype)
    return combined


def swatch_layout(n, size):
    """
    Circle centers/radii (pixel coords of a ``size`` x ``size`` canvas) for the
    overlapping color-combination swatch.

    1 -> single disc; 2 -> side-by-side pair; 3 -> classic three-circle
    overlap; N>=4 -> a rosette of mutually overlapping circles.
    Centers stay within a small frame margin so nothing is clipped.
    """
    cx0 = cy0 = size / 2.0
    if n <= 0:
        return []
    if n == 1:
        return [(cx0, cy0, 0.42 * size)]
    r = {2: 0.30, 3: 0.30}.get(n, max(0.16, 0.34 - 0.02 * n)) * size
    # d < r keeps every circle over the centre (a true all-channel overlap);
    # ~1/sqrt(3) gives the tidy symmetric three-circle layout.
    d = (0.58 if n == 3 else 0.60) * r
    maxreach = 0.47 * size
    if d + r > maxreach:
        scale = maxreach / (d + r)
        r *= scale
        d *= scale
    if n == 2:
        return [(cx0 - d, cy0, r), (cx0 + d, cy0, r)]
    layout = []
    for k in range(n):
        ang = np.pi / 2.0 - 2.0 * np.pi * k / n   # first circle at top, clockwise
        layout.append((cx0 + d * np.cos(ang), cy0 + d * np.sin(ang), r))
    return layout


def combo_swatch(colors, mode='rgb', blend='screen', gamma=2.2, background='black',
                 inverse=False, size=320, labels=None):
    """
    Render the composite-mode-aware color-combination swatch.

    Parameters
    ----------
    colors : list of str
        Channel colors (hex strings), one per circle.
    mode, blend, gamma, background, inverse : see ``combine_colorized_layers``.
    size : int
        Pixel size of the (square) swatch image.
    labels : list of str or None
        Optional per-circle labels (stored in the returned geometry).

    Returns
    -------
    dict or None
        ``{'rgba': (size,size,4) float array in [0,1], 'circles': [...],
        'size': size}`` where each circle is ``{'x','y','r','color','label'}``
        in swatch pixel coords (origin lower-left).  None if ``colors`` is empty.
    """
    colors = list(colors)
    n = len(colors)
    if n == 0:
        return None
    if labels is None:
        labels = ['' for _ in colors]
    layout = swatch_layout(n, size)

    yy, xx = np.mgrid[0:size, 0:size].astype(float)
    masks = []
    for (cx, cy, r) in layout:
        dist = np.hypot(xx - cx, yy - cy)
        masks.append(np.clip(r - dist + 0.5, 0.0, 1.0))  # ~1px anti-aliased edge
    coverage = masks[0] if n == 1 else np.max(masks, axis=0)

    # Colorize each circle mask exactly like a real panel layer, then combine.
    bg_key = str(background).strip().lower() if background is not None else 'transparent'
    eff_inverse = inverse if bg_key in ('black', '') else False
    colorized = []
    for mask, color in zip(masks, colors):
        grey_rgb = np.stack([adjust_gamma(mask, gamma)] * 3, axis=-1)
        c = hex_complement(color) if eff_inverse else color
        colorized.append(colorize_image(grey_rgb, c, colorintype='hex',
                                         gammacorr_color=gamma))
    combined = combine_colorized_layers(colorized, mode=mode, blend=blend, gamma=gamma,
                                        background=background, inverse=inverse)
    rgb = combined[..., :3] if combined.shape[-1] == 4 else combined
    rgb = np.clip(np.nan_to_num(rgb), 0.0, 1.0)
    rgba = np.dstack([rgb, np.clip(coverage, 0.0, 1.0)])

    circles = [{'x': cx, 'y': cy, 'r': r, 'color': colors[k], 'label': labels[k]}
               for k, (cx, cy, r) in enumerate(layout)]
    return {'rgba': rgba, 'circles': circles, 'size': size}
