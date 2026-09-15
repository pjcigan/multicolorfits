"""
Curated color palettes and palette tools for multi-layer composites.

Choosing well-separated colors is one of the harder parts of building a
multicolor image.  This module provides:

* Curated palettes suited to 2–5 layer astronomical composites (including
  sets used in published multicolorfits examples), plus general hue-geometry
  and color-vision-friendly collections.
* ``suggest_colors(n)`` — *n* perceptually spaced colors on a constant
  lightness / chroma ring in CIE LCh.
* Paletton-style hue patterns (complement, triad, square, split-complement,
  analogous) via ``colors_for_hue_pattern``.
* ``rotate_colors_from_base`` — rotate a palette around the hue wheel while
  keeping relative spacing between layers.
* ``check_palette_colorblind`` — check distinguishability under color-vision
  deficiency simulation.
"""

import numpy as np
from skimage import color as ski_color

from .core.colorize import hex_to_rgb, rgb_to_hex
from .core.colormix import simulate_colorblindness

__all__ = [
    'PALETTES', 'HUE_PATTERNS', 'PALETTE_MENU', 'CVD_KINDS',
    'list_palettes', 'list_hue_patterns', 'list_palette_menu',
    'get_palette', 'suggest_colors', 'colors_from_hue_angles',
    'colors_from_hsv', 'colors_for_hue_pattern', 'hex_to_lch',
    'rotate_colors_from_base',
    'check_palette_colorblind', 'palette_colorblind_report',
    'resolve_palette_colors', 'is_palette_name',
]

# CVD types checked by palette_colorblind_report / session reports, in order.
CVD_KINDS = ('deuteranopia', 'protanopia', 'tritanopia')


# ------------------------------------------------------------------ static palettes

PALETTES = {
    # Classic additive primaries
    'rgb': ['#FF0000', '#00FF00', '#0000FF'],
    # Red-yellow-blue (as in the kepler_RYB example image)
    'ryb': ['#C11B17', '#EAC117', '#2B65EC'],
    # Red-yellow-blue-purple (as in the m101_RYBP example image)
    'rybp': ['#C11B17', '#EAC117', '#4677F0', '#B048B5'],
    # Purple-orange-blue (published multicolorfits examples)
    'pob': ['#BE599E', '#DEA215', '#77C0F9'],
    # Warm sequence for 2-4 layers of similar emission
    'warm': ['#7D0541', '#C11B17', '#E8862C', '#EAC117'],
    # Cool sequence
    'cool': ['#4B0082', '#2B65EC', '#43BFC7', '#8AFB17'],
    # Teal/orange complement pair -- high contrast for 2 layers
    'tealorange': ['#008080', '#E8862C'],
    # Magenta/green complement pair
    'magentagreen': ['#C71585', '#3EA055'],
    # Colorblind-safe (subset of Wong 2011, Nature Methods)
    'cb_safe': ['#D55E00', '#F0E442', '#009E73', '#56B4E9', '#CC79A7'],
    # Paul Tol "bright" qualitative scheme
    'tol': ['#4477AA', '#EE6677', '#228833', '#CCBB44', '#66CCEE', '#AA3377'],
    # IBM Design Library accessible palette
    'ibm': ['#648FFF', '#785EF0', '#DC267F', '#FE6100', '#FFB000'],
}

# Hue-geometry pattern names (resolved dynamically for n layers).
HUE_PATTERNS = ('complement', 'triad', 'split', 'square', 'analogous')

# Menu structure for GUIs: (group label, [(label, key), ...]).
# ``perceptual`` and hue-pattern keys are resolved in ``palette_colors_for``.
PALETTE_MENU = [
    ('Hue geometry', [
        ('Even (auto-N)', 'perceptual'),
        ('Complement (180°)', 'complement'),
        ('Triad (120°)', 'triad'),
        ('Split complement', 'split'),
        ('Square / tetrad (90°)', 'square'),
        ('Analogous (30° steps)', 'analogous'),
    ]),
    ('General moods', [
        ('Pastel ring', 'pastel'),
        ('Vivid ring', 'vivid'),
    ]),
    ('CVD-friendly', [
        ('Wong (cb_safe)', 'cb_safe'),
        ('Tol bright', 'tol'),
        ('IBM accessible', 'ibm'),
    ]),
    ('From examples', [
        ('RGB primaries', 'rgb'),
        ('RYB paint', 'ryb'),
        ('RYBP', 'rybp'),
        ('POB', 'pob'),
        ('Warm sequence', 'warm'),
        ('Cool sequence', 'cool'),
        ('Teal / orange', 'tealorange'),
        ('Magenta / green', 'magentagreen'),
    ]),
]

# Fixed angle offsets (degrees) for Paletton-style patterns (first hue at rotation).
_PATTERN_ANGLES = {
    'complement': [0., 180.],
    'triad': [0., 120., 240.],
    'split': [0., 150., 210.],
    'square': [0., 90., 180., 270.],
    'analogous': [0., 30., 60., 90., 120.],
}


def _init_generated_palettes():
    PALETTES['pastel'] = suggest_colors(5, lightness=78., chroma=38., hue_start=20.)
    PALETTES['vivid'] = suggest_colors(5, lightness=55., chroma=72., hue_start=10.)


def list_palettes():
    """
    Return the names of curated (static) palettes.

    See also :func:`list_palette_menu` for the GUI-oriented listing that
    includes dynamic hue patterns.
    """
    return sorted(PALETTES)


def list_hue_patterns():
    """Names of dynamic hue-geometry patterns."""
    return list(HUE_PATTERNS)


def list_palette_menu():
    """Grouped palette menu entries for GUIs."""
    return list(PALETTE_MENU)


def is_palette_name(name):
    """True if ``name`` is a known palette or hue pattern."""
    return (name in ('perceptual', 'suggest', 'auto')
            or name in PALETTES or name in _PATTERN_ANGLES)


def get_palette(name, n=None):
    """
    Get a curated palette by name.

    Parameters
    ----------
    name : str
        One of list_palettes() or a hue pattern from list_hue_patterns().
    n : int or None
        If given, return only the first n colors (must be <= palette length
        for static palettes).

    Returns
    -------
    list of str
        Hex color strings.
    """
    if name in _PATTERN_ANGLES:
        if n is None:
            n = len(_PATTERN_ANGLES[name])
        return colors_for_hue_pattern(name, n)
    if name not in PALETTES:
        raise KeyError('Unknown palette %r; available: %s' % (name, list_palettes()))
    colors = list(PALETTES[name])
    if n is not None:
        if n > len(colors):
            raise ValueError('Palette %r has only %d colors' % (name, len(colors)))
        colors = colors[:n]
    return colors


def hex_to_lch(hex_color):
    """CIE LCh (degrees) for a hex color."""
    rgb = np.array(hex_to_rgb(hex_color), dtype=float) / 255.
    lab = ski_color.rgb2lab(rgb.reshape(1, 1, 3))[0, 0]
    L, a, b = lab
    C = float(np.hypot(a, b))
    h = float(np.degrees(np.arctan2(b, a)) % 360.)
    return L, C, h


def _lch_to_hex(L, C, h_deg):
    h = np.deg2rad(float(h_deg))
    lab = np.array([[[float(L), C * np.cos(h), C * np.sin(h)]]])
    rgb = np.clip(ski_color.lab2rgb(lab), 0, 1)[0, 0]
    return rgb_to_hex(tuple(int(round(v * 255)) for v in rgb)).upper()


def colors_from_hue_angles(angles_deg, rotation=0., lightness=65., chroma=55.):
    """
    Build hex colors from hue angles (degrees) on a constant L/C ring in CIE LCh.

    Parameters
    ----------
    angles_deg : sequence of float
        Hue offsets in degrees (added to ``rotation``).
    rotation, lightness, chroma : float
        Ring parameters.

    Returns
    -------
    list of str
        Hex color strings.
    """
    angles = np.asarray(angles_deg, dtype=float)
    if angles.size < 1:
        return []
    hues = np.deg2rad(rotation + angles)
    n = len(angles)
    lab = np.stack([np.full(n, float(lightness)),
                    chroma * np.cos(hues),
                    chroma * np.sin(hues)], axis=-1).reshape(1, n, 3)
    rgb = np.clip(ski_color.lab2rgb(lab), 0, 1)[0]
    return [rgb_to_hex(tuple(int(round(v * 255)) for v in c)).upper() for c in rgb]


def suggest_colors(n, lightness=65., chroma=55., hue_start=25.):
    """
    Generate n perceptually evenly-spaced colors: equally spaced hue angles
    on a constant-lightness, constant-chroma ring in CIE LCh space.

    This is what the GUI **Palette → Even (auto-N) / Perceptual** entry uses.
    Because the spacing is perceptual (unlike naive HSV hue stepping), the
    colors remain distinct for any n, and equal lightness means no layer
    visually dominates.

    Parameters
    ----------
    n : int
        Number of colors.
    lightness : float
        CIE L* of the ring [0..100].  Higher = pastel, lower = deep.
    chroma : float
        CIE chroma of the ring.  Higher = more saturated (large values may
        clip at the sRGB gamut boundary for some hues).
    hue_start : float
        Hue angle in degrees of the first color.

    Returns
    -------
    list of str
        Hex color strings.

    See Also
    --------
    colors_from_hue_angles, colors_for_hue_pattern, colors_from_hsv,
    preview_palette
    """
    if n < 1:
        raise ValueError('n must be >= 1')
    angles = hue_start + np.arange(n) * 360. / n
    return colors_from_hue_angles(angles, rotation=0., lightness=lightness, chroma=chroma)


def colors_from_hsv(n, saturation=0.85, value=0.90, hue_start=0.):
    """
    Generate *n* colors with evenly spaced **HSV** hues (classical color wheel).

    Prefer :func:`suggest_colors` (CIE LCh) for most scientific composites —
    HSV steps are not perceptually even (yellows look brighter than blues at
    the same V). Use this when you explicitly want a textbook HSV ramp.

    Parameters
    ----------
    n : int
        Number of colors.
    saturation, value : float
        HSV S and V in ``[0, 1]``.
    hue_start : float
        First hue in degrees (0–360).

    Returns
    -------
    list of str
        Hex color strings.

    Notes
    -----
    This is **not** the same as ``compose.combine_mode = 'hsv'``, which
    controls how already-chosen layer colors are *mixed* in the composite.
    """
    if n < 1:
        raise ValueError('n must be >= 1')
    from .core.colorize import hsv_to_rgb as _hsv_to_rgb
    hues = (float(hue_start) + np.arange(n) * (360.0 / n)) % 360.0
    hsv = np.zeros((n, 3), dtype=float)
    hsv[:, 0] = hues / 360.0
    hsv[:, 1] = float(saturation)
    hsv[:, 2] = float(value)
    rgb = np.clip(_hsv_to_rgb(hsv.reshape(1, n, 3))[0], 0, 1)
    return [rgb_to_hex(tuple(int(round(v * 255)) for v in c)).upper() for c in rgb]


def colors_for_hue_pattern(pattern, n, rotation=0., lightness=65., chroma=55.):
    """
    Colors for a Paletton-style hue relationship, adapted to ``n`` layers.

    When more layers are loaded than the pattern defines, extra hues are
    inserted at even intervals around the wheel (after the pattern anchors).
    """
    if n < 1:
        return []
    if pattern in ('perceptual', 'suggest', 'auto', 'even'):
        return suggest_colors(n, lightness=lightness, chroma=chroma, hue_start=rotation)
    if pattern not in _PATTERN_ANGLES:
        raise KeyError('Unknown hue pattern %r' % pattern)
    base = list(_PATTERN_ANGLES[pattern])
    if n <= len(base):
        return colors_from_hue_angles(base[:n], rotation=rotation,
                                      lightness=lightness, chroma=chroma)
    extras = list(np.linspace(0., 360., n - len(base) + 1, endpoint=False))[1:]
    merged = (base + extras)[:n]
    return colors_from_hue_angles(merged, rotation=rotation,
                                  lightness=lightness, chroma=chroma)


def rotate_colors_from_base(base_hex_colors, delta_deg):
    """
    Rotate layer colors around the hue wheel, preserving each color's L/C and
    its angular offset from the first layer.

    Parameters
    ----------
    base_hex_colors : list of str
        Starting hex colors (one per layer).
    delta_deg : float
        Rotation to apply in degrees.

    Returns
    -------
    list of str
        Rotated hex colors.
    """
    if not base_hex_colors:
        return []
    lchs = [hex_to_lch(h) for h in base_hex_colors]
    h0 = lchs[0][2]
    offsets = [(h - h0) % 360. for (_, _, h) in lchs]
    new_h0 = (h0 + float(delta_deg)) % 360.
    return [_lch_to_hex(L, C, (new_h0 + off) % 360.) for (L, C, _), off in zip(lchs, offsets)]


def resolve_palette_colors(name_or_colors, n=None):
    """
    Resolve a palette name or hex list to ``n`` hex colors.

    Same rules as :meth:`~multicolorfits.McfSession.palette_colors_for`:
    ``'perceptual'`` / ``'suggest'`` / ``'auto'`` use :func:`suggest_colors`;
    hue-pattern names use :func:`colors_for_hue_pattern`; curated names use
    :func:`get_palette`, padded with perceptual suggestions when shorter than
    *n*. A sequence of hex strings is returned as-is (truncated or rejected
    if *n* is set and the length does not match, except a single color
    broadcasts).

    Parameters
    ----------
    name_or_colors : str or sequence of str
        Palette / pattern name, or an explicit hex list.
    n : int or None
        Desired length. Required when *name_or_colors* is a name. For an
        explicit list, defaults to ``len(colors)``.

    Returns
    -------
    list of str
        Hex color strings.
    """
    if isinstance(name_or_colors, (str, bytes)):
        name = str(name_or_colors)
        if n is None:
            if name in ('perceptual', 'suggest', 'auto', 'even'):
                raise ValueError('n is required for perceptual / suggest / auto / even')
            if name in _PATTERN_ANGLES:
                n = len(_PATTERN_ANGLES[name])
            else:
                n = len(get_palette(name))
        n = int(n)
        if n < 1:
            return []
        if name in ('perceptual', 'suggest', 'auto'):
            return suggest_colors(n)
        if name in ('complement', 'triad', 'split', 'square', 'analogous', 'even'):
            return colors_for_hue_pattern(name, n)
        avail = get_palette(name)
        if len(avail) >= n:
            return list(avail[:n])
        return list(avail) + suggest_colors(n)[len(avail):]

    colors = [str(c) for c in list(name_or_colors)]
    if not colors:
        return []
    if n is None:
        return colors
    n = int(n)
    if n < 1:
        return []
    if len(colors) == 1 and n != 1:
        return colors * n
    if len(colors) != n:
        raise ValueError('colors has %d entries; expected %d or 1'
                         % (len(colors), n))
    return colors


def check_palette_colorblind(hex_colors, kind='deuteranopia', min_distance=25.):
    """
    Check whether a palette's colors remain distinguishable under a color
    vision deficiency, by simulating the CVD and measuring pairwise CIELAB
    distances.

    Parameters
    ----------
    hex_colors : list of str
        Palette to check.
    kind : str
        'protanopia', 'deuteranopia' (default), or 'tritanopia'.
    min_distance : float
        Minimum acceptable pairwise Delta-E (CIE76).  ~25 is comfortably
        distinguishable for large image regions.

    Returns
    -------
    (bool, list)
        (all pairs pass, list of (i, j, distance) for failing pairs)
    """
    rgb = np.array([np.array(hex_to_rgb(h)) / 255. for h in hex_colors]).reshape(1, -1, 3)
    sim = simulate_colorblindness(rgb, kind=kind)
    lab = ski_color.rgb2lab(sim)[0]
    failures = []
    for i in range(len(hex_colors)):
        for j in range(i + 1, len(hex_colors)):
            dist = float(np.linalg.norm(lab[i] - lab[j]))
            if dist < min_distance:
                failures.append((i, j, dist))
    return (len(failures) == 0), failures


def palette_colorblind_report(hex_colors, kinds=None, min_distance=25.):
    """
    Multi-kind CVD report for a hex list (no session required).

    Parameters
    ----------
    hex_colors : sequence of str
    kinds : sequence of str or None
        Defaults to :data:`CVD_KINDS` (deuteranopia, protanopia, tritanopia).
    min_distance : float
        Passed to :func:`check_palette_colorblind`.

    Returns
    -------
    dict
        ``ok`` (bool), ``colors`` (list), ``kinds`` mapping each kind to
        ``{'ok': bool, 'failures': [(i, j, dE), ...]}`` with 0-based indices.
    """
    colors = [str(c) for c in list(hex_colors)]
    kind_list = tuple(kinds) if kinds is not None else CVD_KINDS
    report = {'ok': True, 'colors': colors, 'kinds': {}}
    if len(colors) < 2:
        for kind in kind_list:
            report['kinds'][kind] = {'ok': True, 'failures': []}
        return report
    for kind in kind_list:
        ok, failures = check_palette_colorblind(
            colors, kind=kind, min_distance=min_distance)
        mapped = [(int(a), int(b), round(float(d), 1)) for a, b, d in failures]
        report['kinds'][kind] = {'ok': ok, 'failures': mapped}
        if not ok:
            report['ok'] = False
    return report


_init_generated_palettes()
