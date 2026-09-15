"""
Intensity scaling / stretching utilities.

These wrap the astropy.visualization stretch classes and provide NaN-safe
helpers used throughout multicolorfits.
"""

import sys

import matplotlib.colors as mcolors
import numpy as np
from astropy.visualization import (
    LinearStretch,
    SqrtStretch,
    SquaredStretch,
    LogStretch,
    PowerDistStretch,
    SinhStretch,
    AsinhStretch,
    ManualInterval,
    ZScaleInterval,
)
from scipy.stats import percentileofscore

__all__ = [
    'stretch_functions',
    'SIGNED_STRETCHES',
    'adjust_gamma',
    'draw_progress_bar',
    'nan_percentile_of_score',
    'make_norm',
    'rescale_image',
    'zscale_limits',
    'suggest_levels',
    'describe_image',
    'describe_images',
]

class _StretchFunctions(dict):
    """Name → astropy stretch class for unsigned ``rescale_image`` / ``make_norm``.

    Keys: ``linear``, ``sqrt``, ``squared``, ``log``, ``power``, ``sinh``,
    ``asinh``.  Signed stretches use :data:`~multicolorfits.SIGNED_STRETCHES`
    instead.
    """


stretch_functions = _StretchFunctions({
    'linear': LinearStretch,
    'sqrt': SqrtStretch,
    'squared': SquaredStretch,
    'log': LogStretch,
    'power': PowerDistStretch,
    'sinh': SinhStretch,
    'asinh': AsinhStretch,
})

#: Signed-data stretch names routed through matplotlib / optional pysymlog norms.
SIGNED_STRETCHES = ('symlog', 'symmetric_log')


def adjust_gamma(array_in, gamma):
    """
    Replacement function for skimage.exposure.adjust_gamma, so that NaNs don't throw errors

    Parameters
    ----------
    array_in : array
        Input image array
    gamma : float
        Gamma correction value

    Returns
    -------
    array
        Gamma-adjusted image values
    """
    return array_in**(float(gamma))


def draw_progress_bar(percent, barlength=20, prefix='', suffix=''):
    """
    Dependency-free ANSI progress bar written to stdout.

    This is the built-in helper used by
    :func:`~multicolorfits.reproject_cube` when ``print_progress=True``.
    It needs no third-party packages and stays in the public API
    (including the legacy alias ``drawProgressBar``) for scripts and
    back-compat.  For richer bars in your own loops, use
    `tqdm <https://tqdm.github.io/>`_ directly — multicolorfits does
    not wrap or depend on it.

    Parameters
    ----------
    percent : float
        Completion fraction in ``[0, 1]`` (values outside that range are
        still rendered; callers usually clamp).
    barlength : int, optional
        Number of character cells in the bar body (default 20).
    prefix, suffix : str, optional
        Text printed immediately before / after the ``[====--] xx.xx%``
        segment (e.g. a label and a ``\"i of n\"`` counter).

    Notes
    -----
    Each call overwrites the current line with ``\\r``.  Pass a trailing
    newline in ``suffix`` (as ``reproject_cube`` does on the final update)
    to leave the finished bar on screen.

    Examples
    --------
    Built-in style (percent updates)::

        import multicolorfits as mcf
        n = 10
        for i in range(n):
            # ... do work ...
            mcf.draw_progress_bar((i + 1) / n, prefix='Working ',
                                  suffix='  %d of %d' % (i + 1, n))
        mcf.draw_progress_bar(1.0, prefix='Working ',
                              suffix='  %d of %d\\n' % (n, n))

    Spectral cube reprojection (uses this helper internally)::

        cube_out = mcf.reproject_cube(
            cube, hdr_from, hdr_to, print_progress=True)

    Same cube loop with tqdm instead (optional third-party)::

        from tqdm import tqdm
        planes = []
        for z in tqdm(range(cube.shape[-3]), desc='Reprojecting'):
            planes.append(
                mcf.reproject_image(cube[z], hdr_from, hdr_to))
        cube_out = np.array(planes)

    See Also
    --------
    reproject_cube : set ``print_progress=True`` to drive this bar
        per spectral plane.
    """
    sys.stdout.write("\r")
    ANSIred = '\033[31m'; ANSIgreen = '\033[32m'; ANSIyellow = '\033[33m'; ANSIblue = '\033[34m'; ANSIreset = '\033[0m'
    progress = '' + ANSIgreen  # Start with ANSI escape code for Green
    for i in range(barlength):
        if i < int(barlength * percent):
            progress += "="  # print (green) double dashes
        else:
            progress += ANSIblue + '-'  # Change to ANSI escape code for blue, print single dashes
    progress += ANSIreset  # Reset style
    sys.stdout.write(prefix + "[ %s ] %s%.2f%%%s" % (progress, ANSIyellow, percent * 100, ANSIreset) + suffix)
    sys.stdout.flush()


def nan_percentile_of_score(array_in, score, **kwargs):
    """
    Usage is identical to scipy.stats.percentileofscore(array_in,score), this just corrects for NaNs

    Parameters
    ----------
    array_in : array
        input dataset
    score : float
        score for calculation
    **kwargs
        keyword arguments to pass to scipy.stats.percentileofscore()

    Returns
    -------
    float, or array the length of score
        percentile of the score
    """
    return percentileofscore(np.ma.masked_invalid(array_in).compressed(), score, **kwargs)


def make_norm(stretch='linear', vmin=None, vmax=None, a=None):
    """
    Create a matplotlib ``Normalize`` for display or ``rescale_image``.

    Parameters
    ----------
    stretch : str
        One of the names in ``stretch_functions`` plus ``'symlog'`` and
        ``'symmetric_log'``.  ``'symlog'`` uses matplotlib's piecewise
        ``SymLogNorm``.  ``'symmetric_log'`` uses the C¹-continuous
        ``pysymlog.SymmetricLogarithmNorm`` when the optional ``pysymlog``
        package is installed (``pip install multicolorfits[pysymlog]``).
    vmin, vmax : float or None
        Display limits.  Required for signed stretches.
    a : float or None
        Stretch parameter: for ``symlog`` this is ``linthresh`` in data units
        (default 1% of the display range); for ``symmetric_log`` this is
        pysymlog's ``shift`` (default 0.1% of the range).

    Returns
    -------
    matplotlib.colors.Normalize
    """
    s = stretch.lower()
    if s in SIGNED_STRETCHES:
        if vmin is None or vmax is None:
            raise ValueError("vmin and vmax are required for stretch=%r" % (stretch,))
        vlo, vhi = float(vmin), float(vmax)
        if s == 'symmetric_log':
            try:
                import pysymlog
            except ImportError as exc:
                raise ImportError(
                    "stretch='symmetric_log' requires the optional `pysymlog` "
                    "package. Install with `pip install pysymlog` or "
                    "`pip install multicolorfits[pysymlog]`. For matplotlib's "
                    "piecewise variant, use stretch='symlog' instead."
                ) from exc
            pysymlog.register_mpl()
            shift = a if a is not None else max(abs(vhi - vlo) * 0.001, 1e-10)
            return pysymlog.SymmetricLogarithmNorm(shift=shift, vmin=vlo, vmax=vhi)
        linthresh = a if a is not None else max(abs(vhi - vlo) * 0.01, 1e-10)
        return mcolors.SymLogNorm(linthresh=linthresh, vmin=vlo, vmax=vhi)

    if s not in stretch_functions:
        raise ValueError(
            "Unknown stretch %r. Available: %s"
            % (stretch, ', '.join(sorted(list(stretch_functions) + list(SIGNED_STRETCHES))))
        )
    return stretch_functions[s]() + ManualInterval(vmin=vmin, vmax=vmax)


def rescale_image(datin, rescalefn='linear', vmin=None, vmax=None, a=None):
    """
    Rescale image data to the [0,1] range with the specified stretch function
    and manual interval.  This is the core scaling operation used by
    to_grey_rgb() and the GUI previews.

    Parameters
    ----------
    datin : array
        Input 2D image data array
    rescalefn : str
        One of 'linear','sqrt','squared','log','power','sinh','asinh',
        'symlog', or 'symmetric_log' (the last requires ``pysymlog``).
    vmin : float or None
        Minimum data value for the interval.  None uses the data (nan)min.
    vmax : float or None
        Maximum data value for the interval.  None uses the data (nan)max.
    a : float or None
        Optional stretch parameter for ``symlog`` / ``symmetric_log`` (see
        ``make_norm``).

    Returns
    -------
    array
        Data rescaled to [0,1]
    """
    if vmin is None:
        vmin = float(np.nanmin(datin))
    if vmax is None:
        vmax = float(np.nanmax(datin))
    if rescalefn.lower() in SIGNED_STRETCHES:
        norm = make_norm(rescalefn, vmin=vmin, vmax=vmax, a=a)
        arr = np.nan_to_num(np.asarray(datin, dtype=float), nan=0.0)
        return np.clip(norm(arr), 0.0, 1.0)
    return (stretch_functions[rescalefn]() + ManualInterval(vmin=vmin, vmax=vmax))(datin)


# Thresholds for suggest_levels. The docstring is the user-facing statement of
# these rules; keep the numbers here and in that docstring in lockstep.
_LEVEL_MIN_SAMPLES = 16
_LEVEL_SAMPLE_MAX = 250_000
_LEVEL_SAMPLE_SEED = 0
_LEVEL_NEG_FRAC = 0.01
_LEVEL_LINEAR_DECADES = 1.0
_LEVEL_SQRT_DECADES = 2.0
_LEVEL_WIDE_DECADES = 4.0
_LEVEL_BRIGHT_TAIL = 2.0
# Lower half is "piled up" when the median is less than this times p5.
_LEVEL_FLOOR_BULK = 3.0
# ...and the top is this many times the median above that pile.
_LEVEL_FLOOR_CONTRAST = 20.0
_LEVEL_PERCENTILES = (1, 5, 16, 50, 84, 95, 99, 99.5, 99.9)


def _sample_finite(data, ignore_zeros):
    """1D finite sample, optional exact-zero drop, and a subsample cap."""
    arr = np.asarray(data, dtype=float).ravel()
    finite = arr[np.isfinite(arr)]
    n_zeros = 0
    if ignore_zeros and finite.size:
        zero = finite == 0
        n_zeros = int(np.count_nonzero(zero))
        finite = finite[~zero]
    n_finite = int(finite.size)
    if n_finite > _LEVEL_SAMPLE_MAX:
        rng = np.random.default_rng(_LEVEL_SAMPLE_SEED)
        finite = finite[rng.choice(n_finite, _LEVEL_SAMPLE_MAX, replace=False)]
    return finite, n_finite, n_zeros


def _levels_result(stretch, vmin, vmax, reason, span, crosses, on_floor,
                   bright_tail, n_zeros, n_finite, n_sample, pct):
    if not np.isfinite(vmin):
        vmin = 0.0
    if not np.isfinite(vmax) or vmax <= vmin:
        eps = max(abs(float(vmin)) * 1e-6, 1e-12)
        vmax = float(vmin) + eps
        if 'constant' not in reason.lower():
            reason = reason.rstrip('.') + '. Image is constant.'
    return {
        'stretch': stretch,
        'vmin': float(vmin),
        'vmax': float(vmax),
        'reason': reason,
        'span_decades': None if span is None else float(span),
        'crosses_zero': bool(crosses),
        'on_floor': bool(on_floor),
        'bright_tail': bool(bright_tail),
        'zeros_excluded': int(n_zeros),
        'n_finite': int(n_finite),
        'n_sample': int(n_sample),
        'percentiles': pct,
    }


def _print_levels(info, name=''):
    label = (' (%s)' % name) if name else ''
    print('Levels suggestion%s' % label)
    print('  Stretch:   %s' % info['stretch'])
    print('  Limits:    %.6g .. %.6g' % (info['vmin'], info['vmax']))
    if info['span_decades'] is None:
        print('  Span:      n/a')
    else:
        print('  Span:      %.2f decades' % info['span_decades'])
    print('  Reason:    %s' % info['reason'])


def suggest_levels(data, ignore_zeros=True, verbose=False):
    """Recommend a stretch and absolute vmin/vmax as a display starting point.

    This does not produce a finished display. A single stretch cannot show
    faint structure and bright peaks across many decades at once (a supernova
    remnant is the usual example). The returned ``reason`` says so when the
    span is very wide, and the bright core may be allowed to saturate so the
    faint end is visible. Treat the result as the first values to put in a
    panel, then adjust by eye. Splitting one layer into faint, medium, and
    bright displays is out of scope here.

    Parameters
    ----------
    data : array
        Image (any shape). Non-finite pixels are ignored.
    ignore_zeros : bool, optional
        Drop exact zeros before measuring (default True) so chip gaps and
        blanked borders do not set the floor. Pass False when zero is a real
        measurement. An image that is all exact zeros raises ``ValueError``
        when this is True, because nothing remains to measure.
    verbose : bool, optional
        Print the same summary :func:`describe_image` prints for the levels.

    Returns
    -------
    dict
        ``stretch`` (``'linear'``, ``'sqrt'``, or ``'asinh'``), ``vmin``,
        ``vmax`` (absolute data values, not percentiles), ``reason`` (one
        line), ``span_decades`` (float or None), ``crosses_zero``,
        ``on_floor``, ``bright_tail``, ``zeros_excluded``, ``n_finite``,
        ``n_sample``, and ``percentiles`` (keys 1, 5, 16, 50, 84, 95, 99,
        99.5, 99.9).

    Raises
    ------
    ValueError
        If no usable pixels remain after dropping non-finite values and,
        when requested, exact zeros.

    Notes
    -----
    At most 250000 pixels are sampled (fixed seed) so a large image stays
    cheap. The stretch is never ``log``: log collapses at or below zero.
    ``symlog`` is mentioned in the reason when the image crosses zero, but
    it is not selected (scripting-only in the GUI). The optional
    ``symmetric_log`` stretch is never selected.

    Rules, applied in order:

    1. Fewer than 16 usable pixels: ``linear`` between the sample min and max.
    2. At least 1% of the sample is negative, or the background sits on the
       floor: ``asinh``. Floor means the 5th percentile is at or below zero,
       or the lower half is piled up (the median is less than 3 times the
       5th percentile) while the 99.5th percentile is at least 20 times the
       median. A distribution that fills many decades is not a floor — the
       lower half is spread out — and uses rule 3 instead.
    3. Otherwise the span is ``log10(p99.5 / p5)`` on that positive sample.
       Under 1 decade: ``linear``. Under 2: ``sqrt``. Wider: ``asinh``.
    4. Limits are a mild clip: vmin is the 1st percentile, vmax is the 99.9th.
       If the 99.9th percentile is at least twice the 99.5th, the bright tail
       is thin, vmax stays at the 99.5th, and the reason says the core will
       saturate.
    5. If the span is more than 4 decades, the stretch from the rules above
       is unchanged and the reason notes that one stretch will not show faint
       structure and bright peaks together.

    A constant image returns ``linear`` with a tiny range so a stretch does
    not divide by zero.

    See Also
    --------
    describe_image, zscale_limits, rescale_image
    """
    sample, n_finite, n_zeros = _sample_finite(data, ignore_zeros)
    n_sample = int(sample.size)
    if n_sample < _LEVEL_MIN_SAMPLES:
        if n_sample == 0:
            raise ValueError(
                'suggest_levels: no finite pixels'
                + (' after excluding exact zeros' if ignore_zeros and n_zeros else ''))
        vmin = float(np.min(sample))
        vmax = float(np.max(sample))
        info = _levels_result(
            'linear', vmin, vmax,
            'too few finite pixels (%d); using linear min/max' % n_sample,
            None, False, False, False, n_zeros, n_finite, n_sample, {})
        if verbose:
            _print_levels(info)
        return info

    qs = list(_LEVEL_PERCENTILES)
    vals = np.percentile(sample, qs)
    pct = {q: float(v) for q, v in zip(qs, vals)}
    p1, p5, p50 = pct[1], pct[5], pct[50]
    p995, p999 = pct[99.5], pct[99.9]
    frac_neg = float(np.count_nonzero(sample < 0)) / n_sample
    crosses = frac_neg >= _LEVEL_NEG_FRAC
    piled = (p5 > 0 and p50 > 0 and (p50 / p5) < _LEVEL_FLOOR_BULK
             and (p995 / p50) >= _LEVEL_FLOOR_CONTRAST)
    on_floor = (not crosses) and (p5 <= 0 or piled)
    span = None
    if p5 > 0 and p995 > p5:
        span = float(np.log10(p995 / p5))
    bright = bool(p995 > 0 and (p999 / p995) >= _LEVEL_BRIGHT_TAIL)
    vmin = p1
    vmax = p995 if bright else p999

    bits = []
    if n_zeros:
        bits.append('Exact zeros excluded (%d).' % n_zeros)
    if bright:
        bits.append('Bright tail above the 99.5th percentile will saturate.')
    if span is not None and span > _LEVEL_WIDE_DECADES:
        bits.append(
            'A single stretch will not show faint structure and bright peaks together.')
    extra = (' ' + ' '.join(bits)) if bits else ''

    if crosses or on_floor:
        stretch = 'asinh'
        if crosses:
            pct_txt = '%.0f%%' % (100.0 * frac_neg) if frac_neg >= 0.01 else 'some'
            reason = (
                '%s of pixels are negative; asinh (log collapses through zero; '
                'scripts may prefer symlog).%s' % (pct_txt, extra))
        else:
            reason = 'background sits on the floor; asinh rather than log.%s' % extra
    elif span is None or span < _LEVEL_LINEAR_DECADES:
        stretch = 'linear'
        if span is None:
            reason = 'narrow positive range; linear with a mild percentile clip.%s' % extra
        else:
            reason = ('%.1f decades; linear with a mild percentile clip.%s'
                      % (span, extra))
    elif span < _LEVEL_SQRT_DECADES:
        stretch = 'sqrt'
        reason = '%.1f decades; sqrt.%s' % (span, extra)
    else:
        stretch = 'asinh'
        reason = '%.1f decades; asinh.%s' % (span, extra)

    info = _levels_result(
        stretch, vmin, vmax, reason.strip(), span, crosses, on_floor, bright,
        n_zeros, n_finite, n_sample, pct)
    if verbose:
        _print_levels(info)
    return info


def describe_image(data, hdr=None, name='', ignore_zeros=True, verbose=True):
    """Summarize a header (if given) and suggest a display stretch.

    The header half is :func:`describe_header` (frame, scale, flip, beam).
    The levels half is :func:`suggest_levels`: a starting stretch and
    vmin/vmax, plus a one-line reason. This does not change the array or
    the header, and it does not claim the suggestion is a finished display.

    For several layers at once (plus suggested colors), use
    :func:`describe_images`.

    Parameters
    ----------
    data : array
        Image whose pixel distribution is measured.
    hdr : header or None, optional
        If given, the celestial summary is included (and printed when
        *verbose*).
    name : str, optional
        Label next to both halves.
    ignore_zeros : bool, optional
        Passed to :func:`suggest_levels` (default True).
    verbose : bool, optional
        Print to stdout (default True). Pass False to build the return
        dict without printing — useful when collecting suggestions in a
        loop.

    Returns
    -------
    dict
        ``header`` is the :func:`describe_header` dict, or None. ``levels``
        is the :func:`suggest_levels` dict.

    Raises
    ------
    ValueError
        If no usable pixels remain (see :func:`suggest_levels`).

    See Also
    --------
    describe_images, suggest_levels, describe_header
    """
    header_info = None
    if hdr is not None:
        from .wcs_tools import describe_header
        header_info = describe_header(hdr, name=name, verbose=verbose)
    levels = suggest_levels(data, ignore_zeros=ignore_zeros, verbose=False)
    if verbose:
        _print_levels(levels, name=name)
    return {'header': header_info, 'levels': levels}


def _coerce_describe_layer(item):
    """Normalize one layer to ``(data, header_or_None)``.

    Accepts a 2D/3D array, a ``(data, header)`` pair, or a FITS path.
    """
    if isinstance(item, (str, bytes)) or hasattr(item, '__fspath__'):
        from .io_output import read_fits
        data, hdr = read_fits(item)
        return np.asarray(data), hdr
    if isinstance(item, (list, tuple)) and len(item) == 2:
        data, hdr = item
        if hdr is None or hasattr(hdr, 'keys') or hasattr(hdr, 'cards'):
            return np.asarray(data), hdr
    return np.asarray(item), None


def describe_images(layers, names=None, colors=None, palette=None,
                    ignore_zeros=True, verbose=False):
    """Characterize several images for a multicolor starting point.

    Runs :func:`describe_image` on each layer and adds a color suggestion
    for the set. Use this instead of looping when you want a dict ready
    for :meth:`~multicolorfits.McfSession.load_files` /
    :meth:`~multicolorfits.McfSession.set_display`.

    Parameters
    ----------
    layers : mapping or sequence
        * dict: ``name -> array`` or ``name -> (data, header)`` (or a FITS
          path). Order follows the dict's insertion order.
        * sequence: arrays, ``(data, header)`` pairs, or FITS paths. Pair
          with *names* when you have labels.
    names : sequence of str or None, optional
        Labels when *layers* is a sequence. Default ``Image1``, ``Image2``,
        …. Ignored when *layers* is a mapping.
    colors : sequence of str or None, optional
        Hex colors, one per layer. If None and *palette* is None, uses
        :func:`~multicolorfits.suggest_colors`.
    palette : str or None, optional
        Named palette (see :func:`~multicolorfits.list_palettes`). When
        set, overrides *colors*.
    ignore_zeros : bool, optional
        Passed to :func:`suggest_levels` (default True).
    verbose : bool, optional
        Print each layer's summary (default False — this helper is meant
        for collecting params). Pass True for an interactive dump.

    Returns
    -------
    dict
        ``names``, ``colors``, ``stretches``, ``vmins``, ``vmaxs`` —
        parallel lists for session helpers — and ``layers``: a mapping
        ``name -> {'header', 'levels', 'color'}``.

    Examples
    --------
    >>> report = mcf.describe_images(aligned)  # doctest: +SKIP
    >>> # aligned is name -> (data, header)
    >>> s.load_files(paths, colors=report['colors'], labels=report['names'],
    ...              stretches=report['stretches'],
    ...              vmins=report['vmins'], vmaxs=report['vmaxs'])  # doctest: +SKIP

    See Also
    --------
    describe_image, suggest_levels, suggest_colors
    """
    if isinstance(layers, dict):
        items = list(layers.items())
        layer_names = [str(k) for k, _ in items]
        layer_vals = [v for _, v in items]
    else:
        layer_vals = list(layers)
        if names is not None:
            layer_names = [str(n) for n in names]
            if len(layer_names) != len(layer_vals):
                raise ValueError('names has %d entries; expected %d'
                                 % (len(layer_names), len(layer_vals)))
        else:
            layer_names = ['Image%d' % (i + 1) for i in range(len(layer_vals))]

    n = len(layer_vals)
    if n == 0:
        raise ValueError('describe_images: no layers given')

    if palette is not None:
        from ..palettes import get_palette
        color_list = list(get_palette(palette, n=n))
    elif colors is not None:
        color_list = list(colors)
        if len(color_list) == 1 and n != 1:
            color_list = color_list * n
        if len(color_list) != n:
            raise ValueError('colors has %d entries; expected %d or 1'
                             % (len(color_list), n))
    else:
        from ..palettes import suggest_colors
        color_list = list(suggest_colors(n))

    layer_map = {}
    stretches, vmins, vmaxs = [], [], []
    for name, raw, color in zip(layer_names, layer_vals, color_list):
        data, hdr = _coerce_describe_layer(raw)
        info = describe_image(data, hdr, name=name, ignore_zeros=ignore_zeros,
                              verbose=verbose)
        info = dict(info)
        info['color'] = color
        layer_map[name] = info
        levels = info['levels']
        stretches.append(levels['stretch'])
        vmins.append(levels['vmin'])
        vmaxs.append(levels['vmax'])

    return {
        'names': layer_names,
        'colors': color_list,
        'stretches': stretches,
        'vmins': vmins,
        'vmaxs': vmaxs,
        'layers': layer_map,
    }


def zscale_limits(datin):
    """
    Compute display limits using the IRAF zscale algorithm
    (astropy.visualization.ZScaleInterval).

    Parameters
    ----------
    datin : array
        Input image data array

    Returns
    -------
    tuple
        (vmin, vmax)
    """
    lims = ZScaleInterval().get_limits(datin)
    return float(lims[0]), float(lims[1])
