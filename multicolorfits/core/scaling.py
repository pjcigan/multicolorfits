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
