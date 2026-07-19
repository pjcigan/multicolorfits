"""
Optional skyplothelper integration for WCS overlays and advanced axes.

mcf keeps compositing and basic WCS plotting in-core; compass roses, beam
markers, scale bars, offset coordinate frames, and other specialized
annotation live behind the optional ``multicolorfits[overlays]`` extra
(``skyplothelper``).  Nothing in this module is imported at package import
time — callers pay for skyplothelper only when they use these helpers.

Simple case (opinionated wrappers)::

    import multicolorfits as mcf

    fig, ax = mcf.make_combined_figure(session), ...
    mcf.overlays.add_compass(ax)
    mcf.overlays.add_beam(ax, session.common_header)
    mcf.overlays.add_scale_bar(ax, session.common_header)

Bring-your-own skyplothelper (full API)::

    sph = mcf.overlays.skyplothelper()
    wcs = sph.offset_coord_WCS(hdr, center)
    sph.add_coord_overlay(ax, wcs, ...)
"""

from __future__ import annotations

import importlib
import importlib.util
import math
import warnings
from typing import Any, List, Optional

__all__ = [
    'overlays_available',
    'require_overlays',
    'skyplothelper',
    'SkyplotHelper',
    'add_compass',
    'add_beam',
    'add_scale_bar',
    'apply_overlays',
    'apply_compose_overlays',
    'offset_coord_wcs',
    'make_offset_figure',
]

_INSTALL_HINT = 'pip install multicolorfits[overlays]'


def overlays_available() -> bool:
    """True if the optional ``skyplothelper`` package is importable."""
    return importlib.util.find_spec('skyplothelper') is not None


def require_overlays():
    """
    Import and return the ``skyplothelper`` module.

    Raises
    ------
    ImportError
        With install instructions when the optional extra is missing.
    """
    try:
        return importlib.import_module('skyplothelper')
    except ImportError as exc:
        raise ImportError(
            "multicolorfits overlays require the optional 'skyplothelper' "
            "package. Install with: %s" % _INSTALL_HINT
        ) from exc


class SkyplotHelper:
    """
    Lazy attribute proxy to the full ``skyplothelper`` namespace.

    For advanced plotting (offset overlays, graticules, globe frames, etc.)
    prefer this over re-vendoring edge-case logic into mcf.
    """

    def __getattr__(self, name: str) -> Any:
        return getattr(require_overlays(), name)

    def __repr__(self) -> str:
        return 'SkyplotHelper(%s)' % ('available' if overlays_available() else 'not installed')


_sph_proxy = SkyplotHelper()


def skyplothelper():
    """Return the ``skyplothelper`` module (lazy). Alias for :func:`require_overlays`."""
    return require_overlays()


def _header_astropy(header: Any) -> Any:
    return getattr(header, 'astropy', header)


def _auto_scale_bar_asec(header: Any) -> float:
    """Pick a readable scale-bar length (~15% of the field width)."""
    from .core.wcs_tools import arcsec_per_pixel

    hdr = _header_astropy(header)
    try:
        pp = abs(float(arcsec_per_pixel(hdr)))
        width = float(hdr['NAXIS1'])
        raw = pp * width * 0.15
    except Exception:
        raw = 10.0
    if raw <= 0 or not math.isfinite(raw):
        return 10.0
    exp = 10.0 ** math.floor(math.log10(raw))
    for mult in (1.0, 2.0, 5.0, 10.0):
        candidate = mult * exp
        if candidate >= raw * 0.8:
            return round(candidate, 6)
    return round(10.0 * exp, 6)


def _scale_bar_label(length_asec: float) -> str:
    if length_asec >= 60.0:
        minutes = length_asec / 60.0
        if abs(minutes - round(minutes)) < 1e-6:
            return "%g'" % round(minutes)
        return "%.1f'" % minutes
    if abs(length_asec - round(length_asec)) < 1e-6:
        return '%g"' % round(length_asec)
    return '%.1f"' % length_asec


def _header_has_beam(header: Any) -> bool:
    hdr = _header_astropy(header)
    try:
        return all(k in hdr for k in ('BMAJ', 'BMIN'))
    except Exception:
        return False


def add_compass(ax: Any, **kwargs: Any) -> Any:
    """Add a WCS-aware north/east compass (delegates to skyplothelper)."""
    return require_overlays().add_compass(ax, **kwargs)


def add_beam(ax: Any, header: Any = None, *, loc: str = 'lower left',
             style: str = 'crosshair', anchored: bool = True, **kwargs: Any) -> Any:
    """
    Add a FITS-beam marker from ``BMAJ`` / ``BMIN`` / ``BPA`` in *header*.

    Parameters
    ----------
    ax : matplotlib Axes
        Target WCS axes.
    header : Header-like
        FITS header or :class:`~multicolorfits.core.McfHeader`.
    loc : str
        Corner for anchored placement (``'lower left'``, ...).
    style : str
        Beam style passed to skyplothelper (``'crosshair'`` is typical for radio).
    anchored : bool
        If True (default), keep the beam in a fixed corner; otherwise draw at
        the beam centre in data coordinates.
    """
    if header is None:
        raise ValueError('add_beam requires a FITS header with BMAJ/BMIN cards')
    sph = require_overlays()
    beam = sph.Beam.from_header(_header_astropy(header), ax=ax, style=style, **kwargs)
    if anchored:
        return beam.add_anchored(ax, loc=loc)
    return beam.add_to(ax)


def add_scale_bar(ax: Any, header: Any = None, *, length_asec: Optional[float] = None,
                  label: Optional[str] = None, **kwargs: Any) -> Any:
    """
    Add an arcsecond scale bar (delegates to skyplothelper ``add_sizebar_asec``).

    When *length_asec* is None or 0, a length of ~15% of the field width is
    chosen automatically from the header pixel scale.
    """
    if header is None:
        raise ValueError('add_scale_bar requires a FITS header for pixel-scale conversion')
    length = float(length_asec) if length_asec else 0.0
    if length <= 0:
        length = _auto_scale_bar_asec(header)
    if label is None:
        label = _scale_bar_label(length)
    return require_overlays().add_sizebar_asec(ax, _header_astropy(header),
                                               length, label, **kwargs)


def apply_overlays(ax: Any, header: Any = None, *, compass: bool = False,
                   beam: bool = False, scale_bar: bool = False,
                   compass_loc: str = 'lower left', beam_loc: str = 'lower left',
                   beam_style: str = 'crosshair', scale_bar_asec: float = 0.0,
                   scale_bar_loc: int = 4, warn_missing: bool = True,
                   **kwargs: Any) -> List[Any]:
    """
    Apply selected skyplothelper overlays to an existing axes.

    Returns a list of overlay artists / handles created.  When
    ``skyplothelper`` is not installed and *warn_missing* is True, emits a
    ``UserWarning`` and returns an empty list.
    """
    if not (compass or beam or scale_bar):
        return []
    if not overlays_available():
        if warn_missing:
            warnings.warn(
                'skyplothelper overlays were requested but the package is not '
                'installed. Install with: %s' % _INSTALL_HINT,
                UserWarning,
                stacklevel=2,
            )
        return []
    artists: List[Any] = []
    compass_kw = dict(kwargs.get('compass_kwargs') or {})
    beam_kw = dict(kwargs.get('beam_kwargs') or {})
    scale_kw = dict(kwargs.get('scale_bar_kwargs') or {})
    if compass:
        artists.append(add_compass(ax, loc=compass_loc, **compass_kw))
    if beam:
        if header is None:
            warnings.warn('Beam overlay requested but no header was provided',
                          UserWarning, stacklevel=2)
        elif not _header_has_beam(header):
            warnings.warn('Beam overlay requested but header has no BMAJ/BMIN cards',
                          UserWarning, stacklevel=2)
        else:
            artists.append(add_beam(ax, header, loc=beam_loc, style=beam_style, **beam_kw))
    if scale_bar:
        if header is None:
            warnings.warn('Scale-bar overlay requested but no header was provided',
                          UserWarning, stacklevel=2)
        else:
            length = scale_bar_asec if scale_bar_asec and scale_bar_asec > 0 else None
            artists.append(add_scale_bar(ax, header, length_asec=length,
                                         loc=scale_bar_loc, **scale_kw))
    return artists


def apply_compose_overlays(ax: Any, header: Any, compose: Any,
                           *, warn_missing: bool = True) -> List[Any]:
    """
    Apply overlay toggles from a :class:`~multicolorfits.session.ComposeState`.

    Intended for :func:`~multicolorfits.figures.setup_combined_axes`; GUI
    controls can flip the compose flags without importing skyplothelper until
    a full WCS render is requested.
    """
    scale_kw = {}
    if getattr(compose, 'scale_bar_color', None):
        scale_kw['color'] = compose.scale_bar_color
    if getattr(compose, 'scale_bar_stroke_color', None) is not None:
        scale_kw['stroke_color'] = compose.scale_bar_stroke_color
    if getattr(compose, 'scale_bar_stroke_lw', None) is not None:
        scale_kw['stroke_lw'] = float(compose.scale_bar_stroke_lw)
    return apply_overlays(
        ax,
        header,
        compass=getattr(compose, 'show_compass', False),
        beam=getattr(compose, 'show_beam', False),
        scale_bar=getattr(compose, 'show_scale_bar', False),
        compass_loc=getattr(compose, 'compass_loc', 'lower left'),
        beam_loc=getattr(compose, 'beam_loc', 'lower left'),
        beam_style=getattr(compose, 'beam_style', 'crosshair'),
        scale_bar_asec=float(getattr(compose, 'scale_bar_asec', 0.0) or 0.0),
        scale_bar_loc=int(getattr(compose, 'scale_bar_loc', 4)),
        scale_bar_kwargs=scale_kw,
        warn_missing=warn_missing,
    )


def offset_coord_wcs(header: Any, center: Any, **kwargs: Any) -> Any:
    """
    Build a locally linear offset WCS centred on *center* (skyplothelper).

    For offset tick labels, graticule overlays, and cutout-style axes, use the
    returned WCS with skyplothelper's ``add_coord_overlay`` / ``add_overlay_ticks``
    rather than reimplementing those paths in mcf.
    """
    return require_overlays().offset_coord_WCS(_header_astropy(header), center, **kwargs)


def make_offset_figure(center: Any, **kwargs: Any):
    """One-call offset field figure (delegates to skyplothelper ``offset_figure``)."""
    return require_overlays().offset_figure(center, **kwargs)
