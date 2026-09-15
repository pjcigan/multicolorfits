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


def add_compass(ax: Any, *, loc: str = 'lower left', length: float = 0.08,
                pad: float = 0.02, color: Any = None, fontsize: float = 10,
                lw: float = 1.5, stroke_color: Any = None, stroke_lw: float = 2,
                label_offset: float = 1.3, north_label: str = 'N',
                east_label: str = 'E', head_width: float = 0.015,
                head_length: float = 0.012, zorder: int = 10,
                **kwargs: Any) -> Any:
    """
    Add a WCS-aware north/east compass.

    The arrows follow the header's north and east, not the pixel axes.
    Delegates to skyplothelper; extra keywords are ``arrowprops`` overrides
    (``linestyle``, ``mutation_scale``), not general artist properties.

    Parameters
    ----------
    ax : matplotlib Axes
        WCS axes (the mosaic hero qualifies).
    loc : str or (x, y), optional
        Corner name (``'lower left'``, ``'upper right'``, …) or an
        axes-fraction pair. Default ``'lower left'``. An ``(x, y)`` pair
        ignores *pad*.
    length : float, optional
        Arrow length as a fraction of the axes. Default 0.08. This also
        pushes the tail inward: a named corner places the tail at
        ``pad + length``, so shortening the arrows pulls the compass toward
        the corner as well as making it smaller.
    pad : float, optional
        Extra inset of the arrow tails, in axes fraction, before *length*
        is added. Default 0.02. skyplothelper's own default is 0.05, which
        with ``length=0.08`` puts the joint 13% into the frame and the
        letters near 23%. Pass ``pad=0`` to sit the tails only one arrow
        length in. Raise it if an outward letter clips the edge.
    color : color, optional
        Arrow and letter ink. Default is the matplotlib text color.
    fontsize : float, optional
        Letter size in points. Default 10.
    lw : float, optional
        Arrow line width in points. Default 1.5.
    stroke_color, stroke_lw : color, float, optional
        Outline behind the arrows and letters. ``stroke_lw`` defaults to 2.
        ``stroke_color`` defaults to the axes facecolor (a contrast stroke,
        not "no stroke"). ``stroke_color='black'`` is a dark halo on a
        bright field.
    label_offset : float, optional
        How far the N/E letters sit past the arrowheads, in arrow lengths.
        Default 1.3. Lower it to keep the letters from reaching into the
        image.
    north_label, east_label : str, optional
        Letter text. ``''`` hides that letter.
    head_width, head_length : float, optional
        Arrow-head size as a fraction of the axes.
    zorder : int, optional
        Default 10.

    Returns
    -------
    list
        Annotation artists. Restyle with ``path_effects`` on those artists
        (and on ``artist.arrow_patch`` for the shafts).

    Examples
    --------
    ::

        import multicolorfits as mcf
        fig = mcf.make_combined_figure(s)
        ax = fig.axes[0]
        mcf.overlays.add_compass(
            ax, loc='upper right', pad=0.02, color='white',
            stroke_color='black')
    """
    return require_overlays().add_compass(
        ax, loc=loc, length=length, pad=pad, color=color, fontsize=fontsize,
        lw=lw, stroke_color=stroke_color, stroke_lw=stroke_lw,
        label_offset=label_offset, north_label=north_label,
        east_label=east_label, head_width=head_width, head_length=head_length,
        zorder=zorder, **kwargs)


def add_beam(ax: Any, header: Any = None, *, loc: str = 'lower left',
             style: str = 'crosshair', anchored: bool = True,
             borderpad: float = 0.4, **kwargs: Any) -> Any:
    """
    Add a FITS-beam marker from BMAJ / BMIN / BPA cards in *header*.

    Parameters
    ----------
    ax : matplotlib Axes
        Target WCS axes.
    header : Header-like
        FITS header or :class:`~multicolorfits.core.McfHeader`.
    loc : str, optional
        Corner when *anchored* is True (``lower left``, ``upper right``,
        etc.). Ignored when *anchored* is False.
    style : str, optional
        ``ellipse``, ``crosshair`` (default; typical for radio), or
        ``crosshairgrid``.
    anchored : bool, optional
        If True (default), pin the marker to *loc*. If False, draw at the
        beam centre in data coordinates.
    borderpad : float, optional
        Padding between an anchored beam and the axes edge, in fractions of
        the font size (matplotlib offset-box units, not axes fraction).
        Default 0.4. Larger values push the beam further into the frame.
        Ignored when *anchored* is False.
    ec, lw, fc, alpha
        Ellipse edge, width, face, and opacity. *fc* defaults to none.
        Pass these as keywords.
    stroke_color, stroke_lw : color, float
        Outline behind the ellipse and crosshair. Default *stroke_color*
        of None draws no stroke. *stroke_lw* defaults to 3 when a stroke
        color is set.
    crosshair_color, crosshair_lw, crosshair_ls
        Crosshair ink, width, and linestyle. Color defaults to the ellipse
        edge.
    grid_marker, grid_density
        Only used when *style* is ``crosshairgrid``.

    Other keywords go to the skyplothelper Beam constructor (patch
    properties), not to the corner box. The return value is the anchored
    offset box when *anchored* is True, otherwise the Beam itself.
    Call Beam.set_stroke(color, lw) on an unanchored beam to restyle it
    afterwards.
    """
    if header is None:
        raise ValueError('add_beam requires a FITS header with BMAJ/BMIN cards')
    sph = require_overlays()
    beam = sph.Beam.from_header(_header_astropy(header), ax=ax, style=style, **kwargs)
    if anchored:
        return beam.add_anchored(ax, loc=loc, borderpad=borderpad)
    return beam.add_to(ax)


def add_scale_bar(ax: Any, header: Any = None, *, length_asec: Optional[float] = None,
                  label: Optional[str] = None, **kwargs: Any) -> Any:
    """
    Add an arcsecond scale bar (delegates to skyplothelper ``add_sizebar_asec``).

    When *length_asec* is None or 0, a length of ~15% of the field width is
    chosen automatically from the header pixel scale. Remaining keywords are
    forwarded to matplotlib's ``AnchoredSizeBar``.

    Parameters
    ----------
    ax : matplotlib Axes
    header : Header-like
        FITS header or :class:`~multicolorfits.core.McfHeader` with a pixel
        scale.
    length_asec : float, optional
        Bar length in arcseconds. None or 0 picks a round length near 15%
        of the field.
    label : str, optional
        Text beside the bar. Default is derived from the length (``30"``,
        ``2'``).
    loc : int, optional
        Matplotlib location code, not a corner name: 1 upper right, 2 upper
        left, 3 lower left, 4 lower right (default).
    color : color, optional
        Bar and label ink. Default is the matplotlib text color.
    fontproperties, prop : FontProperties or dict, optional
        Label font. There is no ``fontsize`` keyword. A dict works::

            add_scale_bar(..., fontproperties={'size': 8})

        or ``FontProperties(size=8)``. ``prop`` is the matplotlib alias.
        After the fact, ``bar.txt_label._text.set_fontsize(8)`` on the
        returned artist does the same thing.
    stroke_color, stroke_lw : color, float, optional
        Outline under the bar and the label. Defaults are ``'k'`` and 1.75.
        ``stroke_color=None`` turns the stroke off.
    path_effects : list, optional
        Replaces ``stroke_color`` / ``stroke_lw``.
    borderpad : float, optional
        Padding between the bar and the axes edge, in fractions of the
        label font size (not axes fraction). skyplothelper's default is 0.8.
        Raise it to pull the bar off the frame; it will not move by a
        percent of the axes the way compass *pad* does.
    sep : float, optional
        Gap between the bar and its label, in points. Default 5.
    pad : float, optional
        Padding inside the scale-bar box, in fractions of the font size.
        Not the compass *pad*.
    frameon : bool, optional
        Draw a box around the bar. Default False.
    size_vertical : float, optional
        Bar thickness in data coordinates. Default 0 (a hairline).

    Returns
    -------
    matplotlib.offsetbox.AnchoredSizeBar
        Stroke is already applied. To change it later, set ``path_effects``
        on the bar line and the label text.

    Examples
    --------
    ::

        import multicolorfits as mcf
        fig = mcf.make_combined_figure(s)
        mcf.overlays.add_scale_bar(
            fig.axes[0], s.common_header, length_asec=30,
            loc=4, color='white', stroke_color='black',
            fontproperties={'size': 8})
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

    Styling is not a flat keyword on this function (``color=`` would be
    ambiguous). Pass a dict per overlay, using the same names as
    :func:`add_compass`, :func:`add_beam`, and :func:`add_scale_bar`.
    Compass *pad* is the tighter multicolorfits default unless the dict
    overrides it::

        apply_overlays(
            ax, hdr, compass=True, scale_bar=True,
            compass_kwargs=dict(pad=0.01, color='white', stroke_color='black'),
            scale_bar_kwargs=dict(color='white', stroke_color='black',
                                  stroke_lw=2, fontproperties={'size': 8},
                                  borderpad=0.4))

    Parameters
    ----------
    compass_loc, beam_loc : str, optional
        Corner name. Default ``'lower left'``.
    beam_style : str, optional
        See :func:`add_beam`. Default ``'crosshair'``.
    scale_bar_asec : float, optional
        Arcseconds. 0 picks a length automatically.
    scale_bar_loc : int, optional
        Matplotlib location code (4 = lower right). Not a corner name.
    compass_kwargs, beam_kwargs, scale_bar_kwargs : dict, optional
        Forwarded to the matching ``add_*`` helper.
    warn_missing : bool, optional
        If skyplothelper is not installed, warn and return ``[]`` instead
        of raising. Default True.

    Returns
    -------
    list
        Overlay artists created.

    Examples
    --------
    ::

        import multicolorfits as mcf
        fig = mcf.make_combined_figure(s)
        mcf.overlays.apply_overlays(
            fig.axes[0], s.common_header,
            compass=True, scale_bar=True, scale_bar_asec=30,
            compass_kwargs=dict(color='white', stroke_color='black'),
            scale_bar_kwargs=dict(color='white', stroke_color='black',
                                  fontproperties={'size': 8}))
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
