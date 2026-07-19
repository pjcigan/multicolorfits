"""
Output helpers: save RGB FITS cubes and make/save WCS-projected matplotlib plots.

matplotlib.pyplot is imported lazily inside the plotting functions so that
importing the core package stays cheap and headless-safe.
"""

import numpy as np
import astropy.io.fits as pyfits
from astropy.wcs import WCS

__all__ = [
    'plot_combined_rgb',
    'compare_multicolor_vs_rgb',
    'save_rgb_fits',
    'annotate_provenance_header',
]


def plot_combined_rgb(multicolorin, hdrin, axtitle, savepath, xaxislabel='RA', yaxislabel='DEC',
                            tickcolor='w', labelcolor='k', facecolor='w', minorticks=True, dpi=150,
                            legend=None, legend_loc='upper right', legend_inverse=False,
                            band_labels=None, band_labels_loc='upper left',
                            swatch=None, swatch_loc='lower right', swatch_labels=False,
                            swatch_label_offset=0.62, swatch_inset_scale=0.24,
                            show_compass=False, compass_loc='lower left',
                            show_beam=False, beam_loc='lower left', beam_style='crosshair',
                            show_scale_bar=False, scale_bar_asec=0.0, scale_bar_loc=4,
                            bare_plot=False):
    """
    Plot a multicolor RGB image with WCS axes and save it to disk.

    Tick and label colors accept any matplotlib color (``'k'``, ``'black'``,
    ``'#000000'``, ``'0.0'``, …).

    Parameters
    ----------
    multicolorin : array
        Multicolor RGB image -- as would be output from combine_multicolor()
    hdrin : astropy.io.fits.header
        Header to use for WCS information
    axtitle : str
        String to use for plot title.  e.g. "My crazy 5-color image", or empty quotes "" for nothing.
    savepath : str
        Path to save file to.  e.g., "./plots/mycrazyimage.pdf"
    xaxislabel : str
        Label to use for x-axis. Default='RA'
    yaxislabel : str
        Label to use for y-axis. Default='DEC'
    tickcolor : str
        Color to use for ticks in plot, default = 'w'.
    labelcolor : str
        Color to use for ticklabels, default = 'k'
    facecolor : str
        Color to use for figure facecolor (border around plot), default ='w'
    minorticks : bool
        True to use minor ticks
    dpi : int
        Default=150. Dots per inch value to use for saved plot.
    legend : list of (color, label) or None
        If given, draw a channel legend (which color = which image) over the
        image using these ``(matplotlib_color, str)`` entries.  Default None.
    legend_loc : str
        Matplotlib legend location, e.g. 'upper right'. Default='upper right'.
    legend_inverse : bool
        Set True when the combined image was made with ``inverse=True`` (light
        background) so legend text/frame colors stay readable. Default False.
    band_labels : list of (color, label) or None
        If given, draw colored per-band text labels in a corner (each label in
        its channel color).  Default None.
    band_labels_loc : str
        Corner for band labels: 'upper left'/'upper right'/'lower left'/'lower
        right'. Default 'upper left'.
    swatch : dict or None
        Pre-rendered color-combination swatch from ``mcf.combo_swatch(...)``
        (keys 'rgba', 'circles', 'size').  If given, drawn as a floating corner
        inset showing the honest overlapping mixed colors.  Default None.
    swatch_loc : str
        Corner for the swatch inset. Default 'lower right'.
    swatch_labels : bool
        Label each swatch circle in its color. Default False.
    swatch_label_offset : float
        Radial push for swatch circle labels (fraction of radius). Default 0.62.
    swatch_inset_scale : float
        Corner inset size as a fraction of the axes (0.08–0.45). Default 0.24.
    show_compass, show_beam, show_scale_bar : bool
        Optional skyplothelper overlays (require ``multicolorfits[overlays]``).
    compass_loc, beam_loc : str
        Corner placement for compass / beam overlays.
    beam_style : str
        Beam style passed to skyplothelper (e.g. ``'crosshair'``).
    scale_bar_asec : float
        Scale-bar length in arcsec; 0 selects an automatic length.
    scale_bar_loc : int
        Matplotlib anchored-artist location code (1–4) for the scale bar.
    """
    from matplotlib import pyplot as plt
    plt.rcParams.update({'font.family': 'serif', 'xtick.major.size': 6, 'ytick.major.size': 6,
                         'xtick.major.width': 1., 'ytick.major.width': 1.,
                         'xtick.direction': 'in', 'ytick.direction': 'in',
                         'xtick.top': True, 'ytick.right': True})
    fig1 = plt.figure(1)
    fig1.set_facecolor(facecolor)
    wcs = WCS(hdrin)
    ax1 = fig1.add_subplot(111, projection=wcs)
    ax1.imshow(multicolorin, origin='lower', interpolation='nearest')
    if bare_plot:
        from ..figures import apply_bare_axes
        apply_bare_axes(ax1, fig1)
    else:
        ax1.set_title(axtitle, color=labelcolor, size=12)
        ax1.coords.frame.set_color(tickcolor)  # This way is particular to astropy wcs frames
        rapars = ax1.coords[0]
        decpars = ax1.coords[1]
        rapars.display_minor_ticks(minorticks); rapars.set_minor_frequency(5)
        decpars.display_minor_ticks(minorticks)
        rapars.set_major_formatter('hh:mm:ss')
        decpars.set_major_formatter('dd:mm:ss')
        rapars.set_separator(('$^\\mathrm{H}$', "'", '"'))
        rapars.set_ticks(number=6, size=8, color=tickcolor)
        decpars.set_ticks(number=6, size=8, color=tickcolor)
        rapars.set_ticklabel(size=10, color=labelcolor); decpars.set_ticklabel(size=10, color=labelcolor)
        ax1.set_xlabel(xaxislabel, color=labelcolor)
        ax1.set_ylabel(yaxislabel, color=labelcolor)
        if legend:
            from matplotlib.patches import Patch
            swatch_edge = 'black' if legend_inverse else 'white'
            text_color = 'black' if legend_inverse else 'white'
            frame_bg = 'white' if legend_inverse else 'black'
            handles = [Patch(facecolor=c, edgecolor=swatch_edge, linewidth=0.5, label=lbl)
                       for c, lbl in legend]
            leg = ax1.legend(handles=handles, loc=legend_loc, framealpha=0.6,
                             facecolor=frame_bg, edgecolor=swatch_edge, fontsize='small')
            for txt in leg.get_texts():
                txt.set_color(text_color)
        if band_labels:
            from ..figures import add_band_labels
            add_band_labels(ax1, band_labels, loc=band_labels_loc)
        if swatch:
            from ..figures import add_swatch_inset
            add_swatch_inset(ax1, swatch, loc=swatch_loc, show_labels=swatch_labels,
                             label_offset=swatch_label_offset,
                             inset_scale=swatch_inset_scale)
        if show_compass or show_beam or show_scale_bar:
            from ..overlays import apply_overlays
            apply_overlays(ax1, hdrin, compass=show_compass, beam=show_beam,
                           scale_bar=show_scale_bar, compass_loc=compass_loc,
                           beam_loc=beam_loc, beam_style=beam_style,
                           scale_bar_asec=scale_bar_asec, scale_bar_loc=scale_bar_loc)
    plt.savefig(savepath, bbox_inches='tight', dpi=dpi, facecolor=fig1.get_facecolor())
    plt.clf(); plt.close('all')


def compare_multicolor_vs_rgb(rgbin, multicolorin, hdrin, ax2title, suptitle, savepath,
                                 xaxislabel='RA', yaxislabel='DEC', tickcolor='w', labelcolor='k',
                                 facecolor='w', supy=.8, dpi=150, minorticks=True):
    """
    Side-by-side comparison of a conventional RGB stack and a multicolor
    composite; saves the figure to disk.

    Tick and label colors accept any matplotlib color (``'k'``, ``'black'``, …).

    Parameters
    ----------
    rgbin : array
        Input pure RGB image  -->  e.g.,   np.dstack( [redframe, gframe, bframe] )
    multicolorin : array
        Multicolor RGB image -- as would be output from combine_multicolor()
    hdrin : astropy.io.fits.header
        Header to use for WCS information
    ax2title : str
        String to use for 2nd axis plot title.  e.g. "My crazy 5-color image", or empty quotes "" for nothing.
    suptitle : str
        String for super title (centered at top of plot, super title applies to whole figure)
    savepath : str
        Path to save file to.  e.g., "./plots/mycrazyimage.pdf"
    xaxislabel : str
        Label to use for x-axes. Default='RA'
    yaxislabel : str
        Label to use for y-axes. Default='DEC'
    tickcolor : str
        Color to use for ticks in plot, default = 'w'.
    labelcolor : str
        Color to use for ticklabels, default = 'k'
    facecolor : str
        Color to use for figure facecolor (border around plot), default ='w'
    minorticks : bool
        True to use minor ticks
    dpi : int
        Default=150. Dots per inch value to use for saved plot.
    supy : float
        Position for suptitle (default = 0.8)
    """
    from matplotlib import pyplot as plt
    plt.rcParams.update({'font.family': 'serif', 'xtick.major.size': 6, 'ytick.major.size': 6,
                         'xtick.major.width': 1., 'ytick.major.width': 1.,
                         'xtick.direction': 'in', 'ytick.direction': 'in',
                         'xtick.top': True, 'ytick.right': True})
    fig1 = plt.figure(1, figsize=(10, 8))
    fig1.set_facecolor(facecolor)
    wcs = WCS(hdrin)
    ax1 = fig1.add_subplot(121, projection=wcs)
    ax1.imshow(rgbin, origin='lower', interpolation='nearest')
    ax1.set_title('Simple RGB', color=labelcolor)
    ax2 = fig1.add_subplot(122, projection=wcs)
    ax2.imshow(multicolorin, origin='lower', interpolation='nearest')
    for ax in [ax1, ax2]:
        rapars = ax.coords[0]
        decpars = ax.coords[1]
        rapars.set_ticks(number=6, color=tickcolor); decpars.set_ticks(number=6, color=tickcolor)
        rapars.set_ticklabel(size=8, color=labelcolor)
        decpars.set_ticklabel(size=8, color=labelcolor)
        rapars.display_minor_ticks(minorticks)
        decpars.display_minor_ticks(minorticks)
        rapars.set_major_formatter('hh:mm:ss')
        decpars.set_major_formatter('dd:mm:ss')
        rapars.set_separator(('$^\\mathrm{H}$', "'", '"'))
        ax.set_xlabel(xaxislabel, color=labelcolor)
        ax.set_ylabel(yaxislabel, color=labelcolor)
    plt.subplots_adjust(wspace=0.3)
    ax2.set_title(ax2title, color=labelcolor)
    for ax in [ax1, ax2]:
        ax.coords.frame.set_color(tickcolor)  # This way is particular to astropy wcs frames
    plt.suptitle(suptitle, y=supy, color=labelcolor)
    plt.savefig(savepath, bbox_inches='tight', dpi=dpi, facecolor=fig1.get_facecolor())
    plt.clf(); plt.close('all')


def annotate_provenance_header(hdr, info=None, version=None):
    """
    Add HISTORY / COMMENT cards noting that a FITS file was created with
    multicolorfits (version, timestamp, optional per-layer settings).

    Parameters
    ----------
    hdr : astropy.io.fits.header
        Header modified in place.
    info : dict or None
        Optional state with ``'panels'`` and ``'compose'`` keys (as from
        :meth:`McfSession.to_dict`).
    version : str or None
        Package version; defaults to ``multicolorfits.__version__``.
    """
    from datetime import datetime, timezone
    try:
        from .. import __version__ as pkg_version
    except ImportError:
        pkg_version = 'unknown'
    ver = version or pkg_version
    now = datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')
    hdr.add_history('Created with multicolorfits v%s on %s' % (ver, now))
    if not info:
        return hdr
    compose = info.get('compose', {})
    if compose.get('gamma') is not None:
        hdr.add_history('  gamma = %.2f' % float(compose['gamma']))
    if compose.get('combine_mode'):
        hdr.add_history('  combine_mode = %s (blend=%s)' % (
            compose.get('combine_mode', 'rgb'), compose.get('combine_blend', 'screen')))
    for i, p in enumerate(info.get('panels', []), start=1):
        if not p.get('loaded'):
            continue
        hdr.add_history('  layer %d: scale=%s vmin=%.4g vmax=%.4g color=%s' % (
            i, p.get('stretch', '?'), p.get('vmin', 0), p.get('vmax', 0),
            p.get('color', '?')))
        if p.get('smooth'):
            hdr.add_history('    smooth sigma=%.2f px' % float(p.get('smooth_sigma', 0)))
    return hdr


def save_rgb_fits(savepath, multicolorRGBdat, commonhdr, overwrite=True, annotate=True):
    """
    Write a multicolor RGB cube to a FITS file with the shared WCS header.

    Parameters
    ----------
    savepath : str
        Output path (e.g. ``'./fits/composite.fits'``).
    multicolorRGBdat : array
        RGB image as returned by :func:`combine_multicolor` (shape ny×nx×3).
    commonhdr : astropy.io.fits.header
        Header supplying the WCS for the output.
    overwrite : bool
        Passed to ``astropy.io.fits.writeto``.
    annotate : bool
        If True (default), ensure the header records multicolorfits provenance
        when it does not already.
    """
    hdr = commonhdr.copy()
    if annotate:
        _ensure_mcf_provenance(hdr)
    pyfits.writeto(savepath, np.swapaxes(np.swapaxes(multicolorRGBdat, 0, 2), 2, 1),
                   hdr, overwrite=overwrite)


def _ensure_mcf_provenance(hdr):
    """Add a COMMENT if the header does not already mention multicolorfits."""
    try:
        parts = []
        for card in hdr.cards:
            if card[0] in ('COMMENT', 'HISTORY'):
                parts.append(str(card[1]))
        blob = ' '.join(parts).lower()
    except Exception:
        blob = ''
    if 'multicolorfits' in blob:
        return hdr
    try:
        from .. import __version__ as ver
    except ImportError:
        ver = 'unknown'
    hdr.add_comment('Created with multicolorfits v%s' % ver)
    return hdr
