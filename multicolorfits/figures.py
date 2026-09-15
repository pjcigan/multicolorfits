"""
Shared matplotlib figure construction for the combined WCS plot.

Used by the web GUI (rendered server-side to PNG / saved to disk) and the Qt
GUI (drawn onto an embedded canvas), so both produce identical output.
"""

import numpy as np

__all__ = ['setup_combined_axes', 'make_combined_figure', 'contrast_color',
           'add_swatch_inset', 'add_band_labels', 'default_axis_labels',
           'axis_labels_for_compose', 'apply_bare_axes', 'apply_bare_plot_style',
           'refresh_wcs_ticklabels', 'component_slot_grid', 'make_component_mosaic']

_TRANSPARENT = ('none', 'transparent', '')


def _is_transparent(facecolor):
    return str(facecolor).strip().lower() in _TRANSPARENT


def contrast_color(facecolor):
    """
    Pick a readable margin-text color (black or white) for a given canvas
    background.  Transparent backgrounds default to black (publication-safe).
    """
    if _is_transparent(facecolor):
        return 'black'
    try:
        from matplotlib.colors import to_rgb
        r, g, b = to_rgb(facecolor)
    except Exception:
        return 'black'
    luminance = 0.2126 * r + 0.7152 * g + 0.0722 * b
    return 'black' if luminance > 0.5 else 'white'


def default_axis_labels(hdr):
    """Default axis titles from a FITS header (CTYPE), or ``x`` / ``y``."""
    xlabel, ylabel = 'x', 'y'
    try:
        if hdr is not None:
            if '-' in hdr['CTYPE1']:
                xlabel = hdr['CTYPE1'].split('-')[0]
            if '-' in hdr['CTYPE2']:
                ylabel = hdr['CTYPE2'].split('-')[0]
    except Exception:
        pass
    return xlabel, ylabel


def axis_labels_for_compose(hdr, compose):
    """Resolve axis labels, honoring user overrides on ``compose`` when set."""
    xlabel, ylabel = default_axis_labels(hdr)
    if (compose.xlabel or '').strip():
        xlabel = compose.xlabel.strip()
    if (compose.ylabel or '').strip():
        ylabel = compose.ylabel.strip()
    return xlabel, ylabel


def _tick_style(compose):
    return {
        'major_size': float(getattr(compose, 'tick_major_size', 8.0) or 8.0),
        'minor_size': float(getattr(compose, 'tick_minor_size', 4.0) or 4.0),
        'major_width': float(getattr(compose, 'tick_major_width', 1.0) or 1.0),
        'minor_width': float(getattr(compose, 'tick_minor_width', 0.5) or 0.5),
        'direction': str(getattr(compose, 'tick_direction', 'in') or 'in'),
    }


def _wcs_tick_direction(direction):
    """
    WCSAxes tick helpers only accept 'in' or 'out'.  Approximate 'inout' as 'in'
    (the classic FITS style); plain matplotlib axes still get true inout.
    """
    d = str(direction or 'in').lower()
    return 'in' if d == 'inout' else d


def _apply_wcs_axis_titles(ax, xlabel, ylabel, text_color):
    """Set RA/Dec axis titles on WCS axes (robust for rotated fields).

    Tick *values* may appear on all four edges (``bltr``) for rotated WCS
    such as NGC 602.  Axis *titles* stay pinned to bottom / left so ``RA``
    and ``DEC`` do not stack on the same edge.
    """
    rapars = ax.coords[0]
    decpars = ax.coords[1]
    rapars.set_auto_axislabel(False)
    decpars.set_auto_axislabel(False)
    rapars.set_axislabel(xlabel, color=text_color, size=11)
    decpars.set_axislabel(ylabel, color=text_color, size=11)
    try:
        rapars.set_axislabel_position('b')
        decpars.set_axislabel_position('l')
    except (AttributeError, TypeError, ValueError):
        pass
    rapars.set_axislabel_visibility_rule('always')
    decpars.set_axislabel_visibility_rule('always')
    try:
        rapars.axislabels.set_clip_on(False)
        decpars.axislabels.set_clip_on(False)
        rapars.ticklabels.set_clip_on(False)
        decpars.ticklabels.set_clip_on(False)
    except AttributeError:
        pass


def _style_wcs_axes(ax, compose, text_color='black'):
    """Apply tick color / coordinate formatting from a ComposeState to WCS axes."""
    tickcolor = compose.tickcolor
    ts = _tick_style(compose)
    tick_kw = dict(number=6, size=ts['major_size'], width=ts['major_width'],
                   color=tickcolor, direction=_wcs_tick_direction(ts['direction']))
    try:
        ax.coords.frame.set_color(tickcolor)
    except AttributeError:
        for s in ax.spines.values():
            s.set_color(tickcolor)
        ax.tick_params(axis='both', which='major', color=tickcolor, labelcolor=text_color,
                       length=ts['major_size'], width=ts['major_width'],
                       direction=ts['direction'])
        ax.tick_params(axis='both', which='minor', color=tickcolor,
                       length=ts['minor_size'], width=ts['minor_width'])
        return
    rapars = ax.coords[0]
    decpars = ax.coords[1]
    rapars.set_ticks(**tick_kw)
    decpars.set_ticks(**tick_kw)
    # Tick-number labels sit in the margin; contrast them with the canvas.
    rapars.set_ticklabel(color=text_color, size=10)
    decpars.set_ticklabel(color=text_color, size=10)
    rapars.set_ticklabel_visible(True)
    decpars.set_ticklabel_visible(True)
    rapars.set_ticks_visible(True)
    decpars.set_ticks_visible(True)
    # Show tick values on every frame edge where WCS places them.  For rotated
    # images (e.g. NGC 602) RA/Dec labels may appear on left/right rather than
    # bottom/left; pinning to 'b'/'l' hides them entirely.
    rapars.set_ticklabel_position('bltr')
    decpars.set_ticklabel_position('bltr')
    rapars.set_ticks_position('bltr')
    decpars.set_ticks_position('bltr')
    try:
        rapars.ticklabels.set_clip_on(False)
        decpars.ticklabels.set_clip_on(False)
    except AttributeError:
        pass
    rapars.set_minor_frequency(5)
    decpars.set_minor_frequency(5)
    rapars.display_minor_ticks(compose.minorticks)
    decpars.display_minor_ticks(compose.minorticks)
    if compose.minorticks:
        rapars.tick_params(which='minor', length=ts['minor_size'])
        decpars.tick_params(which='minor', length=ts['minor_size'])
    if compose.coord_style.lower().startswith('sex'):
        rapars.set_major_formatter(compose.x_format or 'hh:mm:ss.ss')
        decpars.set_major_formatter(compose.y_format or 'dd:mm:ss.ss')
        rapars.set_separator(('$^\\mathrm{H}$ ', "' ", '" '))
    else:
        rapars.set_major_formatter(compose.x_format if ':' not in (compose.x_format or '') and compose.x_format else 'd.dddddd')
        decpars.set_major_formatter(compose.y_format if ':' not in (compose.y_format or '') and compose.y_format else 'd.dddddd')


def refresh_wcs_ticklabels(ax, compose, xlabel=None, ylabel=None):
    """
    Re-apply WCS coordinate tick styling after a canvas draw.

    Some interactive backends (notably QtAgg) drop or clip tick *values* unless
    the coordinate helpers are configured again post-draw.
    """
    if not hasattr(ax, 'coords'):
        return
    facecolor = getattr(compose, 'facecolor', 'white') or 'white'
    text_color = contrast_color(facecolor)
    _style_wcs_axes(ax, compose, text_color=text_color)
    if xlabel is not None and ylabel is not None:
        _apply_wcs_axis_titles(ax, xlabel, ylabel, text_color)


def apply_bare_axes(ax, fig=None, *, fill=True):
    """
    Hide ticks, coordinate labels, title, and the WCS frame on *ax*.

    Used for publication / image-only figures that should show only the
    composite pixels (optionally with zero margin when *fill* is True).
    """
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_title('')
    try:
        ax.coords.frame.set_color('none')
        ax.coords.frame.set_linewidth(0)
        for i in range(ax.coords.naxis):
            ax.coords[i].set_ticks_visible(False)
            ax.coords[i].set_ticklabel_visible(False)
            ax.coords[i].set_axislabel('')
            ax.coords[i].set_auto_axislabel(False)
    except AttributeError:
        pass
    ax.axis('off')
    if fill and fig is not None:
        ax.set_position([0, 0, 1, 1])
        fig.subplots_adjust(left=0, right=1, top=1, bottom=0)


def apply_bare_plot_style(compose, *, transparent=False):
    """
    Configure *compose* for publication / image-only output.

    Hides ticks, axis labels, title, legend, swatch, band labels, and
    skyplothelper overlays.  When *transparent* is True, also sets
    ``facecolor='none'`` so saved PNGs have no margin matting.

    Returns *compose* (same object, mutated in place) for chaining.

    Examples
    --------
    >>> import multicolorfits as mcf
    >>> s = mcf.McfSession()
    >>> mcf.apply_bare_plot_style(s.compose, transparent=True)
    """
    compose.bare_plot = True
    if transparent:
        compose.facecolor = 'none'
    compose.show_legend = False
    compose.show_combo_swatch = False
    compose.show_band_labels = False
    compose.show_compass = False
    compose.show_beam = False
    compose.show_scale_bar = False
    return compose


def setup_combined_axes(fig, session, combined=None):
    """
    Add a styled WCS-projected axes displaying the combined image to a figure.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Target figure (cleared first).
    session : McfSession
        Session providing panels and compose state.
    combined : array or None
        Pre-rendered combined RGB image; computed from the session if None.

    Returns
    -------
    matplotlib axes
    """
    if combined is None:
        combined = session.render_combined()
    hdr = session.common_header  # McfHeader (or None if no active panels)
    compose = session.compose
    facecolor = getattr(compose, 'facecolor', 'white') or 'white'
    transparent = _is_transparent(facecolor)
    text_color = contrast_color(facecolor)
    fig.clf()
    # Canvas (figure) background: 'none' => transparent (shows through on save).
    fig.set_facecolor('none' if transparent else facecolor)
    ax = None
    if hdr is not None:
        try:
            ax = fig.add_subplot(111, aspect=1, projection=hdr.wcs)
        except Exception:
            ax = None
    if ax is None:
        ax = fig.add_subplot(111, aspect=1)
    # The image fills the axes, but keep the axes patch transparent too so a
    # 'none' canvas keeps transparent margins on save.
    ax.patch.set_alpha(0.0 if transparent else 1.0)
    ax.imshow(np.clip(np.nan_to_num(combined), 0, 1), origin='lower', interpolation='nearest')

    if getattr(compose, 'bare_plot', False):
        apply_bare_axes(ax, fig)
        return ax

    # Axis labels from the header CTYPEs (overridable in compose state).
    hdr_ast = hdr.astropy if hdr is not None else None
    xlabel, ylabel = axis_labels_for_compose(hdr_ast, compose)
    if hasattr(ax, 'coords'):
        _apply_wcs_axis_titles(ax, xlabel, ylabel, text_color)
    else:
        ax.set_xlabel(xlabel, color=text_color)
        ax.set_ylabel(ylabel, color=text_color)
    if compose.title:
        ax.set_title(compose.title, color=text_color)
    _style_wcs_axes(ax, compose, text_color=text_color)
    if getattr(compose, 'show_legend', False):
        _add_channel_legend(ax, session, compose)
    if getattr(compose, 'show_combo_swatch', False):
        _add_combo_swatch(ax, session, compose)
    if getattr(compose, 'show_band_labels', False):
        _add_band_labels(ax, session, compose)
    if hdr is not None and (
        getattr(compose, 'show_compass', False)
        or getattr(compose, 'show_beam', False)
        or getattr(compose, 'show_scale_bar', False)
    ):
        from .overlays import apply_compose_overlays
        apply_compose_overlays(ax, hdr, compose)
    if not getattr(compose, 'bare_plot', False):
        # Leave margin room for WCS tick labels (especially on Qt canvases).
        fig.subplots_adjust(left=0.14, bottom=0.14, right=0.97, top=0.95)
    return ax


def _corner_xy(loc, pad):
    """(x, y, va, ha) in axes fraction for a named corner ('upper left', ...)."""
    loc = str(loc or 'upper left').lower()
    x = pad if 'left' in loc else 1.0 - pad
    y = 1.0 - pad if 'upper' in loc else pad
    va = 'top' if 'upper' in loc else 'bottom'
    ha = 'left' if 'left' in loc else 'right'
    return x, y, va, ha


def add_swatch_inset(ax, data, loc='lower right', show_labels=False, label_offset=0.62,
                     inset_scale=0.24):
    """
    Draw a pre-rendered color-combination swatch (from ``session.render_combo_swatch``
    / ``mcf.combo_swatch``) into a corner inset of ``ax`` whose frame is invisible,
    so it floats over the figure.  Optionally label each circle in its color.

    label_offset : float
        How far to push each circle's label radially outward from the cluster
        centre, as a fraction of the circle radius.  Circle centres coincide with
        the multi-circle overlap nodes, so the labels are shifted into each
        circle's single-color outer lobe to label it unambiguously.  0 keeps the
        label at the circle centre; ~0.6 (default) sits in the outer lobe.
    inset_scale : float
        Width and height of the inset as a fraction of the parent axes (0.08–0.45).
        Default 0.24 (~24% of the plot).
    """
    if not data:
        return
    rgba = data['rgba']
    sz = data['size']
    w = h = float(max(0.08, min(0.45, inset_scale or 0.24)))
    pad = 0.015
    loc = str(loc or 'lower right').lower()
    x0 = pad if 'left' in loc else 1.0 - w - pad
    y0 = 1.0 - h - pad if 'upper' in loc else pad
    iax = ax.inset_axes([x0, y0, w, h])
    iax.imshow(rgba, origin='lower', interpolation='bilinear', extent=[0, sz, 0, sz])
    iax.set_xlim(0, sz)
    iax.set_ylim(0, sz)
    iax.set_aspect('equal')
    iax.axis('off')
    iax.patch.set_alpha(0.0)
    if show_labels:
        import matplotlib.patheffects as pe
        cx0 = cy0 = sz / 2.0
        for c in data['circles']:
            dx, dy = c['x'] - cx0, c['y'] - cy0
            nrm = (dx * dx + dy * dy) ** 0.5
            # Single circle (n==1): no cluster direction, so leave at centre.
            ux, uy = (dx / nrm, dy / nrm) if nrm > 1e-9 else (0.0, 0.0)
            tx = c['x'] + ux * label_offset * c['r']
            ty = c['y'] + uy * label_offset * c['r']
            iax.text(tx, ty, c['label'], color=c['color'],
                     fontsize=max(5.0, 7.0 * w / 0.24),
                     ha='center', va='center', clip_on=False,
                     path_effects=[pe.withStroke(linewidth=1.4,
                                                 foreground=contrast_color(c['color']))])


def add_band_labels(ax, entries, loc='upper left'):
    """
    Draw per-band colored text labels chained in a corner (skyplothelper
    add_bandlabels style): each ``(color, label)`` entry is drawn in its color, so
    viewers can read which band maps to which color.  A thin contrasting outline
    keeps the text legible over any part of the image.
    """
    from matplotlib.text import OffsetFrom
    import matplotlib.patheffects as pe
    if not entries:
        return
    colors = [c for c, _ in entries]
    labels = [lbl for _, lbl in entries]
    pad = 0.035
    x, y, va, ha = _corner_xy(loc, pad)
    yfrac = {'top': 1.0, 'center': 0.5, 'bottom': 0.0, 'baseline': 0.0}.get(va, 1.0)
    textpad = 0.08

    def stroke(color):
        return [pe.withStroke(linewidth=1.8, foreground=contrast_color(color))]

    texts = []
    an = ax.annotate(labels[0], xy=(x, y), xycoords=ax.transAxes, color=colors[0],
                     va=va, ha=ha, fontsize=11, fontweight='bold',
                     path_effects=stroke(colors[0]))
    texts.append(an)
    for i in range(1, len(labels)):
        if ha == 'right':
            off = OffsetFrom(texts[i - 1], (0.0 - textpad, yfrac))
            iha = 'right'
        else:
            off = OffsetFrom(texts[i - 1], (1.0 + textpad, yfrac))
            iha = 'left'
        an = ax.annotate(labels[i], xy=(0, 0), xycoords=off, color=colors[i],
                         va=va, ha=iha, fontsize=11, fontweight='bold',
                         path_effects=stroke(colors[i]))
        texts.append(an)


def _add_combo_swatch(ax, session, compose):
    """Session-aware wrapper: render + place the composite-mode-aware swatch."""
    add_swatch_inset(ax, session.render_combo_swatch(),
                     loc=getattr(compose, 'combo_swatch_loc', 'lower right'),
                     show_labels=getattr(compose, 'combo_swatch_labels', False),
                     label_offset=float(getattr(compose, 'combo_swatch_label_offset', 0.62)),
                     inset_scale=float(getattr(compose, 'combo_swatch_inset_scale', 0.24)))


def _add_band_labels(ax, session, compose):
    """Session-aware wrapper: draw colored per-band labels in a corner."""
    add_band_labels(ax, session.legend_entries(),
                    loc=getattr(compose, 'band_labels_loc', 'upper left'))


def _add_channel_legend(ax, session, compose):
    """Add a 'which color = which image' legend over the combined image."""
    from matplotlib.patches import Patch
    entries = session.legend_entries()
    if not entries:
        return
    # Legend sits over the image data: pick text/frame to contrast the
    # compositing background (white/color/transparent) or legacy inverse.
    bg_key = str(getattr(compose, 'combine_background', 'black') or 'black').strip().lower()
    if bg_key in ('transparent', 'none'):
        eff_bg = getattr(compose, 'facecolor', 'white') or 'white'
    elif bg_key not in ('black', ''):
        eff_bg = getattr(compose, 'combine_background')
    else:
        eff_bg = 'white' if getattr(compose, 'inverse', False) else 'black'
    text_color = contrast_color(eff_bg)
    swatch_edge = text_color
    frame_bg = eff_bg if str(eff_bg).strip().lower() not in ('none', 'transparent', '') else 'white'
    handles = [Patch(facecolor=color, edgecolor=swatch_edge, linewidth=0.5, label=label)
               for color, label in entries]
    leg = ax.legend(handles=handles, loc=getattr(compose, 'legend_loc', 'upper right'),
                    framealpha=0.6, facecolor=frame_bg, edgecolor=swatch_edge,
                    fontsize='small')
    for txt in leg.get_texts():
        txt.set_color(text_color)


def make_combined_figure(session, combined=None, figsize=(7, 7), facecolor='w'):
    """
    Build a standalone matplotlib Figure of the combined WCS plot
    (independent of any GUI canvas; safe with the Agg backend).

    Returns the figure only. The combined axes is ``fig.axes[0]`` (grab it
    before adding insets). This is not a session panel — ``session.panels``
    are the input layers. Save with ``fig.savefig``, not ``plt.savefig``.
    To have the axes handed back, use :func:`setup_combined_axes`.
    """
    from matplotlib.figure import Figure
    fig = Figure(figsize=figsize, facecolor=facecolor)
    setup_combined_axes(fig, session, combined=combined)
    return fig


# --------------------------------------------------------------------------- component mosaic


def component_slot_grid(n, max_per_line=3, side='top'):
    """
    Pack *n* component indices into a strip grid next to a hero panel.

    Lines closest to the hero fill first; leftover cells in the outer line are
    ``None`` (empty slots for annotations).  Returned rows/columns are ordered
    for figure layout (outer line first when the strip is above or left of the
    hero; closest line first when below or right).

    Parameters
    ----------
    n : int
        Number of component panels.
    max_per_line : int
        Max panels per strip row (``side`` in top/bottom) or column (left/right).
    side : {'top', 'bottom', 'left', 'right'}
        Where the component strip sits relative to the combined hero.

    Returns
    -------
    grid : list of list
        For top/bottom: list of rows, each length *max_per_line*, entries are
        panel indices or ``None``.
        For left/right: list of columns (same structure; caller treats as cols).
    """
    side = str(side or 'top').strip().lower()
    if side not in ('top', 'bottom', 'left', 'right'):
        raise ValueError("side must be 'top', 'bottom', 'left', or 'right'")
    max_per_line = max(1, int(max_per_line))
    n = max(0, int(n))
    # Fill closest-to-hero first.
    closest_first = []
    remaining = list(range(n))
    while remaining:
        chunk = remaining[:max_per_line]
        remaining = remaining[max_per_line:]
        closest_first.append(chunk + [None] * (max_per_line - len(chunk)))
    if not closest_first:
        return []
    # Figure order: top/left put the outer (furthest) line first in the grid.
    if side in ('top', 'left'):
        return list(reversed(closest_first))
    return closest_first


def _mosaic_apply_ticks(ax, compose, mode, *, is_hero, text_color, xlabel, ylabel):
    """Apply plain / minimal / full tick policy to one mosaic axes."""
    mode = str(mode or 'plain').strip().lower()
    if mode in ('plain', 'none', 'off', 'bare'):
        apply_bare_axes(ax, fig=None, fill=False)
        return
    if not hasattr(ax, 'coords'):
        if is_hero and mode == 'full':
            ax.set_xlabel(xlabel, color=text_color)
            ax.set_ylabel(ylabel, color=text_color)
        elif is_hero and mode == 'minimal':
            ax.set_xlabel(xlabel, color=text_color)
            ax.set_ylabel(ylabel, color=text_color)
        else:
            ax.set_xticks([])
            ax.set_yticks([])
            if not is_hero:
                ax.set_xlabel('')
                ax.set_ylabel('')
        return

    if is_hero or mode == 'full':
        _style_wcs_axes(ax, compose, text_color=text_color)
        if is_hero or mode == 'full':
            _apply_wcs_axis_titles(ax, xlabel, ylabel, text_color)
        return

    # minimal + component: same ticks as hero styling, hide labels / titles.
    _style_wcs_axes(ax, compose, text_color=text_color)
    try:
        nax = len(ax.coords)
    except TypeError:
        nax = 2
    for i in range(nax):
        ax.coords[i].set_ticklabel_visible(False)
        ax.coords[i].set_axislabel('')
        ax.coords[i].set_auto_axislabel(False)


def _mosaic_add_component_label(ax, label, color, loc='upper left'):
    """
    Colored band label for a component panel.

    *loc* is an interior corner (``'upper left'`` default), ``'title'`` for an
    external axes title, or ``'none'`` to skip.
    """
    import matplotlib.patheffects as pe
    loc = str(loc or 'upper left').strip().lower()
    if not label or loc in ('none', 'off', ''):
        return
    stroke = [pe.withStroke(linewidth=1.8, foreground=contrast_color(color))]
    if loc == 'title':
        title = ax.set_title(str(label), color=color, fontsize=10, pad=3)
        title.set_path_effects(stroke)
        return
    x, y, va, ha = _corner_xy(loc, pad=0.045)
    txt = ax.text(x, y, str(label), transform=ax.transAxes, color=color,
                  fontsize=9, fontweight='bold', va=va, ha=ha, zorder=10,
                  clip_on=False)
    txt.set_path_effects(stroke)


def _mosaic_block_rects(side, n_lines, ncol, img_aspect, cell_in, gap_in, pad_in):
    """
    Figure size and strip/hero rectangles (figure fraction) for an ImageGrid strip.

    Strip outer box matches ``n_lines × ncol`` square cells of side *cell_in*
    with *gap_in* padding — the geometry :class:`~mpl_toolkits.axes_grid1.ImageGrid`
    expects when ``aspect=True``.
    """
    n_lines = max(int(n_lines), 1)
    ncol = max(int(ncol), 1)
    gap_in = max(float(gap_in), 0.0)
    pad_in = max(float(pad_in), 0.0)
    cell_in = max(float(cell_in), 0.5)
    img_aspect = float(img_aspect) if img_aspect and img_aspect > 0 else 1.0

    if side in ('top', 'bottom'):
        strip_w = ncol * cell_in + (ncol - 1) * gap_in
        strip_h = n_lines * cell_in + (n_lines - 1) * gap_in
        hero_w = strip_w
        hero_h = hero_w * img_aspect
        total_w = strip_w + 2 * pad_in
        total_h = strip_h + gap_in + hero_h + 2 * pad_in
        left = pad_in / total_w
        width = strip_w / total_w
        if side == 'top':
            strip_bottom = (pad_in + hero_h + gap_in) / total_h
            hero_bottom = pad_in / total_h
        else:
            hero_bottom = (pad_in + strip_h + gap_in) / total_h
            strip_bottom = pad_in / total_h
        strip_rect = [left, strip_bottom, width, strip_h / total_h]
        hero_rect = [left, hero_bottom, width, hero_h / total_h]
    else:
        strip_h = ncol * cell_in + (ncol - 1) * gap_in
        strip_w = n_lines * cell_in + (n_lines - 1) * gap_in
        hero_h = strip_h
        hero_w = hero_h / img_aspect
        total_h = strip_h + 2 * pad_in
        total_w = strip_w + gap_in + hero_w + 2 * pad_in
        bottom = pad_in / total_h
        height = strip_h / total_h
        if side == 'left':
            strip_left = pad_in / total_w
            hero_left = (pad_in + strip_w + gap_in) / total_w
        else:
            hero_left = pad_in / total_w
            strip_left = (pad_in + hero_w + gap_in) / total_w
        strip_rect = [strip_left, bottom, strip_w / total_w, height]
        hero_rect = [hero_left, bottom, hero_w / total_w, height]

    return (total_w, total_h), strip_rect, hero_rect


def _mosaic_axes_class(wcs):
    """ImageGrid ``axes_class`` for shared WCS (same pattern as sph channel_map)."""
    if wcs is None:
        return None
    try:
        from astropy.visualization.wcsaxes import WCSAxes
    except ImportError:
        return None
    return (WCSAxes, {'wcs': wcs})


def make_component_mosaic(session, combined=None, *, components='top',
                          max_per_line=3, ticks='plain', show=None,
                          gap=0.06, panel_pad=0.2, cell_size=1.7,
                          figsize=None, facecolor=None, overlays='hero',
                          title=None, label_loc='upper left'):
    """
    Combined image plus colorized component frames in a shared-WCS mosaic.

    Default layout: a full-width **hero** (combined) panel with a strip of
    smaller single-layer renderings on one side.  Each strip panel is the
    display-ready layer that was fed to the combiner (that panel's stretch,
    vmin/vmax, color, and display gamma — the same image
    :meth:`~multicolorfits.session.PanelState.render_display` returns), not a
    fresh min/max stretch and not the darker pre-display buffer the mixer
    stores internally.  The strip is an
    :class:`~mpl_toolkits.axes_grid1.ImageGrid` (same approach as
    skyplothelper ``channel_map``): equal square cells, ``axes_pad`` in inches,
    optional shared WCS.  Component lines closest to the hero fill first;
    leftover empty cells stay blank for custom labels.

    Parameters
    ----------
    session : McfSession
    combined : array or None
        Pre-rendered combined RGB; computed if None.
    components : {'top', 'bottom', 'left', 'right'}
        Side of the hero for the component strip.
    max_per_line : int
        Max panels per strip row (top/bottom) or column (left/right). Default 3.
    ticks : {'plain', 'minimal', 'full'}
        ``plain`` — no ticks/labels (default); ``minimal`` — full ticks on the
        hero, matching ticks without labels on components; ``full`` — ticks and
        labels on every panel.
    show : sequence of int or None
        Active-panel indices to include (and order).  ``None`` = all active.
    gap : float
        Padding between adjacent strip panels in **inches** (ImageGrid
        ``axes_pad``).  Default 0.06.
    panel_pad : float
        Outer figure margin in **inches**.  Default 0.2.
    cell_size : float
        Side length in inches of each strip cell (sets default figsize).
        Default 1.7.
    figsize : (w, h) or None
        If None, sized tightly from *cell_size* / *gap* / image aspect.
    facecolor : color or None
        Canvas background; defaults to ``session.compose.facecolor``.
    overlays : {'hero', 'none'}
        Apply compose beam/compass/scale-bar overlays on the hero only, or skip.
    title : str or None
        Optional figure-level title (not per-axes).
    label_loc : str
        Component band-label placement: interior corner
        (``'upper left'`` default, also ``'upper right'`` / ``'lower left'`` /
        ``'lower right'``), ``'title'`` for an external axes title, or
        ``'none'`` to omit.

    Returns
    -------
    fig : matplotlib.figure.Figure
        Not registered as the pyplot current figure. Save with
        ``fig.savefig(...)``, not ``plt.savefig`` (that writes a blank canvas).
    axes : dict
        ``combined`` — hero axes (the composited image). This is not a
        session panel; ``session.panels`` are the input layers. Add a scale
        bar, compass, or swatch to ``axes['combined']`` after the call.
        ``components`` — list of component axes in *show*/active order.
        ``slots`` — 2D grid of axes or ``None`` matching
        :func:`component_slot_grid`.
    """
    from matplotlib.figure import Figure
    from mpl_toolkits.axes_grid1 import ImageGrid

    active = session.active_panels()
    if not active:
        raise ValueError('No images loaded in any panel')

    side = str(components or 'top').strip().lower()
    if side not in ('top', 'bottom', 'left', 'right'):
        raise ValueError("components must be 'top', 'bottom', 'left', or 'right'")
    tick_mode = str(ticks or 'plain').strip().lower()
    if tick_mode not in ('plain', 'minimal', 'full', 'none', 'off', 'bare'):
        raise ValueError("ticks must be 'plain', 'minimal', or 'full'")
    max_per_line = max(1, int(max_per_line))
    overlays_mode = str(overlays or 'hero').strip().lower()
    label_loc = str(label_loc if label_loc is not None else 'upper left')

    if show is None:
        panel_indices = list(session.active_indices())
    else:
        active_idx = list(session.active_indices())
        panel_indices = []
        for idx in show:
            i = int(idx)
            if i < 0 or i >= len(active_idx):
                raise IndexError('show index %s out of range for %d active panels'
                                 % (idx, len(active_idx)))
            panel_indices.append(active_idx[i])
    panels = [session.panels[i] for i in panel_indices]
    n = len(panels)
    slot_grid = component_slot_grid(n, max_per_line=max_per_line, side=side)
    n_lines = len(slot_grid)

    compose = session.compose
    fc = facecolor if facecolor is not None else (getattr(compose, 'facecolor', 'white') or 'white')
    transparent = _is_transparent(fc)
    text_color = contrast_color(fc)
    hdr = session.common_header
    hdr_ast = hdr.astropy if hdr is not None else None
    wcs = hdr.wcs if hdr is not None else None
    xlabel, ylabel = axis_labels_for_compose(hdr_ast, compose)

    if combined is None:
        combined = session.render_combined()
    combined = np.clip(np.nan_to_num(np.asarray(combined)), 0, 1)
    _, _, eff_inverse = session._resolve_background()
    gamma = compose.gamma
    # Display-ready single layers, not the pre-display colorize buffers.
    # Those buffers are gamma-encoded for the mixer and imshow as if the
    # panel were still on its load-time linear min/max stretch.
    component_rgbs = [
        p.render_display(gamma=gamma, inverse=eff_inverse)
        for p in panels
    ]

    img_aspect = combined.shape[0] / float(combined.shape[1])
    gap_in = float(gap)
    pad_in = float(panel_pad)
    if 0 < gap_in < 0.02:
        # Values under ~0.02" look like zero with ImageGrid pads; bump to a
        # readable default so "small gap" still separates panels.
        gap_in = 0.06
    if 0 < pad_in < 0.05:
        pad_in = 0.2
    cell_in = float(cell_size)

    def _add_hero(rect):
        if wcs is not None:
            try:
                return fig.add_axes(rect, projection=wcs)
            except Exception:
                pass
        ax = fig.add_axes(rect)
        ax.set_aspect('equal')
        return ax

    def _show(ax, rgb, *, is_hero, label=None, color=None):
        ax.patch.set_alpha(0.0 if transparent else 1.0)
        ax.imshow(rgb, origin='lower', interpolation='nearest')
        _mosaic_apply_ticks(ax, compose, tick_mode, is_hero=is_hero,
                            text_color=text_color, xlabel=xlabel, ylabel=ylabel)
        if (not is_hero) and label:
            _mosaic_add_component_label(ax, label, color or text_color, loc=label_loc)

    if n == 0:
        fig = Figure(figsize=figsize or (6, 6), facecolor=('none' if transparent else fc))
        ax_hero = _add_hero([0.08, 0.08, 0.84, 0.84])
        _show(ax_hero, combined, is_hero=True)
        if title:
            fig.suptitle(title, color=text_color)
        return fig, {'combined': ax_hero, 'components': [], 'slots': []}

    auto_size, strip_rect, hero_rect = _mosaic_block_rects(
        side, n_lines, max_per_line, img_aspect, cell_in, gap_in, pad_in)
    fig = Figure(figsize=figsize or auto_size, facecolor=('none' if transparent else fc))

    if side in ('top', 'bottom'):
        nrows, ncols = n_lines, max_per_line
    else:
        nrows, ncols = max_per_line, n_lines

    grid_kw = dict(
        nrows_ncols=(nrows, ncols),
        axes_pad=gap_in,
        share_all=True,
        aspect=True,
        label_mode='L',
        cbar_mode=None,
    )
    axes_class = _mosaic_axes_class(wcs)
    if axes_class is not None:
        grid_kw['axes_class'] = axes_class
    strip_grid = ImageGrid(fig, strip_rect, **grid_kw)

    slot_axes = [[None] * max_per_line for _ in range(n_lines)]
    component_axes = [None] * n
    ax_hero = _add_hero(hero_rect)

    if side in ('top', 'bottom'):
        for r, row in enumerate(slot_grid):
            for c, pidx in enumerate(row):
                ax = strip_grid[r * ncols + c]
                if pidx is None:
                    ax.set_visible(False)
                    slot_axes[r][c] = None
                    continue
                p = panels[pidx]
                lbl = p.label or ('Image %d' % (panel_indices[pidx] + 1))
                _show(ax, component_rgbs[pidx], is_hero=False,
                      label=lbl, color=p.color)
                slot_axes[r][c] = ax
                component_axes[pidx] = ax
    else:
        for c, col in enumerate(slot_grid):
            for r, pidx in enumerate(col):
                ax = strip_grid[r * ncols + c]
                if pidx is None:
                    ax.set_visible(False)
                    slot_axes[c][r] = None
                    continue
                p = panels[pidx]
                lbl = p.label or ('Image %d' % (panel_indices[pidx] + 1))
                _show(ax, component_rgbs[pidx], is_hero=False,
                      label=lbl, color=p.color)
                slot_axes[c][r] = ax
                component_axes[pidx] = ax

    _show(ax_hero, combined, is_hero=True)
    if compose.title and tick_mode not in ('plain', 'none', 'off', 'bare'):
        ax_hero.set_title(compose.title, color=text_color)
    if overlays_mode == 'hero' and hdr is not None and (
        getattr(compose, 'show_compass', False)
        or getattr(compose, 'show_beam', False)
        or getattr(compose, 'show_scale_bar', False)
    ):
        from .overlays import apply_compose_overlays
        apply_compose_overlays(ax_hero, hdr, compose)
    if getattr(compose, 'show_legend', False) and tick_mode not in ('plain',):
        _add_channel_legend(ax_hero, session, compose)
    if getattr(compose, 'show_combo_swatch', False):
        _add_combo_swatch(ax_hero, session, compose)
    if getattr(compose, 'show_band_labels', False):
        _add_band_labels(ax_hero, session, compose)

    if title:
        fig.suptitle(title, color=text_color, y=0.99)

    return fig, {
        'combined': ax_hero,
        'components': component_axes,
        'slots': slot_axes,
    }
