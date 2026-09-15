"""
Scripting-side palette preview: tiles, combo swatch, and CVD checks.

Keeps matplotlib out of :mod:`multicolorfits.palettes` at import time.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, List, Optional, Sequence, Union

import numpy as np

from .core.colorize import hex_to_rgb
from .core.colormix import simulate_colorblindness
from .core.swatch import combo_swatch
from .figures import contrast_color
from .palettes import (
    CVD_KINDS,
    palette_colorblind_report,
    resolve_palette_colors,
)

__all__ = ['PalettePreview', 'preview_palette']


@dataclass
class PalettePreview:
    """Result of :func:`preview_palette`.

    Attributes
    ----------
    fig : matplotlib.figure.Figure
        A **pyplot-managed** figure (so ``plt.show()`` works after the call).
        Save with ``fig.savefig(...)`` or :meth:`show`.
    colors : list of str
        Hex colors shown.
    labels : list of str
        Per-layer labels (may be empty strings).
    cvd : dict
        :func:`~multicolorfits.palettes.palette_colorblind_report` output.
    swatch : dict or None
        ``combo_swatch`` dict (``rgba``, ``circles``, ``size``).
    """
    fig: Any
    colors: List[str]
    labels: List[str]
    cvd: dict
    swatch: Optional[dict] = None
    extra: dict = field(default_factory=dict)

    def show(self, **kwargs):
        """Display this figure (interactive backends / IPython).

        Equivalent to making *fig* current and calling ``plt.show(**kwargs)``.
        In a notebook, ending a cell with ``res.fig`` or ``res.show()`` usually
        embeds the plot; with a GUI backend, ``plt.show()`` after
        ``preview_palette`` also works because the figure is pyplot-owned.
        """
        import matplotlib.pyplot as plt
        try:
            plt.figure(self.fig.number)
        except Exception:
            pass
        return plt.show(**kwargs)


def _hex_to_rgb01(hex_color):
    return np.array(hex_to_rgb(hex_color), dtype=float) / 255.0


def _tile_rgb(colors, height=48, width=64):
    """(n, height, width, 3) solid tiles for the given hex colors."""
    n = len(colors)
    out = np.zeros((n, height, width, 3), dtype=float)
    for i, c in enumerate(colors):
        out[i] = _hex_to_rgb01(c)
    return out


def _draw_tile_row(ax, colors, labels, background, title=None, face_rgbs=None):
    from matplotlib.colors import to_rgb

    n = len(colors)
    if face_rgbs is None:
        tiles = _tile_rgb(colors)
    else:
        h, w = 48, 64
        tiles = np.zeros((n, h, w, 3), dtype=float)
        for i in range(n):
            tiles[i] = np.asarray(face_rgbs[i], dtype=float).reshape(3)
    # Concatenate horizontally with a thin gutter matching the canvas.
    gutter = 2
    h, w = tiles.shape[1], tiles.shape[2]
    bg = np.array(to_rgb(background), dtype=float)
    canvas = np.broadcast_to(bg, (h, n * w + max(n - 1, 0) * gutter, 3)).copy()
    for i in range(n):
        x0 = i * (w + gutter)
        canvas[:, x0:x0 + w] = tiles[i]
    ax.imshow(canvas, origin='upper', interpolation='nearest')
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_facecolor(background)
    text_color = contrast_color(background)
    for i, (c, lab) in enumerate(zip(colors, labels)):
        x0 = i * (w + gutter)
        cx = x0 + w / 2.0
        label = (lab or '').strip()
        line = label + '\n' + c.upper() if label else c.upper()
        ax.text(cx, h * 0.55, line, ha='center', va='center',
                color=text_color, fontsize=8, fontfamily='monospace',
                linespacing=1.2)
    if title:
        ax.set_title(title, color=text_color, fontsize=10, loc='left', pad=4)


def _format_cvd_warning(kind, info, colors, labels):
    if info['ok']:
        return '%s: ok' % kind
    bits = []
    for i, j, d in info['failures']:
        a = (labels[i] or colors[i]).strip() or colors[i]
        b = (labels[j] or colors[j]).strip() or colors[j]
        bits.append('%s~%s (ΔE=%.1f)' % (a, b, d))
    return '%s: confusing — %s' % (kind, '; '.join(bits))


def preview_palette(
    colors,
    n=None,
    labels=None,
    *,
    background='black',
    mode='rgb',
    blend='screen',
    gamma=2.2,
    combine_background='black',
    inverse=False,
    cvd=True,
    min_distance=25.0,
    swatch_size=220,
    figsize=None,
    title=None,
):
    """
    Preview a palette: tiles, mode-aware combo swatch, and CVD checks.

    Parameters
    ----------
    colors : str or sequence of str
        Palette / hue-pattern name, or an explicit hex list. Names use
        :func:`~multicolorfits.palettes.resolve_palette_colors`.
    n : int or None
        Layer count when *colors* is a name. Ignored for an explicit list
        unless you want a single color broadcast.
    labels : sequence of str or None
        Optional per-layer names under each tile.
    background : color
        Canvas behind the tile row (``'black'``, ``'white'``, …).
    mode, blend, gamma, combine_background, inverse
        Passed to :func:`~multicolorfits.combo_swatch` so overlaps match
        the intended composite.
    cvd : bool or sequence of str
        If True (default), show deuteranopia / protanopia / tritanopia
        tile rows and warnings. False skips them. A sequence selects kinds.
    min_distance : float
        CVD ΔE threshold (same as :func:`check_palette_colorblind`).
    swatch_size : int
        Pixel size of the combo swatch render.
    figsize : (w, h) or None
    title : str or None
        Optional figure title.

    Returns
    -------
    PalettePreview
        ``fig``, ``colors``, ``labels``, ``cvd``, ``swatch``. The figure is
        created with ``pyplot.figure`` so it is the current figure —
        ``plt.show()`` works after the call. In notebooks, ``res.fig`` or
        ``res.show()`` also displays it. Save with ``res.fig.savefig(...)``.

    Examples
    --------
    ::

        import matplotlib.pyplot as plt
        import multicolorfits as mcf

        res = mcf.preview_palette('pob', n=3, mode='lab', blend='screen')
        plt.show()                 # interactive backends
        # res.show()               # same idea
        # res.fig.savefig('pal.png', dpi=120, facecolor=res.fig.get_facecolor())
    """
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec

    resolved = resolve_palette_colors(colors, n=n)
    if not resolved:
        raise ValueError('preview_palette: no colors to show')
    n_colors = len(resolved)
    if labels is None:
        label_list = [''] * n_colors
    else:
        label_list = [str(x) for x in list(labels)]
        if len(label_list) == 1 and n_colors != 1:
            label_list = label_list * n_colors
        if len(label_list) != n_colors:
            raise ValueError('labels has %d entries; expected %d'
                             % (len(label_list), n_colors))

    if cvd is True:
        kind_list = list(CVD_KINDS)
    elif cvd is False:
        kind_list = []
    else:
        kind_list = list(cvd)

    cvd_report = palette_colorblind_report(
        resolved, kinds=kind_list or CVD_KINDS, min_distance=min_distance)
    # If CVD panels are off, still return the full three-kind report for scripts.
    if not kind_list:
        cvd_report = palette_colorblind_report(
            resolved, min_distance=min_distance)

    swatch = combo_swatch(
        resolved,
        mode=mode,
        blend=blend,
        gamma=gamma,
        background=combine_background,
        inverse=inverse,
        size=int(swatch_size),
        labels=label_list,
    )

    n_cvd = len(kind_list)
    # rows: title spacer handled by suptitle; tiles; swatch; optional warnings; cvd rows
    n_rows = 2 + n_cvd
    height_ratios = [1.0, 1.4] + [0.85] * n_cvd
    if figsize is None:
        figsize = (max(6.0, 1.4 * n_colors + 2.5), 2.2 + 1.1 * n_cvd + 1.6)
    # pyplot-managed so plt.show() / IPython see the figure (unlike bare Figure()).
    fig = plt.figure(figsize=figsize, facecolor=background)
    gs = GridSpec(
        n_rows, 2, figure=fig,
        height_ratios=height_ratios,
        width_ratios=[1.2, 1.0],
        hspace=0.45, wspace=0.25,
        left=0.06, right=0.98, top=0.90, bottom=0.06,
    )

    text_color = contrast_color(background)
    ax_tiles = fig.add_subplot(gs[0, :])
    _draw_tile_row(ax_tiles, resolved, label_list, background,
                   title=None)
    ax_tiles.set_ylabel('layers', color=text_color, fontsize=9)

    ax_swatch = fig.add_subplot(gs[1, 0])
    ax_swatch.set_facecolor(background)
    if swatch is not None:
        rgba = swatch['rgba']
        # Composite onto the tile background so transparent margins match.
        from matplotlib.colors import to_rgb
        bg = np.array(to_rgb(
            combine_background if str(combine_background).lower()
            not in ('transparent', 'none') else background), dtype=float)
        rgb = rgba[..., :3] * rgba[..., 3:4] + bg * (1.0 - rgba[..., 3:4])
        ax_swatch.imshow(rgb, origin='lower', interpolation='bilinear')
    ax_swatch.set_xticks([])
    ax_swatch.set_yticks([])
    for spine in ax_swatch.spines.values():
        spine.set_visible(False)
    ax_swatch.set_title(
        'combo swatch (%s / %s)' % (mode, blend),
        color=text_color, fontsize=9, loc='left')

    ax_warn = fig.add_subplot(gs[1, 1])
    ax_warn.set_axis_off()
    ax_warn.set_facecolor(background)
    lines = ['CVD check (ΔE ≥ %.0f)' % float(min_distance)]
    for kind in (kind_list or CVD_KINDS):
        info = cvd_report['kinds'][kind]
        lines.append(_format_cvd_warning(kind, info, resolved, label_list))
    ax_warn.text(
        0.0, 1.0, '\n'.join(lines), transform=ax_warn.transAxes,
        va='top', ha='left', color=text_color, fontsize=8,
        fontfamily='monospace', wrap=True)

    for row, kind in enumerate(kind_list):
        ax = fig.add_subplot(gs[2 + row, :])
        rgb = np.array([_hex_to_rgb01(c) for c in resolved]).reshape(1, -1, 3)
        sim = simulate_colorblindness(rgb, kind=kind)[0]
        _draw_tile_row(ax, resolved, label_list, background, face_rgbs=sim)
        status = cvd_report['kinds'][kind]
        mark = 'ok' if status['ok'] else 'warn'
        ax.set_ylabel('%s\n[%s]' % (kind[:5], mark), color=text_color, fontsize=8)

    if title:
        fig.suptitle(title, color=text_color, fontsize=12)
    elif isinstance(colors, str):
        fig.suptitle('palette: %s' % colors, color=text_color, fontsize=11)

    return PalettePreview(
        fig=fig,
        colors=resolved,
        labels=label_list,
        cvd=cvd_report,
        swatch=swatch,
    )
