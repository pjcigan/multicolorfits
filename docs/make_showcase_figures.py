#!/usr/bin/env python3
"""
Render real-data showcase figures for the docs Examples pages.

Outputs land in ``docs/_static/showcase/`` as explicit ``*_light.png`` /
``*_dark.png`` pairs (rendered twice — no chrome remapping heuristics), plus
``m74_stamp.png``, a genuinely transparent RGBA cutout that lets the site
background show through in either theme.

Data sources (each block skips cleanly when its data is absent):

* NGC 602 — Chandra OpenFITS IR / R / B frames (``$MCF_NGC602_DIR`` or a
  local ``ngc602_*.fits`` cache).
* M74 (NGC 628) — SDSS gri cutouts from the skyplothelper query cache
  (``$MCF_M74_CACHE`` / ``$SPH_QUERY_CACHE`` / sibling checkout).

Usage (from repo root)::

    python docs/make_showcase_figures.py
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

import multicolorfits as mcf

ROOT = Path(__file__).resolve().parents[1]
OUT = Path(__file__).resolve().parent / '_static' / 'showcase'

LIGHT_FACE = 'white'
DARK_FACE = '#14161b'
DARK_TEXT = '#d8dce4'

POB_COLORS = ['#BE599E', '#DEA215', '#77C0F9']   # purple / orange / blue
POB_LABELS = ['IR', 'R', 'B']


def _save_pair(render, stem):
    """render(facecolor, textcolor) -> Figure; saved as _light + _dark."""
    for suffix, face, text in (('light', LIGHT_FACE, 'black'),
                               ('dark', DARK_FACE, DARK_TEXT)):
        fig = render(face, text)
        path = OUT / f'{stem}_{suffix}.png'
        fig.savefig(path, dpi=110, facecolor=face, bbox_inches='tight')
        plt.close(fig)
        print(f'  {path.relative_to(ROOT)}')


# --------------------------------------------------------------------------- NGC 602

def _find_ngc602_dir():
    env = os.environ.get('MCF_NGC602_DIR')
    for d in ([Path(env).expanduser()] if env else []) + [ROOT / 'hidden' / 'testing' / 'ngc602']:
        if d.is_dir() and (d / 'ngc602_ir.fits').is_file():
            return d
    return None


def _ngc602_session(datadir, background='black'):
    s = mcf.McfSession(n_panels=3)
    s.load_files(
        [str(datadir / 'ngc602_ir.fits'),
         str(datadir / 'ngc602_optical_R.fits'),
         str(datadir / 'ngc602_optical_B.fits')],
        colors=POB_COLORS, labels=POB_LABELS,
    )
    for p in s.panels:
        p.stretch = 'linear'
    s.compose.combine_mode = 'lab'
    s.compose.combine_blend = 'screen'
    s.compose.gamma = 2.2
    s.compose.combine_background = background
    return s


def ngc602_figures(datadir):
    print('[ngc602]')
    s = _ngc602_session(datadir)
    # No figure title — it crowds the top WCS tick labels; the legend and
    # combo swatch carry the annotation story.
    s.compose.show_legend = True
    s.compose.show_combo_swatch = True

    # Hero — Lab / screen with legend + combo swatch, real WCS axes.
    rgb_lab = s.render_combined()

    def hero(face, _text):
        s.compose.facecolor = face
        return mcf.make_combined_figure(s, combined=rgb_lab,
                                        figsize=(7.5, 7.5), facecolor=face)

    _save_pair(hero, 'ngc602_hero')

    # RGB vs Lab on identical colorized layers (downsampled; no WCS chrome).
    small = [mcf.downsample_for_preview(p.data, max_size=1200) for p in s.panels]
    colorized = [
        mcf.colorize_image(mcf.to_grey_rgb(d, rescalefn='linear'),
                           c, colorintype='hex', gammacorr_color=2.2)
        for d, c in zip(small, POB_COLORS)
    ]
    img_rgb = mcf.combine_multicolor(colorized, gamma=2.2)
    img_lab = mcf.combine_multicolor_colorspace(
        colorized, colorspace='lab', blend='screen', gamma=2.2)

    def rgb_vs_lab(face, text):
        fig, axes = plt.subplots(1, 2, figsize=(11, 5.6), facecolor=face)
        for ax, img, label in zip(
                axes, (img_rgb, img_lab),
                ('rgb (classic channel sum)', 'lab / screen')):
            ax.imshow(img, origin='lower', interpolation='nearest')
            ax.set_title(label, size=12, color=text)
            ax.set_xticks([]); ax.set_yticks([])
        fig.suptitle('NGC 602 — same layers, same colors, different combine',
                     size=13, color=text)
        fig.tight_layout(rect=[0, 0, 1, 0.94])
        return fig

    _save_pair(rgb_vs_lab, 'ngc602_rgb_vs_lab')

    # Component mosaic — hero + IR/R/B strip on the shared WCS.
    s.compose.show_legend = False
    s.compose.show_combo_swatch = False

    def mosaic(face, _text):
        s.compose.facecolor = face
        fig, _axes = mcf.make_component_mosaic(
            s, combined=rgb_lab, components='top', max_per_line=3,
            ticks='plain', facecolor=face,
        )
        return fig

    _save_pair(mosaic, 'ngc602_mosaic')


# --------------------------------------------------------------------------- M74

def _find_m74_cache():
    need = ('m74_sdss_g.fits', 'm74_sdss_r.fits', 'm74_sdss_i.fits')
    dirs = []
    env = os.environ.get('MCF_M74_CACHE') or os.environ.get('SPH_QUERY_CACHE')
    if env:
        dirs.append(Path(env).expanduser())
    dirs.append(ROOT.parent / 'skyplothelper' / 'examples' / 'data' / 'query_cache')
    for d in dirs:
        if d.is_dir() and all((d / n).is_file() for n in need):
            return d
    return None


def _m74_session(cache):
    from astropy.io import fits

    s = mcf.McfSession(n_panels=3)
    bands = ('m74_sdss_i.fits', 'm74_sdss_r.fits', 'm74_sdss_g.fits')
    # gri -> straight R/G/B channels: near-true-color, blue arms / warm core.
    colors = ['#FF0000', '#00FF00', '#0000FF']
    labels = ['i', 'r', 'g']
    for idx, (band, color, label) in enumerate(zip(bands, colors, labels)):
        data = np.nan_to_num(np.asarray(fits.getdata(cache / band), float))
        s.panels[idx].set_data(data, fits.getheader(cache / band))
        s.panels[idx].color = color
        s.panels[idx].label = label
        s.panels[idx].stretch = 'asinh'
        s.panels[idx].set_percentiles(30.0, 99.5)
    s.compose.combine_mode = 'rgb'
    s.compose.combine_background = 'black'
    return s


def m74_figures(cache):
    print('[m74]')
    s = _m74_session(cache)
    rgb = s.render_combined()
    # Smoothed-threshold matte: blur the intensity map (alpha_smooth) so the
    # matte follows whole arms, then treat alpha_hi as the full-opacity
    # threshold with a soft fade to zero below it.  Arms stay ~opaque even on
    # white; outer sky is exactly clear (no hazy box when overlaid).
    MATTE = dict(alpha_smooth=3, alpha_lo=60, alpha_hi=90, alpha_gamma=0.4)
    rgba = mcf.make_transparent_cutout(rgb, crop='auto', pad=0.06, **MATTE)

    # The real deliverable: a transparent stamp, committed as-is.  On the docs
    # site the page background shows through, in either theme.
    stamp = mcf.make_transparent_cutout(rgb, crop='auto', pad=0.06, size=512,
                                        **MATTE)
    path = OUT / 'm74_stamp.png'
    mcf.save_transparent_cutout(stamp, str(path))
    print(f'  {path.relative_to(ROOT)}')

    # Same cutout over checkerboard / dark / light backgrounds.
    def backgrounds(face, text):
        fig = mcf.preview_cutout_on_backgrounds(rgba)
        fig.patch.set_facecolor(face)
        for ax in fig.axes:
            ax.title.set_color(text)
        for t in fig.texts:   # suptitle
            t.set_color(text)
        return fig

    _save_pair(backgrounds, 'm74_cutout_backgrounds')

    # alpha_gamma: linear fade vs bold body, composited over a light deck.
    ny, nx = rgba.shape[:2]
    light_bg = np.full((ny, nx, 3), 0.94)

    def over(bg, a_gamma):
        kw = dict(MATTE, alpha_gamma=a_gamma)
        r = mcf.make_transparent_cutout(rgb, crop='auto', pad=0.06, **kw)
        al = r[..., 3:4]
        return np.clip(r[..., :3] * al + bg * (1 - al), 0, 1)

    def gamma_fig(face, text):
        cases = ((0.4, 'bold body'), (1.0, 'linear fade'), (1.4, 'tight matte'))
        fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.8), facecolor=face)
        for ax, (g, label) in zip(axes, cases):
            ax.imshow(over(light_bg, g), origin='lower', interpolation='nearest')
            ax.set_title(f'alpha_gamma = {g:g}  ({label})', size=11, color=text)
            ax.set_xticks([]); ax.set_yticks([])
        fig.suptitle('M74 over a light slide background — fade shape',
                     size=12.5, color=text)
        fig.tight_layout(rect=[0, 0, 1, 0.93])
        return fig

    _save_pair(gamma_fig, 'm74_alpha_gamma')

    # Deblend round trip: flatten the RGBA stamp onto a solid color, then
    # recover alpha with the closed-form solve and show it on a checkerboard.
    flat_bg = '#dbe4f0'   # a light "slide" blue-grey
    flat = mcf.flatten_rgba(stamp, background=flat_bg)
    flat8 = np.round(flat * 255) / 255.0            # honest 8-bit copy
    rec = mcf.deblend_background(flat8, background=flat_bg, tol=1.5 / 255)
    again = mcf.flatten_rgba(rec, background=flat_bg)
    err255 = float(np.abs(again - flat8).max()) * 255

    def deblend_fig(face, text):
        fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.8), facecolor=face)
        for ax, img, label in zip(
                axes,
                (flat8, rec[..., 3], again),
                (f'flattened onto {flat_bg} (8-bit)',
                 'recovered alpha (closed-form min-\u03b1 solve)',
                 f're-flattened — max err {err255:.2f}/255')):
            if img.ndim == 2:
                ax.imshow(img, origin='lower', cmap='gray', vmin=0, vmax=1,
                          interpolation='nearest')
            else:
                ax.imshow(img, origin='lower', interpolation='nearest')
            ax.set_title(label, size=11, color=text)
            ax.set_xticks([]); ax.set_yticks([])
        fig.suptitle('deblend_background — un-premultiply from a known solid color',
                     size=12.5, color=text)
        fig.tight_layout(rect=[0, 0, 1, 0.93])
        return fig

    _save_pair(deblend_fig, 'm74_deblend')


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    ngc = _find_ngc602_dir()
    if ngc:
        ngc602_figures(ngc)
    else:
        print('[ngc602] SKIP — data not found (set MCF_NGC602_DIR)')

    m74 = _find_m74_cache()
    if m74:
        m74_figures(m74)
    else:
        print('[m74] SKIP — SDSS cache not found (set MCF_M74_CACHE)')

    print('Done.')


if __name__ == '__main__':
    main()
