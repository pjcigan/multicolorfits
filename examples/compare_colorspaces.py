"""
Visual comparison of classic RGB compositing vs Lab / HSV / HSL on real
tutorial data, plus the NGC 602 subtractive-paint (RYB) demo.

Datasets:
* Kepler's SNR (Chandra openFITS): IR (MIPS 24um), soft X-ray (0.3-1.4 keV),
  hard X-ray (4-6 keV) -- same three bands / linear stretches as the classic
  README POT example. Optical (WFPC2) is a tiny corner patch and is omitted.
* NGC 602 (Chandra openFITS): IR / optical R / optical B -- downsampled

Writes light/dark PNG pairs into examples/output/ (science layers computed
once; matplotlib chrome re-rendered with the docs theme canvas colors).
Copy into docs/_static/compositing/ via::

    python examples/compare_colorspaces.py
    python docs/make_docs_figures.py
"""

from __future__ import annotations

import os
import sys

import numpy as np
import astropy.io.fits as pyfits
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

# Larger titles + DejaVu so letter counters stay open after docs downscaling.
plt.rcParams.update({
    'font.family': 'DejaVu Sans',
    'axes.titlesize': 13,
    'figure.titlesize': 15,
    'axes.titlepad': 8,
    'text.antialiased': True,
})
SAVE_DPI = 160

# Match docs/_static/custom.css charcoal + body text (same as make_showcase_figures).
LIGHT_FACE = 'white'
LIGHT_TEXT = 'black'
DARK_FACE = '#14161b'
DARK_TEXT = '#d8dce4'

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import multicolorfits as mcf

def _data_dir(env_key, *candidates):
    env = os.environ.get(env_key)
    if env and os.path.isdir(env):
        return env
    for path in candidates:
        if path and os.path.isdir(path):
            return path
    return candidates[-1] if candidates else ''


ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
OUTDIR = os.path.join(os.path.dirname(__file__), 'output')
os.makedirs(OUTDIR, exist_ok=True)

_home = os.path.expanduser('~')
KEPLER_DIR = _data_dir(
    'MCF_KEPLER_DIR',
    os.path.join(_home, 'Desktop', 'misc', 'testcode', 'astro_3color', 'kepler'),
    os.path.join(_home, 'data', 'astro_3color', 'kepler'),
)
NGC602_DIR = _data_dir(
    'MCF_NGC602_DIR',
    os.path.join(ROOT, 'hidden', 'testing', 'ngc602'),
)


MODES = ['rgb', 'lab', 'hsv', 'hsl']
THEMES = (
    ('light', LIGHT_FACE, LIGHT_TEXT),
    ('dark', DARK_FACE, DARK_TEXT),
)


def composite(colorized, mode, blend='screen', gamma=2.2):
    if mode == 'rgb':
        return mcf.combine_multicolor(colorized, gamma=gamma)
    return mcf.combine_multicolor_colorspace(colorized, colorspace=mode,
                                             blend=blend, gamma=gamma)


def _apply_chrome(fig, axes, face, text, title):
    fig.patch.set_facecolor(face)
    for ax in np.atleast_1d(axes).ravel():
        ax.set_facecolor(face)
        ax.title.set_color(text)
        for spine in ax.spines.values():
            spine.set_edgecolor(text)
            spine.set_alpha(0.35)
    fig.suptitle(title, color=text)


def _save_theme_pair(build_fig, stem):
    """build_fig(face, text) -> fig; write stem_{light,dark}.png (+ stem.png = light)."""
    light_path = None
    for suffix, face, text in THEMES:
        fig = build_fig(face, text)
        path = os.path.join(OUTDIR, f'{stem}_{suffix}.png')
        fig.savefig(path, dpi=SAVE_DPI, facecolor=face)
        plt.close(fig)
        print('wrote', path)
        if suffix == 'light':
            light_path = path
    # Unsuffixed alias for scripts / make_docs_figures that still look for stem.png.
    if light_path is not None:
        alias = os.path.join(OUTDIR, f'{stem}.png')
        if os.path.lexists(alias):
            os.remove(alias)
        # Hard-link when possible so we do not double the disk footprint.
        try:
            os.link(light_path, alias)
        except OSError:
            import shutil
            shutil.copy2(light_path, alias)


def render_grid(colorized, title, stem, blend='screen'):
    # Composite once; only figure chrome differs per theme.
    panels = [(mode, composite(colorized, mode, blend=blend)) for mode in MODES]

    def build(face, text):
        fig, axes = plt.subplots(1, len(MODES), figsize=(4 * len(MODES), 4.4),
                                 facecolor=face)
        for ax, (mode, img) in zip(axes, panels):
            ax.imshow(img, origin='lower', interpolation='nearest')
            label = 'rgb (classic)' if mode == 'rgb' else '%s / %s' % (mode, blend)
            ax.set_title(label)
            ax.set_xticks([]); ax.set_yticks([])
        _apply_chrome(fig, axes, face, text, title)
        fig.tight_layout(rect=[0, 0, 1, 0.95])
        return fig

    _save_theme_pair(build, stem)


def render_blend_grid(colorized, mode, title, stem):
    blends = ['screen', 'sum', 'max', 'mean']
    panels = [(blend, composite(colorized, mode, blend=blend)) for blend in blends]

    def build(face, text):
        fig, axes = plt.subplots(1, len(blends), figsize=(4 * len(blends), 4.4),
                                 facecolor=face)
        for ax, (blend, img) in zip(axes, panels):
            ax.imshow(img, origin='lower', interpolation='nearest')
            ax.set_title('%s / %s' % (mode, blend))
            ax.set_xticks([]); ax.set_yticks([])
        _apply_chrome(fig, axes, face, text, title)
        fig.tight_layout(rect=[0, 0, 1, 0.95])
        return fig

    _save_theme_pair(build, stem)


def load_grey(path, stretch='linear', scaletype='abs', min_max=(None, None)):
    data = pyfits.getdata(path)
    return mcf.to_grey_rgb(data, rescalefn=stretch, scaletype=scaletype,
                           min_max=list(min_max), gamma=2.2)


# ---------------------------------------------------------------------------
# Kepler's SNR -- 3 bands (IR + soft/hard X-ray), classic linear stretches
# ---------------------------------------------------------------------------
print('--- Kepler SNR ---')
kep_ir = load_grey(os.path.join(KEPLER_DIR, 'kepler_ir.fits'), stretch='linear')
kep_xlo = load_grey(os.path.join(KEPLER_DIR, 'kepler_xray_le.fits'), stretch='linear')
kep_xhi = load_grey(os.path.join(KEPLER_DIR, 'kepler_xray_he.fits'), stretch='linear')

# Classic POT hues from the README / generate_plots recipe (HSV).
kep_colorized = [
    mcf.colorize_image(kep_ir, [282. / 360, 1.0, 0.3], colorintype='hsv',
                       gammacorr_color=2.2),
    mcf.colorize_image(kep_xlo, [21. / 360, 0.90, 0.5], colorintype='hsv',
                       gammacorr_color=2.2),
    mcf.colorize_image(kep_xhi, [180. / 360, 0.77, 0.4], colorintype='hsv',
                       gammacorr_color=2.2),
]

render_grid(
    kep_colorized,
    "Kepler SNR -- IR purple, soft X orange, hard X teal (linear stretch, blend=screen)",
    'kepler_colorspaces')
render_blend_grid(
    kep_colorized, 'lab',
    'Kepler SNR -- Lab compositing, lightness blend comparison',
    'kepler_lab_blends')
render_blend_grid(
    kep_colorized, 'hsv',
    'Kepler SNR -- HSV compositing, lightness blend comparison',
    'kepler_hsv_blends')

# Complementary stress test: two overlapping X-ray-bright layers (no optical).
kep_2lay = [
    mcf.colorize_image(kep_xlo, '#FF8800', colorintype='hex', gammacorr_color=2.2),
    mcf.colorize_image(kep_xhi, '#0088FF', colorintype='hex', gammacorr_color=2.2),
]
render_grid(
    kep_2lay,
    'Kepler SNR -- complementary orange + blue X-ray (overlaps neutralize in additive light)',
    'kepler_complementary')

# ---------------------------------------------------------------------------
# NGC 602 -- tutorial POB palette (purple / orange / blue), linear stretch
# ---------------------------------------------------------------------------
print('--- NGC 602 ---')
n602_ir = mcf.downsample_for_preview(
    pyfits.getdata(os.path.join(NGC602_DIR, 'ngc602_ir.fits')), max_size=1200)
n602_R = mcf.downsample_for_preview(
    pyfits.getdata(os.path.join(NGC602_DIR, 'ngc602_optical_R.fits')), max_size=1200)
n602_B = mcf.downsample_for_preview(
    pyfits.getdata(os.path.join(NGC602_DIR, 'ngc602_optical_B.fits')), max_size=1200)

n602_greys = [
    mcf.to_grey_rgb(n602_ir, rescalefn='linear', gamma=2.2),
    mcf.to_grey_rgb(n602_R, rescalefn='linear', gamma=2.2),
    mcf.to_grey_rgb(n602_B, rescalefn='linear', gamma=2.2),
]

n602_colorized = [
    mcf.colorize_image(n602_greys[0], '#BE599E', colorintype='hex', gammacorr_color=2.2),
    mcf.colorize_image(n602_greys[1], '#DEA215', colorintype='hex', gammacorr_color=2.2),
    mcf.colorize_image(n602_greys[2], '#77C0F9', colorintype='hex', gammacorr_color=2.2),
]

render_grid(
    n602_colorized,
    'NGC 602 -- tutorial POB palette: IR purple, R orange, B blue (blend=screen)',
    'ngc602_colorspaces')
render_blend_grid(
    n602_colorized, 'lab',
    'NGC 602 -- Lab compositing, lightness blend comparison',
    'ngc602_lab_blends')

# ---------------------------------------------------------------------------
# Subtractive paint demo -- same linear greys as above, R/Y/B paint colors
# ---------------------------------------------------------------------------
print('--- paintmix demo ---')
GAMMA = 2.2
paint_layers = [
    mcf.colorize_image(n602_greys[1], '#FF3030', colorintype='hex', gammacorr_color=GAMMA),  # R
    mcf.colorize_image(n602_greys[0], '#FFD000', colorintype='hex', gammacorr_color=GAMMA),  # IR
    mcf.colorize_image(n602_greys[2], '#2E6BFF', colorintype='hex', gammacorr_color=GAMMA),  # B
]

additive = mcf.combine_multicolor(paint_layers, gamma=GAMMA)
old_inverse = mcf.combine_multicolor(paint_layers, gamma=GAMMA, inverse=True)
white_bg = mcf.composite_over_background(additive, background='white')
slide_bg = mcf.combine_multicolor_alpha(
    paint_layers, background='#12243A', mode='rgb', gamma=GAMMA)
ryb_paint = mcf.combine_multicolor_alpha(
    paint_layers, background='white', mode='ryb', gamma=GAMMA)
transparent = mcf.combine_multicolor_alpha(
    paint_layers, background=None, mode='rgb', gamma=GAMMA)

paint_panels = [
    ('Classic additive (black bg)', additive, None),
    ('Classic inverse (naive 1-x)', old_inverse, None),
    ('Alpha over white', white_bg, None),
    ('Alpha over slide color', slide_bg, None),
    ('RYB paint over white', ryb_paint, None),
    ('Transparent RGBA on checker', transparent, 'checker'),
]


def _build_paintmix(face, text):
    fig, axes = plt.subplots(2, 3, figsize=(15, 10), facecolor=face)
    for ax, (title, img, special) in zip(axes.ravel(), paint_panels):
        if special == 'checker':
            ny, nx = img.shape[:2]
            c = (((np.arange(ny)[:, None] // 32) + (np.arange(nx)[None, :] // 32)) % 2)
            # Slightly darker checker on the dark canvas so the plate reads clearly.
            lo, hi = (0.22, 0.38) if face != LIGHT_FACE else (0.45, 0.75)
            checker = np.where(c[..., None] == 0, hi, lo) * np.ones((1, 1, 3))
            ax.imshow(checker, origin='lower', interpolation='nearest')
        ax.imshow(np.clip(np.asarray(img)[..., :3], 0, 1),
                  origin='lower', interpolation='nearest')
        ax.set_title(title)
        ax.set_xticks([]); ax.set_yticks([])
    _apply_chrome(
        fig, axes, face, text,
        'multicolorfits compositing modes -- NGC 602 (R=red, IR=yellow, B=blue; linear stretch)')
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    return fig


_save_theme_pair(_build_paintmix, 'paintmix_demo')

print('done.')
