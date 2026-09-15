#!/usr/bin/env python3
"""
Build light/dark documentation figure pairs under docs/_static/.

Prefers pre-rendered ``*_light.png`` / ``*_dark.png`` pairs from
``examples/output/`` (e.g. from ``compare_colorspaces.py``, which draws
twice with the docs theme canvas colors). Falls back to remapping near-white
chrome for assets that only exist as a single light PNG (GUI screenshots,
legacy gallery stills without a hand-tuned dark).

Also builds theme pairs for the classic Gallery page
(``docs/_static/examples/*_{light,dark}.png``), and dual-renders the
guide mosaic layout demos under ``docs/_static/mosaics/`` (synthetic
cont/Ha/OIII/SII session — never chrome-remapped).

Usage (from repo root or docs/)::

    python docs/make_docs_figures.py
    python docs/make_docs_figures.py --mosaics-only
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from PIL import Image

try:
    from scipy import ndimage
except ImportError:  # pragma: no cover
    ndimage = None

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / 'examples' / 'output'
STATIC = Path(__file__).resolve().parent / '_static'
IMAGES = ROOT / 'images'

# docs/_static/<subdir>/<stem>_light.png + _dark.png
# Names are the light source stem (with or without _light); see resolve_pair().
ASSETS = {
    'compositing': [
        'kepler_colorspaces.png',
        'kepler_lab_blends.png',
        'kepler_hsv_blends.png',
        'kepler_complementary.png',
        'ngc602_colorspaces.png',
        'ngc602_lab_blends.png',
        'paintmix_demo.png',
        'swatch_legend_demo.png',
    ],
    # mosaics: dual-rendered by render_mosaic_layout_figures() (not remapped)
    'gui': [
        'gui_colorspace_modes.png',
        'gui_new_compositing.png',
        'qt_gui_colorspace.png',
    ],
}

# Classic gallery page (docs/examples/gallery.md). Sources live under
# docs/_static/examples/ (and matching hand-tuned darks under images/).
GALLERY = [
    'n602_POB.jpg',
    'WLM_testplot.jpg',
    'kepler_POT.png',
    'm51_RGBL.png',
    'm51_RGB_inverse.png',
    'm106_pureRGB.png',
    'mcf_gui_Crab.png',
]

DARK_BG = np.array([0x14, 0x16, 0x1B], dtype=np.float64)   # docs charcoal
LIGHT_FG = np.array([0xD8, 0xDC, 0xE4], dtype=np.float64)  # docs body text


def chrome_to_dark(rgb: np.ndarray) -> np.ndarray:
    """Remap light figure chrome to a dark-theme canvas; keep science color.

    Fallback for screenshots / legacy stills that were never dual-rendered.
    Prefer re-rendering with the docs face/text colors when a script can
    produce both themes (see examples/compare_colorspaces.py,
    docs/make_showcase_figures.py).
    """
    arr = rgb.astype(np.float64)
    r, g, b = arr[..., 0], arr[..., 1], arr[..., 2]
    mx = np.maximum(np.maximum(r, g), b)
    mn = np.minimum(np.minimum(r, g), b)
    sat = np.divide(mx - mn, mx, out=np.zeros_like(mx, dtype=np.float64), where=mx > 1.0)
    lum = mx / 255.0

    paper_cand = (sat < 0.12) & (lum > 0.94)
    if ndimage is not None:
        labeled, _n = ndimage.label(paper_cand)
        border = np.unique(np.concatenate([
            labeled[0, :], labeled[-1, :], labeled[:, 0], labeled[:, -1],
        ]))
        border = border[border != 0]
        paper = np.isin(labeled, border) if border.size else np.zeros_like(paper_cand)
        chrome = ndimage.binary_dilation(paper, iterations=4)
    else:
        paper = paper_cand.copy()
        paper[1:-1, 1:-1] = False
        chrome = paper.copy()
        for _ in range(4):
            pad = np.pad(chrome, 1, mode='edge')
            chrome = (
                pad[1:-1, 1:-1] | pad[:-2, 1:-1] | pad[2:, 1:-1]
                | pad[1:-1, :-2] | pad[1:-1, 2:]
            ) & ((sat < 0.18) | paper_cand)

    out = arr.copy()
    greyscale = sat < 0.14
    H = arr.shape[0]
    title_band = np.zeros_like(paper_cand)
    title_band[: max(1, int(0.22 * H)), :] = True
    if ndimage is not None:
        near_chrome = ndimage.binary_dilation(chrome, iterations=3)
    else:
        near_chrome = chrome
    interior = paper_cand & ~paper & (title_band | near_chrome)
    if ndimage is not None and np.any(interior):
        labeled_i, n_i = ndimage.label(interior)
        keep = np.zeros(n_i + 1, dtype=bool)
        for lab in range(1, n_i + 1):
            if int((labeled_i == lab).sum()) <= 400:
                keep[lab] = True
        letter_counters = keep[labeled_i]
    else:
        letter_counters = interior
    remap = (chrome | letter_counters) & greyscale
    mapped = LIGHT_FG + (DARK_BG - LIGHT_FG) * lum[..., None]
    out[remap] = mapped[remap]
    return np.clip(out, 0, 255).astype(np.uint8)


def resolve_pair(src: Path) -> tuple[Path, Path | None]:
    """Return (light_src, dark_src_or_None) preferring dual-rendered siblings."""
    stem = src.stem
    parent = src.parent
    # Already a *_light.png path, or stem.png with stem_{light,dark}.png beside it.
    if stem.endswith('_light'):
        base = stem[: -len('_light')]
        light = src if src.is_file() else parent / f'{base}_light.png'
        dark = parent / f'{base}_dark.png'
        return light, (dark if dark.is_file() else None)
    light_sib = parent / f'{stem}_light.png'
    dark_sib = parent / f'{stem}_dark.png'
    if light_sib.is_file() and dark_sib.is_file():
        return light_sib, dark_sib
    if src.is_file() and dark_sib.is_file():
        return src, dark_sib
    if light_sib.is_file():
        return light_sib, (dark_sib if dark_sib.is_file() else None)
    return src, (dark_sib if dark_sib.is_file() else None)


def write_pair(src: Path, dest_stem: Path, dark_src: Path | None = None) -> None:
    dest_stem.parent.mkdir(parents=True, exist_ok=True)
    light = dest_stem.with_name(dest_stem.name + '_light.png')
    dark = dest_stem.with_name(dest_stem.name + '_dark.png')

    light_src, auto_dark = resolve_pair(src)
    if dark_src is None:
        dark_src = auto_dark
    if not light_src.is_file():
        print(f'  SKIP missing {src}')
        return

    Image.open(light_src).convert('RGB').save(light, optimize=True)
    if dark_src is not None and dark_src.is_file():
        Image.open(dark_src).convert('RGB').save(dark, optimize=True)
        how = 'pre-rendered'
    else:
        Image.fromarray(
            chrome_to_dark(np.asarray(Image.open(light_src).convert('RGB'))),
            mode='RGB',
        ).save(dark, optimize=True)
        how = 'remapped'
    print(f'  {light.relative_to(STATIC.parent)}')
    print(f'  {dark.relative_to(STATIC.parent)}  [{how}]')


def write_gallery_pairs() -> None:
    """Theme-aware pairs for docs/examples/gallery.md."""
    src_dir = STATIC / 'examples'
    print(f'\n[examples gallery]')
    for name in GALLERY:
        src = src_dir / name
        if not src.is_file():
            alt = IMAGES / name
            if alt.is_file():
                src = alt
            else:
                print(f'  SKIP missing {name}')
                continue
        stem = Path(name).stem
        # Prefer a hand-rendered dark canvas from images/ when present
        # (README science figures + GUI screenshot are authored twice).
        dark_src = IMAGES / f'{stem}_dark.png'
        if not dark_src.is_file():
            dark_src = None
        write_pair(src, src_dir / stem, dark_src=dark_src)


def render_mosaic_layout_figures() -> None:
    """Dual-render guide mosaic layout demos (never chrome-remap).

    Synthetic 4-band session matching the cont / Ha / OIII / SII layout
    figures in ``docs/guide/figures_and_mosaics.md``.  Writes
    ``docs/_static/mosaics/*_{light,dark}.png``.
    """
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    import multicolorfits as mcf

    light_face, dark_face = 'white', '#14161b'
    dark_text = '#d8dce4'
    out = STATIC / 'mosaics'
    out.mkdir(parents=True, exist_ok=True)

    def _blob(nx, ny, seed, scale=8.0):
        rng = np.random.default_rng(seed)
        y, x = np.mgrid[0:ny, 0:nx]
        blob = np.exp(-(((x - nx / 2.0) / (nx / scale)) ** 2
                        + ((y - ny / 2.0) / (ny / scale)) ** 2))
        return blob + rng.normal(0, 0.04, size=(ny, nx))

    print('\n[mosaics — dual-render]')
    s = mcf.McfSession(n_panels=4)
    bands = [
        ('cont', '#DEA215', 0),
        ('Ha', '#E24A33', 1),
        ('OIII', '#2CA02C', 2),
        ('SII', '#1F77B4', 3),
    ]
    nx = ny = 96
    hdr = None
    try:
        import astropy.io.fits as pyfits
        hdr = pyfits.Header()
        hdr['NAXIS'] = 2
        hdr['NAXIS1'] = nx
        hdr['NAXIS2'] = ny
        hdr['CTYPE1'] = 'RA---TAN'
        hdr['CTYPE2'] = 'DEC--TAN'
        hdr['CRPIX1'] = nx / 2.0
        hdr['CRPIX2'] = ny / 2.0
        hdr['CRVAL1'] = 150.0
        hdr['CRVAL2'] = -30.0
        hdr['CDELT1'] = -0.0002777778
        hdr['CDELT2'] = 0.0002777778
        hdr['RADESYS'] = 'FK5'
        hdr['EQUINOX'] = 2000.0
    except Exception:
        hdr = None
    for i, (label, color, seed) in enumerate(bands):
        s.panels[i].set_data(_blob(nx, ny, seed), hdr)
        s.panels[i].color = color
        s.panels[i].label = label
        s.panels[i].stretch = 'linear'
    s.compose.combine_mode = 'lab'
    s.compose.combine_blend = 'screen'
    s.compose.gamma = 2.2
    s.compose.combine_background = 'black'
    s.compose.show_legend = False
    s.compose.show_combo_swatch = False
    rgb = s.render_combined()

    layouts = (
        ('mosaic_top_max3', dict(components='top', max_per_line=3)),
        ('mosaic_top_max2', dict(components='top', max_per_line=2)),
        ('mosaic_left_max3', dict(components='left', max_per_line=3)),
    )
    for stem, kw in layouts:
        for suffix, face, tick in (
            ('light', light_face, '0.9'),
            ('dark', dark_face, dark_text),
        ):
            s.compose.facecolor = face
            s.compose.tickcolor = tick
            fig, _axes = mcf.make_component_mosaic(
                s, combined=rgb, ticks='minimal', facecolor=face,
                label_loc='upper left', **kw,
            )
            path = out / f'{stem}_{suffix}.png'
            fig.savefig(path, dpi=120, facecolor=face, bbox_inches='tight')
            plt.close(fig)
            print(f'  {path.relative_to(ROOT)}')


def refresh_autogen_darks() -> None:
    """Re-derive remapped ``*_dark.png`` from ``*_light.png`` (legacy only).

    Skips gallery stems with a hand-tuned ``images/<stem>_dark.png``, skips
    dual-rendered dirs (``showcase/``, ``mosaics/``), and skips compositing
    stems that already have a dual-rendered dark in ``examples/output/``
    (so this does not overwrite true theme renders).

    Remapped darks look speckled on WCS spines / thick white frames — prefer
    dual-rendering whenever a generator can produce both themes.
    """
    hand_tuned = {p.name.replace('_dark.png', '') for p in IMAGES.glob('*_dark.png')}
    skip_dirs = {'showcase', 'mosaics'}  # dual-rendered; never chrome-remap
    print('\n[refresh auto darks from *_light.png]')
    for light in sorted(STATIC.rglob('*_light.png')):
        if light.parent.name in skip_dirs:
            print(f'  SKIP dual-render dir {light.relative_to(STATIC)}')
            continue
        stem = light.name.replace('_light.png', '')
        if stem in hand_tuned and light.parent.name == 'examples':
            print(f'  SKIP hand-tuned {light.relative_to(STATIC)}')
            continue
        # Prefer dual-rendered dark from examples/output when present.
        out_dark = SRC / f'{stem}_dark.png'
        if not out_dark.is_file():
            matches = list(SRC.rglob(f'{stem}_dark.png'))
            out_dark = matches[0] if matches else out_dark
        dark = light.with_name(stem + '_dark.png')
        if out_dark.is_file():
            Image.open(out_dark).convert('RGB').save(dark, optimize=True)
            print(f'  {dark.relative_to(STATIC.parent)}  [pre-rendered]')
            continue
        arr = np.asarray(Image.open(light).convert('RGB'))
        Image.fromarray(chrome_to_dark(arr), mode='RGB').save(dark, optimize=True)
        print(f'  {dark.relative_to(STATIC.parent)}  [remapped]')


def main() -> None:
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--refresh-darks', action='store_true',
        help='Refresh *_dark.png from pre-rendered output or remap from *_light.png')
    parser.add_argument(
        '--mosaics-only', action='store_true',
        help='Only dual-render docs/_static/mosaics/ layout figures')
    args = parser.parse_args()
    if args.refresh_darks:
        refresh_autogen_darks()
        print('\nDone.')
        return
    if args.mosaics_only:
        render_mosaic_layout_figures()
        print('\nDone.')
        return
    if not SRC.is_dir():
        raise SystemExit(f'Missing examples output dir: {SRC}')
    print(f'Source: {SRC}')
    print(f'Target: {STATIC}')
    for subdir, names in ASSETS.items():
        print(f'\n[{subdir}]')
        for name in names:
            src = SRC / name
            stem = Path(name).stem
            write_pair(src, STATIC / subdir / stem)
    render_mosaic_layout_figures()
    write_gallery_pairs()
    print('\nDone.')


if __name__ == '__main__':
    main()
