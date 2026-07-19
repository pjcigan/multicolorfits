#!/usr/bin/env python3
"""
Build light/dark documentation figure pairs under docs/_static/.

Prefers pre-rendered ``*_light.png`` / ``*_dark.png`` pairs from
``examples/output/`` (e.g. from ``compare_colorspaces.py``, which draws
twice with the docs theme canvas colors). Falls back to remapping near-white
chrome for assets that only exist as a single light PNG (GUI screenshots,
legacy gallery stills without a hand-tuned dark).

Also builds theme pairs for the classic Gallery page
(``docs/_static/examples/*_{light,dark}.png``).

Usage (from repo root or docs/)::

    python docs/make_docs_figures.py
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
    'mosaics': [
        'component_mosaic/mosaic_top_max3.png',
        'component_mosaic/mosaic_top_max2.png',
        'component_mosaic/mosaic_left_max3.png',
    ],
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


def refresh_autogen_darks() -> None:
    """Re-derive remapped ``*_dark.png`` from ``*_light.png`` (legacy only).

    Skips gallery stems with a hand-tuned ``images/<stem>_dark.png``, and
    skips compositing stems that already have a dual-rendered dark in
    ``examples/output/`` (so this does not overwrite true theme renders).
    """
    hand_tuned = {p.name.replace('_dark.png', '') for p in IMAGES.glob('*_dark.png')}
    print('\n[refresh auto darks from *_light.png]')
    for light in sorted(STATIC.rglob('*_light.png')):
        stem = light.name.replace('_light.png', '')
        if stem in hand_tuned and light.parent.name == 'examples':
            print(f'  SKIP hand-tuned {light.relative_to(STATIC)}')
            continue
        # Prefer dual-rendered dark from examples/output when present.
        out_dark = SRC / f'{stem}_dark.png'
        if not out_dark.is_file():
            # mosaics live under a subdir in examples/output
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
    args = parser.parse_args()
    if args.refresh_darks:
        refresh_autogen_darks()
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
    write_gallery_pairs()
    print('\nDone.')


if __name__ == '__main__':
    main()
