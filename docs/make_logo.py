#!/usr/bin/env python3
"""Build MultiColorFits docs branding marks — committed, NOT run at build time.

Two navbar-mark options are always generated so they can be swapped later:

* ``logo_navbar-color.png`` — filled 3-circle ``combo_swatch`` (ryb, γ=2.2)
  pre-composited over black so lobe colours stay rich; outside circles is
  transparent.  Wordmark hues (C / l / r).  **Active default.**
* ``logo_navbar-mono-{light,dark}.png`` — blank Venn-outline rings in theme
  ink.  Swap the navbar ``<img>`` to these for the quieter sph-like look.

Favicon uses the same ``combo_swatch`` pipeline and circle arrangement as the
color navbar mark, but with bright RYB primaries so it still reads at 16 px:

* ``logo_favicon.png`` (64×64) — also written as ``logo_favicon-32.png``.

Hero / social:

* ``logo_hero_crab.png`` — Crab Nebula cutout (press-release look) for the
  docs landing page.
* ``logo_social.png`` — 1200×630 OG / Twitter card with the Crab + wordmark.

Colour / layout conventions (must stay matched between banner & favicon)::

    after flipud:  idx0 = TOP,  idx1 = BOTTOM-RIGHT,  idx2 = BOTTOM-LEFT
    banner  :  #46618A / #C29B3C / #8A4540   (blue / gold / brick)
    favicon :  #2E6BFF / #FFD000 / #FF3030   (same hue families, bright)

Regenerate::

    python docs/make_logo.py
"""
from __future__ import annotations

import base64
import io
import os
import sys

import numpy as np
from PIL import Image, ImageDraw, ImageFont

_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

from multicolorfits.core.swatch import combo_swatch  # noqa: E402

OUT = os.path.join(_REPO, 'docs', '_static', 'logo')
CRAB_SRC = os.path.join(_REPO, 'hidden', 'crab_icons', 'crab_red_orange_cutout.png')
# Prefer cutout; fall back to full / icon if cutout missing.
for _alt in ('crab_red_orange_full.png', 'crab_red_orange_icon.png'):
    if not os.path.isfile(CRAB_SRC):
        CRAB_SRC = os.path.join(_REPO, 'hidden', 'crab_icons', _alt)

BANNER = ['#46618A', '#C29B3C', '#8A4540']   # top, bottom-right, bottom-left
FAVICON = ['#2E6BFF', '#FFD000', '#FF3030']  # same positions, bright
URAN = ['#46618A', '#B97C52', '#C29B3C', '#5E8C7E', '#8A4540']  # C o l o r

FONT_PATH = os.path.join(
    os.path.dirname(np.__file__), '..', 'matplotlib', 'mpl-data', 'fonts', 'ttf',
    'DejaVuSans-Bold.ttf',
)
_CANDIDATE_FONTS = [FONT_PATH]


def _font(size):
    try:
        from matplotlib import font_manager
        path = font_manager.findfont('DejaVu Sans:style=Bold', fallback_to_default=True)
        return ImageFont.truetype(path, size)
    except Exception:
        for p in _CANDIDATE_FONTS:
            if os.path.isfile(p):
                return ImageFont.truetype(p, size)
        return ImageFont.load_default()


def hexrgb(h):
    h = h.lstrip('#')
    return tuple(int(h[i:i + 2], 16) for i in (0, 2, 4))


def combo_mark(colors, size, *, supersample=4):
    """combo_swatch over black (rich) + transparent outside; one-top / two-bottom."""
    d = combo_swatch(colors, mode='ryb', gamma=2.2, background='black',
                     size=size * supersample)
    arr = np.flipud(np.asarray(d['rgba']))
    arr = (np.clip(arr, 0, 1) * 255).astype(np.uint8)
    return Image.fromarray(arr, 'RGBA').resize((size, size), Image.LANCZOS)


def mono_mark(size, ink, *, supersample=4):
    """Blank Venn-outline rings (one top, two bottom) in a single ink colour."""
    n = size * supersample
    img = Image.new('RGBA', (n, n), (0, 0, 0, 0))
    dr = ImageDraw.Draw(img)
    r = n * 0.30
    cx, cy = n / 2, n / 2
    off = r * 0.62
    centers = [(cx, cy - off),
               (cx - off * 0.87, cy + off * 0.5),
               (cx + off * 0.87, cy + off * 0.5)]
    w = max(2, int(n * 0.045))
    for (x, y) in centers:
        dr.ellipse([x - r, y - r, x + r, y + r], outline=ink + (255,), width=w)
    return img.resize((size, size), Image.LANCZOS)


def save_png(img, name):
    path = os.path.join(OUT, name)
    img.save(path, format='PNG')
    print('wrote', path)
    return path


def wrap_svg(png_img, name, view=32):
    """Tiny SVG wrapper around a PNG (keeps current html_favicon=*.svg option)."""
    buf = io.BytesIO()
    png_img.save(buf, format='PNG')
    b64 = base64.b64encode(buf.getvalue()).decode('ascii')
    svg = (
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        f'<!-- MultiColorFits mark. Regenerate: python docs/make_logo.py -->\n'
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'xmlns:xlink="http://www.w3.org/1999/xlink" '
        f'viewBox="0 0 {view} {view}" role="img" aria-label="MultiColorFits">\n'
        f'  <image width="{view}" height="{view}" '
        f'xlink:href="data:image/png;base64,{b64}"/>\n'
        f'</svg>\n'
    )
    path = os.path.join(OUT, name)
    with open(path, 'w', encoding='utf-8') as f:
        f.write(svg)
    print('wrote', path)
    return path


def make_navbar_and_favicon():
    os.makedirs(OUT, exist_ok=True)
    # Active color mark (works on light & dark — transparent bg, rich lobes).
    color = combo_mark(BANNER, 128)
    save_png(color, 'logo_navbar-color.png')
    # Mono options for future swap.
    save_png(mono_mark(128, (30, 30, 30)), 'logo_navbar-mono-light.png')
    save_png(mono_mark(128, (235, 235, 235)), 'logo_navbar-mono-dark.png')
    # Favicon — bright RYB, same arrangement.
    fav64 = combo_mark(FAVICON, 64)
    fav32 = combo_mark(FAVICON, 32)
    save_png(fav64, 'logo_favicon.png')
    save_png(fav32, 'logo_favicon-32.png')
    wrap_svg(fav32, 'logo_favicon.svg', view=32)
    # Also drop a copy next to the old path so existing refs keep working
    # until conf.py is updated; the SVG is the favicon.
    legacy = os.path.join(_REPO, 'docs', '_static', 'mcf-mark.svg')
    with open(os.path.join(OUT, 'logo_favicon.svg'), encoding='utf-8') as src:
        data = src.read()
    with open(legacy, 'w', encoding='utf-8') as dst:
        dst.write(data.replace('aria-label="MultiColorFits"',
                               'aria-label="MultiColorFits favicon"'))
    print('wrote', legacy, '(synced from logo_favicon.svg)')
    # Browser GUI uses the same bright-RYB favicon.
    web_static = os.path.join(_REPO, 'multicolorfits', 'gui_web', 'static')
    if os.path.isdir(web_static):
        import shutil
        for name in ('logo_favicon.png', 'logo_favicon-32.png', 'logo_favicon.svg'):
            src = os.path.join(OUT, name)
            dst = os.path.join(web_static, name)
            shutil.copy2(src, dst)
            print('wrote', dst)


def make_hero_and_social():
    if not os.path.isfile(CRAB_SRC):
        print('skip hero/social — crab source not found:', CRAB_SRC)
        return
    crab = Image.open(CRAB_SRC).convert('RGBA')
    # Hero: square-ish, keep alpha, pad to neat size for the landing page.
    hero_side = 720
    # Fit crab into hero_side box preserving aspect.
    crab_fit = crab.copy()
    crab_fit.thumbnail((hero_side, hero_side), Image.LANCZOS)
    hero = Image.new('RGBA', (hero_side, hero_side), (0, 0, 0, 255))
    hx = (hero_side - crab_fit.width) // 2
    hy = (hero_side - crab_fit.height) // 2
    hero.alpha_composite(crab_fit, (hx, hy))
    save_png(hero, 'logo_hero_crab.png')

    # Social / OG card 1200×630: Crab left, wordmark + tagline right on black.
    W, H = 1200, 630
    card = Image.new('RGBA', (W, H), (10, 12, 16, 255))
    # Crab panel ~ square on the left.
    panel = 560
    crab_panel = crab.copy()
    crab_panel.thumbnail((panel - 40, panel - 40), Image.LANCZOS)
    cx = 40 + (panel - 40 - crab_panel.width) // 2
    cy = (H - crab_panel.height) // 2
    card.alpha_composite(crab_panel, (cx, cy))

    # Wordmark + small color mark on the right.
    mark = combo_mark(BANNER, 72)
    text_x = panel + 60
    card.alpha_composite(mark, (text_x, H // 2 - 90))
    d = ImageDraw.Draw(card)
    font_lg = _font(54)
    font_sm = _font(26)
    # "Multi" + colored "Color" + "Fits"
    tx, ty = text_x + 90, H // 2 - 80
    for text, col in [('Multi', (232, 232, 232))]:
        d.text((tx, ty), text, font=font_lg, fill=col + (255,))
        tx += d.textlength(text, font=font_lg)
    for ch, hexc in zip('Color', URAN):
        d.text((tx, ty), ch, font=font_lg, fill=hexrgb(hexc) + (255,))
        tx += d.textlength(ch, font=font_lg)
    d.text((tx, ty), 'Fits', font=font_lg, fill=(232, 232, 232, 255))
    d.text((text_x + 90, H // 2 + 10),
           'Colorize & combine FITS images',
           font=font_sm, fill=(170, 175, 185, 255))
    save_png(card.convert('RGBA'), 'logo_social.png')


def main():
    make_navbar_and_favicon()
    make_hero_and_social()
    print('done. Active navbar mark: logo_navbar-color.png')
    print('Swap option: logo_navbar-mono-{light,dark}.png')


if __name__ == '__main__':
    main()
