#!/usr/bin/env python3
"""Regenerate the browser GUI brand mark from combo_swatch.

Default matches ComposeState: rgb mode, screen blend (for lab/hsv), gamma 2.2,
black background, POB palette.

Writes a real PNG (``mcf-mark.png``) for ``<img src=...>``. Nested PNG-in-SVG
is not reliable as an ``<img>`` target, and SVG comments must not contain ``--``.

  python render_mcf_mark.py
  python render_mcf_mark.py --mode lab --size 256 -o /tmp/mcf-mark.png
"""
from __future__ import annotations

import argparse
import os

import numpy as np
from PIL import Image

from multicolorfits.core.swatch import combo_swatch

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_COLORS = ['#BE599E', '#DEA215', '#77C0F9']  # POB: purple, orange, blue


def render_mark(colors, *, mode='rgb', blend='screen', gamma=2.2, background='black',
                size=128, out_path=None):
    data = combo_swatch(colors, mode=mode, blend=blend, gamma=gamma,
                        background=background, size=size)
    if data is None:
        raise SystemExit('combo_swatch returned None (empty colors?)')
    out_path = out_path or os.path.join(HERE, 'mcf-mark.png')
    Image.fromarray((np.clip(data['rgba'], 0, 1) * 255).astype(np.uint8), 'RGBA').save(
        out_path, format='PNG')
    print('wrote', out_path, '(%dx%d, mode=%s, blend=%s)' % (size, size, mode, blend))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--colors', nargs='+', default=DEFAULT_COLORS)
    p.add_argument('--mode', default='rgb')
    p.add_argument('--blend', default='screen')
    p.add_argument('--gamma', type=float, default=2.2)
    p.add_argument('--background', default='black')
    p.add_argument('--size', type=int, default=128)
    p.add_argument('-o', '--output', default=None)
    args = p.parse_args()
    render_mark(args.colors, mode=args.mode, blend=args.blend, gamma=args.gamma,
                background=args.background, size=args.size, out_path=args.output)


if __name__ == '__main__':
    main()
