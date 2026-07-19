#!/usr/bin/env python3
"""Regenerate mcf-mark.svg from combo_swatch (honest overlap colors).

Default matches ComposeState: rgb mode, screen blend (for lab/hsv), gamma 2.2,
black background, POB tutorial palette.

  python render_mcf_mark.py
  python render_mcf_mark.py --mode lab --size 256 -o /tmp/mcf-mark.svg
"""
from __future__ import annotations

import argparse
import base64
import io
import os

import numpy as np
from PIL import Image

from multicolorfits.core.swatch import combo_swatch

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_COLORS = ['#BE599E', '#DEA215', '#77C0F9']  # POB: purple, orange, blue


def render_svg(colors, *, mode='rgb', blend='screen', gamma=2.2, background='black',
               size=128, out_path=None):
    data = combo_swatch(colors, mode=mode, blend=blend, gamma=gamma,
                        background=background, size=size)
    if data is None:
        raise SystemExit('combo_swatch returned None (empty colors?)')
    buf = io.BytesIO()
    Image.fromarray((np.clip(data['rgba'], 0, 1) * 255).astype(np.uint8), 'RGBA').save(
        buf, format='PNG')
    b64 = base64.b64encode(buf.getvalue()).decode('ascii')
    svg = f'''<?xml version="1.0" encoding="UTF-8"?>
<!-- MultiColorFits mark: combo_swatch ({mode!r}, blend={blend!r}, gamma={gamma}, bg={background!r}).
     Regenerate: python {os.path.basename(__file__)} [--mode lab] [--size 256] -->
<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink"
     viewBox="0 0 32 32" role="img" aria-label="MultiColorFits">
  <image width="32" height="32" xlink:href="data:image/png;base64,{b64}"/>
</svg>
'''
    out_path = out_path or os.path.join(HERE, 'mcf-mark.svg')
    with open(out_path, 'w', encoding='utf-8') as f:
        f.write(svg)
    print('wrote', out_path)


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
    render_svg(args.colors, mode=args.mode, blend=args.blend, gamma=args.gamma,
               background=args.background, size=args.size, out_path=args.output)


if __name__ == '__main__':
    main()
