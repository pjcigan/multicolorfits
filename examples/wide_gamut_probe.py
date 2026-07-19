"""
Wide-gamut research probe.

Quantifies how often Lab compositing leaves the linear sRGB cube when layer
colors are ordinary sRGB hex (today's workflow):

* **overflow** — any linear channel > 1  (brightness / screen-blend headroom)
* **chromatic OOG** — any linear channel < 0  (true out-of-sRGB chroma)

Saves optional soft-proof masks into ``examples/output/``.

Usage (repo root)::

    python examples/wide_gamut_probe.py
    python examples/wide_gamut_probe.py --ngc602   # if local FITS are present
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np
from matplotlib import pyplot as plt
from skimage import color as ski_color

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import multicolorfits as mcf
from multicolorfits.core.colormix import _blend_stack

OUTDIR = os.path.join(os.path.dirname(__file__), 'output')
# Optional real-data probe (override with MCF_NGC602_DIR).
NGC602_DIR = os.environ.get(
    'MCF_NGC602_DIR',
    os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'hidden', 'testing', 'ngc602')))

# XYZ (D65) → linear sRGB, IEC 61966-2-1 (no encoding, no clip)
_XYZ_TO_LINEAR_SRGB = np.array([
    [3.2404542, -1.5371385, -0.4985314],
    [-0.9692660, 1.8760108, 0.0415560],
    [0.0556434, -0.2040259, 1.0572252],
])

BLENDS = ('screen', 'max', 'sum')


def lab_combine_linear_srgb(im_list_colorized, blend='screen', gamma=2.2,
                            weights=None):
    """
    Same Lab mix as ``combine_multicolor_colorspace(..., colorspace='lab')``,
    but return **linear** sRGB without clip (encoding / display gamut check).
    """
    n = len(im_list_colorized)
    if weights is None:
        weights = np.ones(n, dtype=float)
    else:
        weights = np.asarray(weights, dtype=float)
        if weights.shape != (n,):
            raise ValueError('weights must have one value per image')

    disp = []
    for im in im_list_colorized:
        arr = np.clip(np.nan_to_num(np.asarray(im, dtype=float)), 0, 1)
        disp.append(arr ** (1. / gamma) if gamma != 1 else arr)

    labs = [ski_color.rgb2lab(d) for d in disp]
    L = np.stack([lab[..., 0] / 100. for lab in labs])
    a = np.stack([lab[..., 1] for lab in labs])
    b = np.stack([lab[..., 2] for lab in labs])
    w = L * weights[:, None, None]
    wsum = w.sum(axis=0)
    eps = 1e-9
    safe = wsum > eps
    a_out = np.zeros_like(wsum)
    b_out = np.zeros_like(wsum)
    a_out[safe] = (w * a).sum(axis=0)[safe] / wsum[safe]
    b_out[safe] = (w * b).sum(axis=0)[safe] / wsum[safe]
    L_out = _blend_stack(L, blend)
    lab = np.stack([L_out * 100., a_out, b_out], axis=-1)
    xyz = ski_color.lab2xyz(lab)
    return xyz @ _XYZ_TO_LINEAR_SRGB.T


def oog_stats(lin, tol=1e-4):
    """Fraction of pixels with overflow (>1) vs chromatic OOG (<0)."""
    over = np.any(lin > 1.0 + tol, axis=-1)
    neg = np.any(lin < -tol, axis=-1)
    return {
        'overflow_frac': float(over.mean()),
        'chromatic_frac': float(neg.mean()),
        'lin_min': float(np.nanmin(lin)),
        'lin_max': float(np.nanmax(lin)),
        'over_mask': over,
        'neg_mask': neg,
    }


def colorize_planes(planes, colors, gamma=2.2):
    return [
        mcf.colorize_image(p, c, colorintype='hex', gammacorr_color=gamma)
        for p, c in zip(planes, colors)
    ]


def gaussian_blob(shape, cx, cy, radius):
    ny, nx = shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    p = np.clip(1.0 - np.hypot(xx - cx, yy - cy) / float(radius), 0, 1)
    return np.stack([p, p, p], axis=-1)


def save_softproof(name, lin, stats):
    """RGB soft-proof: green = overflow, magenta = chromatic OOG, grey = ok."""
    os.makedirs(OUTDIR, exist_ok=True)
    base = np.clip(lin, 0, 1)  # naive encode for backdrop
    # encode ~display for context
    display = np.power(base, 1.0 / 2.2)
    rgb = np.stack([display[..., 0], display[..., 1], display[..., 2]], axis=-1)
    over = stats['over_mask']
    neg = stats['neg_mask']
    rgb[over] = np.array([0.15, 0.95, 0.25])
    rgb[neg] = np.array([0.95, 0.15, 0.85])
    path = os.path.join(OUTDIR, 'wide_gamut_%s.png' % name)
    fig, ax = plt.subplots(figsize=(4.2, 4.2))
    ax.imshow(rgb, origin='lower', interpolation='nearest')
    ax.set_title(
        '%s\nover=%.1f%%  chrom=%.1f%%' % (
            name,
            100 * stats['overflow_frac'],
            100 * stats['chromatic_frac']),
        size=10)
    ax.set_xticks([]); ax.set_yticks([])
    fig.tight_layout()
    fig.savefig(path, dpi=110, facecolor='w')
    plt.close(fig)
    print('wrote', path)


def report_row(label, blend, stats):
    print('  %-28s  blend=%-6s  over=%5.1f%%  chrom=%5.1f%%  '
          'lin=[%+.3f, %+.3f]' % (
              label, blend,
              100 * stats['overflow_frac'],
              100 * stats['chromatic_frac'],
              stats['lin_min'], stats['lin_max']))


def run_case(label, layers, blends=BLENDS, save_masks=True):
    print('\n== %s ==' % label)
    for blend in blends:
        lin = lab_combine_linear_srgb(layers, blend=blend)
        st = oog_stats(lin)
        report_row(label, blend, st)
        if save_masks and blend == 'screen':
            safe = ''.join(c if c.isalnum() or c in '-_' else '_'
                           for c in label.lower())
            save_softproof('%s_%s' % (safe, blend), lin, st)


def synthetic_cases(save_masks=True):
    ny = nx = 128
    shape = (ny, nx)

    pob_planes = [
        gaussian_blob(shape, 40, 64, 45),
        gaussian_blob(shape, 64, 64, 45),
        gaussian_blob(shape, 88, 64, 45),
    ]
    run_case(
        'POB NGC-like',
        colorize_planes(pob_planes, ['#a020f0', '#ff8c00', '#1e90ff']),
        save_masks=save_masks)

    run_case(
        'complementary O+B',
        colorize_planes(
            [gaussian_blob(shape, 40, 64, 45),
             gaussian_blob(shape, 88, 64, 45)],
            ['#FF8C00', '#0000FF']),
        save_masks=save_masks)

    wheel_cols = ['#FF0000', '#FFFF00', '#00FF00',
                  '#00FFFF', '#0000FF', '#FF00FF']
    wheel_planes = [
        gaussian_blob(shape, 30 + 15 * i, 40 + 10 * (i % 2), 42)
        for i in range(6)
    ]
    run_case(
        'max-sat 6-wheel',
        colorize_planes(wheel_planes, wheel_cols),
        save_masks=save_masks)


def ngc602_case(save_masks=True, max_size=800):
    import astropy.io.fits as pyfits

    paths = {
        'ir': os.path.join(NGC602_DIR, 'ngc602_ir.fits'),
        'R': os.path.join(NGC602_DIR, 'ngc602_optical_R.fits'),
        'B': os.path.join(NGC602_DIR, 'ngc602_optical_B.fits'),
    }
    missing = [p for p in paths.values() if not os.path.isfile(p)]
    if missing:
        print('\n[skip NGC 602] missing files under %s' % NGC602_DIR)
        return

    print('\nLoading NGC 602 (downsample max_size=%d) ...' % max_size)
    datas = [
        mcf.downsample_for_preview(
            np.asarray(pyfits.getdata(paths[k]), dtype=np.float64),
            max_size=max_size)
        for k in ('ir', 'R', 'B')
    ]
    # Tutorial-ish POB hexes used in compare_colorspaces.py
    colors = ['#BE599E', '#DEA215', '#77C0F9']
    layers = [
        mcf.colorize_image(
            mcf.to_grey_rgb(d, rescalefn='linear', gamma=2.2),
            c, colorintype='hex', gammacorr_color=2.2)
        for d, c in zip(datas, colors)
    ]
    run_case('NGC602 tutorial POB', layers, save_masks=save_masks)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--ngc602', action='store_true',
                   help='Also probe real NGC 602 FITS (if present)')
    p.add_argument('--no-masks', action='store_true',
                   help='Skip writing soft-proof PNGs')
    args = p.parse_args(argv)

    print('Wide-gamut probe — Lab mix → unclipped linear sRGB')
    print('overflow = channel > 1; chromatic = channel < 0')
    synthetic_cases(save_masks=not args.no_masks)
    if args.ngc602:
        ngc602_case(save_masks=not args.no_masks)
    print('\nDone. Soft-proof masks (if any) are in', OUTDIR)
    print('Conclusion template: with sRGB hex inputs, Lab OOG is almost always')
    print('screen-blend brightness overflow — not chromatic Display-P3 territory.')


if __name__ == '__main__':
    main()
