# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     kernelspec:
#       display_name: Python 3
#       name: python3
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.4
# ---

# %% [markdown]
# # Prepare and align a rotated image
#
# Walk through the processing path in
# {doc}`../guide/preparing_images`: tidy the header, put north up in the
# image's own frame, oversample, and (if a second layer is present) crop to
# the overlap.
#
# Needs a local FITS crop with a celestial WCS. Point ``MCF_PREP_FITS`` at a
# file, or place the Chandra NGC 602 IR frame where the NGC 602 tutorial
# looks (``MCF_NGC602_DIR`` / ``hidden/testing/ngc602``). The cell skips
# cleanly when no file is available, so this page still builds on Read the
# Docs.

# %%
import os
from pathlib import Path

import numpy as np

import multicolorfits as mcf

# %%
def _find_fits():
    env = os.environ.get('MCF_PREP_FITS')
    if env and Path(env).expanduser().is_file():
        return Path(env).expanduser()
    root = Path(mcf.__file__).resolve().parents[1]
    ngc = os.environ.get('MCF_NGC602_DIR')
    candidates = []
    if ngc:
        candidates.append(Path(ngc).expanduser() / 'ngc602_ir.fits')
    candidates.append(root / 'hidden' / 'testing' / 'ngc602' / 'ngc602_ir.fits')
    for path in candidates:
        if path.is_file():
            return path
    return None

fits_path = _find_fits()
print(fits_path or 'no local FITS — remaining cells will skip')

# %% [markdown]
# ## Read and describe
#
# ``read_fits`` runs ``tidy_header`` so conflicting PC/CD/CROTA cards are
# dropped and a mirrored WCS is announced. ``describe_header`` /
# ``describe_image`` print by default; pass ``verbose=False`` to keep the
# dict. For a set of layers plus suggested colors, see
# ``describe_images`` in {doc}`../guide/intensity_scaling`.

# %%
if fits_path is None:
    data = hdr = None
    print('skip')
else:
    data, hdr = mcf.read_fits(fits_path)
    info = mcf.describe_header(hdr, name=fits_path.name)
    print('flipped:', info['flipped'], 'scale arcsec:', info['pixscale_arcsec'])
    levels = mcf.describe_image(data, hdr, name=fits_path.name, verbose=False)['levels']
    print('levels:', levels['stretch'], levels['vmin'], levels['vmax'])

# %% [markdown]
# ## North-up, same frame
#
# ``rotation_deg=0`` is north-up in the header frame (Galactic stays
# Galactic). ``oversample=2`` halves the pixel scale so the resample is less
# blocky. This does not invent resolution.

# %%
if data is None:
    print('skip')
else:
    # Downsample first so the tutorial stays light on a full HST/Chandra frame.
    small = mcf.downsample_for_preview(data, max_size=400)
    # Scale the header to match the preview grid.
    factor = data.shape[1] / float(small.shape[1])
    hdr_small = hdr.copy()
    hdr_small['NAXIS1'] = small.shape[1]
    hdr_small['NAXIS2'] = small.shape[0]
    hdr_small['CRPIX1'] = (hdr['CRPIX1'] - 1) / factor + 1
    hdr_small['CRPIX2'] = (hdr['CRPIX2'] - 1) / factor + 1
    cd1, cd2 = mcf.get_cdelts(hdr)
    hdr_small['CDELT1'] = cd1 * factor
    hdr_small['CDELT2'] = cd2 * factor
    for key in ('PC1_1', 'PC1_2', 'PC2_1', 'PC2_2', 'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2'):
        if key in hdr_small and key.startswith('CD'):
            hdr_small[key] = float(hdr[key]) * factor

    north, hdr_n = mcf.reproject_north_up(small, hdr_small, oversample=2, order=1)
    print('input', small.shape, 'north-up', north.shape, hdr_n['CTYPE1'], hdr_n['CTYPE2'])
    print('CDELT', float(hdr_n['CDELT1']), float(hdr_n['CDELT2']))

# %% [markdown]
# ## Optional second layer and overlap crop
#
# If a sibling optical frame sits next to the file, ``prep_layers`` aligns
# both and drops the NaN margin.

# %%
if data is None:
    print('skip')
else:
    sibling = fits_path.with_name('ngc602_optical_R.fits')
    if not sibling.is_file():
        print('no second layer at', sibling, '— north-up result stands alone')
    else:
        other, ohdr = mcf.read_fits(sibling)
        prepared = mcf.prep_layers(
            [(small, hdr_small), (mcf.downsample_for_preview(other, max_size=400), ohdr)],
            north_up=True, oversample=1, crop='overlap')
        print('prepared shapes', [p[0].shape for p in prepared])
