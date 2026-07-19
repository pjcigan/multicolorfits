"""
Agent-facing capability map: keep ``mcf.overview()`` / ``mcf.recipes()``,
``llms.txt`` / ``llms-full.txt``, and the real API in lockstep.

* Every recipe in ``multicolorfits._overview.RECIPES`` executes against
  synthetic data (skipped only when an OPTIONAL dependency is missing).
* The committed ``llms.txt`` / ``llms-full.txt`` must equal a fresh render.
* Teaching errors at the high-confusion API boundaries stay actionable.
"""

from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pytest
from astropy.io import fits
from astropy.wcs import WCS

import multicolorfits as mcf
from multicolorfits._overview import RECIPES, render_llms

ROOT = Path(__file__).resolve().parent.parent


def _header(ny, nx):
    w = WCS(naxis=2)
    w.wcs.crpix = [nx / 2.0, ny / 2.0]
    w.wcs.cdelt = [-0.0008, 0.0008]
    w.wcs.crval = [150.0, -30.0]
    w.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    hdr = w.to_header()
    hdr['NAXIS'] = 2
    hdr['NAXIS1'] = nx
    hdr['NAXIS2'] = ny
    return hdr


def _blob(xx, yy, cx, cy, s):
    return np.exp(-(((xx - cx) ** 2 + (yy - cy) ** 2) / (2.0 * s ** 2)))


def _bindings(tmp_path):
    ny, nx = 48, 40
    yy, xx = np.mgrid[0:ny, 0:nx]
    data = _blob(xx, yy, 20, 24, 8) * 10.0
    data_r = _blob(xx, yy, 15, 24, 7) * 10.0
    data_g = _blob(xx, yy, 20, 24, 7) * 10.0
    data_b = _blob(xx, yy, 25, 24, 7) * 10.0
    header = _header(ny, nx)

    paths = {}
    for key, arr in (('fits_r', data_r), ('fits_g', data_g), ('fits_b', data_b)):
        p = tmp_path / (key + '.fits')
        fits.writeto(p, arr.astype('float32'), header, overwrite=True)
        paths[key] = str(p)

    # A ready session (several recipes reference `s`).  Load from the temp
    # FITS files so save_state/load_state round-trips (state stores paths).
    s = mcf.McfSession(n_panels=3)
    s.load_files([paths['fits_r'], paths['fits_g'], paths['fits_b']],
                 colors=['#E4002B', '#33CC33', '#0088FF'],
                 labels=['R', 'G', 'B'])
    for panel in s.panels:
        panel.stretch = 'asinh'
        panel.set_percentiles(1, 99.5)
    s.compose.gamma = 2.2

    layers = [
        mcf.colorize_image(mcf.to_grey_rgb(data_r, rescalefn='asinh'), '#E4002B'),
        mcf.colorize_image(mcf.to_grey_rgb(data_b, rescalefn='asinh'), '#0088FF'),
    ]
    combined = mcf.combine_multicolor(layers, gamma=2.2)
    rgba = mcf.make_transparent_cutout(combined, crop='auto', pad=0.05)

    ns = dict(
        np=np, plt=plt, mcf=mcf,
        data=data, data_r=data_r, data_g=data_g, data_b=data_b,
        header=header, s=s, layers=layers, combined=combined, rgba=rgba,
        **paths,
    )
    return ns


@pytest.mark.parametrize('recipe', RECIPES, ids=lambda r: r.task[:45])
def test_recipe_runs(recipe, tmp_path, monkeypatch):
    """Each catalog recipe executes cleanly (skipped only if an OPTIONAL
    dependency it needs isn't installed)."""
    monkeypatch.chdir(tmp_path)
    ns = _bindings(tmp_path)
    code = '\n'.join(ln for ln in recipe.code.splitlines()
                     if not ln.strip().endswith('.show()'))
    try:
        exec(code, dict(ns))
    except ImportError as exc:
        pytest.skip(f'optional dependency missing: {exc}')
    finally:
        plt.close('all')


@pytest.mark.parametrize('name,full', [('llms.txt', False),
                                       ('llms-full.txt', True)])
def test_llms_txt_in_sync(name, full):
    """The committed llms.txt / llms-full.txt must match a fresh render of the
    catalog — otherwise run `python scripts/make_llms_txt.py`."""
    committed = (ROOT / name).read_text()
    assert committed == render_llms(full=full), (
        f'{name} is stale — regenerate with `python scripts/make_llms_txt.py`')


def test_llms_txt_carries_recipe_code():
    """llms.txt must contain runnable CODE, not just a function index (an agent
    reading only the concise file otherwise guesses arguments)."""
    text = render_llms(full=False)
    assert '```python' in text
    assert 'to_grey_rgb(' in text and 'combine_multicolor(' in text
    assert "colorintype='hex'" in text


def test_overview_as_dict_shape():
    cat = mcf.overview(as_dict=True)
    assert set(cat) == {'layer_first', 'conventions', 'recipes'}
    assert cat['recipes'] and all(
        {'task', 'category', 'functions', 'code', 'notes'} <= set(r)
        for r in cat['recipes'])


def test_recipes_query_matches(capsys):
    mcf.recipes('cutout')
    out = capsys.readouterr().out
    assert 'make_transparent_cutout' in out


def test_recipes_no_match_is_helpful(capsys):
    mcf.recipes('definitely-not-a-topic')
    out = capsys.readouterr().out
    assert 'No recipe matched' in out


class TestTeachingErrors:
    def test_colorize_image_rejects_2d_data(self):
        with pytest.raises(TypeError, match='greyscale RGB cube'):
            mcf.colorize_image(np.zeros((8, 8)), '#FF0000')

    def test_combine_multicolor_rejects_single_array(self):
        with pytest.raises(TypeError, match='LIST of colorized layers'):
            mcf.combine_multicolor(np.zeros((8, 8, 3)))

    def test_combine_multicolor_rejects_2d_layers(self):
        with pytest.raises(TypeError, match=r'shape \(ny, nx, 3\)'):
            mcf.combine_multicolor([np.zeros((8, 8)), np.zeros((8, 8))])
