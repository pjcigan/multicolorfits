"""
Real-data check: M74 (NGC 628) transparent cutout, matching the
skyplothelper markers tutorial recipe.

Data source (first hit wins):
  1. ``$MCF_M74_CACHE`` or ``$SPH_QUERY_CACHE`` directory
  2. Sibling checkout ``../skyplothelper/examples/data/query_cache``
     (next to this multicolorfits repo)

Requires the SDSS gri FITS cutouts written by the sph tutorial's
``cached_band`` helper (``m74_sdss_{g,r,i}.fits``).  Skips cleanly when
the cache is absent so CI without the files stays green.
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pytest

import multicolorfits as mcf


def _candidate_cache_dirs():
    env = os.environ.get('MCF_M74_CACHE') or os.environ.get('SPH_QUERY_CACHE')
    if env:
        yield Path(env).expanduser()
    repo = Path(__file__).resolve().parents[1]  # .../multicolorfits
    yield repo.parent / 'skyplothelper' / 'examples' / 'data' / 'query_cache'
    yield Path.home() / 'programs' / 'github' / 'skyplothelper' / 'examples' / 'data' / 'query_cache'


def _find_m74_cache():
    need = ('m74_sdss_g.fits', 'm74_sdss_r.fits', 'm74_sdss_i.fits')
    for d in _candidate_cache_dirs():
        if d.is_dir() and all((d / name).is_file() for name in need):
            return d
    return None


M74_CACHE = _find_m74_cache()
requires_m74 = pytest.mark.skipif(
    M74_CACHE is None,
    reason='M74 SDSS cache not found (set MCF_M74_CACHE or clone skyplothelper next to this repo)',
)


def _stretch(chan, lo=30, hi=99.5, soft=0.12):
    """Tutorial arcsinh stretch — background-subtract then soft arcsinh to [0, 1]."""
    base = np.percentile(chan, lo)
    x = np.clip(chan - base, 0, None)
    top = np.percentile(x, hi)
    if top <= 0:
        return np.zeros_like(x)
    return np.arcsinh(x / top / soft) / np.arcsinh(1 / soft)


def _composite(red, grn, blu):
    return np.clip(np.dstack([_stretch(red), _stretch(grn), _stretch(blu)]), 0, 1)


def _sph_to_transparent(rgb, lo=55, hi=99.3, gamma=0.5):
    """Exact numpy recipe from the skyplothelper markers tutorial."""
    lum = 0.2126 * rgb[..., 0] + 0.7152 * rgb[..., 1] + 0.0722 * rgb[..., 2]
    a_lo, a_hi = np.percentile(lum, lo), np.percentile(lum, hi)
    alpha = np.clip((lum - a_lo) / (a_hi - a_lo + 1e-9), 0, 1) ** gamma
    return np.dstack([rgb, alpha])


def _load_m74_sdss_rgb():
    from astropy.io import fits
    g = np.nan_to_num(np.asarray(fits.getdata(M74_CACHE / 'm74_sdss_g.fits'), float))
    r = np.nan_to_num(np.asarray(fits.getdata(M74_CACHE / 'm74_sdss_r.fits'), float))
    i_ = np.nan_to_num(np.asarray(fits.getdata(M74_CACHE / 'm74_sdss_i.fits'), float))
    # R<-i, G<-r, B<-g  (tutorial order)
    return _composite(i_, r, g)


@requires_m74
class TestM74TransparentCutout:
    def test_matches_sph_numpy_recipe(self):
        rgb = _load_m74_sdss_rgb()
        assert rgb.shape[-1] == 3 and min(rgb.shape[:2]) >= 200

        sph = _sph_to_transparent(rgb, lo=55, hi=99.3, gamma=0.5)
        mcf_rgba = mcf.make_transparent_cutout(
            rgb, alpha_lo=55, alpha_hi=99.3, alpha_gamma=0.5,
            crop='none', size=None,
        )

        # Same RGB, nearly identical alpha (mcf uses finite-only percentiles
        # and a tiny floor when hi==lo — irrelevant for this image).
        assert np.allclose(mcf_rgba[..., :3], sph[..., :3], atol=1e-6)
        assert np.allclose(mcf_rgba[..., 3], sph[..., 3], atol=2e-3)

        # Galaxy body opaque; outer sky nearly clear.
        cy, cx = rgb.shape[0] // 2, rgb.shape[1] // 2
        assert mcf_rgba[cy, cx, 3] > 0.85
        assert mcf_rgba[2, 2, 3] < 0.05
        assert mcf_rgba[-3, -3, 3] < 0.05

    def test_gamma_boldens_disk_vs_linear(self):
        rgb = _load_m74_sdss_rgb()
        a05 = mcf.make_transparent_cutout(
            rgb, alpha_gamma=0.5, crop='none', size=None)[..., 3]
        a10 = mcf.make_transparent_cutout(
            rgb, alpha_gamma=1.0, crop='none', size=None)[..., 3]
        mid = (a10 > 0.15) & (a10 < 0.85)
        assert mid.any()
        assert float(np.mean(a05[mid])) > float(np.mean(a10[mid]))

    def test_auto_crop_and_stamp_size(self, tmp_path):
        rgb = _load_m74_sdss_rgb()
        stamp = mcf.make_transparent_cutout(
            rgb, size=128, alpha_gamma=0.5, crop='auto', pad=0.05,
        )
        assert stamp.shape[-1] == 4
        assert max(stamp.shape[:2]) == 128
        # Crop should have removed empty sky vs the full 400 px frame.
        assert stamp.shape[0] * stamp.shape[1] < rgb.shape[0] * rgb.shape[1]

        path = tmp_path / 'm74_stamp.png'
        mcf.save_transparent_cutout(stamp, str(path))
        assert path.is_file() and path.stat().st_size > 500

    def test_session_export_from_fits(self, tmp_path):
        """Compose M74 with mcf's own pipeline, then export a transparent stamp."""
        from astropy.io import fits

        s = mcf.McfSession()
        colors = ['#FFAA55', '#55FF88', '#5588FF']  # warm / mid / blue-ish
        for idx, (band, color) in enumerate(zip(
                ('m74_sdss_i.fits', 'm74_sdss_r.fits', 'm74_sdss_g.fits'), colors)):
            data = np.nan_to_num(np.asarray(
                fits.getdata(M74_CACHE / band), float))
            hdr = fits.getheader(M74_CACHE / band)
            s.panels[idx].set_data(data, hdr)
            s.panels[idx].color = color
            s.panels[idx].stretch = 'asinh'
            s.panels[idx].set_percentiles(30.0, 99.5)

        s.compose.combine_mode = 'rgb'
        s.compose.combine_background = 'black'
        out = tmp_path / 'm74_session_stamp.png'
        rgba = s.export_transparent_cutout(
            size=160, savepath=str(out), alpha_gamma=0.5, crop='auto',
        )
        assert rgba.shape[-1] == 4
        assert max(rgba.shape[:2]) == 160
        assert out.is_file()
        assert s.compose.combine_background == 'black'
        cy, cx = rgba.shape[0] // 2, rgba.shape[1] // 2
        assert rgba[cy, cx, 3] > 0.5
        assert rgba[0, 0, 3] < 0.15
