"""Tests for transparent cutout helpers (slides / stamps)."""
import os

import numpy as np
import pytest

import multicolorfits as mcf
from multicolorfits.core.cutout import (
    alpha_from_intensity,
    batch_transparent_cutouts,
    deblend_background,
    flatten_rgba,
    luminance,
    make_checkerboard,
    make_transparent_cutout,
    preview_cutout_on_backgrounds,
    save_transparent_cutout,
)
from conftest import make_test_data, make_test_header


def _blob_rgb(ny=80, nx=100, soft=True):
    """Bright Gaussian blob on a black field (RGB)."""
    yy, xx = np.mgrid[0:ny, 0:nx].astype(float)
    cy, cx = ny / 2.0, nx / 2.0
    r2 = ((yy - cy) / (ny * 0.18)) ** 2 + ((xx - cx) / (nx * 0.18)) ** 2
    inten = np.exp(-0.5 * r2)
    if soft:
        inten = np.clip(inten, 0, 1)
    rgb = np.zeros((ny, nx, 3), dtype=float)
    rgb[..., 0] = inten
    rgb[..., 1] = 0.6 * inten
    rgb[..., 2] = 0.2 * inten
    return rgb


class TestLuminanceAndAlpha:
    def test_luminance_2d_passthrough(self):
        a = np.linspace(0, 1, 12).reshape(3, 4)
        assert np.allclose(luminance(a), a)

    def test_luminance_rgb(self):
        rgb = np.zeros((2, 2, 3))
        rgb[..., 1] = 1.0  # pure green
        assert luminance(rgb)[0, 0] == pytest.approx(0.7152)

    def test_alpha_percentile_and_gamma(self):
        inten = _blob_rgb()[..., 0]
        a_lin = alpha_from_intensity(inten, lo=50, hi=99, gamma=1.0)
        a_bold = alpha_from_intensity(inten, lo=50, hi=99, gamma=0.5)
        # At mid-bright pixels, gamma < 1 raises alpha (bolder body).
        mid = (a_lin > 0.1) & (a_lin < 0.9)
        assert mid.any()
        assert np.nanmean(a_bold[mid]) > np.nanmean(a_lin[mid])

    def test_alpha_gamma_above_one_tightens(self):
        inten = _blob_rgb()[..., 0]
        a_lin = alpha_from_intensity(inten, lo=50, hi=99, gamma=1.0)
        a_tight = alpha_from_intensity(inten, lo=50, hi=99, gamma=1.4)
        mid = (a_lin > 0.1) & (a_lin < 0.9)
        assert mid.any()
        # gamma > 1 lowers alpha in the fade zone (tighter matte)...
        assert np.nanmean(a_tight[mid]) < np.nanmean(a_lin[mid])
        # ...but pixels at/above the hi percentile stay fully opaque.
        assert a_tight[a_lin >= 1.0].min() == pytest.approx(1.0)

    def test_alpha_invert(self):
        inten = np.linspace(0, 1, 100).reshape(10, 10)
        a = alpha_from_intensity(inten, lo=10, hi=90, gamma=1.0, invert=True)
        assert a[0, 0] > a[-1, -1]


class TestMakeTransparentCutout:
    def test_shape_and_range(self):
        rgba = make_transparent_cutout(_blob_rgb(), crop='none')
        assert rgba.shape == (80, 100, 4)
        assert rgba.min() >= 0.0 and rgba.max() <= 1.0
        # Corners should be mostly transparent; center opaque.
        assert rgba[0, 0, 3] < 0.05
        assert rgba[40, 50, 3] > 0.5

    def test_auto_crop_shrinks(self):
        rgba = make_transparent_cutout(_blob_rgb(120, 160), crop='auto', pad=0.05)
        assert rgba.shape[0] < 120 and rgba.shape[1] < 160
        assert rgba.shape[-1] == 4

    def test_size_downsample(self):
        rgba = make_transparent_cutout(_blob_rgb(200, 200), size=64, crop='none')
        assert max(rgba.shape[:2]) == 64

    def test_circle_matte(self):
        rgba = make_transparent_cutout(_blob_rgb(), matte='circle', crop='none',
                                       soft_edge=0.0)
        # Far corner outside inscribed circle → near-zero alpha.
        assert rgba[0, 0, 3] < 0.02

    def test_soft_edge_requires_scipy(self):
        pytest.importorskip('scipy')
        rgba = make_transparent_cutout(_blob_rgb(), soft_edge=2.0, crop='none')
        assert rgba.shape[-1] == 4

    def test_alpha_source_max_mean_layer(self):
        rgb = _blob_rgb()
        a_max = make_transparent_cutout(rgb, alpha_source='max', crop='none')
        a_mean = make_transparent_cutout(rgb, alpha_source='mean', crop='none')
        a_lay = make_transparent_cutout(rgb, alpha_source='layer', layer=0, crop='none')
        assert a_max.shape == a_mean.shape == a_lay.shape

    def test_alpha_smooth_fills_pinholes(self):
        pytest.importorskip('scipy')
        rng = np.random.default_rng(42)
        rgb = _blob_rgb()
        # Salt the bright body with dark pinholes (per-pixel noise).
        holes = rng.random(rgb.shape[:2]) < 0.05
        noisy = rgb.copy()
        noisy[holes] = 0.0
        kw = dict(alpha_lo=50, alpha_hi=90, alpha_gamma=1.0, crop='none')
        a_raw = make_transparent_cutout(noisy, alpha_smooth=0, **kw)[..., 3]
        a_sm = make_transparent_cutout(noisy, alpha_smooth=2.5, **kw)[..., 3]
        body = rgb[..., 0] > 0.6
        # Smoothing keeps the noisy body coherent/opaque...
        assert a_sm[body].mean() > a_raw[body].mean()
        assert a_sm[body].min() > 0.5
        # ...while the far corners stay fully transparent.
        assert a_sm[0, 0] < 0.02 and a_sm[-1, -1] < 0.02

    def test_alpha_source_sat_and_dist(self):
        # Colorful blob on a white field: luminance can't separate it, but
        # saturation and distance-from-white can.
        rgb = np.ones((60, 60, 3))
        rgb[20:40, 20:40] = [0.9, 0.2, 0.1]
        kw = dict(alpha_lo=50, alpha_hi=95, alpha_gamma=1.0, crop='none')
        a_sat = make_transparent_cutout(rgb, alpha_source='sat', **kw)[..., 3]
        a_dist = make_transparent_cutout(rgb, alpha_source='dist',
                                         sky_color='white', **kw)[..., 3]
        for a in (a_sat, a_dist):
            assert a[30, 30] > 0.9    # colored blob survives
            assert a[5, 5] < 0.05     # white "sky" is transparent
        # 'dist' also keeps near-black features on white; 'sat' does not.
        rgb2 = np.ones((60, 60, 3))
        rgb2[20:40, 20:40] = 0.05     # dark grey blob
        a2 = make_transparent_cutout(rgb2, alpha_source='dist',
                                     sky_color='white', **kw)[..., 3]
        assert a2[30, 30] > 0.9 and a2[5, 5] < 0.05

    def test_existing_alpha_modes(self):
        rgb = _blob_rgb()
        base = make_transparent_cutout(rgb, crop='none')
        keep = make_transparent_cutout(base, existing_alpha='keep', crop='none')
        assert np.allclose(keep[..., 3], base[..., 3])
        mult = make_transparent_cutout(base, existing_alpha='multiply',
                                      alpha_lo=0, alpha_hi=100, alpha_gamma=1.0,
                                      crop='none')
        assert mult[..., 3].max() <= base[..., 3].max() + 1e-6

    def test_batch(self):
        outs = batch_transparent_cutouts([_blob_rgb(), _blob_rgb()], size=32, crop='none')
        assert len(outs) == 2
        assert outs[0].shape[0] <= 32

    def test_sky_box_crop(self):
        from astropy.wcs import WCS
        hdr = make_test_header(nx=100, ny=80)
        wcs = WCS(hdr)
        rgb = _blob_rgb(80, 100)
        # Small box around CRVAL should shrink the frame.
        sky_box = (149.995, -30.005, 150.005, -29.995)
        rgba = make_transparent_cutout(rgb, wcs=wcs, sky_box=sky_box,
                                       crop='none', size=None)
        assert rgba.shape[0] < 80 and rgba.shape[1] < 100


class TestSaveAndPreview:
    def test_save_png(self, tmp_path):
        rgba = make_transparent_cutout(_blob_rgb(), size=48, crop='auto')
        path = tmp_path / 'stamp.png'
        out = save_transparent_cutout(rgba, str(path))
        assert os.path.isfile(out)
        assert os.path.getsize(out) > 0

    def test_save_flips_for_image_viewers(self, tmp_path):
        """PNG row 0 is the top of the sky (FITS row -1), matching GUI exports."""
        pytest.importorskip('PIL')
        from PIL import Image
        rgba = np.zeros((40, 30, 4), dtype=float)
        rgba[0, :, :] = (1.0, 0.0, 0.0, 1.0)   # FITS bottom row = red
        rgba[-1, :, :] = (0.0, 0.0, 1.0, 1.0)  # FITS top row = blue
        path = tmp_path / 'orient.png'
        save_transparent_cutout(rgba, str(path))
        saved = np.asarray(Image.open(path))
        assert saved[0, 15, 2] > 200   # file top == blue (sky north)
        assert saved[-1, 15, 0] > 200  # file bottom == red

    def test_save_tiff(self, tmp_path):
        pytest.importorskip('PIL')
        rgba = make_transparent_cutout(_blob_rgb(), size=32, crop='none')
        path = tmp_path / 'stamp.tif'
        out = save_transparent_cutout(rgba, str(path))
        assert os.path.isfile(out)

    def test_checkerboard_shape(self):
        bg = make_checkerboard(40, 50, cell=8)
        assert bg.shape == (40, 50, 3)

    def test_preview_figure(self):
        rgba = make_transparent_cutout(_blob_rgb(), size=40, crop='none')
        fig = preview_cutout_on_backgrounds(rgba)
        assert fig is not None
        import matplotlib.pyplot as plt
        plt.close(fig)


class TestSessionExport:
    def test_export_transparent_cutout(self, tmp_path):
        s = mcf.McfSession()
        s.panels[0].set_data(make_test_data(seed=1), make_test_header())
        s.panels[0].color = '#FF4400'
        s.compose.combine_background = 'black'  # must not stay transparent after
        path = tmp_path / 'sess_stamp.png'
        rgba = s.export_transparent_cutout(size=48, savepath=str(path), crop='auto')
        assert rgba.shape[-1] == 4
        assert os.path.isfile(path)
        assert s.compose.combine_background == 'black'

    def test_export_preserves_black_rgb_colors(self):
        """Classic black RGB cutouts must match render_combined (not the alpha path)."""
        s = mcf.McfSession()
        s.panels[0].set_data(make_test_data(seed=2), make_test_header())
        s.panels[0].color = '#A40000'
        s.panels[1].set_data(make_test_data(seed=3), make_test_header())
        s.panels[1].color = '#E9B96E'
        s.compose.combine_background = 'black'
        s.compose.gamma = 2.2
        truth = s.render_combined()
        # Soft matte: body stays fully opaque so RGB is unchanged there.
        rgba = s.export_transparent_cutout(
            crop='none', alpha_lo=0, alpha_hi=10, alpha_gamma=1.0,
            alpha_smooth=0)
        opaque = rgba[..., 3] > 0.99
        assert opaque.any()
        assert np.allclose(rgba[opaque, :3], truth[opaque], atol=1e-6)
        # Forcing transparent render must diverge for multi-layer RGB.
        blown = s.export_transparent_cutout(
            crop='none', alpha_lo=0, alpha_hi=10, alpha_gamma=1.0,
            alpha_smooth=0, render_background='transparent')
        assert not np.allclose(blown[opaque, :3], truth[opaque], atol=1e-3)

    def test_public_exports(self):
        assert callable(mcf.make_transparent_cutout)
        assert callable(mcf.save_transparent_cutout)
        assert callable(mcf.preview_cutout_on_backgrounds)
        assert callable(mcf.deblend_background)
        assert callable(mcf.flatten_rgba)


def _blob_rgba(ny=60, nx=70):
    """Colorful blob with a smoothly varying alpha skirt (RGBA test image)."""
    yy, xx = np.mgrid[0:ny, 0:nx].astype(float)
    r2 = ((yy - ny / 2) / (ny * 0.22)) ** 2 + ((xx - nx / 2) / (nx * 0.22)) ** 2
    alpha = np.clip(np.exp(-0.5 * r2), 0, 1)
    alpha[alpha < 0.02] = 0.0
    rgb = np.zeros((ny, nx, 3))
    rgb[..., 0] = 0.9
    rgb[..., 1] = 0.3 + 0.4 * (xx / nx)
    rgb[..., 2] = 0.15
    return np.dstack([rgb, alpha])


class TestDeblendBackground:
    @pytest.mark.parametrize('bg', ['white', 'black', '#1a6b3c'])
    def test_roundtrip_reconstruction(self, bg):
        """Flatten -> deblend -> re-flatten must reproduce the composite."""
        flat = flatten_rgba(_blob_rgba(), background=bg)
        rec = deblend_background(flat, background=bg, tol=0.0)
        assert rec.shape[-1] == 4
        again = flatten_rgba(rec, background=bg)
        assert np.allclose(again, flat, atol=1e-9)

    def test_minimal_alpha_special_cases(self):
        rng = np.random.default_rng(7)
        rgb = rng.random((20, 20, 3))
        a_black = deblend_background(rgb, background='black', tol=0.0)[..., 3]
        assert np.allclose(a_black, rgb.max(axis=-1))
        a_white = deblend_background(rgb, background='white', tol=0.0)[..., 3]
        assert np.allclose(a_white, (1.0 - rgb).max(axis=-1))

    def test_alpha_not_more_than_original(self):
        """Minimal-alpha solution can't exceed the alpha that made the blend."""
        src = _blob_rgba()
        flat = flatten_rgba(src, background='#406080')
        rec = deblend_background(flat, background='#406080', tol=0.0)
        assert (rec[..., 3] <= src[..., 3] + 1e-9).all()

    def test_dark_features_survive_white_background(self):
        """Near-black foreground on white must come back opaque (luma fails here)."""
        src = np.zeros((10, 10, 4))
        src[3:7, 3:7] = [0.05, 0.05, 0.05, 1.0]   # dark opaque square
        flat = flatten_rgba(src, background='white')
        rec = deblend_background(flat, background='white')
        assert rec[5, 5, 3] > 0.9
        assert rec[0, 0, 3] == 0.0

    def test_tol_snaps_quantized_background_to_zero(self):
        """8-bit quantization noise on the background must not leave haze."""
        flat = flatten_rgba(_blob_rgba(), background='white')
        flat8 = np.round(flat * 255) / 255.0                       # quantize
        rec = deblend_background(flat8, background='white', tol=1.5 / 255)
        corner = rec[:5, :5, 3]
        assert corner.max() == 0.0

    def test_gamma_roundtrip(self):
        src = _blob_rgba()
        flat = flatten_rgba(src, background='#d0d0d0', gamma=2.2)
        rec = deblend_background(flat, background='#d0d0d0', gamma=2.2, tol=0.0)
        again = flatten_rgba(rec, background='#d0d0d0', gamma=2.2)
        assert np.allclose(again, flat, atol=1e-7)

    def test_alpha_smooth_runs_and_stays_bounded(self):
        pytest.importorskip('scipy')
        flat = flatten_rgba(_blob_rgba(), background='white')
        rec = deblend_background(flat, background='white', alpha_smooth=1.5)
        assert rec[..., 3].min() >= 0.0 and rec[..., 3].max() <= 1.0

    def test_fill_color_under_transparent(self):
        flat = flatten_rgba(_blob_rgba(), background='white')
        rec = deblend_background(flat, background='white', fill='black')
        transparent = rec[..., 3] == 0.0
        assert transparent.any()
        assert np.allclose(rec[..., :3][transparent], 0.0)

    def test_pipes_into_make_transparent_cutout(self):
        """Deblended RGBA can be trimmed further via existing_alpha='multiply'."""
        flat = flatten_rgba(_blob_rgba(), background='white')
        rec = deblend_background(flat, background='white')
        out = make_transparent_cutout(rec, existing_alpha='multiply',
                                      alpha_source='dist', sky_color='white',
                                      alpha_lo=20, alpha_hi=90, crop='auto')
        assert out.shape[-1] == 4
        assert out[..., 3].max() > 0.5

    def test_shared_matte_knobs(self):
        """Same threshold/gamma/crop vocabulary as make_transparent_cutout."""
        flat = flatten_rgba(_blob_rgba(), background='white')
        plain = deblend_background(flat, background='white')
        # alpha_gamma alone: plain power on the physical alpha.
        tight = deblend_background(flat, background='white', alpha_gamma=1.4)
        mid = (plain[..., 3] > 0.1) & (plain[..., 3] < 0.9)
        assert mid.any()
        assert tight[..., 3][mid].mean() < plain[..., 3][mid].mean()
        # Percentile ramp: alpha_hi threshold saturates the body to opaque.
        ramped = deblend_background(flat, background='white',
                                    alpha_lo=40, alpha_hi=90, alpha_gamma=0.5)
        assert (ramped[..., 3] == 1.0).sum() > (plain[..., 3] == 1.0).sum()
        # Crop/size passthrough shares the cutout implementation.
        small = deblend_background(flat, background='white', crop='auto',
                                   pad=0.02, size=32)
        assert max(small.shape[:2]) <= 32 and small.shape[-1] == 4
