"""Tests for CMYK color space + signed-data stretches."""

import numpy as np
import pytest

import multicolorfits as mcf
from multicolorfits.core.paintmix import (
    rgb_to_cmyk, cmyk_to_rgb, hex_to_cmyk, cmyk_to_hex,
    combine_multicolor_alpha,
)
from multicolorfits.core.scaling import make_norm, rescale_image, SIGNED_STRETCHES


# Standard RGB cube corners -> CMYK (device-independent, no UCR)
RGB_TO_CMYK_CORNERS = {
    'red': ([1, 0, 0], [0, 1, 1, 0]),
    'green': ([0, 1, 0], [1, 0, 1, 0]),
    'blue': ([0, 0, 1], [1, 1, 0, 0]),
    'yellow': ([1, 1, 0], [0, 0, 1, 0]),
    'cyan': ([0, 1, 1], [1, 0, 0, 0]),
    'magenta': ([1, 0, 1], [0, 1, 0, 0]),
    'white': ([1, 1, 1], [0, 0, 0, 0]),
    'black': ([0, 0, 0], [0, 0, 0, 1]),
}


class TestCMYKConversion:
    @pytest.mark.parametrize('name', list(RGB_TO_CMYK_CORNERS))
    def test_rgb_to_cmyk_corners(self, name):
        rgb, exp = RGB_TO_CMYK_CORNERS[name]
        got = rgb_to_cmyk(np.array(rgb, dtype=float))
        assert np.allclose(got, exp, atol=1e-6), '%s -> %s (want %s)' % (name, got, exp)

    def test_roundtrip_is_reversible(self):
        rng = np.random.default_rng(2)
        X = rng.random((2000, 3))
        back = cmyk_to_rgb(rgb_to_cmyk(X))
        assert np.max(np.abs(back - X)) < 1e-6

    def test_vectorized_shape_preserved(self):
        img = np.random.default_rng(0).random((7, 5, 3))
        assert rgb_to_cmyk(img).shape == (7, 5, 4)
        cmyk = np.random.default_rng(1).random((7, 5, 4))
        assert cmyk_to_rgb(cmyk).shape == (7, 5, 3)

    def test_hex_helpers(self):
        assert np.allclose(hex_to_cmyk('#FFFF00'), [0, 0, 1, 0], atol=1e-6)
        assert cmyk_to_hex([0, 0, 1, 0]).upper() == '#FFFF00'

    def test_cmyk_yellow_plus_cyan_is_green(self):
        y = rgb_to_cmyk(np.array([1., 1, 0]))
        c = rgb_to_cmyk(np.array([0., 1, 1]))
        mixed = cmyk_to_rgb(np.clip(y + c, 0, 1))
        assert np.allclose(mixed, [0, 1, 0], atol=1e-6)

    def test_cmyk_yellow_plus_blue_is_black(self):
        """CMYK differs from RYB: full yellow + blue inks -> black."""
        y = rgb_to_cmyk(np.array([1., 1, 0]))
        b = rgb_to_cmyk(np.array([0., 0, 1]))
        mixed = cmyk_to_rgb(np.clip(y + b, 0, 1))
        assert np.allclose(mixed, [0, 0, 0], atol=1e-6)


class TestCombineMulticolorAlphaCMYK:
    def test_cmyk_yellow_cyan_is_green(self):
        h = w = 2
        yel = np.zeros((h, w, 3)); yel[..., 0] = 1; yel[..., 1] = 1
        cyn = np.zeros((h, w, 3)); cyn[..., 1] = 1; cyn[..., 2] = 1
        out = combine_multicolor_alpha([yel, cyn], background='white',
                                       mode='cmyk', gamma=1)
        assert np.allclose(out[0, 0], [0, 1, 0], atol=0.02)

    def test_exported_at_top_level(self):
        for name in ('rgb_to_cmyk', 'cmyk_to_rgb',
                     'hex_to_cmyk', 'cmyk_to_hex'):
            assert hasattr(mcf, name)


class TestSignedStretches:
    def test_signed_stretches_tuple(self):
        assert 'symlog' in SIGNED_STRETCHES
        assert 'symmetric_log' in SIGNED_STRETCHES

    def test_rescale_image_symlog_signed(self):
        data = np.linspace(-100, 100, 201)
        out = rescale_image(data, rescalefn='symlog', vmin=-100, vmax=100, a=5)
        assert out.shape == data.shape
        assert out[100] == pytest.approx(0.5, abs=0.05)
        assert out[0] < out[100] < out[-1]

    def test_rescale_image_symlog_nan_safe(self):
        img = np.array([[1., np.nan, -1.], [0., 2., -2.]])
        out = rescale_image(img, rescalefn='symlog', vmin=-10, vmax=10)
        assert np.all(np.isfinite(out))

    def test_make_norm_symlog_returns_symlognorm(self):
        from matplotlib.colors import SymLogNorm
        norm = make_norm('symlog', vmin=-10, vmax=10, a=0.1)
        assert isinstance(norm, SymLogNorm)

    def test_make_norm_symmetric_log_optional(self):
        try:
            out = rescale_image(np.linspace(-10, 10, 50),
                                rescalefn='symmetric_log', vmin=-10, vmax=10)
            assert np.all((out >= 0) & (out <= 1))
        except ImportError as exc:
            assert 'pysymlog' in str(exc)
            pytest.skip('pysymlog not installed')

    def test_greyRGBize_symlog(self, test_data):
        grey = mcf.to_grey_rgb(test_data, rescalefn='symlog',
                                    scaletype='abs', min_max=[-1., 1.])
        assert grey.shape[-1] == 3
        assert np.nanmax(grey) <= 1.0

    def test_make_norm_exported(self):
        assert hasattr(mcf, 'make_norm')
        assert hasattr(mcf, 'SIGNED_STRETCHES')
