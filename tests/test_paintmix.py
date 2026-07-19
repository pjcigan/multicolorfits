"""Tests for RYB color space + alpha/background compositing."""

import numpy as np
import pytest

import multicolorfits as mcf
from multicolorfits.core.paintmix import (
    rgb_to_ryb, ryb_to_rgb, hex_to_ryb, ryb_to_hex,
    layer_coverage, composite_over_background, combine_multicolor_alpha,
    ALPHA_SPACES,
)


# Sugita & Takahashi 2017, Table 1 (RGB cube corners -> RYB)
RGB_TO_RYB_CORNERS = {
    'red': ([1, 0, 0], [1, 0, 0]),
    'green': ([0, 1, 0], [0, 1, 1]),
    'blue': ([0, 0, 1], [0, 0, 1]),
    'yellow': ([1, 1, 0], [0, 1, 0]),
    'cyan': ([0, 1, 1], [0, 0.5, 1]),
    'magenta': ([1, 0, 1], [1, 0, 0.5]),
    'black': ([0, 0, 0], [1, 1, 1]),
    'white': ([1, 1, 1], [0, 0, 0]),
}


class TestRYBConversion:
    @pytest.mark.parametrize('name', list(RGB_TO_RYB_CORNERS))
    def test_rgb_to_ryb_corners_match_paper(self, name):
        rgb, exp = RGB_TO_RYB_CORNERS[name]
        got = rgb_to_ryb(np.array(rgb, dtype=float))
        assert np.allclose(got, exp, atol=1e-6), '%s -> %s (want %s)' % (name, got, exp)

    def test_roundtrip_is_reversible(self):
        rng = np.random.default_rng(1)
        X = rng.random((2000, 3))
        back = ryb_to_rgb(rgb_to_ryb(X))
        assert np.max(np.abs(back - X)) < 1e-6

    def test_vectorized_shape_preserved(self):
        img = np.random.default_rng(0).random((7, 5, 3))
        assert rgb_to_ryb(img).shape == (7, 5, 3)
        assert ryb_to_rgb(img).shape == (7, 5, 3)

    def test_hex_helpers(self):
        assert np.allclose(hex_to_ryb('#FFFF00'), [0, 1, 0], atol=1e-6)  # yellow
        assert ryb_to_hex([0, 1, 0]).upper() == '#FFFF00'                # -> yellow rgb

    def test_paint_yellow_plus_blue_is_green(self):
        """The defining feature of RYB: yellow + blue mix to green, not white."""
        ry = rgb_to_ryb(np.array([1., 1, 0]))
        rb = rgb_to_ryb(np.array([0., 0, 1]))
        mixed = ryb_to_rgb(np.clip(ry + rb, 0, 1))
        assert np.allclose(mixed, [0, 1, 0], atol=1e-6)


class TestImageValue:
    def test_value_is_max_channel(self):
        img = np.array([[[0.2, 0.5, 0.1]]])
        assert layer_coverage(img)[0, 0] == pytest.approx(0.5)

    def test_black_is_zero_white_is_one(self):
        assert layer_coverage(np.zeros((2, 2, 3)))[0, 0] == 0.0
        assert layer_coverage(np.ones((2, 2, 3)))[0, 0] == 1.0


class TestCompositeOverBackground:
    def test_zero_signal_becomes_background(self):
        black = np.zeros((3, 3, 3))
        out = composite_over_background(black, background='white')
        assert np.allclose(out, 1.0)  # all white

    def test_bright_colors_preserved_on_any_background(self):
        # A fully-bright pure color (value=1) should be unchanged over any bg.
        red = np.zeros((2, 2, 3)); red[..., 0] = 1.0
        for bg in ('white', 'black', '#204060'):
            out = composite_over_background(red, background=bg)
            assert np.allclose(out, red, atol=1e-6)

    def test_faint_signal_fades_toward_background(self):
        red_half = np.zeros((1, 1, 3)); red_half[..., 0] = 0.5
        out = composite_over_background(red_half, background='white')
        # value=0.5 => 50% toward white: (1, 0.5, 0.5)
        assert np.allclose(out[0, 0], [1.0, 0.5, 0.5], atol=1e-6)

    def test_transparent_output_has_alpha(self):
        red_half = np.zeros((1, 1, 3)); red_half[..., 0] = 0.5
        rgba = composite_over_background(red_half, background=None)
        assert rgba.shape == (1, 1, 4)
        assert rgba[0, 0, 3] == pytest.approx(0.5)          # alpha = value
        assert np.allclose(rgba[0, 0, :3], [1, 0, 0], atol=1e-6)  # straight color = pure red


class TestCombineMulticolorAlpha:
    def _layers(self):
        h = w = 4
        red = np.zeros((h, w, 3)); red[..., 0] = np.linspace(0, 1, w)[None, :]
        grn = np.zeros((h, w, 3)); grn[..., 1] = np.linspace(0, 1, w)[None, :]
        return [red, grn]

    def test_white_background_zero_is_white(self):
        out = combine_multicolor_alpha(self._layers(), background='white',
                                       mode='rgb', gamma=1)
        assert np.allclose(out[0, 0], 1.0)  # col 0 has no signal

    def test_black_background_matches_additive(self):
        layers = self._layers()
        out = combine_multicolor_alpha(layers, background='black', mode='rgb', gamma=1)
        assert np.allclose(out, np.clip(sum(layers), 0, 1), atol=1e-6)

    def test_transparent_output_shape_and_alpha(self):
        rgba = combine_multicolor_alpha(self._layers(), background=None,
                                        mode='rgb', gamma=1)
        assert rgba.shape[-1] == 4
        assert rgba[0, 0, 3] == pytest.approx(0.0)   # no signal -> transparent

    def test_ryb_paint_yellow_blue_is_green(self):
        h = w = 2
        yel = np.zeros((h, w, 3)); yel[..., 0] = 1; yel[..., 1] = 1
        blu = np.zeros((h, w, 3)); blu[..., 2] = 1
        paint = combine_multicolor_alpha([yel, blu], background='white',
                                         mode='ryb', gamma=1)
        assert np.allclose(paint[0, 0], [0, 1, 0], atol=0.02)

    def test_ryb_gamma_darkens_midtones(self):
        # Session-like gamma-space layers (I**g * C**g).  Display-referred mix
        # alone cancels gamma; relative tone must still respond.
        g_lo, g_hi = 1.5, 3.0
        yel = np.ones((4, 4, 3)) * (0.6 ** g_lo) * (np.array([1.0, 0.8, 0.0]) ** g_lo)
        blu = np.ones((4, 4, 3)) * (0.6 ** g_lo) * (np.array([0.2, 0.4, 1.0]) ** g_lo)
        # Rebuild at each gamma as the session does
        def layers(g):
            return [
                np.ones((4, 4, 3)) * (0.6 ** g) * (np.array([1.0, 0.8, 0.0]) ** g),
                np.ones((4, 4, 3)) * (0.6 ** g) * (np.array([0.2, 0.4, 1.0]) ** g),
            ]
        bright = combine_multicolor_alpha(layers(g_lo), background='black',
                                          mode='ryb', gamma=g_lo)
        dark = combine_multicolor_alpha(layers(g_hi), background='black',
                                        mode='ryb', gamma=g_hi)
        assert dark.mean() < bright.mean() - 0.02

    def test_rgb_additive_yellow_blue_is_white(self):
        h = w = 2
        yel = np.zeros((h, w, 3)); yel[..., 0] = 1; yel[..., 1] = 1
        blu = np.zeros((h, w, 3)); blu[..., 2] = 1
        add = combine_multicolor_alpha([yel, blu], background='white',
                                       mode='rgb', gamma=1)
        assert np.allclose(add[0, 0], [1, 1, 1], atol=1e-6)

    def test_bad_mode_raises(self):
        with pytest.raises(ValueError):
            combine_multicolor_alpha(self._layers(), mode='bogus')

    def test_empty_raises(self):
        with pytest.raises(ValueError):
            combine_multicolor_alpha([], background='white')

    def test_exported_at_top_level(self):
        for name in ('rgb_to_ryb', 'ryb_to_rgb',
                     'combine_multicolor_alpha', 'composite_over_background'):
            assert hasattr(mcf, name)

    def test_alpha_spaces_includes_cmyk(self):
        assert 'cmyk' in ALPHA_SPACES
