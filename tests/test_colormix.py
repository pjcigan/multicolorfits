"""Tests for experimental color-space compositing (colormix)."""

import colorsys

import numpy as np
import pytest

import multicolorfits as mcf


class TestHslConversions:
    def test_hsl_roundtrip(self):
        rng = np.random.default_rng(0)
        rgb = rng.uniform(0.05, 0.95, size=(16, 16, 3)).astype(np.float64)
        back = mcf.hsl_to_rgb(mcf.rgb_to_hsl(rgb))
        assert np.allclose(rgb, back, atol=1e-3)

    def test_hsl_matches_colorsys_sample(self):
        for r, g, b in [(1, 0, 0), (0, 1, 0), (0.5, 0.5, 0.5)]:
            rgb = np.array([[[r, g, b]]], dtype=np.float64)
            h, s, l = mcf.rgb_to_hsl(rgb)[0, 0]
            # colorsys.rgb_to_hls returns (hue, lightness, saturation)
            hr, lr, sr = colorsys.rgb_to_hls(r, g, b)
            assert np.isclose(h, hr, atol=1e-3)
            assert np.isclose(s, sr, atol=1e-3)
            assert np.isclose(l, lr, atol=1e-3)


class TestLabConversions:
    def test_lab_roundtrip_grey(self):
        grey = np.full((8, 8, 3), 0.5)
        lab = mcf.rgb_to_lab(grey)
        back = mcf.lab_to_rgb(lab)
        assert np.allclose(grey, back, atol=0.05)


class TestCombineColorspace:
    @pytest.fixture
    def colorized_pair(self, test_data):
        g = mcf.to_grey_rgb(test_data, rescalefn='linear')
        return [
            mcf.colorize_image(g, '#FF0000'),
            mcf.colorize_image(g, '#0000FF'),
        ]

    @pytest.mark.parametrize('colorspace', ['lab', 'hsv', 'hsl'])
    @pytest.mark.parametrize('blend', ['screen', 'sum', 'max', 'mean'])
    def test_runs_and_bounded(self, colorized_pair, colorspace, blend):
        out = mcf.combine_multicolor_colorspace(
            colorized_pair, colorspace=colorspace, blend=blend, gamma=2.2)
        assert out.shape == colorized_pair[0].shape
        assert out.min() >= 0 and out.max() <= 1.0 + 1e-6

    def test_lab_differs_from_rgb_sum(self, colorized_pair):
        """Lab compositing should not be identical to classic RGB addition."""
        rgb = mcf.combine_multicolor(colorized_pair, gamma=2.2)
        lab = mcf.combine_multicolor_colorspace(colorized_pair, colorspace='lab', gamma=2.2)
        assert not np.allclose(rgb, lab, atol=0.05)

    def test_red_yellow_hsv_not_green(self):
        """Hue vector sum of red+yellow should land in orange, not green."""
        mixed = mcf.mix_colors_hex(['#FF0000', '#FFFF00'], colorspace='hsv', blend='mean')
        r, g, b = mcf.hex_to_rgb(mixed)
        assert r > g > b  # orange: high R, moderate G, low B

    def test_complementary_hues_neutralize_in_hsv(self):
        """Yellow + blue at equal brightness -> grey (additive light physics)."""
        mixed = mcf.mix_colors_hex(['#FFFF00', '#0000FF'], colorspace='hsv', blend='mean')
        r, g, b = mcf.hex_to_rgb(mixed)
        assert abs(r - g) < 30 and abs(g - b) < 30  # roughly neutral grey


class TestAccessibility:
    def test_greyscale_lab(self, test_data):
        g = mcf.to_grey_rgb(test_data, rescalefn='linear')
        c = mcf.colorize_image(g, '#C11B17')
        grey = mcf.greyscale_image(c, method='lab')
        assert grey.shape == c.shape
        assert np.allclose(grey[..., 0], grey[..., 1])

    def test_colorblind_simulation_bounded(self, test_data):
        g = mcf.to_grey_rgb(test_data, rescalefn='linear')
        c = mcf.colorize_image(g, '#C11B17')
        sim = mcf.simulate_colorblindness(c, kind='deuteranopia')
        assert sim.shape == c.shape
        assert sim.min() >= 0 and sim.max() <= 1.0

    def test_invalid_cvd_kind_raises(self):
        with pytest.raises(ValueError):
            mcf.simulate_colorblindness(np.zeros((2, 2, 3)), kind='bogus')
