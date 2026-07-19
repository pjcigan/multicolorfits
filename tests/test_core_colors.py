"""Color math: hex/rgb/hsv conversions and vectorized round trips."""

import colorsys

import numpy as np
import pytest

import multicolorfits as mcf


class TestHexRgb:
    def test_hex_to_rgb_white(self):
        assert mcf.hex_to_rgb('#FFFFFF') == (255, 255, 255)

    def test_hex_to_rgb_no_hash(self):
        assert mcf.hex_to_rgb('FF0000') == (255, 0, 0)

    def test_rgb_to_hex(self):
        assert mcf.rgb_to_hex((255, 255, 255)).lower() == '#ffffff'
        assert mcf.rgb_to_hex((193, 27, 23)).lower() == '#c11b17'

    def test_roundtrip(self):
        for hexstr in ['#C11B17', '#00FF7F', '#123456', '#000000', '#FFFFFF']:
            assert mcf.rgb_to_hex(mcf.hex_to_rgb(hexstr)).lower() == hexstr.lower()

    def test_hexinv(self):
        assert mcf.hex_complement('#FF0000').upper() == '#00FFFF'
        assert mcf.hex_complement('#000000').upper() == '#FFFFFF'
        # Involution: inverse of inverse is identity
        assert mcf.hex_complement(mcf.hex_complement('#C11B17')).upper() == '#C11B17'


class TestHexHsv:
    def test_hex_to_hsv_matches_colorsys(self):
        for hexstr in ['#C11B17', '#00FF7F', '#123456', '#FFFFFF']:
            r, g, b = np.array(mcf.hex_to_rgb(hexstr)) / 255.
            expected = colorsys.rgb_to_hsv(r, g, b)
            result = mcf.hex_to_hsv(hexstr)
            assert np.allclose(result, expected, atol=1e-8)


class TestVectorizedHsv:
    def test_roundtrip_random(self):
        rng = np.random.default_rng(0)
        rgb = rng.uniform(0, 1, size=(32, 32, 3)).astype(np.float32)
        rgb_back = mcf.hsv_to_rgb(mcf.rgb_to_hsv(rgb))
        assert np.allclose(rgb, rgb_back, atol=1e-3)

    def test_matches_colorsys_pixelwise(self):
        rng = np.random.default_rng(1)
        rgb = rng.uniform(0, 1, size=(4, 4, 3)).astype(np.float64)
        hsv = mcf.rgb_to_hsv(rgb)
        for i in range(4):
            for j in range(4):
                expected = colorsys.rgb_to_hsv(*rgb[i, j])
                assert np.allclose(hsv[i, j], expected, atol=1e-3)

    def test_grey_has_zero_saturation(self):
        grey = np.ones((8, 8, 3)) * 0.5
        hsv = mcf.rgb_to_hsv(grey)
        assert np.allclose(hsv[..., 1], 0)
        assert np.allclose(hsv[..., 2], 0.5)


class TestToHex:
    def test_named_color(self):
        assert mcf.to_hex('white').lower() == '#ffffff'
        assert mcf.to_hex('k').lower() == '#000000'
