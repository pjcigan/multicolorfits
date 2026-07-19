"""Tests for the composite-mode-aware color-combination swatch and band labels."""
import numpy as np
import pytest

import multicolorfits as mcf
from multicolorfits.core.swatch import (
    combine_colorized_layers, swatch_layout, combo_swatch,
)
from conftest import make_test_data, make_test_header


def _blob_session(colors, labels=None):
    s = mcf.McfSession()
    for i, col in enumerate(colors):
        s.panels[i].set_data(make_test_data(seed=i), make_test_header())
        s.panels[i].color = col
        if labels:
            s.panels[i].label = labels[i]
    return s


class TestSwatchLayout:
    @pytest.mark.parametrize('n', [1, 2, 3, 4, 5, 6])
    def test_count(self, n):
        assert len(swatch_layout(n, 300)) == n

    def test_empty(self):
        assert swatch_layout(0, 300) == []

    def test_circles_stay_in_frame(self):
        size = 300
        for n in range(1, 7):
            for (cx, cy, r) in swatch_layout(n, size):
                assert -1 <= cx - r and cx + r <= size + 1
                assert -1 <= cy - r and cy + r <= size + 1


class TestComboSwatch:
    def test_none_when_empty(self):
        assert combo_swatch([]) is None

    def test_shape_and_alpha(self):
        d = combo_swatch(['#FF0000', '#00FF00', '#0000FF'], size=128)
        assert d['rgba'].shape == (128, 128, 4)
        assert len(d['circles']) == 3
        # Corner is outside every circle -> transparent; centre is covered.
        assert d['rgba'][0, 0, 3] == pytest.approx(0.0, abs=1e-6)
        assert d['rgba'][64, 64, 3] > 0.5

    def test_labels_carried(self):
        d = combo_swatch(['#FF0000', '#0000FF'], labels=['Ha', 'CO'], size=64)
        assert [c['label'] for c in d['circles']] == ['Ha', 'CO']

    def test_ryb_yellow_plus_blue_is_green(self):
        # The whole point: paint mixing must be honest in the overlap region.
        d = combo_swatch(['#FFFF00', '#0000FF'], mode='ryb', background='white', size=200)
        r, g, b, a = d['rgba'][100, 100]
        assert a > 0.5
        assert g > r and g > b        # green dominates the yellow+blue overlap
        assert g > 0.3

    def test_rgb_yellow_plus_blue_is_bright(self):
        # Additive: yellow + blue -> near white in the overlap.
        d = combo_swatch(['#FFFF00', '#0000FF'], mode='rgb', background='black', size=200)
        r, g, b, a = d['rgba'][100, 100]
        assert a > 0.5
        assert min(r, g, b) > 0.5     # all channels lit -> whitish, not green

    def test_transparent_corner_regardless_of_background(self):
        for bg in ('black', 'white', 'transparent', '#204060'):
            d = combo_swatch(['#FF0000', '#00FF00', '#0000FF'], background=bg, size=64)
            assert d['rgba'][0, 0, 3] == pytest.approx(0.0, abs=1e-6)


class TestCombineColorizedLayers:
    def test_matches_render_combined_rgb(self):
        # combine_colorized_layers is the shared dispatch; the session must agree.
        s = _blob_session(['#FF0000', '#00FF00'])
        s.compose.combine_mode = 'rgb'
        s.compose.combine_background = 'black'
        img = s.render_combined()
        assert img.ndim == 3 and img.shape[-1] == 3

    def test_transparent_gives_rgba(self):
        out = combine_colorized_layers(
            [np.ones((4, 4, 3)) * 0.5], mode='rgb', background='transparent')
        assert out.shape[-1] == 4

    def test_white_background_flattens(self):
        out = combine_colorized_layers(
            [np.zeros((4, 4, 3))], mode='rgb', background='white')
        assert out.shape[-1] == 3
        assert np.allclose(out, 1.0)   # zero signal over white -> white


class TestSessionSwatch:
    def test_render_combo_swatch(self):
        s = _blob_session(['#FF0000', '#00FF00', '#0000FF'], labels=['a', 'b', 'c'])
        d = s.render_combo_swatch(size=100)
        assert d['rgba'].shape == (100, 100, 4)
        assert [c['label'] for c in d['circles']] == ['a', 'b', 'c']

    def test_render_combo_swatch_none_when_empty(self):
        assert mcf.McfSession().render_combo_swatch() is None

    def test_swatch_tracks_mode(self):
        s = _blob_session(['#FFFF00', '#0000FF'])
        s.compose.combine_background = 'white'
        s.compose.combine_mode = 'ryb'
        ryb = s.render_combo_swatch(size=200)['rgba'][100, 100]
        s.compose.combine_mode = 'rgb'
        s.compose.combine_background = 'black'
        rgb = s.render_combo_swatch(size=200)['rgba'][100, 100]
        # RYB overlap is green-dominant; RGB overlap is whitish.
        assert ryb[1] > ryb[2]
        assert min(rgb[:3]) > 0.5


class TestStateSerialization:
    def test_roundtrip(self):
        s = _blob_session(['#FF0000'])
        s.compose.show_combo_swatch = True
        s.compose.combo_swatch_loc = 'upper left'
        s.compose.combo_swatch_labels = True
        s.compose.show_band_labels = True
        s.compose.band_labels_loc = 'lower right'
        d = s.compose.to_dict()
        assert d['show_combo_swatch'] is True
        assert d['combo_swatch_loc'] == 'upper left'
        assert d['combo_swatch_labels'] is True
        assert d['show_band_labels'] is True
        assert d['band_labels_loc'] == 'lower right'


class TestExportScript:
    def test_includes_swatch_and_band_labels(self, tmp_path):
        import astropy.io.fits as pyfits
        fp = tmp_path / 'a.fits'
        pyfits.writeto(str(fp), make_test_data(), make_test_header())
        s = mcf.McfSession()
        s.panels[0].load_fits(str(fp))
        s.panels[0].color = '#FF0000'
        s.panels[0].label = 'Red'
        s.compose.show_combo_swatch = True
        s.compose.show_band_labels = True
        script = s.export_script()
        assert 'mcf.combo_swatch(' in script
        assert 'band_labels=' in script
        assert 'swatch=_swatch' in script
        compile(script, '<export>', 'exec')  # valid python

    def test_no_annotations_when_off(self, tmp_path):
        import astropy.io.fits as pyfits
        fp = tmp_path / 'a.fits'
        pyfits.writeto(str(fp), make_test_data(), make_test_header())
        s = mcf.McfSession()
        s.panels[0].load_fits(str(fp))
        script = s.export_script()
        assert 'mcf.combo_swatch(' not in script
        assert 'band_labels=' not in script
