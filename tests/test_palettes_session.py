"""Tests for palettes, session palette picker, and hue rotation."""

import pytest

import multicolorfits as mcf
from multicolorfits.palettes import (
    get_palette, list_palettes, colors_for_hue_pattern,
    rotate_colors_from_base, is_palette_name,
)

from conftest import make_test_header, make_test_data


def _session(n_loaded=3):
    s = mcf.McfSession()
    for i in range(n_loaded):
        s.panels[i].set_data(make_test_data(32, 32), make_test_header(32, 32))
    return s


class TestHuePatterns:
    def test_triad_three_layers(self):
        c = colors_for_hue_pattern('triad', 3)
        assert len(c) == 3
        assert len(set(c)) == 3

    def test_complement_two_layers(self):
        c = colors_for_hue_pattern('complement', 2)
        assert len(c) == 2
        assert c[0] != c[1]

    def test_rotate_full_circle(self):
        base = colors_for_hue_pattern('triad', 3)
        full = rotate_colors_from_base(base, 360.)
        assert full == base

    def test_is_palette_name(self):
        assert is_palette_name('perceptual')
        assert is_palette_name('triad')
        assert is_palette_name('pob')
        assert not is_palette_name('nope')


class TestApplyPalette:
    def test_apply_curated_palette(self):
        s = _session(3)
        res = s.apply_palette('pob')
        assert res['colors'] == get_palette('pob')
        assert s.panel_colors() == get_palette('pob')

    def test_apply_only_loaded_panels(self):
        s = _session(2)
        res = s.apply_palette('rgb')
        assert len(res['colors']) == 2
        assert res['indices'] == [0, 1]
        # unloaded panels keep their default color
        assert s.panels[3].color == PanelDefault().color

    def test_perceptual_palette_length(self):
        s = _session(4)
        res = s.apply_palette('perceptual')
        assert len(res['colors']) == 4
        assert all(c.startswith('#') and len(c) == 7 for c in res['colors'])

    def test_short_palette_is_padded(self):
        s = _session(3)  # tealorange has only 2 colors
        res = s.apply_palette('tealorange')
        assert len(res['colors']) == 3
        assert res['colors'][:2] == get_palette('tealorange')

    def test_palette_colors_for_no_panels(self):
        s = mcf.McfSession()
        assert s.apply_palette('rgb')['colors'] == []

    def test_apply_hue_pattern(self):
        s = _session(3)
        res = s.apply_palette('triad')
        assert len(res['colors']) == 3

    def test_rotate_panel_hues(self):
        s = _session(2)
        s.apply_palette('complement')
        base = s.panel_colors()
        s.rotate_panel_hues(30., base_colors=base)
        assert s.panel_colors() != base
        s.rotate_panel_hues(0., base_colors=base)
        assert s.panel_colors() == base

    def test_all_curated_palettes_apply(self):
        for name in list_palettes():
            s = _session(2)
            res = s.apply_palette(name)
            assert len(res['colors']) == 2


class TestColorblindReport:
    def test_single_panel_ok(self):
        s = _session(1)
        rep = s.colorblind_report()
        assert rep['ok'] is True
        assert set(rep['kinds']) == {'deuteranopia', 'protanopia', 'tritanopia'}

    def test_well_separated_palette_passes(self):
        s = _session(3)
        s.apply_palette('ryb')  # red/yellow/blue: distinct under all CVD types
        rep = s.colorblind_report()
        assert rep['ok'] is True

    def test_looser_threshold_passes_more(self):
        s = _session(3)
        s.apply_palette('pob')  # borderline under tritanopia at the default 25
        strict = s.colorblind_report(min_distance=25.)
        loose = s.colorblind_report(min_distance=15.)
        assert strict['ok'] is False
        assert loose['ok'] is True

    def test_confusable_palette_flagged(self):
        s = _session(2)
        # two near-identical reds -> confusable under CVD
        s.panels[0].color = '#FF0000'
        s.panels[1].color = '#FA0505'
        rep = s.colorblind_report()
        assert rep['ok'] is False
        # failing pair reported as 1-based panel numbers
        failed = any(rep['kinds'][k]['failures'] for k in rep['kinds'])
        assert failed
        for k in rep['kinds']:
            for a, b, d in rep['kinds'][k]['failures']:
                assert (a, b) == (1, 2)

    def test_indices_are_active_panels(self):
        s = mcf.McfSession()
        s.panels[1].set_data(make_test_data(16, 16), make_test_header(16, 16))
        s.panels[2].set_data(make_test_data(16, 16), make_test_header(16, 16))
        rep = s.colorblind_report()
        assert rep['indices'] == [1, 2]


def PanelDefault():
    from multicolorfits.session import PanelState
    return PanelState()
