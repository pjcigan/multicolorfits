"""Session JSON, cursor sampling, CMYK mode, and GUI preload helpers."""

import json
import numpy as np
import pytest

import multicolorfits as mcf
from multicolorfits.session import COMBINE_MODES


class TestCombineModes:
    def test_cmyk_in_combine_modes(self):
        assert 'cmyk' in COMBINE_MODES
        assert hasattr(mcf, 'COMBINE_MODES')


class TestSessionJson:
    def test_to_state_dict_roundtrip(self, test_fits_file, tmp_path):
        s = mcf.McfSession()
        s.panels[0].load_fits(test_fits_file)
        s.compose.combine_mode = 'lab'
        s.compose.title = 'test'
        path = tmp_path / 'sess.json'
        s.save_state(str(path))
        s2 = mcf.McfSession()
        warnings = s2.load_state(str(path))
        assert warnings == []
        assert s2.compose.combine_mode == 'lab'
        assert s2.compose.title == 'test'
        assert s2.panels[0].in_use

    def test_load_state_from_dict(self, test_fits_file):
        s = mcf.McfSession()
        s.panels[0].load_fits(test_fits_file)
        s.compose.combine_mode = 'cmyk'
        d = s.to_state_dict()
        s2 = mcf.McfSession()
        s2.load_state(d)
        assert s2.compose.combine_mode == 'cmyk'

    def test_load_files(self, test_fits_file):
        s = mcf.McfSession()
        w = s.load_files([test_fits_file], colors=['#FF0000'], labels=['A'])
        assert w == []
        assert s.panels[0].in_use
        assert s.panels[0].color == '#FF0000'
        assert s.panels[0].label == 'A'
        assert s.panels[0].stretch == 'linear'

    def test_load_files_display(self, test_fits_file):
        s = mcf.McfSession()
        w = s.load_files([test_fits_file], stretches='asinh', vmins=0.2, vmaxs=3)
        assert w == []
        assert s.panels[0].stretch == 'asinh'
        assert s.panels[0].vmin == pytest.approx(0.2)
        assert s.panels[0].vmax == pytest.approx(3)

    def test_add_remove_panel(self):
        s = mcf.McfSession(n_panels=2)
        assert len(s.panels) == 2
        idx = s.add_panel(color='#00FF00', label='extra')
        assert idx == 2
        assert len(s.panels) == 3
        assert s.panels[2].color == '#00FF00'
        assert s.panels[2].label == 'extra'
        s.remove_panel(1)
        assert len(s.panels) == 2
        assert s.panels[1].label == 'extra'
        with pytest.raises(ValueError):
            s.set_n_panels(0)
        s.set_n_panels(1)
        assert len(s.panels) == 1
        with pytest.raises(ValueError):
            s.remove_panel(0)

    def test_load_state_grows_panels(self, test_fits_file, tmp_path):
        s = mcf.McfSession(n_panels=5)
        for i in range(5):
            s.panels[i].load_fits(test_fits_file)
            s.panels[i].color = '#%02X0000' % (40 + 40 * i)
        path = tmp_path / 'five.json'
        s.save_state(str(path))
        s2 = mcf.McfSession()  # default 4 slots
        assert len(s2.panels) == 4
        warnings = s2.load_state(str(path))
        assert warnings == []
        assert len(s2.panels) == 5
        assert s2.panels[4].in_use

    def test_load_files_grows_panels(self, test_fits_file):
        s = mcf.McfSession(n_panels=2)
        paths = [test_fits_file] * 5
        w = s.load_files(paths)
        assert w == []
        assert len(s.panels) == 5
        assert all(p.in_use for p in s.panels)


class TestSampleAtPixel:
    def test_sample_returns_layer_value(self, test_fits_file):
        s = mcf.McfSession()
        s.panels[0].load_fits(test_fits_file)
        ny, nx = s.panels[0].data.shape
        info = s.sample_at_pixel(nx // 2, ny // 2)
        assert info['x'] == nx // 2
        assert len(info['layers']) == 1
        assert np.isfinite(info['layers'][0]['value'])

    def test_display_to_data_pixel(self):
        s = mcf.McfSession()
        s.panels[0].set_data(np.arange(100).reshape(10, 10))
        x, y = s.display_to_data_pixel(0, 0, 10, 10)
        assert x == 0
        assert y == 9


class TestSessionReset:
    def test_reset_clears_panels(self, test_fits_file):
        s = mcf.McfSession()
        s.load_files([test_fits_file])
        s.compose.title = 'gone'
        assert s.active_panels()
        s.reset()
        assert not s.active_panels()
        assert s.compose.title == ''
        assert s.compose.gamma == 2.2
