"""Optional skyplothelper overlay shim (multicolorfits[overlays])."""

import importlib

import pytest

import multicolorfits as mcf
from multicolorfits import McfSession, PanelState
from multicolorfits.figures import make_combined_figure
from conftest import make_test_data, make_test_header

_sph = pytest.importorskip('skyplothelper')


class TestOverlaysAvailability:
    def test_overlays_available_is_bool(self):
        assert isinstance(mcf.overlays_available(), bool)
        assert mcf.overlays_available() == mcf.overlays.overlays_available()

    def test_require_overlays_raises_clear_error(self, monkeypatch):
        def _fail_import(name, package=None):
            raise ImportError('nope')

        monkeypatch.setattr(importlib, 'import_module', _fail_import)
        with pytest.raises(ImportError, match=r'\[overlays\]'):
            mcf.overlays.require_overlays()


class TestOverlayHelpers:
    def test_auto_scale_bar_asec(self, test_header):
        length = mcf.overlays._auto_scale_bar_asec(test_header)
        assert 0.5 <= length <= 30.0

    def test_scale_bar_label_formats(self):
        assert mcf.overlays._scale_bar_label(5.0) == '5"'
        assert mcf.overlays._scale_bar_label(120.0) == "2'"

    def test_apply_overlays_warns_when_missing(self, monkeypatch):
        monkeypatch.setattr(mcf.overlays, 'overlays_available', lambda: False)
        import matplotlib
        matplotlib.use('Agg')
        from matplotlib.figure import Figure
        fig = Figure()
        ax = fig.add_subplot(111)
        with pytest.warns(UserWarning, match='not installed'):
            out = mcf.overlays.apply_overlays(ax, compass=True)
        assert out == []

    def test_compose_state_serializes_overlay_fields(self):
        c = mcf.ComposeState(show_compass=True, show_scale_bar=True, scale_bar_asec=10.0)
        d = c.to_dict()
        assert d['show_compass'] is True
        assert d['show_scale_bar'] is True
        assert d['scale_bar_asec'] == 10.0


class TestSkyplothelperIntegration:
    def test_require_overlays_returns_module(self):
        sph = mcf.overlays.require_overlays()
        assert sph.__name__ == 'skyplothelper'

    def test_skyplothelper_proxy_delegates(self):
        proxy = mcf.overlays.SkyplotHelper()
        assert proxy.add_compass is mcf.overlays.require_overlays().add_compass

    def test_add_compass_on_wcs_axes(self, test_header):
        import matplotlib
        matplotlib.use('Agg')
        import astropy.wcs as pywcs
        from matplotlib.figure import Figure
        fig = Figure(figsize=(4, 4))
        ax = fig.add_subplot(111, projection=pywcs.WCS(test_header).celestial)
        artists = mcf.overlays.add_compass(ax, loc='lower left')
        assert artists

    def test_add_beam_from_header(self, test_header):
        import matplotlib
        matplotlib.use('Agg')
        import astropy.wcs as pywcs
        from matplotlib.figure import Figure
        hdr = test_header.copy()
        hdr['BMAJ'] = 5.0 / 3600.0
        hdr['BMIN'] = 3.0 / 3600.0
        hdr['BPA'] = 30.0
        fig = Figure(figsize=(4, 4))
        ax = fig.add_subplot(111, projection=pywcs.WCS(hdr).celestial)
        beam = mcf.overlays.add_beam(ax, hdr, loc='lower left')
        assert beam is not None

    def test_add_scale_bar_from_header(self, test_header):
        import matplotlib
        matplotlib.use('Agg')
        import astropy.wcs as pywcs
        from matplotlib.figure import Figure
        fig = Figure(figsize=(4, 4))
        ax = fig.add_subplot(111, projection=pywcs.WCS(test_header).celestial)
        bar = mcf.overlays.add_scale_bar(ax, test_header, length_asec=5.0)
        assert bar is not None

    def test_offset_coord_wcs_smoke(self, test_header):
        wcs = mcf.overlays.offset_coord_wcs(test_header, [150.0, -30.0])
        assert wcs.wcs.ctype[0].endswith('OFFSET')

    def test_setup_combined_axes_applies_compose_overlays(self, test_header):
        import matplotlib
        matplotlib.use('Agg')
        s = McfSession()
        s.panels[0].set_data(make_test_data(), test_header)
        s.compose.show_compass = True
        s.compose.show_scale_bar = True
        s.compose.scale_bar_asec = 5.0
        fig = make_combined_figure(s)
        assert len(fig.axes) >= 1

    def test_panel_load_with_beam_header(self, tmp_path, test_data, test_header):
        import astropy.io.fits as pyfits
        hdr = test_header.copy()
        hdr['BMAJ'] = 4.0 / 3600.0
        hdr['BMIN'] = 2.0 / 3600.0
        hdr['BPA'] = 0.0
        path = tmp_path / 'beam.fits'
        pyfits.writeto(str(path), test_data, hdr, overwrite=True)
        p = PanelState()
        p.load_fits(str(path))
        s = McfSession()
        s.panels[0] = p
        s.compose.show_beam = True
        fig = make_combined_figure(s)
        assert len(fig.axes) >= 1
