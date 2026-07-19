"""Tests for grid-mismatch detection and in-session reprojection/alignment."""

import numpy as np
import pytest

import multicolorfits as mcf
from multicolorfits.session import reproject_available

from conftest import make_test_header, make_test_data

requires_reproject = pytest.mark.skipif(
    not reproject_available(), reason='reproject package not installed')


def _mismatched_header(nx=100, ny=80):
    h = make_test_header(nx, ny)
    h['CRVAL1'] = 150.02  # shifted center + different shape => genuine mismatch
    return h


class TestGridReport:
    def test_single_panel_is_aligned(self):
        s = mcf.McfSession()
        s.panels[0].set_data(make_test_data(), make_test_header())
        rep = s.grid_report()
        assert rep['aligned'] is True
        assert rep['reference'] == 0
        assert rep['mismatched'] == []

    def test_matching_grids_aligned(self):
        s = mcf.McfSession()
        s.panels[0].set_data(make_test_data(64, 64), make_test_header(64, 64))
        s.panels[1].set_data(make_test_data(64, 64), make_test_header(64, 64))
        assert s.grid_report()['aligned'] is True

    def test_shape_mismatch_flagged(self):
        s = mcf.McfSession()
        s.panels[0].set_data(make_test_data(64, 64), make_test_header(64, 64))
        s.panels[1].set_data(make_test_data(80, 100), _mismatched_header())
        rep = s.grid_report()
        assert rep['aligned'] is False
        assert 1 in rep['mismatched']

    def test_wcs_mismatch_same_shape_flagged(self):
        s = mcf.McfSession()
        s.panels[0].set_data(make_test_data(64, 64), make_test_header(64, 64))
        shifted = make_test_header(64, 64)
        shifted['CRVAL1'] = 151.0  # same shape, different pointing
        s.panels[1].set_data(make_test_data(64, 64), shifted)
        rep = s.grid_report()
        assert rep['aligned'] is False
        assert 1 in rep['mismatched']

    def test_reference_is_first_loaded(self):
        s = mcf.McfSession()
        s.panels[2].set_data(make_test_data(), make_test_header())
        s.panels[3].set_data(make_test_data(), make_test_header())
        assert s.grid_report()['reference'] == 2


@requires_reproject
class TestAlignPanels:
    def _mismatched_session(self):
        s = mcf.McfSession()
        s.panels[0].set_data(make_test_data(64, 64), make_test_header(64, 64))
        s.panels[0].color = '#FF0000'
        s.panels[1].set_data(make_test_data(80, 100), _mismatched_header())
        s.panels[1].color = '#00FF00'
        return s

    def test_align_to_reference_makes_shapes_match(self):
        s = self._mismatched_session()
        assert not s.grid_report()['aligned']
        res = s.align_panels(target='reference')
        assert res['ok'] is True
        assert set(res['changed']) == {0, 1}
        assert s.grid_report()['aligned'] is True
        shapes = {p.data.shape for p in s.active_panels()}
        assert len(shapes) == 1

    def test_align_enables_combine(self):
        s = self._mismatched_session()
        s.align_panels(target='reference')
        combined = s.render_combined()
        assert combined.ndim == 3 and combined.shape[2] == 3

    def test_align_to_galactic_sets_frame(self):
        s = self._mismatched_session()
        s.align_panels(target='galactic')
        assert s.grid_report()['aligned'] is True
        for p in s.active_panels():
            assert p.header.frame == 'galactic'

    def test_align_to_icrs_sets_frame(self):
        s = self._mismatched_session()
        s.align_panels(target='icrs')
        assert s.grid_report()['aligned'] is True
        for p in s.active_panels():
            assert p.header.frame == 'icrs'

    def test_align_preserves_display_settings(self):
        s = self._mismatched_session()
        s.panels[0].stretch = 'asinh'
        vmin0, vmax0 = s.panels[0].vmin, s.panels[0].vmax
        s.align_panels(target='reference')
        assert s.panels[0].stretch == 'asinh'
        assert s.panels[0].color == '#FF0000'
        assert s.panels[0].vmin == vmin0 and s.panels[0].vmax == vmax0

    def test_align_custom_reference(self):
        s = self._mismatched_session()
        res = s.align_panels(target='reference', reference=1)
        assert res['reference'] == 1
        # all layers now on panel 1's (80x100) grid
        assert s.panels[0].data.shape == (80, 100)

    def test_bad_reference_raises(self):
        s = self._mismatched_session()
        with pytest.raises(ValueError):
            s.align_panels(reference=3)  # panel 3 not loaded

    def test_fewer_than_two_is_noop(self):
        s = mcf.McfSession()
        s.panels[0].set_data(make_test_data(), make_test_header())
        res = s.align_panels()
        assert res['ok'] is True and res['changed'] == []


class TestReprojectedFlag:
    def test_fresh_panel_not_reprojected(self):
        s = mcf.McfSession()
        s.panels[0].set_data(make_test_data(), make_test_header())
        assert s.panels[0].reprojected is False

    @requires_reproject
    def test_flag_set_after_align(self):
        s = mcf.McfSession()
        s.panels[0].set_data(make_test_data(64, 64), make_test_header(64, 64))
        s.panels[1].set_data(make_test_data(80, 100), _mismatched_header())
        s.align_panels(target='reference')
        assert all(p.reprojected for p in s.active_panels())


@requires_reproject
class TestFullFidelityScript:
    def _real_file_session(self, tmp_path):
        import astropy.io.fits as pyfits
        s = mcf.McfSession()
        f0 = tmp_path / 'a.fits'
        f1 = tmp_path / 'b.fits'
        pyfits.writeto(str(f0), make_test_data(64, 64), make_test_header(64, 64))
        pyfits.writeto(str(f1), make_test_data(80, 100), _mismatched_header())
        s.panels[0].load_fits(str(f0)); s.panels[0].color = '#FF0000'
        s.panels[1].load_fits(str(f1)); s.panels[1].color = '#00FF00'
        return s

    def test_export_embeds_common_header_and_reproject(self, tmp_path):
        s = self._real_file_session(tmp_path)
        s.align_panels(target='reference')
        script = s.export_script()
        assert '_common_hdr' in script
        assert 'reproject_image' in script
        compile(script, '<x>', 'exec')

    def test_no_reproject_when_grids_match(self):
        s = mcf.McfSession()
        s.panels[0].set_data(make_test_data(64, 64), make_test_header(64, 64))
        s.panels[1].set_data(make_test_data(64, 64), make_test_header(64, 64))
        script = s.export_script()
        assert '_common_hdr' not in script
        assert 'reproject_image' not in script

    def test_exported_script_reproduces_aligned_combined(self, tmp_path):
        s = self._real_file_session(tmp_path)
        s.align_panels(target='reference')
        expected = s.render_combined()
        script = s.export_script()
        lines = [ln for ln in script.split('\n')
                 if not ln.startswith('mcf.plot_combined_rgb')]
        ns = {}
        exec('\n'.join(lines), ns)
        got = ns['combined']
        assert got.shape == expected.shape
        finite = np.isfinite(got) & np.isfinite(expected)
        assert finite.any()
        assert np.allclose(got[finite], expected[finite], atol=1e-4)
